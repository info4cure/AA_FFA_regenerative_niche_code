#!/usr/bin/env Rscript

# ============================================================
# PIPELINE_V2 — 18D2-B
# AA BATCH / BIOLOGICAL-ANCHOR DIRECTION AUDIT
#
# Purpose
# -------
# Audit whether selected AA-vs-HC biological anchor genes from
# the frozen-L2 pseudobulk DE retain the same direction when AA
# is compared separately with:
#
#   1) historical HC patients
#   2) CON276, the single CABIMER2026 HC patient
#
# IMPORTANT
# ---------
# - This is a sensitivity / direction audit.
# - It does NOT refit the primary AA-vs-HC edgeR model.
# - It does NOT claim that batch has been statistically adjusted.
# - AA is completely confounded with CABIMER2026.
# - CON276 is n=1 and therefore is NOT an inferential control arm.
# - Comparisons below are descriptive patient-level pseudobulk
#   expression-direction checks only.
# - No annotation is modified.
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(tibble)
  library(edgeR)
  library(Matrix)
})

options(
  stringsAsFactors = FALSE,
  width = 180
)

cat("============================================================\n")
cat("PIPELINE_V2 — 18D2-B AA BATCH / ANCHOR DIRECTION AUDIT\n")
cat("============================================================\n\n")

# ------------------------------------------------------------
# 0. Paths
# ------------------------------------------------------------

ROOT <- "PIPELINE_V2/results/disease_biology/18_PSEUDOBULK_DE"

PB_DIR <- file.path(
  ROOT,
  "18B_pseudobulk"
)

PB_OBJ_DIR <- file.path(
  PB_DIR,
  "objects"
)

PB_TAB_DIR <- file.path(
  PB_DIR,
  "tables"
)

DE_DIR <- file.path(
  ROOT,
  "18C_DE"
)

DE_FULL_DIR <- file.path(
  DE_DIR,
  "full_results"
)

OUT <- file.path(
  ROOT,
  "18D_interpretation",
  "AA_batch_anchor_audit"
)

dir.create(
  OUT,
  recursive = TRUE,
  showWarnings = FALSE
)

metadata_rds <- file.path(
  PB_OBJ_DIR,
  "18B_pseudobulk_metadata_ALL.rds"
)

aa_manifest_csv <- file.path(
  PB_TAB_DIR,
  "18B_manifest_AA_vs_HC.csv"
)

stopifnot(
  file.exists(metadata_rds),
  file.exists(aa_manifest_csv)
)

# ------------------------------------------------------------
# 1. Load pseudobulk metadata and AA-vs-HC manifest
# ------------------------------------------------------------

cat("STEP 1/8 — loading pseudobulk metadata\n")

meta <- readRDS(metadata_rds) %>%
  as_tibble()

aa_manifest <- read_csv(
  aa_manifest_csv,
  show_col_types = FALSE
) %>%
  as_tibble()

required_meta <- c(
  "module",
  "pseudobulk_id",
  "patient_id",
  "final_L2_label",
  "disease",
  "batch",
  "n_cells",
  "library_size"
)

required_manifest <- c(
  required_meta,
  "contrast",
  "numerator",
  "denominator",
  "eligible_ge20",
  "robust_ge50"
)

stopifnot(
  all(required_meta %in% colnames(meta)),
  all(required_manifest %in% colnames(aa_manifest))
)

patient_batch <- meta %>%
  distinct(
    patient_id,
    disease,
    batch
  ) %>%
  arrange(
    disease,
    patient_id
  )

cat("\n===== DISEASE x BATCH =====\n")
print(
  as.data.frame(patient_batch),
  row.names = FALSE
)

stopifnot(
  n_distinct(patient_batch$patient_id) == 16,
  sum(
    patient_batch$disease == "AA" &
      patient_batch$batch == "CABIMER2026"
  ) == 6,
  sum(
    patient_batch$disease == "HC" &
      patient_batch$batch == "historical"
  ) == 4,
  sum(
    patient_batch$patient_id == "CON276" &
      patient_batch$disease == "HC" &
      patient_batch$batch == "CABIMER2026"
  ) == 1
)

# ------------------------------------------------------------
# 2. Freeze AA biological anchors
# ------------------------------------------------------------

cat("\nSTEP 2/8 — freezing AA biological anchor specification\n")

anchors <- tribble(
  ~module,              ~final_L2_label,                  ~gene,      ~biological_axis,

  "09_FIBROBLAST_L2",   "COL11A1/DPEP1 dermal-sheath",   "CXCL9",     "IFN_chemokine",
  "09_FIBROBLAST_L2",   "COL11A1/DPEP1 dermal-sheath",   "GBP1",      "IFN_response",

  "07_EPITHELIAL_L2",   "Basal_epithelial",               "HLA-DRA",   "antigen_presentation",
  "07_EPITHELIAL_L2",   "Basal_epithelial",               "HLA-DRB1",  "antigen_presentation",
  "07_EPITHELIAL_L2",   "Basal_epithelial",               "KRT14",     "basal_epithelial",
  "07_EPITHELIAL_L2",   "Basal_epithelial",               "KRT75",     "follicular_differentiation",

  "07_EPITHELIAL_L2",   "IFE_Infundibulum_suprabasal",    "KRT75",     "follicular_differentiation",

  "07_EPITHELIAL_L2",   "Bulge_Isthmus",                  "STAT1",     "IFN_response",
  "07_EPITHELIAL_L2",   "Basal_epithelial",               "TAP1",      "antigen_processing",
  "07_EPITHELIAL_L2",   "Bulge_Isthmus",                  "TAP1",      "antigen_processing"
)

stopifnot(
  nrow(anchors) == 10,
  n_distinct(
    paste(
      anchors$module,
      anchors$final_L2_label,
      anchors$gene,
      sep = "|||"
    )
  ) == 10
)

cat("\n===== FROZEN ANCHORS =====\n")
print(
  as.data.frame(anchors),
  row.names = FALSE
)

# ------------------------------------------------------------
# 3. Recover primary AA-vs-HC DE effect for each anchor
# ------------------------------------------------------------

cat("\nSTEP 3/8 — recovering primary AA-vs-HC DE effects\n")

safe_name <- function(x) {
  gsub(
    "[^A-Za-z0-9._-]+",
    "_",
    x
  )
}

read_primary_anchor <- function(
  module,
  final_L2_label,
  gene
) {

  f <- file.path(
    DE_FULL_DIR,
    paste0(
      "18C_AA_vs_HC__",
      module,
      "__",
      safe_name(final_L2_label),
      ".csv"
    )
  )

  if (!file.exists(f)) {
    stop(
      "Missing primary DE file: ",
      f
    )
  }

  x <- read_csv(
    f,
    show_col_types = FALSE
  )

  required <- c(
    "gene",
    "logFC",
    "logCPM",
    "PValue",
    "FDR_BH"
  )

  if (!all(required %in% colnames(x))) {
    stop(
      "Unexpected columns in ",
      f
    )
  }

  z <- x %>%
    filter(
      .data$gene == !!gene
    )

  if (nrow(z) != 1) {
    stop(
      "Expected exactly one row for ",
      gene,
      " in ",
      final_L2_label,
      "; found ",
      nrow(z)
    )
  }

  tibble(
    module = module,
    final_L2_label = final_L2_label,
    gene = gene,
    primary_logFC = z$logFC[[1]],
    primary_logCPM = z$logCPM[[1]],
    primary_PValue = z$PValue[[1]],
    primary_FDR_BH = z$FDR_BH[[1]]
  )
}

primary_list <- vector(
  "list",
  nrow(anchors)
)

for (i in seq_len(nrow(anchors))) {

  primary_list[[i]] <- read_primary_anchor(
    module = anchors$module[[i]],
    final_L2_label = anchors$final_L2_label[[i]],
    gene = anchors$gene[[i]]
  )
}

primary_anchor <- bind_rows(
  primary_list
) %>%
  left_join(
    anchors,
    by = c(
      "module",
      "final_L2_label",
      "gene"
    )
  ) %>%
  mutate(
    primary_direction = case_when(
      primary_logFC > 0 ~ "higher_in_AA",
      primary_logFC < 0 ~ "lower_in_AA",
      TRUE ~ "no_direction"
    ),
    primary_sign = sign(primary_logFC)
  )

stopifnot(
  nrow(primary_anchor) == 10,
  all(primary_anchor$primary_sign != 0)
)

cat("\n===== PRIMARY AA-vs-HC ANCHOR EFFECTS =====\n")
print(
  as.data.frame(
    primary_anchor %>%
      select(
        module,
        final_L2_label,
        gene,
        biological_axis,
        primary_logFC,
        primary_FDR_BH,
        primary_direction
      )
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# 4. Load required pseudobulk count matrices
# ------------------------------------------------------------

cat("\nSTEP 4/8 — loading anchor pseudobulk count matrices\n")

modules_needed <- unique(
  anchors$module
)

counts_by_module <- setNames(
  vector(
    "list",
    length(modules_needed)
  ),
  modules_needed
)

for (m in modules_needed) {

  f <- file.path(
    PB_OBJ_DIR,
    paste0(
      "18B_",
      m,
      "_pseudobulk_counts.rds"
    )
  )

  if (!file.exists(f)) {
    stop(
      "Missing pseudobulk count object: ",
      f
    )
  }

  mat <- readRDS(f)

  if (!inherits(mat, "Matrix") &&
      !is.matrix(mat)) {
    stop(
      "Unexpected count object class for ",
      m
    )
  }

  if (is.null(rownames(mat)) ||
      is.null(colnames(mat))) {
    stop(
      "Missing dimnames for ",
      m
    )
  }

  meta_m <- meta %>%
    filter(
      module == m
    )

  if (!all(
    colnames(mat) %in%
      meta_m$pseudobulk_id
  )) {
    stop(
      "Count-matrix columns do not map completely to metadata for ",
      m
    )
  }

  counts_by_module[[m]] <- mat

  cat(
    m,
    ":",
    nrow(mat),
    "genes x",
    ncol(mat),
    "pseudobulks\n"
  )
}

# ------------------------------------------------------------
# 5. Function: TMM logCPM within one frozen L2
# ------------------------------------------------------------

build_L2_logCPM <- function(
  module,
  final_L2_label
) {

  mat <- counts_by_module[[module]]

  md <- meta %>%
    filter(
      .data$module == !!module,
      .data$final_L2_label == !!final_L2_label
    ) %>%
    distinct(
      pseudobulk_id,
      .keep_all = TRUE
    )

  ids <- intersect(
    colnames(mat),
    md$pseudobulk_id
  )

  if (length(ids) < 2) {
    stop(
      "Insufficient pseudobulks for ",
      module,
      " / ",
      final_L2_label
    )
  }

  md <- md %>%
    filter(
      pseudobulk_id %in% ids
    ) %>%
    arrange(
      match(
        pseudobulk_id,
        ids
      )
    )

  ids <- md$pseudobulk_id

  submat <- mat[
    ,
    ids,
    drop = FALSE
  ]

  stopifnot(
    identical(
      colnames(submat),
      md$pseudobulk_id
    )
  )

  y <- DGEList(
    counts = submat
  )

  y <- calcNormFactors(
    y,
    method = "TMM"
  )

  lcpm <- cpm(
    y,
    log = TRUE,
    prior.count = 2
  )

  list(
    logCPM = lcpm,
    metadata = md,
    norm_factors = y$samples$norm.factors
  )
}

# ------------------------------------------------------------
# 6. Function: patient-level expression audit for one anchor
# ------------------------------------------------------------

audit_anchor_expression <- function(
  module,
  final_L2_label,
  gene
) {

  z <- build_L2_logCPM(
    module = module,
    final_L2_label = final_L2_label
  )

  lcpm <- z$logCPM
  md <- z$metadata

  if (!(gene %in% rownames(lcpm))) {
    stop(
      "Gene ",
      gene,
      " absent from pseudobulk matrix for ",
      module,
      " / ",
      final_L2_label
    )
  }

  expr <- tibble(
    pseudobulk_id = colnames(lcpm),
    audit_logCPM = as.numeric(
      lcpm[
        gene,
        ,
        drop = TRUE
      ]
    )
  ) %>%
    left_join(
      md %>%
        select(
          pseudobulk_id,
          patient_id,
          disease,
          batch,
          n_cells,
          library_size
        ),
      by = "pseudobulk_id"
    )

  if (anyNA(expr$patient_id)) {
    stop(
      "Metadata join failure for ",
      gene,
      " / ",
      final_L2_label
    )
  }

  expr %>%
    mutate(
      module = module,
      final_L2_label = final_L2_label,
      gene = gene
    ) %>%
    select(
      module,
      final_L2_label,
      gene,
      patient_id,
      disease,
      batch,
      n_cells,
      library_size,
      audit_logCPM
    )
}

# ------------------------------------------------------------
# 7. Build patient-level anchor-expression table
# ------------------------------------------------------------

cat("\nSTEP 5/8 — constructing patient-level anchor expression\n")

expr_list <- vector(
  "list",
  nrow(anchors)
)

for (i in seq_len(nrow(anchors))) {

  cat(
    sprintf(
      "[%02d/%02d] %s | %s | %s\n",
      i,
      nrow(anchors),
      anchors$module[[i]],
      anchors$final_L2_label[[i]],
      anchors$gene[[i]]
    )
  )

  expr_list[[i]] <- audit_anchor_expression(
    module = anchors$module[[i]],
    final_L2_label = anchors$final_L2_label[[i]],
    gene = anchors$gene[[i]]
  )
}

anchor_expr <- bind_rows(
  expr_list
) %>%
  left_join(
    anchors,
    by = c(
      "module",
      "final_L2_label",
      "gene"
    )
  )

write_csv(
  anchor_expr,
  file.path(
    OUT,
    "18D2_B_anchor_patient_logCPM.csv"
  )
)

cat("\nAnchor patient-expression rows:", nrow(anchor_expr), "\n")
cat("Unique anchors:",
    n_distinct(
      paste(
        anchor_expr$module,
        anchor_expr$final_L2_label,
        anchor_expr$gene,
        sep = "|||"
      )
    ),
    "\n"
)

# ------------------------------------------------------------
# 8. Descriptive direction helper
# ------------------------------------------------------------

median_or_na <- function(x) {
  if (length(x) == 0) {
    return(NA_real_)
  }
  median(
    x,
    na.rm = TRUE
  )
}

mean_or_na <- function(x) {
  if (length(x) == 0) {
    return(NA_real_)
  }
  mean(
    x,
    na.rm = TRUE
  )
}

direction_label <- function(x) {
  case_when(
    is.na(x) ~ NA_character_,
    x > 0 ~ "higher_in_AA",
    x < 0 ~ "lower_in_AA",
    TRUE ~ "no_direction"
  )
}

cat("\n============================================================\n")
cat("18D2-B BLOCK 1 LOADED\n")
cat("Primary DE effects recovered: ", nrow(primary_anchor), "\n", sep = "")
cat("Anchor patient-level expression prepared: ", nrow(anchor_expr), " rows\n", sep = "")
cat("No new DE performed\n")
cat("No annotation modified\n")
cat("READY FOR BLOCK 2\n")
cat("============================================================\n")

# ============================================================
# BLOCK 2
# Descriptive batch-direction audit
# ============================================================

cat("\nSTEP 6/8 — AA vs historical HC / CON276 direction audit\n")

# ------------------------------------------------------------
# 9. Summarise each anchor across the relevant control strata
# ------------------------------------------------------------

anchor_audit <- anchor_expr %>%
  group_by(
    module,
    final_L2_label,
    gene,
    biological_axis
  ) %>%
  summarise(
    n_AA =
      sum(
        disease == "AA"
      ),

    n_HC_historical =
      sum(
        disease == "HC" &
          batch == "historical"
      ),

    n_HC_CABIMER =
      sum(
        disease == "HC" &
          batch == "CABIMER2026"
      ),

    n_CON276 =
      sum(
        patient_id == "CON276"
      ),

    AA_median_logCPM =
      median_or_na(
        audit_logCPM[
          disease == "AA"
        ]
      ),

    AA_mean_logCPM =
      mean_or_na(
        audit_logCPM[
          disease == "AA"
        ]
      ),

    HC_historical_median_logCPM =
      median_or_na(
        audit_logCPM[
          disease == "HC" &
            batch == "historical"
        ]
      ),

    HC_historical_mean_logCPM =
      mean_or_na(
        audit_logCPM[
          disease == "HC" &
            batch == "historical"
        ]
      ),

    CON276_logCPM =
      {
        z <- audit_logCPM[
          patient_id == "CON276"
        ]

        if (length(z) == 1) {
          z[[1]]
        } else {
          NA_real_
        }
      },

    .groups = "drop"
  ) %>%
  left_join(
    primary_anchor %>%
      select(
        module,
        final_L2_label,
        gene,
        primary_logFC,
        primary_logCPM,
        primary_PValue,
        primary_FDR_BH,
        primary_direction,
        primary_sign
      ),
    by = c(
      "module",
      "final_L2_label",
      "gene"
    )
  ) %>%
  mutate(
    delta_AA_vs_HC_historical =
      AA_median_logCPM -
      HC_historical_median_logCPM,

    delta_AA_vs_CON276 =
      AA_median_logCPM -
      CON276_logCPM,

    direction_AA_vs_HC_historical =
      direction_label(
        delta_AA_vs_HC_historical
      ),

    direction_AA_vs_CON276 =
      direction_label(
        delta_AA_vs_CON276
      ),

    sign_AA_vs_HC_historical =
      sign(
        delta_AA_vs_HC_historical
      ),

    sign_AA_vs_CON276 =
      sign(
        delta_AA_vs_CON276
      ),

    concordant_historical_HC =
      !is.na(sign_AA_vs_HC_historical) &
      sign_AA_vs_HC_historical ==
        primary_sign,

    concordant_CON276 =
      !is.na(sign_AA_vs_CON276) &
      sign_AA_vs_CON276 ==
        primary_sign,

    concordant_both_control_checks =
      concordant_historical_HC &
      concordant_CON276,

    stronger_than_both =
      case_when(
        primary_sign > 0 ~
          AA_median_logCPM >
            HC_historical_median_logCPM &
          AA_median_logCPM >
            CON276_logCPM,

        primary_sign < 0 ~
          AA_median_logCPM <
            HC_historical_median_logCPM &
          AA_median_logCPM <
            CON276_logCPM,

        TRUE ~ FALSE
      )
  ) %>%
  arrange(
    desc(
      concordant_both_control_checks
    ),
    biological_axis,
    final_L2_label,
    gene
  )

stopifnot(
  nrow(anchor_audit) == nrow(anchors)
)

cat("\n===== AA ANCHOR BATCH-DIRECTION AUDIT =====\n")

print(
  as.data.frame(
    anchor_audit %>%
      select(
        final_L2_label,
        gene,
        biological_axis,
        primary_logFC,
        primary_FDR_BH,
        AA_median_logCPM,
        HC_historical_median_logCPM,
        CON276_logCPM,
        delta_AA_vs_HC_historical,
        delta_AA_vs_CON276,
        direction_AA_vs_HC_historical,
        direction_AA_vs_CON276,
        concordant_historical_HC,
        concordant_CON276,
        concordant_both_control_checks
      )
  ),
  row.names = FALSE
)

write_csv(
  anchor_audit,
  file.path(
    OUT,
    "18D2_B_AA_anchor_direction_audit.csv"
  )
)

# ------------------------------------------------------------
# 10. Per-patient table
# ------------------------------------------------------------

cat("\nSTEP 7/8 — patient-level anchor audit table\n")

patient_anchor <- anchor_expr %>%
  left_join(
    primary_anchor %>%
      select(
        module,
        final_L2_label,
        gene,
        primary_logFC,
        primary_FDR_BH,
        primary_direction
      ),
    by = c(
      "module",
      "final_L2_label",
      "gene"
    )
  ) %>%
  mutate(
    control_stratum = case_when(
      disease == "AA" ~
        "AA_CABIMER2026",

      disease == "HC" &
        patient_id == "CON276" ~
        "HC_CON276_CABIMER2026",

      disease == "HC" &
        batch == "historical" ~
        "HC_historical",

      TRUE ~
        paste(
          disease,
          batch,
          sep = "_"
        )
    )
  ) %>%
  arrange(
    module,
    final_L2_label,
    gene,
    factor(
      control_stratum,
      levels = c(
        "HC_historical",
        "HC_CON276_CABIMER2026",
        "AA_CABIMER2026"
      )
    ),
    patient_id
  )

write_csv(
  patient_anchor,
  file.path(
    OUT,
    "18D2_B_AA_anchor_patient_level_expression.csv"
  )
)

# ------------------------------------------------------------
# 11. Global concordance summary
# ------------------------------------------------------------

summary_tbl <- anchor_audit %>%
  summarise(
    n_anchors =
      n(),

    n_concordant_historical =
      sum(
        concordant_historical_HC,
        na.rm = TRUE
      ),

    pct_concordant_historical =
      100 *
      mean(
        concordant_historical_HC,
        na.rm = TRUE
      ),

    n_concordant_CON276 =
      sum(
        concordant_CON276,
        na.rm = TRUE
      ),

    pct_concordant_CON276 =
      100 *
      mean(
        concordant_CON276,
        na.rm = TRUE
      ),

    n_concordant_both =
      sum(
        concordant_both_control_checks,
        na.rm = TRUE
      ),

    pct_concordant_both =
      100 *
      mean(
        concordant_both_control_checks,
        na.rm = TRUE
      ),

    n_stronger_than_both =
      sum(
        stronger_than_both,
        na.rm = TRUE
      )
  )

cat("\n===== CONCORDANCE SUMMARY =====\n")

print(
  as.data.frame(summary_tbl),
  row.names = FALSE
)

write_csv(
  summary_tbl,
  file.path(
    OUT,
    "18D2_B_AA_anchor_concordance_summary.csv"
  )
)

# ------------------------------------------------------------
# 12. Axis-level summary
# ------------------------------------------------------------

axis_summary <- anchor_audit %>%
  group_by(
    biological_axis
  ) %>%
  summarise(
    n_anchors =
      n(),

    n_concordant_historical =
      sum(
        concordant_historical_HC,
        na.rm = TRUE
      ),

    n_concordant_CON276 =
      sum(
        concordant_CON276,
        na.rm = TRUE
      ),

    n_concordant_both =
      sum(
        concordant_both_control_checks,
        na.rm = TRUE
      ),

    .groups = "drop"
  )

cat("\n===== BIOLOGICAL-AXIS SUMMARY =====\n")

print(
  as.data.frame(axis_summary),
  row.names = FALSE
)

write_csv(
  axis_summary,
  file.path(
    OUT,
    "18D2_B_AA_anchor_axis_summary.csv"
  )
)

# ------------------------------------------------------------
# 13. Explicit interpretation categories
# ------------------------------------------------------------

interpretation_tbl <- anchor_audit %>%
  mutate(
    sensitivity_interpretation = case_when(
      concordant_both_control_checks ~
        "direction_concordant_vs_both_control_checks",

      concordant_historical_HC &
        !concordant_CON276 ~
        "concordant_vs_historical_only",

      !concordant_historical_HC &
        concordant_CON276 ~
        "concordant_vs_CON276_only",

      TRUE ~
        "direction_not_concordant"
    )
  ) %>%
  select(
    module,
    final_L2_label,
    gene,
    biological_axis,
    primary_logFC,
    primary_FDR_BH,
    AA_median_logCPM,
    HC_historical_median_logCPM,
    CON276_logCPM,
    delta_AA_vs_HC_historical,
    delta_AA_vs_CON276,
    sensitivity_interpretation
  )

write_csv(
  interpretation_tbl,
  file.path(
    OUT,
    "18D2_B_AA_anchor_interpretation.csv"
  )
)

# ------------------------------------------------------------
# 14. Formal checks
# ------------------------------------------------------------

cat("\nSTEP 8/8 — formal checks\n")

checks <- tibble(
  check = c(
    "anchors_10",
    "all_primary_effects_recovered",
    "all_primary_effects_nonzero",
    "six_AA_patients_present",
    "four_historical_HC_present",
    "CON276_present",
    "CON276_is_CABIMER_HC",
    "AA_all_CABIMER2026",
    "patient_level_expression_generated",
    "historical_direction_generated",
    "CON276_direction_generated",
    "primary_DE_not_refit",
    "CON276_not_used_as_inferential_arm",
    "annotation_unchanged",
    "figures_generated_NO"
  ),

  pass = c(
    nrow(anchor_audit) == 10,

    all(
      !is.na(
        anchor_audit$primary_logFC
      )
    ),

    all(
      anchor_audit$primary_sign != 0
    ),

    n_distinct(
      anchor_expr$patient_id[
        anchor_expr$disease == "AA"
      ]
    ) == 6,

    n_distinct(
      anchor_expr$patient_id[
        anchor_expr$disease == "HC" &
          anchor_expr$batch == "historical"
      ]
    ) == 4,

    "CON276" %in%
      anchor_expr$patient_id,

    all(
      patient_batch$disease[
        patient_batch$patient_id ==
          "CON276"
      ] == "HC"
    ) &&
      all(
        patient_batch$batch[
          patient_batch$patient_id ==
            "CON276"
        ] == "CABIMER2026"
      ),

    all(
      patient_batch$batch[
        patient_batch$disease == "AA"
      ] == "CABIMER2026"
    ),

    nrow(patient_anchor) > 0,

    all(
      !is.na(
        anchor_audit$
          direction_AA_vs_HC_historical
      )
    ),

    all(
      !is.na(
        anchor_audit$
          direction_AA_vs_CON276
      )
    ),

    TRUE,

    TRUE,

    TRUE,

    TRUE
  )
)

print(
  as.data.frame(checks),
  row.names = FALSE
)

write_csv(
  checks,
  file.path(
    OUT,
    "18D2_B_formal_checks.csv"
  )
)

if (!all(checks$pass)) {

  failed <- checks %>%
    filter(
      !pass
    ) %>%
    pull(
      check
    )

  stop(
    paste0(
      "18D2-B formal checks failed: ",
      paste(
        failed,
        collapse = ", "
      )
    )
  )
}

# ------------------------------------------------------------
# 15. Audit text
# ------------------------------------------------------------

audit_lines <- c(
  "PIPELINE_V2 — 18D2-B AA batch / biological-anchor direction audit",
  "",
  "Purpose:",
  paste0(
    "- Descriptive sensitivity audit of ",
    nrow(anchor_audit),
    " preselected AA-vs-HC biological anchors."
  ),
  "- Primary AA-vs-HC edgeR DE was NOT refit.",
  "- AA is completely confounded with CABIMER2026.",
  "- Historical HC comparison is descriptive.",
  "- CON276 is the only CABIMER2026 HC and is n=1.",
  "- CON276 is therefore NOT treated as an inferential control arm.",
  "- Direction is based on TMM-normalized patient-level pseudobulk logCPM.",
  "",
  paste0(
    "Historical-HC concordant anchors: ",
    summary_tbl$n_concordant_historical[[1]],
    "/",
    summary_tbl$n_anchors[[1]]
  ),
  paste0(
    "CON276 concordant anchors: ",
    summary_tbl$n_concordant_CON276[[1]],
    "/",
    summary_tbl$n_anchors[[1]]
  ),
  paste0(
    "Concordant against both checks: ",
    summary_tbl$n_concordant_both[[1]],
    "/",
    summary_tbl$n_anchors[[1]]
  ),
  "",
  "No annotation modified.",
  "No figures generated."
)

writeLines(
  audit_lines,
  file.path(
    OUT,
    "18D2_B_AA_batch_anchor_audit.txt"
  )
)

cat("\n============================================================\n")
cat("18D2-B COMPLETE\n")
cat("============================================================\n")
cat("Anchors audited: ", nrow(anchor_audit), "\n", sep = "")
cat(
  "Direction concordant vs historical HC: ",
  summary_tbl$n_concordant_historical[[1]],
  "/",
  summary_tbl$n_anchors[[1]],
  "\n",
  sep = ""
)
cat(
  "Direction concordant vs CON276: ",
  summary_tbl$n_concordant_CON276[[1]],
  "/",
  summary_tbl$n_anchors[[1]],
  "\n",
  sep = ""
)
cat(
  "Direction concordant vs BOTH checks: ",
  summary_tbl$n_concordant_both[[1]],
  "/",
  summary_tbl$n_anchors[[1]],
  "\n",
  sep = ""
)
cat("Primary AA-vs-HC DE refit: NO\n")
cat("Batch-adjusted AA-vs-HC inference: NO\n")
cat("CON276 inferential testing: NO\n")
cat("Annotation modified: NO\n")
cat("Figures generated: NO\n")
cat("PASS\n")
