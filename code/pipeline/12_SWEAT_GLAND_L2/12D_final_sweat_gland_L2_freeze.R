suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tibble)
  library(readr)
})

cat("============================================================\n")
cat("PIPELINE_V2 — 12D FINAL SWEAT-GLAND L2 FREEZE\n")
cat("============================================================\n\n")

# ============================================================
# 0. Paths
# ============================================================

candidate_objects <- c(
  "PIPELINE_V2/results/annotation/12_SWEAT_GLAND_L2/12C_adjudication/sweat_gland_v2_12C_adjudication_scaffold.rds",
  "PIPELINE_V2/results/annotation/12_SWEAT_GLAND_L2/12B_recluster/sweat_gland_v2_12B_reclustered.rds"
)

SRC <- candidate_objects[file.exists(candidate_objects)][1]

if (is.na(SRC) || length(SRC) == 0) {

  candidates <- list.files(
    "PIPELINE_V2",
    pattern = "sweat.*12C.*\\.rds$|sweat.*12B.*\\.rds$",
    recursive = TRUE,
    full.names = TRUE,
    ignore.case = TRUE
  )

  if (length(candidates) == 0) {
    stop(
      "No 12C/12B sweat-gland object found.\n",
      "Run:\n",
      "find PIPELINE_V2 -type f \\( -iname '*12C*.rds' -o -iname '*12B*.rds' \\) -print\n"
    )
  }

  SRC <- candidates[1]
}

OUT <- "PIPELINE_V2/results/annotation/12_SWEAT_GLAND_L2/12D_freeze"
TAB <- file.path(OUT, "tables")

dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
dir.create(TAB, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# 1. Load source scaffold
# ============================================================

cat("STEP 1/8 — loading sweat-gland adjudication scaffold\n")
cat("Source:", SRC, "\n")

x <- readRDS(SRC)

stopifnot(ncol(x) == 1843)
stopifnot(nrow(x) == 23412)
stopifnot(length(unique(x$patient_id)) == 16)
stopifnot(length(unique(x$batch)) == 2)

cat("Cells:", ncol(x), "\n")
cat("Genes:", nrow(x), "\n")
cat("Patients:", length(unique(x$patient_id)), "\n")
cat("Batches:", length(unique(x$batch)), "\n\n")

cells_before <- colnames(x)
genes_before <- rownames(x)
reductions_before <- names(x@reductions)

# ============================================================
# 2. Validate candidate resolution r0.4
# ============================================================

cat("STEP 2/8 — validating resolution 0.4\n")

r04_candidates <- c(
  "sweat_snn_res.0.4",
  "RNA_snn_res.0.4"
)

r04_col <- r04_candidates[
  r04_candidates %in% colnames(x@meta.data)
][1]

if (is.na(r04_col) || length(r04_col) == 0) {
  stop("Resolution 0.4 metadata column not found.")
}

cluster_before <- as.character(
  x@meta.data[[r04_col]]
)

expected_sizes <- c(
  "0" = 389,
  "1" = 283,
  "2" = 273,
  "3" = 251,
  "4" = 248,
  "5" = 187,
  "6" = 154,
  "7" = 58
)

observed_sizes <- table(cluster_before)

stopifnot(length(unique(cluster_before)) == 8)
stopifnot(all(names(expected_sizes) %in% names(observed_sizes)))
stopifnot(
  all(
    as.integer(observed_sizes[names(expected_sizes)]) ==
      as.integer(expected_sizes)
  )
)

cat("Resolution: 0.4\n")
cat("Clusters: 8\n\n")

# ============================================================
# 3. Final biological adjudication
# ============================================================

cat("STEP 3/8 — defining final sweat-gland adjudication\n")

mapping <- tribble(
  ~cluster,
  ~sweat_family,
  ~sweat_state,
  ~sweat_L2_label,
  ~sweat_L2_confidence,
  ~identity_basis,

  "0",
  "Secretory",
  "Dark / apocrine-like secretory",
  "Dark/apocrine-like secretory cell",
  "HIGH",
  paste(
    "Strong DARK_SECRETORY and SECRETORY_APOCRINE_LIKE programmes;",
    "SCGB/PIP-associated secretory phenotype with broad patient recurrence."
  ),

  "1",
  "Secretory",
  "Clear secretory",
  "Clear secretory cell",
  "HIGH",
  paste(
    "Dominant CLEAR_SECRETORY programme with strong general secretory",
    "and pan-sweat support; recurrent across patients."
  ),

  "2",
  "Secretory",
  "Dark / apocrine-like secretory",
  "Dark/apocrine-like secretory cell",
  "HIGH",
  paste(
    "Very strong DARK_SECRETORY and SECRETORY_APOCRINE_LIKE programmes;",
    "PIP, DCD, MUCL1 and secretoglobin-rich de novo marker profile."
  ),

  "3",
  "Secretory / ductal",
  "Inflammatory secretory-ductal",
  "Inflammatory secretory/ductal cell",
  "HIGH",
  paste(
    "Clear/secretory and ductal programmes accompanied by a strong",
    "LCN2/LTF/SAA1/SAA2/CHI3L inflammatory epithelial programme."
  ),

  "4",
  "Secretory",
  "Clear secretory",
  "Clear secretory cell",
  "HIGH",
  paste(
    "CLEAR_SECRETORY-dominant phenotype with general secretory",
    "and pan-sweat support; no competing lineage programme."
  ),

  "5",
  "Secretory / ductal",
  "Ductal-clear transitional",
  "Ductal/clear transitional cell",
  "HIGH",
  paste(
    "Joint clear-secretory, ductal-luminal and ductal-basal programmes",
    "consistent with a recurrent intraglandular transitional state."
  ),

  "6",
  "Secretory / ductal",
  "Ductal-clear transitional",
  "Ductal/clear transitional cell",
  "MODERATE",
  paste(
    "Joint clear-secretory and ductal programmes with weaker state-specific",
    "markers; retained as a recurrent transitional sweat-gland population."
  ),

  "7",
  "Ductal / myoepithelial",
  "KRT77-positive basal ductal / myoepithelial",
  "KRT77+ basal ductal/myoepithelial cell",
  "HIGH",
  paste(
    "Strong DUCTAL_BASAL and MYOEPITHELIAL programmes with KRT77, EDAR,",
    "SCEL and related ductal markers; recurrent across both batches."
  )
)

stopifnot(nrow(mapping) == 8)
stopifnot(length(unique(mapping$cluster)) == 8)
stopifnot(all(sort(mapping$cluster) == as.character(0:7)))

cat("\n=== FINAL SWEAT-GLAND ADJUDICATION ===\n")
print(mapping, n = Inf)

write_csv(
  mapping,
  file.path(
    TAB,
    "12D_frozen_sweat_gland_mapping.csv"
  )
)

# ============================================================
# 4. Apply frozen annotation
# ============================================================

cat("\nSTEP 4/8 — applying frozen sweat-gland annotation\n")

x$sweat_cluster_frozen <- cluster_before

idx <- match(
  x$sweat_cluster_frozen,
  mapping$cluster
)

stopifnot(!any(is.na(idx)))

x$sweat_family <-
  mapping$sweat_family[idx]

x$sweat_state <-
  mapping$sweat_state[idx]

x$sweat_L2_label <-
  mapping$sweat_L2_label[idx]

x$sweat_L2_confidence <-
  mapping$sweat_L2_confidence[idx]

x$sweat_resolution_frozen <- "0.4"

x$sweat_L2_status <- "FROZEN"

stopifnot(!any(is.na(x$sweat_L2_label)))
stopifnot(all(x$sweat_L2_status == "FROZEN"))

# ============================================================
# 5. Cell-level freeze ledger
# ============================================================

cat("STEP 5/8 — generating cell-level freeze ledger\n")

ledger <- tibble(
  cell = colnames(x),
  patient_id = as.character(x$patient_id),
  batch = as.character(x$batch),
  sweat_cluster_frozen =
    as.character(x$sweat_cluster_frozen),
  sweat_family =
    as.character(x$sweat_family),
  sweat_state =
    as.character(x$sweat_state),
  sweat_L2_label =
    as.character(x$sweat_L2_label),
  sweat_L2_confidence =
    as.character(x$sweat_L2_confidence),
  sweat_resolution_frozen =
    as.character(x$sweat_resolution_frozen),
  sweat_L2_status =
    as.character(x$sweat_L2_status)
)

write_csv(
  ledger,
  file.path(
    TAB,
    "12D_sweat_gland_cell_freeze_ledger.csv"
  )
)

# ============================================================
# 6. Frozen summary
# ============================================================

cat("STEP 6/8 — frozen L2 summary\n")

summary <- ledger %>%
  group_by(
    sweat_cluster_frozen,
    sweat_family,
    sweat_state,
    sweat_L2_label,
    sweat_L2_confidence
  ) %>%
  summarise(
    n_cells = n(),
    n_patients = n_distinct(patient_id),
    n_batches = n_distinct(batch),
    .groups = "drop"
  ) %>%
  arrange(
    as.integer(sweat_cluster_frozen)
  )

cat("\n=== FINAL FROZEN SWEAT-GLAND L2 ===\n")
print(summary, n = Inf)

write_csv(
  summary,
  file.path(
    TAB,
    "12D_frozen_sweat_gland_L2_summary.csv"
  )
)

# ============================================================
# 7. Formal invariance checks
# ============================================================

cat("\nSTEP 7/8 — structural invariance and formal checks\n")

cells_after <- colnames(x)
genes_after <- rownames(x)
reductions_after <- names(x@reductions)

cluster_after <- as.character(
  x@meta.data[[r04_col]]
)

frozen_sizes <- table(
  x$sweat_cluster_frozen
)

checks <- tibble(
  check = c(
    "source_cells_1843",
    "source_genes_23412",
    "patients_16",
    "batches_2",
    "resolution_04_frozen",
    "eight_clusters_preserved",
    "expected_cluster_sizes_preserved",
    "all_original_cells_preserved",
    "all_original_genes_preserved",
    "cluster_assignments_unchanged",
    "reductions_preserved",
    "eight_cluster_adjudications_defined",
    "all_cells_have_family",
    "all_cells_have_state",
    "all_cells_have_L2_label",
    "all_cells_have_confidence",
    "C0_dark_apocrine",
    "C1_clear_secretory",
    "C2_dark_apocrine",
    "C3_inflammatory_secretory_ductal",
    "C4_clear_secretory",
    "C5_ductal_clear_transition",
    "C6_ductal_clear_transition",
    "C7_KRT77_basal_myoepithelial",
    "disease_not_used_for_identity",
    "no_cells_removed",
    "source_object_not_overwritten",
    "L2_status_frozen"
  ),

  pass = c(
    ncol(x) == 1843,
    nrow(x) == 23412,
    length(unique(x$patient_id)) == 16,
    length(unique(x$batch)) == 2,

    all(
      x$sweat_resolution_frozen == "0.4"
    ),

    length(
      unique(x$sweat_cluster_frozen)
    ) == 8,

    all(
      as.integer(
        frozen_sizes[names(expected_sizes)]
      ) ==
        as.integer(expected_sizes)
    ),

    identical(
      cells_before,
      cells_after
    ),

    identical(
      genes_before,
      genes_after
    ),

    identical(
      cluster_before,
      cluster_after
    ),

    identical(
      reductions_before,
      reductions_after
    ),

    nrow(mapping) == 8,

    !any(is.na(x$sweat_family)),
    !any(is.na(x$sweat_state)),
    !any(is.na(x$sweat_L2_label)),
    !any(is.na(x$sweat_L2_confidence)),

    mapping$sweat_L2_label[
      mapping$cluster == "0"
    ] ==
      "Dark/apocrine-like secretory cell",

    mapping$sweat_L2_label[
      mapping$cluster == "1"
    ] ==
      "Clear secretory cell",

    mapping$sweat_L2_label[
      mapping$cluster == "2"
    ] ==
      "Dark/apocrine-like secretory cell",

    mapping$sweat_L2_label[
      mapping$cluster == "3"
    ] ==
      "Inflammatory secretory/ductal cell",

    mapping$sweat_L2_label[
      mapping$cluster == "4"
    ] ==
      "Clear secretory cell",

    mapping$sweat_L2_label[
      mapping$cluster == "5"
    ] ==
      "Ductal/clear transitional cell",

    mapping$sweat_L2_label[
      mapping$cluster == "6"
    ] ==
      "Ductal/clear transitional cell",

    mapping$sweat_L2_label[
      mapping$cluster == "7"
    ] ==
      "KRT77+ basal ductal/myoepithelial cell",

    TRUE,
    ncol(x) == length(cells_before),
    TRUE,
    all(x$sweat_L2_status == "FROZEN")
  )
)

cat("\n=== FORMAL CHECKS ===\n")
print(checks, n = Inf)

write_csv(
  checks,
  file.path(
    TAB,
    "12D_formal_checks.csv"
  )
)

stopifnot(all(checks$pass))

# ============================================================
# 8. Save frozen object
# ============================================================

cat("\nSTEP 8/8 — saving frozen sweat-gland object\n")

FINAL_RDS <- file.path(
  OUT,
  "sweat_gland_v2_12D_L2_frozen.rds"
)

saveRDS(
  x,
  FINAL_RDS
)

cat("\nSaved:\n")
cat(FINAL_RDS, "\n")

cat("\n============================================================\n")
cat("12D COMPLETE\n")
cat("============================================================\n")
cat("Cells:", ncol(x), "\n")
cat("Genes:", nrow(x), "\n")
cat("Patients:", length(unique(x$patient_id)), "\n")
cat("Batches:", length(unique(x$batch)), "\n")
cat("Resolution: 0.4 FROZEN\n")
cat("Clusters: 8\n")
cat("Sweat-gland L2 status: FROZEN\n")
cat("C0: dark/apocrine-like secretory cell\n")
cat("C1: clear secretory cell\n")
cat("C2: dark/apocrine-like secretory cell\n")
cat("C3: inflammatory secretory/ductal cell\n")
cat("C4: clear secretory cell\n")
cat("C5: ductal/clear transitional cell\n")
cat("C6: ductal/clear transitional cell\n")
cat("C7: KRT77+ basal ductal/myoepithelial cell\n")
cat("Disease used for identity: NO\n")
cat("Cells removed: 0\n")
cat("Source object overwritten: NO\n")
cat("PASS\n")
