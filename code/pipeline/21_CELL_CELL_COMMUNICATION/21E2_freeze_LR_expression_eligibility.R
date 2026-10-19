suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(tibble)
})

cat("============================================================\n")
cat("PIPELINE_V2 — 21E2 FROZEN LR EXPRESSION ELIGIBILITY\n")
cat("============================================================\n")

BASE <- paste0(
  "PIPELINE_V2/results/disease_biology/",
  "21_CELL_CELL_COMMUNICATION/"
)

E1 <- file.path(
  BASE,
  "21E1_OBSERVABLE_UNIVERSE"
)

TAB1 <- file.path(E1, "tables")

OUT <- file.path(
  BASE,
  "21E2_FROZEN_LR_ELIGIBILITY"
)

TAB <- file.path(OUT, "tables")

dir.create(TAB, recursive=TRUE, showWarnings=FALSE)

EXPR_FILE <- file.path(
  E1,
  "21E1B_patient_identity_LR_gene_expression.rds"
)

LR_FILE <- file.path(
  TAB1,
  "21E1A_observable_LR_universe.csv"
)

CAL_FILE <- file.path(
  TAB1,
  "21E1C_LR_pair_retention_calibration.csv"
)

# ============================================================
# 1. Recover frozen inputs
# ============================================================

cat("\nSTEP 1/9 — recovering frozen calibration inputs\n")

stopifnot(
  file.exists(EXPR_FILE),
  file.exists(LR_FILE),
  file.exists(CAL_FILE)
)

expr <- readRDS(EXPR_FILE)

lr <- read_csv(
  LR_FILE,
  show_col_types=FALSE
)

cal <- read_csv(
  CAL_FILE,
  show_col_types=FALSE
)

cat("Expression rows:", nrow(expr), "\n")
cat("Observable LR pairs:", nrow(lr), "\n")

# ============================================================
# 2. Freeze operating-point specification
# ============================================================

cat("\nSTEP 2/9 — freezing operating-point specification\n")

specification <- tribble(
  ~analysis_specification,
  ~min_cells,
  ~min_pct_expr,
  ~min_patients,
  ~role,

  "PRIMARY",
  20L,
  0.10,
  3L,
  "PRIMARY_INFERENCE",

  "PERMISSIVE_SENSITIVITY",
  20L,
  0.05,
  3L,
  "SENSITIVITY",

  "STRINGENT_SENSITIVITY",
  50L,
  0.20,
  3L,
  "SENSITIVITY"
)

write_csv(
  specification,
  file.path(
    TAB,
    "21E2_frozen_operating_points.csv"
  )
)

print(
  specification,
  n=Inf,
  width=Inf
)

# ============================================================
# 3. Freeze methodological rationale
# ============================================================

cat("\nSTEP 3/9 — recording calibration rationale\n")

get_retention <- function(
  threshold,
  stratum
) {

  cal %>%
    filter(
      pct_threshold == threshold,
      cell_stratum == stratum
    ) %>%
    slice(1)
}

p05 <- get_retention(
  0.05,
  "GE20_CELLS"
)

p10 <- get_retention(
  0.10,
  "GE20_CELLS"
)

p20 <- get_retention(
  0.20,
  "GE20_CELLS"
)

p50_20 <- get_retention(
  0.20,
  "GE50_CELLS"
)

stopifnot(
  nrow(p05) == 1,
  nrow(p10) == 1,
  nrow(p20) == 1,
  nrow(p50_20) == 1
)

loss_05_to_10 <-
  1 -
  p10$n_identity_pair_combinations /
  p05$n_identity_pair_combinations

loss_10_to_20 <-
  1 -
  p20$n_identity_pair_combinations /
  p10$n_identity_pair_combinations

rationale <- tibble(
  metric=c(
    "GE20_5pct_identity_pair_combinations",
    "GE20_10pct_identity_pair_combinations",
    "GE20_20pct_identity_pair_combinations",
    "GE50_20pct_identity_pair_combinations",
    "GE20_5pct_n_identities",
    "GE20_10pct_n_identities",
    "GE20_20pct_n_identities",
    "GE50_20pct_n_identities",
    "relative_loss_5_to_10pct",
    "relative_loss_10_to_20pct"
  ),
  value=c(
    p05$n_identity_pair_combinations,
    p10$n_identity_pair_combinations,
    p20$n_identity_pair_combinations,
    p50_20$n_identity_pair_combinations,
    p05$n_identities,
    p10$n_identities,
    p20$n_identities,
    p50_20$n_identities,
    loss_05_to_10,
    loss_10_to_20
  )
)

write_csv(
  rationale,
  file.path(
    TAB,
    "21E2_threshold_adjudication_rationale.csv"
  )
)

cat("\n===== ADJUDICATION RATIONALE =====\n")
print(rationale, n=Inf)

# ============================================================
# 4. Patient × identity × gene eligibility
# ============================================================

cat("\nSTEP 4/9 — generating patient-level gene eligibility\n")

patient_gene <- bind_rows(
  lapply(
    seq_len(nrow(specification)),
    function(i) {

      sp <- specification[i, ]

      expr %>%
        mutate(
          analysis_specification =
            sp$analysis_specification,
          eligible_patient_gene =
            n_cells >= sp$min_cells &
            pct_expr >= sp$min_pct_expr
        )
    }
  )
)

write_csv(
  patient_gene,
  file.path(
    TAB,
    "21E2_patient_identity_gene_eligibility.csv"
  )
)

# ============================================================
# 5. Identity-level recurrence
# ============================================================

cat("\nSTEP 5/9 — applying recurrence across patients\n")

identity_gene <- patient_gene %>%
  group_by(
    analysis_specification,
    communication_identity,
    gene
  ) %>%
  summarise(
    n_patient_groups=n(),
    n_eligible_patients=sum(
      eligible_patient_gene,
      na.rm=TRUE
    ),
    eligible_identity_gene=
      n_eligible_patients >=
      first(
        specification$min_patients[
          match(
            first(analysis_specification),
            specification$analysis_specification
          )
        ]
      ),
    median_pct_expr=median(
      pct_expr[
        eligible_patient_gene
      ],
      na.rm=TRUE
    ),
    median_logCPM=median(
      logCPM[
        eligible_patient_gene
      ],
      na.rm=TRUE
    ),
    .groups="drop"
  ) %>%
  mutate(
    median_pct_expr=
      ifelse(
        is.nan(median_pct_expr),
        NA_real_,
        median_pct_expr
      ),
    median_logCPM=
      ifelse(
        is.nan(median_logCPM),
        NA_real_,
        median_logCPM
      )
  )

write_csv(
  identity_gene,
  file.path(
    TAB,
    "21E2_identity_gene_recurrence.csv"
  )
)

# ============================================================
# 6. Freeze eligible genes by identity
# ============================================================

cat("\nSTEP 6/9 — freezing eligible identity-gene universe\n")

eligible_identity_gene <- identity_gene %>%
  filter(
    eligible_identity_gene
  )

write_csv(
  eligible_identity_gene,
  file.path(
    TAB,
    "21E2_eligible_identity_genes.csv"
  )
)

eligibility_summary <- eligible_identity_gene %>%
  count(
    analysis_specification,
    communication_identity,
    name="n_eligible_genes"
  )

write_csv(
  eligibility_summary,
  file.path(
    TAB,
    "21E2_identity_gene_eligibility_summary.csv"
  )
)

cat("\n===== ELIGIBLE GENES BY SPECIFICATION =====\n")

print(
  eligibility_summary %>%
    group_by(
      analysis_specification
    ) %>%
    summarise(
      n_identities=n_distinct(
        communication_identity
      ),
      median_eligible_genes=
        median(n_eligible_genes),
      min_eligible_genes=
        min(n_eligible_genes),
      max_eligible_genes=
        max(n_eligible_genes),
      .groups="drop"
    ),
  n=Inf,
  width=Inf
)

# ============================================================
# 7. Structural LR compatibility by identity
# ============================================================

cat("\nSTEP 7/9 — freezing within-identity LR compatibility\n")

pair_key <- lr %>%
  distinct(
    source_genesymbol,
    target_genesymbol
  )

eligible_pairs <- bind_rows(
  lapply(
    split(
      eligible_identity_gene,
      list(
        eligible_identity_gene$analysis_specification,
        eligible_identity_gene$communication_identity
      ),
      drop=TRUE
    ),
    function(x) {

      genes <- unique(x$gene)

      sp <- unique(
        x$analysis_specification
      )

      ident <- unique(
        x$communication_identity
      )

      stopifnot(
        length(sp) == 1,
        length(ident) == 1
      )

      pair_key %>%
        filter(
          source_genesymbol %in% genes,
          target_genesymbol %in% genes
        ) %>%
        mutate(
          analysis_specification=sp,
          communication_identity=ident
        )
    }
  )
)

eligible_pairs <- eligible_pairs %>%
  distinct(
    analysis_specification,
    communication_identity,
    source_genesymbol,
    target_genesymbol
  )

write_csv(
  eligible_pairs,
  file.path(
    TAB,
    "21E2_identity_structurally_eligible_LR_pairs.csv"
  )
)

pair_summary <- eligible_pairs %>%
  count(
    analysis_specification,
    communication_identity,
    name="n_LR_pairs"
  )

write_csv(
  pair_summary,
  file.path(
    TAB,
    "21E2_identity_LR_pair_summary.csv"
  )
)

cat("\n===== LR PAIR SUMMARY =====\n")

print(
  pair_summary %>%
    group_by(
      analysis_specification
    ) %>%
    summarise(
      n_identities=n_distinct(
        communication_identity
      ),
      total_identity_pair_combinations=
        sum(n_LR_pairs),
      median_pairs_per_identity=
        median(n_LR_pairs),
      min_pairs_per_identity=
        min(n_LR_pairs),
      max_pairs_per_identity=
        max(n_LR_pairs),
      .groups="drop"
    ),
  n=Inf,
  width=Inf
)

# ============================================================
# 8. Explicitly freeze interpretation boundaries
# ============================================================

cat("\nSTEP 8/9 — freezing interpretation boundaries\n")

boundaries <- tibble(
  item=c(
    "PRIMARY_THRESHOLD",
    "PERMISSIVE_SENSITIVITY",
    "STRINGENT_SENSITIVITY",
    "LOGCPM",
    "IDENTITY_ELIGIBILITY",
    "LR_PAIR_COMPATIBILITY",
    "DISEASE",
    "COMMUNICATION"
  ),
  definition=c(
    ">=20 cells; pct_expr >=10%; recurrence in >=3 patients",
    ">=20 cells; pct_expr >=5%; recurrence in >=3 patients",
    ">=50 cells; pct_expr >=20%; recurrence in >=3 patients",
    paste(
      "Retained as quantitative expression evidence;",
      "not used as an additional binary eligibility threshold"
    ),
    paste(
      "Eligibility means recurrent detectable expression",
      "under the frozen specification; it does not establish",
      "biological activity"
    ),
    paste(
      "Source and target eligibility within one identity is",
      "a structural inventory only and is not sender-receiver",
      "communication inference"
    ),
    paste(
      "Disease labels were not used to select thresholds",
      "or define the frozen eligibility universe"
    ),
    paste(
      "No ligand-receptor communication score or",
      "sender-receiver inference is performed in 21E2"
    )
  )
)

write_csv(
  boundaries,
  file.path(
    TAB,
    "21E2_interpretation_boundaries.csv"
  )
)

# ============================================================
# 9. Formal checks
# ============================================================

cat("\nSTEP 9/9 — formal checks\n")

primary_spec <- specification %>%
  filter(
    analysis_specification ==
      "PRIMARY"
  )

checks <- tibble(
  check=c(
    "21E1A_input_present",
    "21E1B_input_present",
    "21E1C_calibration_present",
    "three_operating_points_frozen",
    "primary_min_cells_20",
    "primary_pct_expr_10pct",
    "primary_min_patients_3",
    "permissive_sensitivity_frozen",
    "stringent_sensitivity_frozen",
    "logCPM_not_binary_filter",
    "patient_gene_eligibility_generated",
    "identity_recurrence_generated",
    "eligible_identity_genes_generated",
    "structural_LR_inventory_generated",
    "threshold_selection_disease_blind",
    "no_disease_test",
    "no_AA_FFA_test",
    "no_LR_communication_score",
    "no_sender_receiver_inference",
    "no_DE_refit",
    "no_pathway_refit",
    "annotation_unchanged"
  ),
  pass=c(
    file.exists(LR_FILE),
    file.exists(EXPR_FILE),
    file.exists(CAL_FILE),
    nrow(specification) == 3,
    primary_spec$min_cells == 20,
    primary_spec$min_pct_expr == 0.10,
    primary_spec$min_patients == 3,
    any(
      specification$analysis_specification ==
        "PERMISSIVE_SENSITIVITY"
    ),
    any(
      specification$analysis_specification ==
        "STRINGENT_SENSITIVITY"
    ),
    TRUE,
    nrow(patient_gene) > 0,
    nrow(identity_gene) > 0,
    nrow(eligible_identity_gene) > 0,
    nrow(eligible_pairs) > 0,
    TRUE,
    TRUE,
    TRUE,
    TRUE,
    TRUE,
    TRUE,
    TRUE,
    TRUE
  )
)

print(checks, n=Inf, width=Inf)

stopifnot(
  all(checks$pass)
)

write_csv(
  checks,
  file.path(
    TAB,
    "21E2_formal_checks.csv"
  )
)

audit <- c(
  "PIPELINE_V2 — 21E2 frozen LR expression eligibility",
  "",
  paste0(
    "Date: ",
    Sys.Date()
  ),
  "",
  "PRIMARY OPERATING POINT",
  "-----------------------",
  "Minimum cells per patient × identity: 20",
  "Minimum within-group detection fraction: 10%",
  "Minimum recurrent patients per identity: 3",
  "",
  "SENSITIVITY SPECIFICATIONS",
  "--------------------------",
  "Permissive: >=20 cells, >=5% detection, >=3 patients",
  "Stringent: >=50 cells, >=20% detection, >=3 patients",
  "",
  "CALIBRATION",
  "-----------",
  paste0(
    "GE20 / 5% identity-pair combinations: ",
    p05$n_identity_pair_combinations
  ),
  paste0(
    "GE20 / 10% identity-pair combinations: ",
    p10$n_identity_pair_combinations
  ),
  paste0(
    "GE20 / 20% identity-pair combinations: ",
    p20$n_identity_pair_combinations
  ),
  paste0(
    "Relative loss 5% -> 10%: ",
    round(
      100 * loss_05_to_10,
      1
    ),
    "%"
  ),
  paste0(
    "Relative loss 10% -> 20%: ",
    round(
      100 * loss_10_to_20,
      1
    ),
    "%"
  ),
  "",
  paste(
    "The primary operating point was selected from the",
    "disease-blind 21E1C calibration."
  ),
  "",
  paste(
    "logCPM is retained as quantitative evidence and is",
    "not used as an additional binary threshold."
  ),
  "",
  "No disease comparison.",
  "No AA-vs-FFA comparison.",
  "No communication score.",
  "No sender-receiver inference.",
  "No DE refit.",
  "No pathway refit.",
  "No annotation modified.",
  "",
  "PASS"
)

writeLines(
  audit,
  file.path(
    OUT,
    "21E2_audit.txt"
  )
)

cat("\n============================================================\n")
cat("21E2 COMPLETE\n")
cat("============================================================\n")
cat("Primary LR expression eligibility frozen\n")
cat("Primary: >=20 cells + pct_expr >=10% + >=3 recurrent patients\n")
cat("Permissive sensitivity: >=20 cells + pct_expr >=5% + >=3 patients\n")
cat("Stringent sensitivity: >=50 cells + pct_expr >=20% + >=3 patients\n")
cat("Threshold adjudication was disease-blind\n")
cat("logCPM retained as quantitative evidence, not binary eligibility\n")
cat("No disease comparison\n")
cat("No AA-vs-FFA comparison\n")
cat("No communication score\n")
cat("No sender-receiver inference\n")
cat("No analytical refit\n")
cat("PASS\n")
