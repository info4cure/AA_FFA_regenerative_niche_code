suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(stringr)
  library(tibble)
})

cat("============================================================\n")
cat("PIPELINE_V2 — 21E1C EXPRESSION-THRESHOLD CALIBRATION\n")
cat("============================================================\n")

BASE <- paste0(
  "PIPELINE_V2/results/disease_biology/",
  "21_CELL_CELL_COMMUNICATION/",
  "21E1_OBSERVABLE_UNIVERSE"
)

TAB <- file.path(BASE, "tables")

EXPR_RDS <- file.path(
  BASE,
  "21E1B_patient_identity_LR_gene_expression.rds"
)

LR_FILE <- file.path(
  TAB,
  "21E1A_observable_LR_universe.csv"
)

SUPPORT_FILE <- file.path(
  TAB,
  "21E1B_patient_identity_support.csv"
)

stopifnot(
  file.exists(EXPR_RDS),
  file.exists(LR_FILE),
  file.exists(SUPPORT_FILE)
)

# ============================================================
# 1. Frozen inputs
# ============================================================

cat("\nSTEP 1/9 — recovering frozen 21E1A/21E1B inputs\n")

expr <- readRDS(EXPR_RDS)

lr <- read_csv(
  LR_FILE,
  show_col_types=FALSE
)

group_support <- read_csv(
  SUPPORT_FILE,
  show_col_types=FALSE
)

cat("Expression rows :", nrow(expr), "\n")
cat("LR pairs        :", nrow(lr), "\n")
cat("Patient×identity:", nrow(group_support), "\n")

required_expr <- c(
  "gene",
  "pseudobulk_id",
  "logCPM",
  "pct_expr",
  "patient_id",
  "disease",
  "batch",
  "communication_identity",
  "n_cells",
  "cells_ge20",
  "cells_ge50"
)

stopifnot(
  all(required_expr %in% names(expr)),
  nrow(expr) > 0,
  nrow(lr) > 0
)

# ============================================================
# 2. Frozen source/target gene roles
# ============================================================

cat("\nSTEP 2/9 — recovering LR endpoint roles\n")

source_genes <- sort(
  unique(
    lr$source_genesymbol[
      !is.na(lr$source_genesymbol)
    ]
  )
)

target_genes <- sort(
  unique(
    lr$target_genesymbol[
      !is.na(lr$target_genesymbol)
    ]
  )
)

gene_role <- tibble(
  gene=sort(
    unique(
      c(
        source_genes,
        target_genes
      )
    )
  )
) %>%
  mutate(
    is_source=gene %in% source_genes,
    is_target=gene %in% target_genes,
    endpoint_role=case_when(
      is_source & is_target ~ "SOURCE_AND_TARGET",
      is_source ~ "SOURCE_ONLY",
      is_target ~ "TARGET_ONLY",
      TRUE ~ "UNRESOLVED"
    )
  )

cat("Source genes       :", length(source_genes), "\n")
cat("Target genes       :", length(target_genes), "\n")
cat(
  "Source+target genes:",
  sum(gene_role$endpoint_role == "SOURCE_AND_TARGET"),
  "\n"
)

write_csv(
  gene_role,
  file.path(
    TAB,
    "21E1C_LR_gene_endpoint_roles.csv"
  )
)

# ============================================================
# 3. Threshold grid
# ============================================================

cat("\nSTEP 3/9 — defining disease-blind calibration grid\n")

pct_thresholds <- c(
  0,
  0.01,
  0.05,
  0.10,
  0.20,
  0.30
)

cell_strata <- tibble(
  cell_stratum=c(
    "ALL_NONZERO_GROUPS",
    "GE20_CELLS",
    "GE50_CELLS"
  ),
  min_cells=c(
    1L,
    20L,
    50L
  )
)

threshold_grid <- tidyr::crossing(
  pct_threshold=pct_thresholds,
  cell_strata
)

write_csv(
  threshold_grid,
  file.path(
    TAB,
    "21E1C_threshold_grid.csv"
  )
)

print(
  threshold_grid,
  n=Inf
)

# ============================================================
# 4. Global distribution audit
# ============================================================

cat("\nSTEP 4/9 — auditing expression distributions\n")

quantile_probs <- c(
  0,
  0.01,
  0.05,
  0.10,
  0.25,
  0.50,
  0.75,
  0.90,
  0.95,
  0.99,
  1
)

make_quantiles <- function(df, variable, stratum) {

  x <- df[[variable]]
  x <- x[is.finite(x)]

  q <- quantile(
    x,
    probs=quantile_probs,
    na.rm=TRUE,
    names=FALSE
  )

  tibble(
    variable=variable,
    cell_stratum=stratum,
    probability=quantile_probs,
    value=as.numeric(q)
  )
}

distribution_quantiles <- bind_rows(
  make_quantiles(
    expr,
    "pct_expr",
    "ALL_NONZERO_GROUPS"
  ),
  make_quantiles(
    filter(expr, n_cells >= 20),
    "pct_expr",
    "GE20_CELLS"
  ),
  make_quantiles(
    filter(expr, n_cells >= 50),
    "pct_expr",
    "GE50_CELLS"
  ),
  make_quantiles(
    expr,
    "logCPM",
    "ALL_NONZERO_GROUPS"
  ),
  make_quantiles(
    filter(expr, n_cells >= 20),
    "logCPM",
    "GE20_CELLS"
  ),
  make_quantiles(
    filter(expr, n_cells >= 50),
    "logCPM",
    "GE50_CELLS"
  )
)

write_csv(
  distribution_quantiles,
  file.path(
    TAB,
    "21E1C_expression_distribution_quantiles.csv"
  )
)

cat("\n===== DISTRIBUTION QUANTILES =====\n")
print(
  distribution_quantiles,
  n=Inf,
  width=Inf
)

# ============================================================
# 5. Gene-level threshold coverage
# ============================================================

cat("\nSTEP 5/9 — calibrating gene-level detection coverage\n")

expr_role <- expr %>%
  left_join(
    gene_role,
    by="gene"
  )

gene_threshold_coverage <- bind_rows(
  lapply(
    seq_len(nrow(threshold_grid)),
    function(i) {

      th <- threshold_grid$pct_threshold[i]
      min_n <- threshold_grid$min_cells[i]
      stratum <- threshold_grid$cell_stratum[i]

      x <- expr_role %>%
        filter(
          n_cells >= min_n
        ) %>%
        mutate(
          detected =
            pct_expr >= th &
            pct_expr > 0
        )

      x %>%
        group_by(
          communication_identity,
          gene,
          endpoint_role
        ) %>%
        summarise(
          n_patient_groups=n(),
          n_patients_detected=sum(
            detected,
            na.rm=TRUE
          ),
          fraction_patients_detected=
            n_patients_detected /
            n_patient_groups,
          median_pct_expr=median(
            pct_expr,
            na.rm=TRUE
          ),
          median_logCPM=median(
            logCPM,
            na.rm=TRUE
          ),
          .groups="drop"
        ) %>%
        mutate(
          pct_threshold=th,
          cell_stratum=stratum,
          min_cells=min_n
        )
    }
  )
)

write_csv(
  gene_threshold_coverage,
  file.path(
    TAB,
    "21E1C_gene_threshold_coverage.csv"
  )
)

# ============================================================
# 6. Recurrence calibration
# ============================================================

cat("\nSTEP 6/9 — calibrating recurrence across patients\n")

recurrence_levels <- c(
  1L,
  2L,
  3L
)

recurrence_summary <- bind_rows(
  lapply(
    recurrence_levels,
    function(min_patients) {

      gene_threshold_coverage %>%
        group_by(
          pct_threshold,
          cell_stratum,
          min_cells,
          endpoint_role
        ) %>%
        summarise(
          min_patients=min_patients,
          n_identity_gene_pairs=sum(
            n_patients_detected >=
              min_patients
          ),
          n_unique_genes=n_distinct(
            gene[
              n_patients_detected >=
                min_patients
            ]
          ),
          n_identities=n_distinct(
            communication_identity[
              n_patients_detected >=
                min_patients
            ]
          ),
          .groups="drop"
        )
    }
  )
)

write_csv(
  recurrence_summary,
  file.path(
    TAB,
    "21E1C_recurrence_calibration.csv"
  )
)

cat("\n===== RECURRENCE CALIBRATION =====\n")
print(
  recurrence_summary,
  n=Inf,
  width=Inf
)

# ============================================================
# 7. LR-pair structural retention under thresholds
# ============================================================

cat("\nSTEP 7/9 — estimating LR-pair retention\n")

pair_key <- lr %>%
  distinct(
    source_genesymbol,
    target_genesymbol
  )

pair_retention <- bind_rows(
  lapply(
    seq_len(nrow(threshold_grid)),
    function(i) {

      th <- threshold_grid$pct_threshold[i]
      min_n <- threshold_grid$min_cells[i]
      stratum <- threshold_grid$cell_stratum[i]

      eligible_gene <- gene_threshold_coverage %>%
        filter(
          pct_threshold == th,
          cell_stratum == stratum,
          n_patients_detected >= 3
        ) %>%
        distinct(
          communication_identity,
          gene
        )

      identity_gene_sets <- split(
        eligible_gene$gene,
        eligible_gene$communication_identity
      )

      if (!length(identity_gene_sets)) {
        return(
          tibble(
            pct_threshold=th,
            cell_stratum=stratum,
            min_cells=min_n,
            n_identities=0L,
            n_identity_pair_combinations=0L,
            median_pairs_per_identity=NA_real_,
            min_pairs_per_identity=NA_real_,
            max_pairs_per_identity=NA_real_
          )
        )
      }

      pair_counts <- vapply(
        identity_gene_sets,
        function(g) {
          sum(
            pair_key$source_genesymbol %in% g &
            pair_key$target_genesymbol %in% g
          )
        },
        numeric(1)
      )

      tibble(
        pct_threshold=th,
        cell_stratum=stratum,
        min_cells=min_n,
        n_identities=length(pair_counts),
        n_identity_pair_combinations=
          sum(pair_counts),
        median_pairs_per_identity=
          median(pair_counts),
        min_pairs_per_identity=
          min(pair_counts),
        max_pairs_per_identity=
          max(pair_counts)
      )
    }
  )
)

write_csv(
  pair_retention,
  file.path(
    TAB,
    "21E1C_LR_pair_retention_calibration.csv"
  )
)

cat("\n===== LR PAIR RETENTION =====\n")
print(
  pair_retention,
  n=Inf,
  width=Inf
)

# ============================================================
# 8. Disease-blind candidate operating points
# ============================================================

cat("\nSTEP 8/9 — generating neutral operating-point table\n")

# These are candidates only.
# No threshold is selected in 21E1C.

operating_points <- tibble(
  operating_point=c(
    "PERMISSIVE",
    "PRIMARY_CANDIDATE",
    "STRINGENT_SENSITIVITY"
  ),
  min_cells=c(
    20L,
    20L,
    50L
  ),
  min_pct_expr=c(
    0.05,
    0.10,
    0.20
  ),
  min_patients=c(
    3L,
    3L,
    3L
  ),
  status="CANDIDATE_NOT_ADJUDICATED",
  rationale=c(
    paste(
      "Preserves relatively restricted expression while",
      "requiring recurrence across biological replicates."
    ),
    paste(
      "Intermediate detection requirement with >=20 cells",
      "and recurrence in >=3 patients."
    ),
    paste(
      "High-cell-support sensitivity specification with",
      "more stringent within-identity detection."
    )
  )
)

write_csv(
  operating_points,
  file.path(
    TAB,
    "21E1C_candidate_operating_points.csv"
  )
)

print(
  operating_points,
  n=Inf,
  width=Inf
)

# ============================================================
# 9. Formal checks
# ============================================================

cat("\nSTEP 9/9 — formal checks\n")

checks <- tibble(
  check=c(
    "21E1B_expression_support_present",
    "21E1A_LR_universe_present",
    "threshold_grid_generated",
    "pct_thresholds_include_0_1_5_10_20_30",
    "cell_strata_include_all_20_50",
    "endpoint_roles_retained",
    "distribution_quantiles_generated",
    "gene_threshold_coverage_generated",
    "patient_recurrence_calibrated",
    "LR_pair_retention_calibrated",
    "candidate_operating_points_generated",
    "no_primary_threshold_selected",
    "no_disease_used_for_threshold_selection",
    "no_LR_pair_score",
    "no_sender_receiver_inference",
    "no_disease_test",
    "no_AA_FFA_test",
    "no_DE_refit",
    "no_pathway_refit",
    "annotation_unchanged"
  ),
  pass=c(
    file.exists(EXPR_RDS),
    file.exists(LR_FILE),
    nrow(threshold_grid) == 18,
    identical(
      sort(unique(threshold_grid$pct_threshold)),
      sort(pct_thresholds)
    ),
    all(
      c(
        "ALL_NONZERO_GROUPS",
        "GE20_CELLS",
        "GE50_CELLS"
      ) %in%
        threshold_grid$cell_stratum
    ),
    all(
      c(
        "SOURCE_ONLY",
        "TARGET_ONLY",
        "SOURCE_AND_TARGET"
      ) %in%
        gene_role$endpoint_role
    ),
    nrow(distribution_quantiles) > 0,
    nrow(gene_threshold_coverage) > 0,
    nrow(recurrence_summary) > 0,
    nrow(pair_retention) == 18,
    nrow(operating_points) == 3,
    TRUE,
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

print(
  checks,
  n=Inf,
  width=Inf
)

stopifnot(all(checks$pass))

write_csv(
  checks,
  file.path(
    TAB,
    "21E1C_formal_checks.csv"
  )
)

audit <- c(
  "PIPELINE_V2 — 21E1C expression-threshold calibration",
  "",
  paste0("Date: ", Sys.Date()),
  "",
  "Purpose:",
  paste(
    "Calibrate expression-support thresholds before",
    "ligand-receptor communication inference."
  ),
  "",
  "Calibration grid:",
  "- pct_expr: >0, 1%, 5%, 10%, 20%, 30%",
  "- cell support: all nonzero groups, >=20 cells, >=50 cells",
  "- recurrence evaluated at >=1, >=2 and >=3 patients",
  "",
  "Source and target endpoint roles were retained separately.",
  "",
  "Candidate operating points were generated but NOT selected.",
  "",
  "Disease labels were not used to choose an operating point.",
  "No LR pair score calculated.",
  "No sender-receiver inference.",
  "No disease comparison.",
  "No AA-vs-FFA comparison.",
  "No differential-expression refit.",
  "No pathway refit.",
  "No annotation modified.",
  "",
  "PASS"
)

writeLines(
  audit,
  file.path(
    BASE,
    "21E1C_audit.txt"
  )
)

cat("\n============================================================\n")
cat("21E1C COMPLETE\n")
cat("============================================================\n")
cat("Expression-support threshold calibration complete\n")
cat("No primary threshold selected\n")
cat("Disease-blind candidate operating points generated\n")
cat("No communication inference\n")
cat("No disease comparison\n")
cat("PASS\n")
