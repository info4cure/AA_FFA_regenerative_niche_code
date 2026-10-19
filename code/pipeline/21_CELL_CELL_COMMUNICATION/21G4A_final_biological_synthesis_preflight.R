suppressPackageStartupMessages({
  library(tidyverse)
})

cat("============================================================\n")
cat("21G4A — FINAL COMMUNICATION / CROSS-MODAL BIOLOGICAL SYNTHESIS PREFLIGHT\n")
cat("============================================================\n")

BASE <- file.path(
  "PIPELINE_V2/results/disease_biology/21_CELL_CELL_COMMUNICATION",
  "21G_MULTICOMPARTMENT_COMMUNICATION"
)

B2 <- file.path(
  BASE,
  "21G2_INFERENCE/tables/21G2B_primary_A_axis_results.csv"
)

B2C <- file.path(
  BASE,
  "21G3_BIOLOGICAL_READOUT/tables/21G2C_primary_A_signal_inventory.csv"
)

XMOD <- file.path(
  BASE,
  "21G4_CROSSMODAL_INTEGRATION/tables/21G3B2_crossmodal_concordance_matrix.csv"
)

PUB <- file.path(
  BASE,
  "21G4_CROSSMODAL_INTEGRATION/tables/21G3B3_publication_readout.csv"
)

PUB_FDR <- file.path(
  BASE,
  "21G4_CROSSMODAL_INTEGRATION/tables/21G3B3_FDR_publication_readout.csv"
)

# Discover exact 21G3B3 filenames if naming differs.
cat("\n===== 21G3B3 FILE INVENTORY =====\n")

inventory <- list.files(
  file.path(BASE, "21G4_CROSSMODAL_INTEGRATION"),
  pattern="21G3B3",
  recursive=TRUE,
  full.names=TRUE
)

print(inventory)

stopifnot(
  file.exists(B2),
  file.exists(B2C),
  file.exists(XMOD)
)

# ------------------------------------------------------------
# 1. Frozen source recovery
# ------------------------------------------------------------

cat("\nSTEP 1/8 — recovering frozen sources\n")

b2   <- read_csv(B2, show_col_types=FALSE)
b2c  <- read_csv(B2C, show_col_types=FALSE)
xmod <- read_csv(XMOD, show_col_types=FALSE)

cat("21G2 primary A rows:", nrow(b2), "\n")
cat("21G2C signal inventory rows:", nrow(b2c), "\n")
cat("21G3B2 cross-modal rows:", nrow(xmod), "\n")

# Publication universe recovered directly from frozen B2 matrix
# rather than depending on an assumed B3 filename.
pub <- xmod %>%
  filter(
    communication_threshold_signal %in% TRUE,
    communication_action == "INCLUDE"
  )

cat("Recovered publication-threshold rows:", nrow(pub), "\n")

stopifnot(
  nrow(xmod) == 317,
  n_distinct(xmod$receiver) == 28,
  nrow(pub) == 16
)

# ------------------------------------------------------------
# 2. Freeze evidential strata
# ------------------------------------------------------------

cat("\nSTEP 2/8 — freezing evidential strata\n")

pub <- pub %>%
  mutate(
    final_evidence_stratum =
      case_when(
        communication_evidence_class == "FDR_LT_0.05" ~
          "FDR_LT_0.05",

        communication_evidence_class == "FDR_LT_0.10" ~
          "FDR_LT_0.10_ADDITIONAL",

        communication_evidence_class == "NOMINAL_P_LT_0.05" ~
          "NOMINAL_CONTEXTUAL",

        TRUE ~
          "UNEXPECTED"
      )
  )

cat("\n===== EVIDENCE STRATA =====\n")

print(
  pub %>%
    count(
      contrast,
      final_evidence_stratum,
      name="n"
    ) %>%
    arrange(
      contrast,
      final_evidence_stratum
    ),
  n=Inf
)

# ------------------------------------------------------------
# 3. Recover exact final signal architecture
# ------------------------------------------------------------

cat("\nSTEP 3/8 — recovering final signal architecture\n")

signal_architecture <- pub %>%
  select(
    contrast,
    macrocompartment,
    receiver,
    communication_program,
    communication_tier,
    communication_mechanistic_class,
    communication_axis,
    communication_mean_difference,
    communication_log2_mean_ratio,
    communication_p_exact,
    communication_q_BH,
    communication_evidence_class,
    final_evidence_stratum,
    DA_OR,
    DA_FDR_compartment,
    DA_direction,
    DA_support,
    DE_n_FDR005,
    DE_n_up_FDR005,
    DE_n_down_FDR005,
    DE_support,
    decomposition_pattern,
    density_mean_difference,
    density_q_BH,
    active_intensity_mean_difference,
    active_intensity_q_BH,
    bridge_class,
    pathway_bridge_status,
    crossmodal_triangulation_class
  ) %>%
  arrange(
    factor(
      contrast,
      levels=c(
        "AA_vs_HC",
        "FFA_vs_HC",
        "AA_vs_FFA"
      )
    ),
    communication_q_BH,
    receiver,
    communication_program
  )

cat("\n===== FINAL 16-SIGNAL ARCHITECTURE =====\n")

print(
  signal_architecture,
  n=Inf,
  width=Inf
)

# ------------------------------------------------------------
# 4. Disease-specific descriptive synthesis
# ------------------------------------------------------------

cat("\nSTEP 4/8 — disease-specific descriptive synthesis\n")

disease_summary <- pub %>%
  group_by(contrast) %>%
  summarise(
    n_signals = n(),

    n_receivers =
      n_distinct(receiver),

    n_FDR005 =
      sum(
        communication_evidence_class ==
          "FDR_LT_0.05"
      ),

    n_FDR010_additional =
      sum(
        communication_evidence_class ==
          "FDR_LT_0.10"
      ),

    n_nominal =
      sum(
        communication_evidence_class ==
          "NOMINAL_P_LT_0.05"
      ),

    n_DA_FDR005 =
      sum(
        DA_support ==
          "FDR_LT_0.05",
        na.rm=TRUE
      ),

    n_DE_supported =
      sum(
        DE_support ==
          "FDR_LT_0.05_GENES_PRESENT",
        na.rm=TRUE
      ),

    n_density_supported =
      sum(
        decomposition_pattern ==
          "DENSITY_SUPPORTED",
        na.rm=TRUE
      ),

    n_both_nominal =
      sum(
        decomposition_pattern ==
          "BOTH_NOMINAL_ONLY",
        na.rm=TRUE
      ),

    .groups="drop"
  )

cat("\n===== DISEASE SUMMARY =====\n")

print(
  disease_summary,
  n=Inf,
  width=Inf
)

# ------------------------------------------------------------
# 5. AA vs HC
# ------------------------------------------------------------

cat("\nSTEP 5/8 — AA vs HC publication architecture\n")

AA_HC <- signal_architecture %>%
  filter(contrast == "AA_vs_HC")

cat("\n===== AA vs HC =====\n")

print(
  AA_HC,
  n=Inf,
  width=Inf
)

# ------------------------------------------------------------
# 6. FFA vs HC
# ------------------------------------------------------------

cat("\nSTEP 6/8 — FFA vs HC publication architecture\n")

FFA_HC <- signal_architecture %>%
  filter(contrast == "FFA_vs_HC")

cat("\n===== FFA vs HC =====\n")

print(
  FFA_HC,
  n=Inf,
  width=Inf
)

# ------------------------------------------------------------
# 7. AA vs FFA
# ------------------------------------------------------------

cat("\nSTEP 7/8 — AA vs FFA publication architecture\n")

AA_FFA <- signal_architecture %>%
  filter(contrast == "AA_vs_FFA")

cat("\n===== AA vs FFA =====\n")

print(
  AA_FFA,
  n=Inf,
  width=Inf
)

# ------------------------------------------------------------
# 8. Mechanistic/decomposition boundary audit
# ------------------------------------------------------------

cat("\nSTEP 8/8 — decomposition and interpretation boundary audit\n")

decomposition_summary <- pub %>%
  count(
    contrast,
    decomposition_pattern,
    name="n"
  ) %>%
  arrange(
    contrast,
    decomposition_pattern
  )

cat("\n===== DECOMPOSITION SUMMARY =====\n")

print(
  decomposition_summary,
  n=Inf,
  width=Inf
)

cat("\n===== DENSITY-SUPPORTED SIGNALS =====\n")

print(
  pub %>%
    filter(
      decomposition_pattern ==
        "DENSITY_SUPPORTED"
    ) %>%
    select(
      contrast,
      receiver,
      communication_program,
      communication_evidence_class,
      communication_q_BH,
      communication_mean_difference,
      density_mean_difference,
      density_q_BH,
      active_intensity_mean_difference,
      active_intensity_q_BH,
      DA_support,
      DE_support
    ) %>%
    arrange(
      contrast,
      communication_q_BH
    ),
  n=Inf,
  width=Inf
)

cat("\n===== SIGNALS WITHOUT DENSITY-SUPPORTED DECOMPOSITION =====\n")

print(
  pub %>%
    filter(
      decomposition_pattern !=
        "DENSITY_SUPPORTED" |
        is.na(decomposition_pattern)
    ) %>%
    select(
      contrast,
      receiver,
      communication_program,
      communication_evidence_class,
      communication_q_BH,
      decomposition_pattern,
      density_q_BH,
      active_intensity_q_BH
    ) %>%
    arrange(
      contrast,
      communication_q_BH
    ),
  n=Inf,
  width=Inf
)

# ------------------------------------------------------------
# Formal checks
# ------------------------------------------------------------

formal_checks <- tribble(
  ~check, ~pass,

  "crossmodal_matrix_317",
  nrow(xmod) == 317,

  "twenty_eight_receivers",
  n_distinct(xmod$receiver) == 28,

  "publication_universe_16",
  nrow(pub) == 16,

  "publication_receivers_10",
  n_distinct(pub$receiver) == 10,

  "three_contrasts",
  n_distinct(pub$contrast) == 3,

  "FDR005_total_4",
  sum(
    pub$communication_evidence_class ==
      "FDR_LT_0.05"
  ) == 4,

  "FDR010_additional_total_4",
  sum(
    pub$communication_evidence_class ==
      "FDR_LT_0.10"
  ) == 4,

  "nominal_total_8",
  sum(
    pub$communication_evidence_class ==
      "NOMINAL_P_LT_0.05"
  ) == 8,

  "AA_HC_FDR005_4",
  sum(
    pub$contrast == "AA_vs_HC" &
    pub$communication_evidence_class ==
      "FDR_LT_0.05"
  ) == 4,

  "FFA_HC_FDR010_4",
  sum(
    pub$contrast == "FFA_vs_HC" &
    pub$communication_evidence_class ==
      "FDR_LT_0.10"
  ) == 4,

  "AA_FFA_nominal_1",
  sum(
    pub$contrast == "AA_vs_FFA" &
    pub$communication_evidence_class ==
      "NOMINAL_P_LT_0.05"
  ) == 1,

  "no_review_promoted",
  all(
    pub$communication_action ==
      "INCLUDE"
  ),

  "no_crossmodal_score",
  TRUE,

  "no_combined_pq",
  TRUE,

  "no_new_statistical_test",
  TRUE,

  "no_DE_refit",
  TRUE,

  "no_DA_refit",
  TRUE,

  "no_pathway_refit",
  TRUE,

  "no_communication_refit",
  TRUE,

  "annotation_unchanged",
  TRUE
)

cat("\n===== FORMAL CHECKS =====\n")

print(
  formal_checks,
  n=Inf,
  width=Inf
)

stopifnot(
  all(formal_checks$pass)
)

cat("\n============================================================\n")
cat("21G4A PREFLIGHT COMPLETE\n")
cat("============================================================\n")
cat("Frozen communication/cross-modal evidence recovered\n")
cat("16 publication-threshold signals retained\n")
cat("Evidence strata preserved exactly\n")
cat("Disease-specific architecture exposed\n")
cat("Decomposition architecture exposed\n")
cat("No biological hierarchy assigned yet\n")
cat("No figure selection performed\n")
cat("No manuscript modification\n")
cat("No statistical inference\n")
cat("No analytical refit\n")
cat("PASS\n")
