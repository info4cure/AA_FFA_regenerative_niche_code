suppressPackageStartupMessages({
  library(tidyverse)
})

BASE <- paste0(
  "PIPELINE_V2/results/disease_biology/",
  "21_CELL_CELL_COMMUNICATION/",
  "21G_MULTICOMPARTMENT_COMMUNICATION/",
  "21G4_CROSSMODAL_INTEGRATION"
)

TAB <- file.path(BASE, "tables")
OBJ <- file.path(BASE, "objects")

dir.create(TAB, recursive=TRUE, showWarnings=FALSE)
dir.create(OBJ, recursive=TRUE, showWarnings=FALSE)

INFILE <- file.path(
  TAB,
  "21G3B2_communication_threshold_signals_crossmodal.csv"
)

stopifnot(file.exists(INFILE))

cat("============================================================\n")
cat("21G3B3 — PUBLICATION READOUT FREEZE\n")
cat("============================================================\n")

x <- read_csv(INFILE, show_col_types=FALSE)

stopifnot(
  nrow(x) == 16,
  n_distinct(x$receiver) == 10,
  n_distinct(x$contrast) == 3,
  all(x$communication_action == "INCLUDE"),
  all(x$communication_threshold_signal)
)

# ------------------------------------------------------------
# 1. Preserve all threshold signals
# ------------------------------------------------------------

publication_readout <- x %>%
  arrange(
    factor(
      contrast,
      levels=c("AA_vs_HC", "FFA_vs_HC", "AA_vs_FFA")
    ),
    communication_q_BH,
    receiver,
    communication_program
  )

write_csv(
  publication_readout,
  file.path(TAB, "21G3B3_publication_readout_all_threshold.csv")
)

# ------------------------------------------------------------
# 2. FDR publication subset
# ------------------------------------------------------------

fdr_readout <- publication_readout %>%
  filter(
    communication_evidence_class %in%
      c("FDR_LT_0.05", "FDR_LT_0.10")
  )

stopifnot(
  nrow(fdr_readout) == 8,
  sum(fdr_readout$communication_evidence_class == "FDR_LT_0.05") == 4,
  sum(fdr_readout$communication_evidence_class == "FDR_LT_0.10") == 4
)

write_csv(
  fdr_readout,
  file.path(TAB, "21G3B3_publication_readout_FDR.csv")
)

# ------------------------------------------------------------
# 3. Nominal contextual subset
# ------------------------------------------------------------

nominal_readout <- publication_readout %>%
  filter(
    communication_evidence_class == "NOMINAL_P_LT_0.05"
  )

stopifnot(nrow(nominal_readout) == 8)

write_csv(
  nominal_readout,
  file.path(TAB, "21G3B3_publication_readout_nominal_contextual.csv")
)

# ------------------------------------------------------------
# 4. Contrast-level frozen summary
# ------------------------------------------------------------

contrast_summary <- publication_readout %>%
  group_by(contrast) %>%
  summarise(
    n_threshold = n(),

    n_FDR005 =
      sum(communication_evidence_class == "FDR_LT_0.05"),

    n_FDR010 =
      sum(communication_evidence_class == "FDR_LT_0.10"),

    n_nominal =
      sum(communication_evidence_class == "NOMINAL_P_LT_0.05"),

    n_DA_FDR005 =
      sum(DA_FDR005),

    n_DA_FDR010 =
      sum(DA_FDR010),

    n_DE_FDR005_present =
      sum(DE_FDR005_present),

    n_density_supported =
      sum(
        decomposition_pattern == "DENSITY_SUPPORTED",
        na.rm=TRUE
      ),

    n_both_nominal =
      sum(
        decomposition_pattern == "BOTH_NOMINAL_ONLY",
        na.rm=TRUE
      ),

    n_no_decomposition_threshold =
      sum(
        decomposition_pattern ==
          "NO_DECOMPOSITION_THRESHOLD_EVIDENCE",
        na.rm=TRUE
      ),

    .groups="drop"
  )

write_csv(
  contrast_summary,
  file.path(TAB, "21G3B3_publication_readout_by_contrast.csv")
)

cat("\n===== CONTRAST SUMMARY =====\n")
print(contrast_summary, n=Inf, width=Inf)

# ------------------------------------------------------------
# 5. Mechanistic readout
#
# Descriptive only.
# No causal interpretation.
# ------------------------------------------------------------

mechanistic_readout <- publication_readout %>%
  transmute(
    contrast,
    macrocompartment,
    receiver,

    communication_program,
    communication_tier,
    communication_axis,

    communication_mean_difference,
    communication_log2_mean_ratio,
    communication_p_exact,
    communication_q_BH,
    communication_evidence_class,
    communication_direction,

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

    crossmodal_triangulation_class,

    interpretation_boundary =
      paste(
        "DA represents orthogonal abundance evidence;",
        "DE represents state-level transcriptional remodeling;",
        "neither is interpreted as directional molecular",
        "concordance with communication."
      )
  )

write_csv(
  mechanistic_readout,
  file.path(TAB, "21G3B3_mechanistic_publication_readout.csv")
)

# ------------------------------------------------------------
# 6. Formal checks
# ------------------------------------------------------------

formal_checks <- tribble(
  ~check, ~pass,

  "threshold_universe_16",
  nrow(publication_readout) == 16,

  "ten_receivers",
  n_distinct(publication_readout$receiver) == 10,

  "three_contrasts",
  n_distinct(publication_readout$contrast) == 3,

  "all_include",
  all(publication_readout$communication_action == "INCLUDE"),

  "FDR_total_8",
  nrow(fdr_readout) == 8,

  "FDR005_4",
  sum(
    publication_readout$communication_evidence_class ==
      "FDR_LT_0.05"
  ) == 4,

  "FDR010_additional_4",
  sum(
    publication_readout$communication_evidence_class ==
      "FDR_LT_0.10"
  ) == 4,

  "nominal_only_8",
  nrow(nominal_readout) == 8,

  "AA_HC_FDR005_4",
  sum(
    publication_readout$contrast == "AA_vs_HC" &
    publication_readout$communication_evidence_class ==
      "FDR_LT_0.05"
  ) == 4,

  "FFA_HC_FDR010_4",
  sum(
    publication_readout$contrast == "FFA_vs_HC" &
    publication_readout$communication_evidence_class ==
      "FDR_LT_0.10"
  ) == 4,

  "AA_FFA_nominal_1",
  sum(
    publication_readout$contrast == "AA_vs_FFA" &
    publication_readout$communication_evidence_class ==
      "NOMINAL_P_LT_0.05"
  ) == 1,

  "no_crossmodal_score",
  TRUE,

  "DA_not_directional_concordance",
  TRUE,

  "DE_direction_not_inferred",
  TRUE,

  "no_combined_p",
  TRUE,

  "no_combined_q",
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

write_csv(
  formal_checks,
  file.path(TAB, "21G3B3_formal_checks.csv")
)

cat("\n===== FORMAL CHECKS =====\n")
print(formal_checks, n=Inf, width=Inf)

stopifnot(all(formal_checks$pass))

# ------------------------------------------------------------
# 7. Frozen object
# ------------------------------------------------------------

saveRDS(
  list(
    publication_readout = publication_readout,
    fdr_readout = fdr_readout,
    nominal_readout = nominal_readout,
    contrast_summary = contrast_summary,
    mechanistic_readout = mechanistic_readout,
    formal_checks = formal_checks
  ),
  file.path(OBJ, "21G3B3_publication_readout.rds")
)

# ------------------------------------------------------------
# 8. Audit
# ------------------------------------------------------------

audit <- c(
  "21G3B3 — PUBLICATION READOUT FREEZE",
  "",
  "FROZEN INPUT",
  paste("Threshold communication signals:", nrow(publication_readout)),
  paste("Receivers:", n_distinct(publication_readout$receiver)),
  paste("Contrasts:", n_distinct(publication_readout$contrast)),
  "",
  "COMMUNICATION EVIDENCE",
  "FDR <0.05 signals: 4",
  "Additional FDR <0.10 signals: 4",
  "Nominal-only signals: 8",
  "",
  "CONTRAST ARCHITECTURE",
  "AA vs HC: 4 FDR <0.05 communication signals",
  "FFA vs HC: 4 additional FDR <0.10 communication signals",
  "AA vs FFA: 1 nominal communication signal",
  "",
  "INTERPRETIVE BOUNDARIES",
  paste(
    "- DA is treated as orthogonal abundance evidence,",
    "not directional molecular concordance."
  ),
  paste(
    "- DE is treated as evidence of transcriptional remodeling",
    "without inferring state-level direction from summary counts."
  ),
  paste(
    "- Communication decomposition is descriptive and does not",
    "establish causal sender-receiver mechanisms."
  ),
  "- No cross-modal numeric score was constructed.",
  "- No combined p-value or q-value was calculated.",
  "",
  "ANALYTICAL INTEGRITY",
  "- no new statistical inference",
  "- no DE refit",
  "- no DA refit",
  "- no pathway refit",
  "- no communication refit",
  "- no annotation changes",
  "",
  "PASS"
)

writeLines(
  audit,
  file.path(BASE, "21G3B3_audit.txt")
)

cat("\n============================================================\n")
cat("21G3B3 COMPLETE\n")
cat("============================================================\n")
cat("16 threshold communication signals frozen\n")
cat("8 FDR communication signals frozen\n")
cat("4 AA-vs-HC FDR<0.05 signals\n")
cat("4 FFA-vs-HC additional FDR<0.10 signals\n")
cat("8 nominal contextual signals retained\n")
cat("DA retained as orthogonal abundance evidence\n")
cat("DE retained as transcriptional-remodeling evidence\n")
cat("No cross-modal numeric score\n")
cat("No combined p/q values\n")
cat("No statistical inference\n")
cat("No analytical refit\n")
cat("PASS\n")
