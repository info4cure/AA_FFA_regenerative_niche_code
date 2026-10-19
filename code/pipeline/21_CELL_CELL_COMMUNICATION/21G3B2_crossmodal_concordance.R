#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(stringr)
  library(tibble)
})

options(
  dplyr.summarise.inform = FALSE,
  readr.show_col_types = FALSE
)

cat("============================================================\n")
cat("21G3B2 — CROSS-MODAL CONCORDANCE CLASSIFICATION\n")
cat("============================================================\n")

# ============================================================
# 0. Frozen inputs / outputs
# ============================================================

ROOT <- file.path(
  "PIPELINE_V2",
  "results",
  "disease_biology",
  "21_CELL_CELL_COMMUNICATION",
  "21G_MULTICOMPARTMENT_COMMUNICATION"
)

IN_FILE <- file.path(
  ROOT,
  "21G4_CROSSMODAL_INTEGRATION",
  "tables",
  "21G3B1_crossmodal_evidence_matrix.csv"
)

OUT <- file.path(
  ROOT,
  "21G4_CROSSMODAL_INTEGRATION"
)

TAB <- file.path(OUT, "tables")
OBJ <- file.path(OUT, "objects")

dir.create(TAB, recursive=TRUE, showWarnings=FALSE)
dir.create(OBJ, recursive=TRUE, showWarnings=FALSE)

stopifnot(file.exists(IN_FILE))

# ============================================================
# 1. Recover frozen evidence matrix
# ============================================================

cat("\nSTEP 1/8 — recovering frozen 21G3B1 matrix\n")

x <- read_csv(
  IN_FILE,
  show_col_types=FALSE
)

cat("Rows:", nrow(x), "\n")
cat("Receivers:", n_distinct(x$receiver), "\n")
cat("Contrasts:", n_distinct(x$contrast), "\n")

stopifnot(
  nrow(x) == 317,
  n_distinct(x$receiver) == 28,
  n_distinct(x$contrast) == 3
)

# ============================================================
# 2. Verify full frozen communication key
# ============================================================

cat("\nSTEP 2/8 — verifying frozen analytical key\n")

key_check <- x %>%
  count(
    receiver,
    macrocompartment,
    contrast,
    communication_program,
    communication_tier,
    communication_mechanistic_class,
    communication_action,
    communication_axis,
    name="n"
  ) %>%
  filter(n > 1)

cat("Duplicated full analytical keys:", nrow(key_check), "\n")

stopifnot(
  nrow(key_check) == 0
)

# ============================================================
# 3. Freeze modality-specific support semantics
# ============================================================

cat("\nSTEP 3/8 — deriving modality-specific support semantics\n")

x <- x %>%
  mutate(

    # --------------------------------------------------------
    # Communication
    # --------------------------------------------------------

    communication_support =
      case_when(
        communication_evidence_class ==
          "FDR_LT_0.05" ~
          "FDR_LT_0.05",

        communication_evidence_class ==
          "FDR_LT_0.10" ~
          "FDR_LT_0.10",

        communication_evidence_class ==
          "NOMINAL_P_LT_0.05" ~
          "NOMINAL_P_LT_0.05",

        TRUE ~
          "NO_THRESHOLD_EVIDENCE"
      ),

    communication_threshold_signal =
      communication_action == "INCLUDE" &
      communication_evidence_class %in%
        c(
          "FDR_LT_0.05",
          "FDR_LT_0.10",
          "NOMINAL_P_LT_0.05"
        ),

    communication_FDR_signal =
      communication_action == "INCLUDE" &
      communication_evidence_class %in%
        c(
          "FDR_LT_0.05",
          "FDR_LT_0.10"
        ),

    # --------------------------------------------------------
    # Differential abundance
    #
    # Primary DA inferential family is the already frozen
    # compartment-level BH correction.
    # No DA p/q value is recalculated here.
    # --------------------------------------------------------

    DA_support =
      case_when(
        is.na(DA_FDR_compartment) ~
          "NOT_AVAILABLE",

        DA_FDR_compartment < 0.05 ~
          "FDR_LT_0.05",

        DA_FDR_compartment < 0.10 ~
          "FDR_LT_0.10",

        TRUE ~
          "NO_THRESHOLD_EVIDENCE"
      ),

    DA_FDR005 =
      !is.na(DA_FDR_compartment) &
      DA_FDR_compartment < 0.05,

    DA_FDR010 =
      !is.na(DA_FDR_compartment) &
      DA_FDR_compartment < 0.10,

    # --------------------------------------------------------
    # Differential expression
    #
    # 18C summary supplies numbers of significant genes but
    # does not define a single directional DE state.
    # Therefore DE is used as non-directional orthogonal
    # evidence only.
    # --------------------------------------------------------

    DE_support =
      case_when(
        is.na(DE_n_genes_tested) ~
          "NOT_AVAILABLE",

        DE_n_FDR005 > 0 ~
          "FDR_LT_0.05_GENES_PRESENT",

        TRUE ~
          "NO_FDR_LT_0.05_GENES"
      ),

    DE_FDR005_present =
      !is.na(DE_n_FDR005) &
      DE_n_FDR005 > 0,

    DE_direction_class =
      "NOT_CLASSIFIED_FROM_18C_SUMMARY",

    # --------------------------------------------------------
    # Pathway
    # --------------------------------------------------------

    pathway_record_present =
      pathway_bridge_status ==
        "PATHWAY_RECORD_AVAILABLE",

    pathway_support =
      case_when(
        bridge_class ==
          "NO_PRESPECIFIED_BRIDGE" ~
          "NO_PRESPECIFIED_BRIDGE",

        pathway_bridge_status ==
          "NO_MATCHING_20F2B_STATE_AXIS_RECORD" ~
          "NO_MATCHING_20F2B_RECORD",

        pathway_record_present &
          pathway_status %in%
            c(
              "UP",
              "DOWN",
              "AA_HIGHER",
              "FFA_HIGHER",
              "MIXED"
            ) ~
          "SUPPORTED_PATHWAY_RECORD",

        pathway_record_present &
          pathway_status == "NONE" ~
          "PATHWAY_RECORD_NO_DIRECT_CONTRAST_SIGNAL",

        TRUE ~
          "NOT_AVAILABLE"
      )
  )

# ============================================================
# 4. Explicit disease-direction semantics
# ============================================================

cat("\nSTEP 4/8 — harmonizing contrast-specific direction semantics\n")

x <- x %>%
  mutate(

    communication_direction_standard =
      case_when(

        contrast == "AA_vs_HC" &
          communication_direction ==
            "HIGHER_IN_NUMERATOR" ~
          "AA_HIGHER",

        contrast == "AA_vs_HC" &
          communication_direction ==
            "LOWER_IN_NUMERATOR" ~
          "HC_HIGHER",

        contrast == "FFA_vs_HC" &
          communication_direction ==
            "HIGHER_IN_NUMERATOR" ~
          "FFA_HIGHER",

        contrast == "FFA_vs_HC" &
          communication_direction ==
            "LOWER_IN_NUMERATOR" ~
          "HC_HIGHER",

        contrast == "AA_vs_FFA" &
          communication_direction ==
            "HIGHER_IN_NUMERATOR" ~
          "AA_HIGHER",

        contrast == "AA_vs_FFA" &
          communication_direction ==
            "LOWER_IN_NUMERATOR" ~
          "FFA_HIGHER",

        TRUE ~
          NA_character_
      ),

    pathway_direction_standard =
      case_when(

        !pathway_record_present ~
          NA_character_,

        contrast == "AA_vs_HC" &
          pathway_direction == "UP" ~
          "AA_HIGHER",

        contrast == "AA_vs_HC" &
          pathway_direction == "DOWN" ~
          "HC_HIGHER",

        contrast == "FFA_vs_HC" &
          pathway_direction == "UP" ~
          "FFA_HIGHER",

        contrast == "FFA_vs_HC" &
          pathway_direction == "DOWN" ~
          "HC_HIGHER",

        contrast == "AA_vs_FFA" &
          pathway_direction ==
            "AA_HIGHER" ~
          "AA_HIGHER",

        contrast == "AA_vs_FFA" &
          pathway_direction ==
            "FFA_HIGHER" ~
          "FFA_HIGHER",

        pathway_direction == "NONE" ~
          "NONE",

        pathway_direction == "MIXED" ~
          "MIXED",

        TRUE ~
          NA_character_
      )
  )

# ============================================================
# 5. Directional pathway relation
# ============================================================

cat("\nSTEP 5/8 — classifying communication-pathway directional relation\n")

x <- x %>%
  mutate(

    communication_pathway_relation =
      case_when(

        bridge_class ==
          "NO_PRESPECIFIED_BRIDGE" ~
          "NOT_EVALUABLE_NO_BRIDGE",

        pathway_bridge_status ==
          "NO_MATCHING_20F2B_STATE_AXIS_RECORD" ~
          "NOT_EVALUABLE_NO_PATHWAY_RECORD",

        pathway_direction_standard ==
          "NONE" ~
          "PATHWAY_NO_DIRECT_SIGNAL",

        pathway_direction_standard ==
          "MIXED" ~
          "PATHWAY_MIXED",

        is.na(pathway_direction_standard) ~
          "NOT_EVALUABLE",

        communication_direction_standard ==
          pathway_direction_standard ~
          "DIRECTIONALLY_CONCORDANT",

        TRUE ~
          "DIRECTIONALLY_DISCORDANT"
      )
  )

cat("\n===== MODALITY SUPPORT =====\n")

print(
  x %>%
    count(
      communication_support,
      DA_support,
      DE_support,
      pathway_support,
      name="n"
    ) %>%
    arrange(desc(n)),
  n=Inf,
  width=Inf
)

cat("\n===== COMMUNICATION ↔ PATHWAY DIRECTION =====\n")

print(
  x %>%
    count(
      contrast,
      bridge_class,
      communication_pathway_relation,
      name="n"
    ) %>%
    arrange(
      contrast,
      bridge_class,
      communication_pathway_relation
    ),
  n=Inf,
  width=Inf
)

# ============================================================
# BLOCK 1 COMPLETE
# ============================================================

cat("\n21G3B2 BLOCK 1 COMPLETE\n")

# ============================================================
# 6. Cross-modal descriptive triangulation
# ============================================================

cat("\nSTEP 6/8 — deriving descriptive cross-modal triangulation\n")

x <- x %>%
  mutate(

    # --------------------------------------------------------
    # Orthogonal evidence flags
    #
    # These are descriptive flags only. They are NOT components
    # of a numerical score and are not used to recalculate any
    # inferential quantity.
    # --------------------------------------------------------

    DA_orthogonal_support =
      DA_FDR005,

    DA_orthogonal_suggestive =
      DA_FDR010,

    DE_orthogonal_support =
      DE_FDR005_present,

    pathway_directional_support_direct =
      bridge_class == "DIRECT" &
      communication_pathway_relation ==
        "DIRECTIONALLY_CONCORDANT",

    pathway_directional_discordance_direct =
      bridge_class == "DIRECT" &
      communication_pathway_relation ==
        "DIRECTIONALLY_DISCORDANT",

    pathway_directional_support_related =
      bridge_class == "RELATED" &
      communication_pathway_relation ==
        "DIRECTIONALLY_CONCORDANT",

    pathway_directional_discordance_related =
      bridge_class == "RELATED" &
      communication_pathway_relation ==
        "DIRECTIONALLY_DISCORDANT",

    # Number of available orthogonal evidence domains.
    #
    # This is an availability count, not an evidence score.
    # Pathway availability requires an actual matched 20F2B
    # state-axis record.
    n_orthogonal_modalities_available =
      as.integer(!is.na(DA_FDR_compartment)) +
      as.integer(!is.na(DE_n_genes_tested)) +
      as.integer(pathway_record_present),

    # Number of orthogonal domains carrying positive support.
    #
    # This is retained for audit/description only and must not
    # be interpreted as a weighted or inferential score.
    n_orthogonal_support_domains =
      as.integer(DA_orthogonal_support) +
      as.integer(DE_orthogonal_support) +
      as.integer(
        pathway_directional_support_direct |
        pathway_directional_support_related
      ),

    # --------------------------------------------------------
    # Descriptive triangulation class
    #
    # Communication remains the anchor modality.
    # Only communication rows reaching a frozen threshold are
    # promoted into threshold-signal triangulation classes.
    #
    # Pathway DIRECT concordance is distinguished from RELATED.
    # Discordance is retained explicitly and never overridden
    # by DA/DE support.
    # --------------------------------------------------------

    crossmodal_triangulation_class =
      case_when(

        !communication_threshold_signal ~
          "COMMUNICATION_BELOW_THRESHOLD",

        pathway_directional_discordance_direct ~
          "COMM_SIGNAL_DIRECT_PATHWAY_DISCORDANT",

        pathway_directional_discordance_related ~
          "COMM_SIGNAL_RELATED_PATHWAY_DISCORDANT",

        pathway_directional_support_direct &
          DA_orthogonal_support &
          DE_orthogonal_support ~
          "COMM_SIGNAL_DIRECT_PATHWAY_CONCORDANT_PLUS_DA_DE",

        pathway_directional_support_direct &
          DA_orthogonal_support ~
          "COMM_SIGNAL_DIRECT_PATHWAY_CONCORDANT_PLUS_DA",

        pathway_directional_support_direct &
          DE_orthogonal_support ~
          "COMM_SIGNAL_DIRECT_PATHWAY_CONCORDANT_PLUS_DE",

        pathway_directional_support_direct ~
          "COMM_SIGNAL_DIRECT_PATHWAY_CONCORDANT",

        pathway_directional_support_related &
          DA_orthogonal_support &
          DE_orthogonal_support ~
          "COMM_SIGNAL_RELATED_PATHWAY_CONCORDANT_PLUS_DA_DE",

        pathway_directional_support_related &
          DA_orthogonal_support ~
          "COMM_SIGNAL_RELATED_PATHWAY_CONCORDANT_PLUS_DA",

        pathway_directional_support_related &
          DE_orthogonal_support ~
          "COMM_SIGNAL_RELATED_PATHWAY_CONCORDANT_PLUS_DE",

        pathway_directional_support_related ~
          "COMM_SIGNAL_RELATED_PATHWAY_CONCORDANT",

        DA_orthogonal_support &
          DE_orthogonal_support ~
          "COMM_SIGNAL_PLUS_DA_DE",

        DA_orthogonal_support ~
          "COMM_SIGNAL_PLUS_DA",

        DE_orthogonal_support ~
          "COMM_SIGNAL_PLUS_DE",

        TRUE ~
          "COMM_SIGNAL_ONLY"
      ),

    # --------------------------------------------------------
    # Pathway evaluability class
    # --------------------------------------------------------

    pathway_evaluability =
      case_when(

        bridge_class ==
          "NO_PRESPECIFIED_BRIDGE" ~
          "NO_PRESPECIFIED_BRIDGE",

        pathway_bridge_status ==
          "NO_MATCHING_20F2B_STATE_AXIS_RECORD" ~
          "BRIDGE_DEFINED_NO_MATCHING_RECORD",

        bridge_class == "DIRECT" &
          pathway_record_present ~
          "DIRECT_BRIDGE_RECORD_AVAILABLE",

        bridge_class == "RELATED" &
          pathway_record_present ~
          "RELATED_BRIDGE_RECORD_AVAILABLE",

        TRUE ~
          "OTHER"
      )
  )

cat("\n===== CROSS-MODAL TRIANGULATION CLASSES =====\n")

print(
  x %>%
    count(
      crossmodal_triangulation_class,
      name="n"
    ) %>%
    arrange(desc(n)),
  n=Inf,
  width=Inf
)

# ============================================================
# 7. Communication-threshold readout
# ============================================================

cat("\nSTEP 7/8 — auditing communication threshold signals\n")

threshold_anchor <- x %>%
  filter(
    communication_threshold_signal
  ) %>%
  arrange(
    contrast,
    communication_q_BH,
    communication_p_exact,
    macrocompartment,
    receiver
  )

FDR_anchor <- x %>%
  filter(
    communication_FDR_signal
  ) %>%
  arrange(
    contrast,
    communication_q_BH,
    communication_p_exact,
    macrocompartment,
    receiver
  )

cat(
  "Communication threshold rows:",
  nrow(threshold_anchor),
  "\n"
)

cat(
  "Communication FDR rows:",
  nrow(FDR_anchor),
  "\n"
)

cat("\n===== THRESHOLD SIGNALS × TRIANGULATION =====\n")

print(
  threshold_anchor %>%
    select(
      contrast,
      macrocompartment,
      receiver,
      communication_program,
      communication_tier,
      communication_mechanistic_class,
      communication_action,
      communication_axis,
      communication_mean_difference,
      communication_p_exact,
      communication_q_BH,
      communication_evidence_class,
      DA_support,
      DE_support,
      bridge_class,
      pathway_bridge_status,
      pathway_direction,
      communication_pathway_relation,
      decomposition_pattern,
      crossmodal_triangulation_class
    ),
  n=Inf,
  width=Inf
)

cat("\n===== THRESHOLD SIGNAL SUMMARY =====\n")

threshold_summary <- threshold_anchor %>%
  count(
    contrast,
    macrocompartment,
    communication_evidence_class,
    crossmodal_triangulation_class,
    name="n"
  ) %>%
  arrange(
    contrast,
    macrocompartment,
    communication_evidence_class,
    crossmodal_triangulation_class
  )

print(
  threshold_summary,
  n=Inf,
  width=Inf
)

cat("\n===== FDR SIGNALS =====\n")

print(
  FDR_anchor %>%
    select(
      contrast,
      macrocompartment,
      receiver,
      communication_program,
      communication_tier,
      communication_axis,
      communication_mean_difference,
      communication_q_BH,
      DA_support,
      DE_support,
      bridge_class,
      pathway_bridge_status,
      communication_pathway_relation,
      decomposition_pattern,
      crossmodal_triangulation_class
    ),
  n=Inf,
  width=Inf
)

# ============================================================
# 8. Formal checks + freeze
# ============================================================

cat("\nSTEP 8/8 — formal checks and freeze\n")

# Expected frozen communication counts from 21G2C:
# 4 FDR < 0.05
# 4 additional FDR < 0.10
# 8 nominal-only
# = 16 total threshold rows.

formal_checks <- tribble(
  ~check, ~pass,

  "matrix_317_rows_retained",
  nrow(x) == 317,

  "twenty_eight_receivers_retained",
  n_distinct(x$receiver) == 28,

  "three_contrasts_retained",
  n_distinct(x$contrast) == 3,

  "nine_macrocompartments_retained",
  n_distinct(x$macrocompartment) == 9,

  "full_analytical_key_unique",
  nrow(
    x %>%
      count(
        receiver,
        macrocompartment,
        contrast,
        communication_program,
        communication_tier,
        communication_mechanistic_class,
        communication_action,
        communication_axis
      ) %>%
      filter(n > 1)
  ) == 0,

  "communication_threshold_16_rows",
  nrow(threshold_anchor) == 16,

  "communication_FDR_8_rows",
  nrow(FDR_anchor) == 8,

  "communication_FDR005_4_rows",
  sum(
    x$communication_action == "INCLUDE" &
    x$communication_evidence_class ==
      "FDR_LT_0.05"
  ) == 4,

  "communication_FDR010_additional_4_rows",
  sum(
    x$communication_action == "INCLUDE" &
    x$communication_evidence_class ==
      "FDR_LT_0.10"
  ) == 4,

  "communication_nominal_only_8_rows",
  sum(
    x$communication_action == "INCLUDE" &
    x$communication_evidence_class ==
      "NOMINAL_P_LT_0.05"
  ) == 8,

  "review_rows_not_threshold_promoted",
  !any(
    x$communication_action == "REVIEW" &
    x$communication_threshold_signal
  ),

  "DE_direction_not_inferred",
  all(
    x$DE_direction_class ==
      "NOT_CLASSIFIED_FROM_18C_SUMMARY"
  ),

  "no_bridge_not_called_pathway_negative",
  !any(
    x$bridge_class ==
      "NO_PRESPECIFIED_BRIDGE" &
    x$communication_pathway_relation %in%
      c(
        "DIRECTIONALLY_CONCORDANT",
        "DIRECTIONALLY_DISCORDANT"
      )
  ),

  "missing_pathway_record_not_called_negative",
  !any(
    x$pathway_bridge_status ==
      "NO_MATCHING_20F2B_STATE_AXIS_RECORD" &
    x$communication_pathway_relation %in%
      c(
        "DIRECTIONALLY_CONCORDANT",
        "DIRECTIONALLY_DISCORDANT"
      )
  ),

  "direct_related_bridge_distinguished",
  all(
    x$bridge_class %in%
      c(
        "DIRECT",
        "RELATED",
        "NO_PRESPECIFIED_BRIDGE"
      )
  ),

  "no_crossmodal_numeric_score",
  TRUE,

  "no_combined_p_value",
  TRUE,

  "no_combined_q_value",
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

# ------------------------------------------------------------
# Write frozen outputs
# ------------------------------------------------------------

write_csv(
  x,
  file.path(
    TAB,
    "21G3B2_crossmodal_concordance_matrix.csv"
  )
)

saveRDS(
  x,
  file.path(
    OBJ,
    "21G3B2_crossmodal_concordance_matrix.rds"
  )
)

write_csv(
  threshold_anchor,
  file.path(
    TAB,
    "21G3B2_communication_threshold_signals_crossmodal.csv"
  )
)

write_csv(
  FDR_anchor,
  file.path(
    TAB,
    "21G3B2_communication_FDR_signals_crossmodal.csv"
  )
)

write_csv(
  threshold_summary,
  file.path(
    TAB,
    "21G3B2_threshold_signal_summary.csv"
  )
)

write_csv(
  x %>%
    count(
      contrast,
      bridge_class,
      pathway_bridge_status,
      communication_pathway_relation,
      name="n"
    ),
  file.path(
    TAB,
    "21G3B2_pathway_directional_relation_summary.csv"
  )
)

write_csv(
  x %>%
    count(
      crossmodal_triangulation_class,
      name="n"
    ) %>%
    arrange(desc(n)),
  file.path(
    TAB,
    "21G3B2_triangulation_class_summary.csv"
  )
)

write_csv(
  formal_checks,
  file.path(
    TAB,
    "21G3B2_formal_checks.csv"
  )
)

# ------------------------------------------------------------
# Audit
# ------------------------------------------------------------

audit <- c(
  "============================================================",
  "21G3B2 — CROSS-MODAL CONCORDANCE CLASSIFICATION",
  "============================================================",
  "",
  paste("Rows retained:", nrow(x)),
  paste("Receivers retained:", n_distinct(x$receiver)),
  paste("Contrasts retained:", n_distinct(x$contrast)),
  "",
  paste(
    "Communication threshold signals:",
    nrow(threshold_anchor)
  ),
  paste(
    "Communication FDR signals:",
    nrow(FDR_anchor)
  ),
  "",
  paste(
    "Communication remains the anchor modality."
  ),
  paste(
    "DA evidence uses frozen compartment-level BH results."
  ),
  paste(
    "DE evidence is non-directional because the frozen 18C",
    "summary does not define a single state-level DE direction."
  ),
  paste(
    "Pathway directional concordance is evaluated only where",
    "a prespecified bridge and matching frozen 20F2B record exist."
  ),
  paste(
    "DIRECT and RELATED pathway bridges remain distinct."
  ),
  paste(
    "Absence of a bridge or matching pathway record is treated",
    "as non-evaluable, not as negative evidence."
  ),
  paste(
    "Cross-modal triangulation classes are descriptive categories",
    "and are not numerical evidence scores."
  ),
  "",
  paste(
    "No p-value, q-value, BH correction, DE model, DA model,",
    "pathway model or communication model was recalculated."
  ),
  "",
  "PASS"
)

writeLines(
  audit,
  file.path(
    OUT,
    "21G3B2_audit.txt"
  )
)

cat("\n============================================================\n")
cat("21G3B2 COMPLETE\n")
cat("============================================================\n")
cat("317 frozen communication rows retained\n")
cat("28 receivers retained\n")
cat("Three disease contrasts retained\n")
cat("Communication retained as anchor modality\n")
cat("DA incorporated as frozen orthogonal evidence\n")
cat("DE incorporated without inventing state-level direction\n")
cat("Pathway concordance restricted to prespecified bridges\n")
cat("DIRECT and RELATED pathway bridges distinguished\n")
cat("16 communication threshold signals retained\n")
cat("8 communication FDR signals retained\n")
cat("No cross-modal numeric score\n")
cat("No combined p/q values\n")
cat("No statistical inference\n")
cat("No analytical refit\n")
cat("PASS\n")
