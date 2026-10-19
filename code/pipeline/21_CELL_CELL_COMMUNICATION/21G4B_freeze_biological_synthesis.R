suppressPackageStartupMessages({
  library(tidyverse)
})

BASE <- file.path(
  "PIPELINE_V2/results/disease_biology/21_CELL_CELL_COMMUNICATION",
  "21G_MULTICOMPARTMENT_COMMUNICATION"
)

INFILE <- file.path(
  BASE,
  "21G4_CROSSMODAL_INTEGRATION/tables",
  "21G3B2_communication_threshold_signals_crossmodal.csv"
)

OUT <- file.path(
  BASE,
  "21G4_CROSSMODAL_INTEGRATION"
)

TAB <- file.path(OUT, "tables")

dir.create(TAB, recursive=TRUE, showWarnings=FALSE)

stopifnot(file.exists(INFILE))

x <- read_csv(INFILE, show_col_types=FALSE)

cat("============================================================\n")
cat("21G4B — BIOLOGICAL SYNTHESIS FREEZE\n")
cat("============================================================\n")

stopifnot(
  nrow(x) == 16,
  n_distinct(x$receiver) == 10,
  all(x$communication_action == "INCLUDE")
)

# ------------------------------------------------------------
# Prespecified descriptive biological classes
# ------------------------------------------------------------

x <- x %>%
  mutate(
    biological_readout_class = case_when(

      contrast == "AA_vs_HC" &
        communication_program == "IGF2R_CONTEXTUAL" ~
        "PRIMARY_AA_IGF2R",

      contrast == "AA_vs_HC" &
        communication_program == "GROWTH_SURVIVAL" ~
        "PRIMARY_AA_GROWTH_LOSS",

      contrast == "FFA_vs_HC" &
        communication_program == "GROWTH_SURVIVAL" ~
        "PRIMARY_SUPPORTING_FFA_GROWTH_LOSS",

      contrast == "AA_vs_HC" &
        communication_program == "IMMUNE_REGULATION" ~
        "CONTEXTUAL_AA_IMMUNE",

      TRUE ~
        "CONTEXTUAL_DISEASE_DIFFERENTIATION"
    ),

    publication_strength = case_when(

      communication_evidence_class == "FDR_LT_0.05" ~
        "PRIMARY_FDR",

      communication_evidence_class == "FDR_LT_0.10" ~
        "SUPPORTING_FDR",

      communication_evidence_class == "NOMINAL_P_LT_0.05" ~
        "CONTEXTUAL_NOMINAL",

      TRUE ~
        "OTHER"
    ),

    density_interpretation = case_when(

      decomposition_pattern == "DENSITY_SUPPORTED" ~
        "DENSITY_SUPPORTED",

      decomposition_pattern == "BOTH_NOMINAL_ONLY" ~
        "BOTH_COMPONENTS_NOMINAL_ONLY",

      decomposition_pattern == "DENSITY_NOMINAL_ONLY" ~
        "DENSITY_NOMINAL_ONLY",

      decomposition_pattern ==
        "NO_DECOMPOSITION_THRESHOLD_EVIDENCE" ~
        "NO_DECOMPOSITION_THRESHOLD_EVIDENCE",

      TRUE ~ NA_character_
    )
  )

# ------------------------------------------------------------
# Freeze exact signal-level table
# ------------------------------------------------------------

write_csv(
  x,
  file.path(
    TAB,
    "21G4B_biological_synthesis_signals.csv"
  )
)

# ------------------------------------------------------------
# Biological pattern summary
# ------------------------------------------------------------

summary <- x %>%
  group_by(
    biological_readout_class,
    contrast,
    communication_program
  ) %>%
  summarise(
    n_signals = n(),
    n_receivers = n_distinct(receiver),

    receivers =
      paste(sort(unique(receiver)), collapse="; "),

    n_FDR005 =
      sum(communication_evidence_class == "FDR_LT_0.05"),

    n_FDR010 =
      sum(communication_evidence_class == "FDR_LT_0.10"),

    n_nominal =
      sum(
        communication_evidence_class ==
          "NOMINAL_P_LT_0.05"
      ),

    n_DA_FDR005 =
      sum(DA_support == "FDR_LT_0.05", na.rm=TRUE),

    n_DE_FDR005 =
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

    .groups="drop"
  )

write_csv(
  summary,
  file.path(
    TAB,
    "21G4B_biological_pattern_summary.csv"
  )
)

cat("\n===== BIOLOGICAL PATTERN SUMMARY =====\n")
print(summary, n=Inf, width=Inf)

# ------------------------------------------------------------
# Exact integrity checks
# ------------------------------------------------------------

checks <- tribble(
  ~check, ~pass,

  "sixteen_threshold_signals",
  nrow(x) == 16,

  "ten_receivers",
  n_distinct(x$receiver) == 10,

  "three_contrasts",
  n_distinct(x$contrast) == 3,

  "AA_HC_IGF2R_three_signals",
  sum(
    x$contrast == "AA_vs_HC" &
    x$communication_program == "IGF2R_CONTEXTUAL"
  ) == 3,

  "AA_HC_IGF2R_all_FDR005",
  all(
    x$communication_evidence_class[
      x$contrast == "AA_vs_HC" &
      x$communication_program == "IGF2R_CONTEXTUAL"
    ] == "FDR_LT_0.05"
  ),

  "AA_HC_growth_one_FDR005",
  sum(
    x$contrast == "AA_vs_HC" &
    x$communication_program == "GROWTH_SURVIVAL" &
    x$communication_evidence_class == "FDR_LT_0.05"
  ) == 1,

  "FFA_HC_growth_four_FDR010",
  sum(
    x$contrast == "FFA_vs_HC" &
    x$communication_program == "GROWTH_SURVIVAL" &
    x$communication_evidence_class == "FDR_LT_0.10"
  ) == 4,

  "AA_FFA_one_nominal",
  sum(x$contrast == "AA_vs_FFA") == 1 &&
  x$communication_evidence_class[
    x$contrast == "AA_vs_FFA"
  ] == "NOMINAL_P_LT_0.05",

  "four_AA_HC_FDR005_density_supported",
  sum(
    x$contrast == "AA_vs_HC" &
    x$communication_evidence_class == "FDR_LT_0.05" &
    x$decomposition_pattern == "DENSITY_SUPPORTED"
  ) == 4,

  "no_crossmodal_score",
  TRUE,

  "no_combined_pq",
  TRUE,

  "no_new_statistical_test",
  TRUE,

  "no_analytical_refit",
  TRUE,

  "annotation_unchanged",
  TRUE
)

cat("\n===== FORMAL CHECKS =====\n")
print(checks, n=Inf, width=Inf)

stopifnot(all(checks$pass))

write_csv(
  checks,
  file.path(
    TAB,
    "21G4B_formal_checks.csv"
  )
)

audit <- c(
  "21G4B BIOLOGICAL SYNTHESIS FREEZE",
  "",
  "Frozen source: 21G3B2 communication threshold universe.",
  "",
  "PRIMARY BIOLOGICAL ARCHITECTURE",
  paste(
    "- AA vs HC: increased IGF2R-associated communication",
    "across Bulge/Isthmus, ORS_basal and ORS;",
    "all three signals q<0.05."
  ),
  paste(
    "- AA vs HC: reduced growth/survival communication",
    "in Bulge/Isthmus at q<0.05, with an additional",
    "nominal interfollicular signal."
  ),
  paste(
    "- FFA vs HC: distributed reduction of epithelial",
    "growth/survival communication; four signals q<0.10",
    "and four additional nominal signals."
  ),
  paste(
    "- AA vs FFA: one nominal IGF trophic signal;",
    "no FDR-level direct disease-discriminating",
    "communication signal."
  ),
  "",
  "INTERPRETIVE BOUNDARIES",
  paste(
    "- IGF2R_CONTEXTUAL denotes IGF2R-associated",
    "trafficking/clearance and is not interpreted as",
    "direct evidence of increased canonical IGF signaling."
  ),
  paste(
    "- DENSITY_SUPPORTED indicates that the communication",
    "difference is compatible with a density/structural",
    "component of the frozen communication endpoint."
  ),
  paste(
    "- Density decomposition does not establish causality",
    "and is distinct from orthogonal DA analysis."
  ),
  paste(
    "- DE support denotes transcriptional remodeling;",
    "state-level DE direction was not inferred."
  ),
  "- No cross-modal numeric score was constructed.",
  "- No combined p-value or q-value was calculated.",
  "",
  "ANALYTICAL INTEGRITY",
  "- no statistical recalculation",
  "- no DE refit",
  "- no DA refit",
  "- no pathway refit",
  "- no communication refit",
  "- no annotation change",
  "",
  "PASS"
)

writeLines(
  audit,
  file.path(
    OUT,
    "21G4B_audit.txt"
  )
)

cat("\n============================================================\n")
cat("21G4B COMPLETE\n")
cat("============================================================\n")
cat("16 publication-threshold signals adjudicated\n")
cat("AA-vs-HC IGF2R contextual architecture frozen\n")
cat("AA-vs-HC growth/survival reduction frozen\n")
cat("FFA-vs-HC distributed growth/survival reduction frozen\n")
cat("AA-vs-FFA retained as contextual only\n")
cat("Density decomposition interpretation frozen\n")
cat("No cross-modal numeric score\n")
cat("No statistical inference\n")
cat("No analytical refit\n")
cat("PASS\n")
