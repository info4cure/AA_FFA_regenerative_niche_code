suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(stringr)
  library(tibble)
})

cat("============================================================\n")
cat("21G4C2 — CANDIDATE FIGURE MANIFEST\n")
cat("============================================================\n")

# ------------------------------------------------------------
# 1. Paths
# ------------------------------------------------------------

BASE <- paste0(
  "PIPELINE_V2/results/disease_biology/",
  "21_CELL_CELL_COMMUNICATION/",
  "21G_MULTICOMPARTMENT_COMMUNICATION/",
  "21G4_CROSSMODAL_INTEGRATION"
)

TAB <- file.path(BASE, "tables")

SRC <- file.path(
  TAB,
  "21G4B_biological_synthesis_signals.csv"
)

stopifnot(file.exists(SRC))

# ------------------------------------------------------------
# 2. Recover canonical frozen synthesis
# ------------------------------------------------------------

cat("\nSTEP 1/6 — recovering frozen 21G4B synthesis\n")

x <- read_csv(
  SRC,
  show_col_types = FALSE
)

cat("Rows:", nrow(x), "\n")
cat("Columns:", ncol(x), "\n")

stopifnot(
  nrow(x) == 16
)

cat("\n===== AVAILABLE COLUMNS =====\n")
print(names(x))

# ------------------------------------------------------------
# 3. Inspect biological hierarchy columns
# ------------------------------------------------------------

cat("\nSTEP 2/6 — inspecting frozen adjudication fields\n")

candidate_hierarchy_cols <- c(
  "biological_pattern",
  "biological_architecture",
  "biological_role",
  "publication_role",
  "synthesis_class",
  "figure_role",
  "evidence",
  "communication_evidence_class",
  "crossmodal_triangulation_class",
  "decomposition_pattern"
)

hierarchy_cols <- intersect(
  candidate_hierarchy_cols,
  names(x)
)

cat("\n===== FROZEN HIERARCHY / EVIDENCE FIELDS =====\n")
print(hierarchy_cols)

for (v in hierarchy_cols) {
  cat("\n---", v, "---\n")
  print(
    x %>%
      count(.data[[v]], name = "n") %>%
      arrange(desc(n))
  )
}

# ------------------------------------------------------------
# 4. Construct lossless candidate manifest
# ------------------------------------------------------------

cat("\nSTEP 3/6 — constructing lossless candidate manifest\n")

manifest <- x %>%
  mutate(
    signal_id = sprintf(
      "21G4C2_%02d",
      row_number()
    ),

    figure_eligibility =
      "CANDIDATE_PENDING_MANUAL_ADJUDICATION",

    figure_selection =
      FALSE,

    proposed_panel =
      NA_character_,

    proposed_display_order =
      NA_integer_,

    selection_rationale =
      NA_character_,

    exclusion_rationale =
      NA_character_
  ) %>%
  relocate(
    signal_id,
    figure_eligibility,
    figure_selection,
    proposed_panel,
    proposed_display_order,
    selection_rationale,
    exclusion_rationale
  )

stopifnot(
  nrow(manifest) == nrow(x)
)

# ------------------------------------------------------------
# 5. Produce descriptive candidate summaries
# ------------------------------------------------------------

cat("\nSTEP 4/6 — descriptive candidate summaries\n")

required_summary_cols <- c(
  "contrast",
  "macrocompartment",
  "receiver"
)

stopifnot(
  all(required_summary_cols %in% names(manifest))
)

contrast_summary <- manifest %>%
  count(
    contrast,
    name = "n_signals"
  ) %>%
  arrange(contrast)

receiver_summary <- manifest %>%
  count(
    contrast,
    macrocompartment,
    receiver,
    name = "n_signals"
  ) %>%
  arrange(
    contrast,
    macrocompartment,
    receiver
  )

cat("\n===== SIGNALS BY CONTRAST =====\n")
print(
  contrast_summary,
  n = Inf,
  width = Inf
)

cat("\n===== SIGNALS BY RECEIVER =====\n")
print(
  receiver_summary,
  n = Inf,
  width = Inf
)

# ------------------------------------------------------------
# 6. Print compact adjudication view
# ------------------------------------------------------------

cat("\nSTEP 5/6 — candidate adjudication view\n")

preferred_cols <- c(
  "signal_id",
  "contrast",
  "macrocompartment",
  "receiver",
  "communication_program",
  "communication_tier",
  "communication_axis",
  "communication_evidence_class",
  "communication_mean_difference",
  "communication_q_BH",
  "DA_OR",
  "DA_FDR_compartment",
  "DA_support",
  "DE_n_FDR005",
  "DE_support",
  "decomposition_pattern",
  "density_mean_difference",
  "density_q_BH",
  "crossmodal_triangulation_class",
  hierarchy_cols,
  "figure_eligibility",
  "figure_selection"
)

preferred_cols <- unique(
  preferred_cols[
    preferred_cols %in% names(manifest)
  ]
)

cat("\n===== CANDIDATE FIGURE MANIFEST =====\n")

print(
  manifest %>%
    select(all_of(preferred_cols)),
  n = Inf,
  width = Inf
)

# ------------------------------------------------------------
# 7. Write manifest — NO FIGURES
# ------------------------------------------------------------

cat("\nSTEP 6/6 — writing candidate manifest\n")

write_csv(
  manifest,
  file.path(
    TAB,
    "21G4C2_candidate_figure_manifest.csv"
  )
)

write_csv(
  contrast_summary,
  file.path(
    TAB,
    "21G4C2_candidate_figure_manifest_by_contrast.csv"
  )
)

write_csv(
  receiver_summary,
  file.path(
    TAB,
    "21G4C2_candidate_figure_manifest_by_receiver.csv"
  )
)

formal_checks <- tribble(
  ~check, ~pass,

  "canonical_21G4B_source_present",
  file.exists(SRC),

  "sixteen_signals_retained",
  nrow(manifest) == 16,

  "three_contrasts_retained",
  n_distinct(manifest$contrast) == 3,

  "no_signal_preselected",
  all(!manifest$figure_selection),

  "all_candidates_pending_manual_adjudication",
  all(
    manifest$figure_eligibility ==
      "CANDIDATE_PENDING_MANUAL_ADJUDICATION"
  ),

  "no_panel_assigned",
  all(is.na(manifest$proposed_panel)),

  "no_display_order_assigned",
  all(is.na(manifest$proposed_display_order)),

  "no_threshold_refit",
  TRUE,

  "no_crossmodal_score",
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

  "no_biological_reclassification",
  TRUE,

  "no_figure_generated",
  TRUE,

  "no_manuscript_modified",
  TRUE
)

write_csv(
  formal_checks,
  file.path(
    TAB,
    "21G4C2_formal_checks.csv"
  )
)

cat("\n===== FORMAL CHECKS =====\n")

print(
  formal_checks,
  n = Inf,
  width = Inf
)

stopifnot(
  all(formal_checks$pass)
)

audit <- c(
  "21G4C2 — CANDIDATE FIGURE MANIFEST",
  "",
  paste("Canonical source:", SRC),
  paste("Signals retained:", nrow(manifest)),
  paste(
    "Receivers:",
    n_distinct(manifest$receiver)
  ),
  paste(
    "Contrasts:",
    n_distinct(manifest$contrast)
  ),
  "",
  "DESIGN",
  "- 21G4B adjudicated synthesis retained as canonical source.",
  "- All 16 publication-threshold signals retained.",
  "- No signal automatically selected.",
  "- No panel automatically assigned.",
  "- No display order automatically assigned.",
  "- Figure selection requires manual biological adjudication.",
  "",
  "ANALYTICAL INTEGRITY",
  "- no threshold refit",
  "- no cross-modal numeric score",
  "- no combined p/q values",
  "- no new statistical inference",
  "- no DE refit",
  "- no DA refit",
  "- no pathway refit",
  "- no communication refit",
  "- no biological reclassification",
  "- no figures generated",
  "- no manuscript modification",
  "",
  "PASS"
)

writeLines(
  audit,
  file.path(
    BASE,
    "21G4C2_audit.txt"
  )
)

cat("\n============================================================\n")
cat("21G4C2 COMPLETE\n")
cat("============================================================\n")
cat("21G4B adjudicated synthesis retained as canonical source\n")
cat("16 signals retained in candidate figure manifest\n")
cat("No signal preselected\n")
cat("No panel assigned\n")
cat("No display order assigned\n")
cat("Manual adjudication required before plotting\n")
cat("No figure generated\n")
cat("No manuscript modification\n")
cat("No statistical inference\n")
cat("No analytical refit\n")
cat("PASS\n")
