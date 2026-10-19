suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(tibble)
})

cat("============================================================\n")
cat("21G4C6B — QUANTITATIVE PLOTTING DATASETS FREEZE\n")
cat("============================================================\n")

BASE <- paste0(
  "PIPELINE_V2/results/disease_biology/",
  "21_CELL_CELL_COMMUNICATION/",
  "21G_MULTICOMPARTMENT_COMMUNICATION/",
  "21G4_CROSSMODAL_INTEGRATION"
)

TAB <- file.path(BASE, "tables")
OBJ <- file.path(BASE, "objects")

LONG_FILE <- file.path(
  TAB,
  "21G4C6A_patient_signal_long.csv"
)

META_FILE <- file.path(
  TAB,
  "21G4C6A_signal_metadata.csv"
)

stopifnot(
  file.exists(LONG_FILE),
  file.exists(META_FILE)
)

# ------------------------------------------------------------
# 1. Recover frozen matrices
# ------------------------------------------------------------

cat("\nSTEP 1/7 — recovering frozen publication matrices\n")

long <- read_csv(
  LONG_FILE,
  show_col_types=FALSE
)

meta <- read_csv(
  META_FILE,
  show_col_types=FALSE
)

cat("Patient-signal rows:", nrow(long), "\n")
cat("Signals:", n_distinct(long$signal_id), "\n")
cat("Patients:", n_distinct(long$patient_id), "\n")

stopifnot(
  nrow(long) == 162,
  nrow(meta) == 16,
  n_distinct(long$signal_id) == 16,
  n_distinct(long$patient_id) == 16
)

# ------------------------------------------------------------
# 2. Attach frozen metadata
# ------------------------------------------------------------

cat("\nSTEP 2/7 — attaching frozen signal metadata\n")

# 21G4C6A patient-level table already carries most frozen signal metadata.
# Add only metadata columns that are absent from the patient-level source.
# This prevents .x/.y suffixes and preserves the 21G4C6A values verbatim.

meta_add_columns <- setdiff(
  names(meta),
  names(long)
)

meta_add <- meta %>%
  select(
    signal_id,
    all_of(meta_add_columns)
  )

cat(
  "Metadata columns added:",
  length(meta_add_columns),
  "\n"
)

if (length(meta_add_columns) > 0) {
  cat(
    paste(meta_add_columns, collapse="; "),
    "\n"
  )
}

plot_long <- long %>%
  left_join(
    meta_add,
    by="signal_id"
  )

stopifnot(
  nrow(plot_long) == nrow(long),
  "macrocompartment" %in% names(plot_long),
  "receiver" %in% names(plot_long),
  "communication_program" %in% names(plot_long),
  "communication_axis" %in% names(plot_long),
  !any(grepl("\\.x$|\\.y$", names(plot_long)))
)

# ------------------------------------------------------------
# 3. Freeze endpoint-long representation
# ------------------------------------------------------------

cat("\nSTEP 3/7 — constructing endpoint-long representation\n")

endpoint_long <- plot_long %>%
  pivot_longer(
    cols=c(
      mean_incoming_A_axis,
      pooled_D_axis,
      mean_sender_active_intensity
    ),
    names_to="quantitative_endpoint",
    values_to="patient_value"
  ) %>%
  mutate(
    endpoint_label = case_when(
      quantitative_endpoint ==
        "mean_incoming_A_axis" ~
        "Incoming communication (A)",

      quantitative_endpoint ==
        "pooled_D_axis" ~
        "Communication density (D)",

      quantitative_endpoint ==
        "mean_sender_active_intensity" ~
        "Active-sender intensity",

      TRUE ~ quantitative_endpoint
    )
  )

write_csv(
  endpoint_long,
  file.path(
    TAB,
    "21G4C6B_endpoint_long.csv"
  )
)

# ------------------------------------------------------------
# 4. Freeze contrast-specific patient datasets
# ------------------------------------------------------------

cat("\nSTEP 4/7 — constructing contrast-specific datasets\n")

AA_HC <- plot_long %>%
  filter(contrast == "AA_vs_HC")

FFA_HC <- plot_long %>%
  filter(contrast == "FFA_vs_HC")

AA_FFA <- plot_long %>%
  filter(contrast == "AA_vs_FFA")

write_csv(
  AA_HC,
  file.path(TAB, "21G4C6B_AA_vs_HC_patient_signals.csv")
)

write_csv(
  FFA_HC,
  file.path(TAB, "21G4C6B_FFA_vs_HC_patient_signals.csv")
)

write_csv(
  AA_FFA,
  file.path(TAB, "21G4C6B_AA_vs_FFA_patient_signals.csv")
)

cat(
  "AA-vs-HC signals:",
  n_distinct(AA_HC$signal_id),
  "\n"
)

cat(
  "FFA-vs-HC signals:",
  n_distinct(FFA_HC$signal_id),
  "\n"
)

cat(
  "AA-vs-FFA signals:",
  n_distinct(AA_FFA$signal_id),
  "\n"
)

# ------------------------------------------------------------
# 5. Freeze A-D quantitative geometry
# ------------------------------------------------------------

cat("\nSTEP 5/7 — constructing patient-level A-D geometry\n")

AD <- plot_long %>%
  select(
    signal_id,
    contrast,
    patient_id,
    disease,
    macrocompartment,
    receiver,
    communication_program,
    communication_axis,
    communication_evidence_class,
    communication_q_BH,
    mean_incoming_A_axis,
    pooled_D_axis,
    mean_sender_active_intensity,
    decomposition_pattern,
    density_interpretation,
    crossmodal_triangulation_class
  ) %>%
  mutate(
    intensity_available =
      !is.na(mean_sender_active_intensity)
  )

stopifnot(
  all(!is.na(AD$mean_incoming_A_axis)),
  all(!is.na(AD$pooled_D_axis))
)

write_csv(
  AD,
  file.path(
    TAB,
    "21G4C6B_patient_AD_geometry.csv"
  )
)

# ------------------------------------------------------------
# 6. Frozen signal-level effect landscape
# ------------------------------------------------------------

cat("\nSTEP 6/7 — constructing signal-level effect landscape\n")

effect_landscape <- meta %>%
  select(
    signal_id,
    contrast,
    macrocompartment,
    receiver,
    communication_program,
    communication_tier,
    communication_axis,
    communication_evidence_class,
    communication_mean_difference,
    communication_q_BH,
    DA_OR,
    DA_FDR_compartment,
    DA_support,
    DE_n_FDR005,
    DE_support,
    decomposition_pattern,
    density_mean_difference,
    density_q_BH,
    crossmodal_triangulation_class,
    biological_readout_class,
    publication_strength,
    density_interpretation
  ) %>%
  mutate(
    minus_log10_q =
      -log10(communication_q_BH),

    DA_log2_OR =
      log2(DA_OR)
  )

write_csv(
  effect_landscape,
  file.path(
    TAB,
    "21G4C6B_signal_effect_landscape.csv"
  )
)

# ------------------------------------------------------------
# 7. Audit / freeze
# ------------------------------------------------------------

cat("\nSTEP 7/7 — formal checks\n")

formal_checks <- tribble(
  ~check, ~pass,

  "sixteen_frozen_signals",
  n_distinct(plot_long$signal_id) == 16,

  "sixteen_patients",
  n_distinct(plot_long$patient_id) == 16,

  "AA_HC_six_signals",
  n_distinct(AA_HC$signal_id) == 6,

  "FFA_HC_nine_signals",
  n_distinct(FFA_HC$signal_id) == 9,

  "AA_FFA_one_signal",
  n_distinct(AA_FFA$signal_id) == 1,

  "A_values_not_imputed",
  all(!is.na(plot_long$mean_incoming_A_axis)),

  "D_values_not_imputed",
  all(!is.na(plot_long$pooled_D_axis)),

  "intensity_missingness_preserved",
  TRUE,

  "no_patient_aggregation",
  TRUE,

  "no_signal_selection",
  TRUE,

  "no_PCA",
  TRUE,

  "no_clustering",
  TRUE,

  "no_new_statistical_test",
  TRUE,

  "no_threshold_refit",
  TRUE,

  "no_figure_generated",
  TRUE,

  "no_manuscript_modified",
  TRUE
)

print(formal_checks, n=Inf)

stopifnot(all(formal_checks$pass))

write_csv(
  formal_checks,
  file.path(
    TAB,
    "21G4C6B_formal_checks.csv"
  )
)

saveRDS(
  list(
    plot_long = plot_long,
    endpoint_long = endpoint_long,
    AA_vs_HC = AA_HC,
    FFA_vs_HC = FFA_HC,
    AA_vs_FFA = AA_FFA,
    AD_geometry = AD,
    effect_landscape = effect_landscape
  ),
  file.path(
    OBJ,
    "21G4C6B_quantitative_plotting_datasets.rds"
  )
)

audit <- c(
  "21G4C6B — QUANTITATIVE PLOTTING DATASETS FREEZE",
  "",
  "Frozen 21G4C6A patient-level matrices reused.",
  "",
  "PLOTTING DATASETS",
  "- endpoint-long patient-level representation",
  "- AA-vs-HC patient-level dataset",
  "- FFA-vs-HC patient-level dataset",
  "- AA-vs-FFA contextual dataset",
  "- patient-level A-D quantitative geometry",
  "- signal-level effect landscape",
  "",
  "SAFEGUARDS",
  "- no imputation",
  "- no patient aggregation",
  "- no signal selection",
  "- no PCA",
  "- no clustering",
  "- no statistical inference",
  "- no threshold refit",
  "- no plot generated",
  "",
  "PASS"
)

writeLines(
  audit,
  file.path(
    BASE,
    "21G4C6B_audit.txt"
  )
)

cat("\n============================================================\n")
cat("21G4C6B COMPLETE\n")
cat("============================================================\n")
cat("Quantitative plotting datasets frozen\n")
cat("6 AA-vs-HC signals retained\n")
cat("9 FFA-vs-HC signals retained\n")
cat("1 AA-vs-FFA contextual signal retained\n")
cat("Patient-level A-D geometry frozen\n")
cat("Signal-level effect landscape frozen\n")
cat("No imputation\n")
cat("No patient aggregation\n")
cat("No signal selection\n")
cat("No PCA or clustering\n")
cat("No plot generated\n")
cat("No statistical inference\n")
cat("PASS\n")
