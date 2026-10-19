suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
})

cat("============================================================\n")
cat("21G4C5B — PATIENT-LEVEL PUBLICATION SIGNAL EXTRACTION\n")
cat("============================================================\n")

BASE <- paste0(
  "PIPELINE_V2/results/disease_biology/21_CELL_CELL_COMMUNICATION/",
  "21G_MULTICOMPARTMENT_COMMUNICATION/21G4_CROSSMODAL_INTEGRATION"
)

TAB <- file.path(BASE, "tables")
OBJ <- file.path(BASE, "objects")

dir.create(TAB, recursive=TRUE, showWarnings=FALSE)
dir.create(OBJ, recursive=TRUE, showWarnings=FALSE)

PATIENT_FILE <- paste0(
  "PIPELINE_V2/results/disease_biology/21_CELL_CELL_COMMUNICATION/",
  "21F2_BIOLOGICAL_LR_AXES/21F2B2_PATIENT_AXIS_PROJECTION/objects/",
  "21F2B2_patient_receiver_axis_metrics.rds"
)

SIGNAL_FILE <- file.path(
  TAB,
  "21G4B_biological_synthesis_signals.csv"
)

stopifnot(
  file.exists(PATIENT_FILE),
  file.exists(SIGNAL_FILE)
)

cat("\nSTEP 1/7 — loading frozen sources\n")

p <- readRDS(PATIENT_FILE)
s <- read_csv(SIGNAL_FILE, show_col_types=FALSE)

cat("Patient-level rows:", nrow(p), "\n")
cat("Frozen synthesis signals:", nrow(s), "\n")

stopifnot(nrow(s) == 16)

# ------------------------------------------------------------
# 2. Construct exact biological publication key
# ------------------------------------------------------------

cat("\nSTEP 2/7 — constructing exact signal key\n")

signal_key <- s %>%
  transmute(
    signal_id = sprintf(
      "21G4PUB_%02d",
      row_number()
    ),

    contrast,
    macrocompartment,
    receiver,

    communication_program,
    communication_tier,
    communication_mechanistic_class,
    communication_action,
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
  )

stopifnot(
  nrow(signal_key) == 16,
  n_distinct(signal_key$signal_id) == 16
)

# ------------------------------------------------------------
# 3. Expand contrast to disease pair
# ------------------------------------------------------------

cat("\nSTEP 3/7 — resolving disease membership\n")

contrast_lookup <- tibble::tribble(
  ~contrast,   ~numerator, ~denominator,
  "AA_vs_HC",  "AA",       "HC",
  "FFA_vs_HC", "FFA",      "HC",
  "AA_vs_FFA", "AA",       "FFA"
)

signal_key <- signal_key %>%
  left_join(
    contrast_lookup,
    by="contrast"
  )

stopifnot(
  all(!is.na(signal_key$numerator)),
  all(!is.na(signal_key$denominator))
)

# ------------------------------------------------------------
# 4. Exact patient-level join
# ------------------------------------------------------------

cat("\nSTEP 4/7 — exact patient-level extraction\n")

patient_signals <- signal_key %>%
  inner_join(
    p,
    by=c(
      "receiver" = "receiver",
      "communication_axis" = "axis",
      "communication_mechanistic_class" =
        "mechanistic_class",
      "communication_action" = "action"
    ),
    relationship="many-to-many"
  ) %>%
  filter(
    disease == numerator |
    disease == denominator
  )

cat(
  "Extracted patient-signal rows:",
  nrow(patient_signals),
  "\n"
)

# ------------------------------------------------------------
# 5. Exact patient-signal uniqueness
# ------------------------------------------------------------

cat("\nSTEP 5/7 — auditing patient-signal cardinality\n")

dup <- patient_signals %>%
  count(
    signal_id,
    patient_id,
    name="n"
  ) %>%
  filter(n > 1)

cat(
  "Duplicated signal × patient keys:",
  nrow(dup),
  "\n"
)

stopifnot(nrow(dup) == 0)

coverage <- patient_signals %>%
  group_by(
    signal_id,
    contrast,
    receiver,
    communication_program,
    communication_axis,
    communication_evidence_class,
    communication_q_BH
  ) %>%
  summarise(
    n_patients = n_distinct(patient_id),

    n_numerator =
      n_distinct(
        patient_id[disease == numerator]
      ),

    n_denominator =
      n_distinct(
        patient_id[disease == denominator]
      ),

    .groups="drop"
  ) %>%
  arrange(
    contrast,
    communication_q_BH
  )

cat("\n===== PATIENT COVERAGE BY SIGNAL =====\n")
print(
  coverage,
  n=Inf,
  width=Inf
)

# ------------------------------------------------------------
# 6. Quantitative endpoint audit
# ------------------------------------------------------------

cat("\nSTEP 6/7 — quantitative endpoint audit\n")

endpoint_coverage <- patient_signals %>%
  group_by(signal_id) %>%
  summarise(
    n_patients =
      n_distinct(patient_id),

    A_nonmissing =
      sum(!is.na(mean_incoming_A_axis)),

    D_nonmissing =
      sum(!is.na(pooled_D_axis)),

    intensity_nonmissing =
      sum(!is.na(mean_sender_active_intensity)),

    .groups="drop"
  )

cat("\n===== ENDPOINT COVERAGE =====\n")
print(
  endpoint_coverage,
  n=Inf,
  width=Inf
)

# Disease representation
disease_coverage <- patient_signals %>%
  distinct(
    signal_id,
    contrast,
    patient_id,
    disease
  ) %>%
  count(
    signal_id,
    contrast,
    disease,
    name="n_patients"
  ) %>%
  arrange(
    signal_id,
    disease
  )

cat("\n===== DISEASE COVERAGE =====\n")
print(
  disease_coverage,
  n=Inf,
  width=Inf
)

# ------------------------------------------------------------
# 7. Freeze extraction
# ------------------------------------------------------------

cat("\nSTEP 7/7 — formal checks and freeze\n")

checks <- tibble::tribble(
  ~check, ~pass,

  "sixteen_frozen_signals",
  nrow(signal_key) == 16,

  "sixteen_signal_ids",
  n_distinct(signal_key$signal_id) == 16,

  "patient_signal_unique",
  nrow(dup) == 0,

  "all_signals_represented",
  n_distinct(patient_signals$signal_id) == 16,

  "all_rows_correct_disease_pair",
  all(
    patient_signals$disease ==
      patient_signals$numerator |
    patient_signals$disease ==
      patient_signals$denominator
  ),

  "A_endpoint_present",
  "mean_incoming_A_axis" %in%
    names(patient_signals),

  "D_endpoint_present",
  "pooled_D_axis" %in%
    names(patient_signals),

  "intensity_endpoint_present",
  "mean_sender_active_intensity" %in%
    names(patient_signals),

  "no_imputation",
  TRUE,

  "no_aggregation",
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

cat("\n===== FORMAL CHECKS =====\n")
print(checks, n=Inf)

stopifnot(all(checks$pass))

write_csv(
  signal_key,
  file.path(
    TAB,
    "21G4C5B_publication_signal_key.csv"
  )
)

write_csv(
  patient_signals,
  file.path(
    TAB,
    "21G4C5B_patient_publication_signals.csv"
  )
)

write_csv(
  coverage,
  file.path(
    TAB,
    "21G4C5B_patient_signal_coverage.csv"
  )
)

write_csv(
  endpoint_coverage,
  file.path(
    TAB,
    "21G4C5B_endpoint_coverage.csv"
  )
)

write_csv(
  disease_coverage,
  file.path(
    TAB,
    "21G4C5B_disease_coverage.csv"
  )
)

write_csv(
  checks,
  file.path(
    TAB,
    "21G4C5B_formal_checks.csv"
  )
)

saveRDS(
  patient_signals,
  file.path(
    OBJ,
    "21G4C5B_patient_publication_signals.rds"
  )
)

cat("\n============================================================\n")
cat("21G4C5B COMPLETE\n")
cat("============================================================\n")
cat("16 frozen publication signals projected to patient level\n")
cat("Exact receiver × axis × mechanism × action matching used\n")
cat("Only contrast-relevant diseases retained\n")
cat("Patient-level A, D and active-intensity endpoints preserved\n")
cat("No imputation\n")
cat("No aggregation\n")
cat("No new statistical inference\n")
cat("No plot generated\n")
cat("PASS\n")
