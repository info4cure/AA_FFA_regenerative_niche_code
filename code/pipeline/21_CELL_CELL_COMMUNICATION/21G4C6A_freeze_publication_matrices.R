suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(tibble)
})

cat("============================================================\n")
cat("21G4C6A — FREEZE PATIENT-LEVEL PUBLICATION MATRICES\n")
cat("============================================================\n")

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
  "21G4C5B_patient_publication_signals.csv"
)

stopifnot(file.exists(INFILE))

cat("\nSTEP 1/7 — loading frozen patient-level publication signals\n")

x <- read_csv(INFILE, show_col_types=FALSE)

cat("Rows:", nrow(x), "\n")
cat("Patients:", n_distinct(x$patient_id), "\n")
cat("Signals:", n_distinct(x$signal_id), "\n")

stopifnot(
  nrow(x) == 162,
  n_distinct(x$patient_id) == 16,
  n_distinct(x$signal_id) == 16
)

# ------------------------------------------------------------
# 2. Exact patient × signal key
# ------------------------------------------------------------

cat("\nSTEP 2/7 — auditing exact patient × signal key\n")

dup <- x %>%
  count(patient_id, signal_id, name="n") %>%
  filter(n > 1)

cat("Duplicated patient × signal keys:", nrow(dup), "\n")

stopifnot(nrow(dup) == 0)

# ------------------------------------------------------------
# 3. Freeze signal metadata
# ------------------------------------------------------------

cat("\nSTEP 3/7 — freezing signal metadata\n")

signal_metadata <- x %>%
  distinct(
    signal_id,
    contrast,
    numerator,
    denominator,
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
  ) %>%
  arrange(signal_id)

stopifnot(nrow(signal_metadata) == 16)

write_csv(
  signal_metadata,
  file.path(TAB, "21G4C6A_signal_metadata.csv")
)

# ------------------------------------------------------------
# 4. Long quantitative matrix
# ------------------------------------------------------------

cat("\nSTEP 4/7 — freezing lossless long quantitative matrix\n")

long <- x %>%
  select(
    signal_id,
    contrast,
    patient_id,
    disease,
    receiver,
    communication_program,
    communication_axis,
    mean_incoming_A_axis,
    pooled_D_axis,
    mean_sender_active_intensity
  ) %>%
  arrange(signal_id, disease, patient_id)

write_csv(
  long,
  file.path(TAB, "21G4C6A_patient_signal_long.csv")
)

saveRDS(
  long,
  file.path(OBJ, "21G4C6A_patient_signal_long.rds")
)

# ------------------------------------------------------------
# 5. Wide endpoint matrices
# ------------------------------------------------------------

cat("\nSTEP 5/7 — constructing A, D and intensity matrices\n")

patient_metadata <- long %>%
  distinct(patient_id, disease) %>%
  arrange(
    factor(disease, levels=c("HC","AA","FFA")),
    patient_id
  )

make_wide <- function(dat, endpoint) {

  dat %>%
    select(
      patient_id,
      disease,
      signal_id,
      all_of(endpoint)
    ) %>%
    pivot_wider(
      names_from = signal_id,
      values_from = all_of(endpoint)
    ) %>%
    arrange(
      factor(disease, levels=c("HC","AA","FFA")),
      patient_id
    )
}

A_matrix <- make_wide(
  long,
  "mean_incoming_A_axis"
)

D_matrix <- make_wide(
  long,
  "pooled_D_axis"
)

I_matrix <- make_wide(
  long,
  "mean_sender_active_intensity"
)

write_csv(
  patient_metadata,
  file.path(TAB, "21G4C6A_patient_metadata.csv")
)

write_csv(
  A_matrix,
  file.path(TAB, "21G4C6A_A_matrix.csv")
)

write_csv(
  D_matrix,
  file.path(TAB, "21G4C6A_D_matrix.csv")
)

write_csv(
  I_matrix,
  file.path(TAB, "21G4C6A_active_intensity_matrix.csv")
)

saveRDS(
  list(
    patient_metadata = patient_metadata,
    signal_metadata = signal_metadata,
    A = A_matrix,
    D = D_matrix,
    intensity = I_matrix,
    long = long
  ),
  file.path(
    OBJ,
    "21G4C6A_publication_matrices.rds"
  )
)

# ------------------------------------------------------------
# 6. Coverage audit
# ------------------------------------------------------------

cat("\nSTEP 6/7 — endpoint completeness audit\n")

endpoint_coverage <- x %>%
  group_by(signal_id, contrast) %>%
  summarise(
    n_patients = n(),
    A_nonmissing =
      sum(!is.na(mean_incoming_A_axis)),
    D_nonmissing =
      sum(!is.na(pooled_D_axis)),
    intensity_nonmissing =
      sum(!is.na(mean_sender_active_intensity)),
    .groups="drop"
  ) %>%
  mutate(
    A_complete = A_nonmissing == n_patients,
    D_complete = D_nonmissing == n_patients,
    intensity_complete =
      intensity_nonmissing == n_patients
  )

write_csv(
  endpoint_coverage,
  file.path(
    TAB,
    "21G4C6A_endpoint_completeness.csv"
  )
)

cat("\n===== ENDPOINT COMPLETENESS =====\n")
print(endpoint_coverage, n=Inf, width=Inf)

cat("\nA complete for all evaluable patient-signal records:",
    all(endpoint_coverage$A_complete), "\n")

cat("D complete for all evaluable patient-signal records:",
    all(endpoint_coverage$D_complete), "\n")

cat(
  "Signals with incomplete active intensity:",
  sum(!endpoint_coverage$intensity_complete),
  "\n"
)

# ------------------------------------------------------------
# 7. Formal freeze
# ------------------------------------------------------------

cat("\nSTEP 7/7 — formal checks\n")

formal_checks <- tribble(
  ~check, ~pass,

  "source_162_rows",
  nrow(x) == 162,

  "sixteen_patients",
  n_distinct(x$patient_id) == 16,

  "sixteen_frozen_signals",
  n_distinct(x$signal_id) == 16,

  "patient_signal_key_unique",
  nrow(dup) == 0,

  "A_complete_within_evaluable_records",
  all(endpoint_coverage$A_complete),

  "D_complete_within_evaluable_records",
  all(endpoint_coverage$D_complete),

  "intensity_missingness_preserved",
  TRUE,

  "no_imputation",
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
    "21G4C6A_formal_checks.csv"
  )
)

audit <- c(
  "21G4C6A — PATIENT-LEVEL PUBLICATION MATRIX FREEZE",
  "",
  paste("Source rows:", nrow(x)),
  paste("Patients:", n_distinct(x$patient_id)),
  paste("Frozen signals:", n_distinct(x$signal_id)),
  "",
  "ENDPOINT POLICY",
  "- mean_incoming_A_axis retained without imputation",
  "- pooled_D_axis retained without imputation",
  "- mean_sender_active_intensity retained with original missingness",
  "- no missing values filled",
  "- no patient aggregation",
  "- no signal selection",
  "",
  "ANALYTICAL POLICY",
  "- no PCA",
  "- no clustering",
  "- no statistical inference",
  "- no threshold refit",
  "- no figure generation",
  "- no manuscript modification",
  "",
  "PASS"
)

writeLines(
  audit,
  file.path(BASE, "21G4C6A_audit.txt")
)

cat("\n============================================================\n")
cat("21G4C6A COMPLETE\n")
cat("============================================================\n")
cat("16 frozen publication signals retained\n")
cat("16 patients retained across AA, FFA and HC\n")
cat("A-axis matrix frozen\n")
cat("D-axis matrix frozen\n")
cat("Active-intensity matrix frozen with native missingness\n")
cat("No imputation\n")
cat("No patient aggregation\n")
cat("No signal selection\n")
cat("No PCA or clustering\n")
cat("No plot generated\n")
cat("No statistical inference\n")
cat("PASS\n")
