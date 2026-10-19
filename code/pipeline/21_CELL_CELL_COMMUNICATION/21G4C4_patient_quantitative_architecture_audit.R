suppressPackageStartupMessages({
  library(tidyverse)
})

cat("============================================================\n")
cat("21G4C4 — PATIENT-LEVEL QUANTITATIVE ARCHITECTURE AUDIT\n")
cat("============================================================\n")

BASE <- "PIPELINE_V2/results/disease_biology/21_CELL_CELL_COMMUNICATION"
G21  <- file.path(BASE, "21G_MULTICOMPARTMENT_COMMUNICATION")
OUT  <- file.path(G21, "21G4_CROSSMODAL_INTEGRATION")
TAB  <- file.path(OUT, "tables")

dir.create(TAB, recursive=TRUE, showWarnings=FALSE)

# ------------------------------------------------------------
# 1. Locate canonical patient-level candidates
# ------------------------------------------------------------

cat("\nSTEP 1/8 — locating canonical patient-level sources\n")

files <- list.files(
  BASE,
  recursive=TRUE,
  full.names=TRUE,
  pattern="\\.(csv|rds)$"
)

targets <- tibble(
  file=files,
  basename=basename(files)
) %>%
  filter(
    basename %in% c(
      "21F2B2_patient_receiver_axis_metrics.csv",
      "21F2B2_patient_receiver_axis_metrics.rds",
      "21F2B2_patient_sender_receiver_axis_metrics.csv",
      "21F2B2_patient_sender_receiver_axis_metrics.rds",
      "21G2B_patient_quantitative_universe.csv",
      "21G2B_patient_quantitative_universe.rds"
    )
  ) %>%
  distinct(file, .keep_all=TRUE)

cat("\n===== TARGET SOURCES =====\n")
print(targets, n=Inf, width=Inf)

stopifnot(nrow(targets) > 0)

# ------------------------------------------------------------
# helper
# ------------------------------------------------------------

read_any <- function(f) {
  if (str_ends(f, ".csv")) {
    read_csv(f, show_col_types=FALSE)
  } else {
    readRDS(f)
  }
}

# ------------------------------------------------------------
# 2. Inspect exact schemas
# ------------------------------------------------------------

cat("\nSTEP 2/8 — exact schema inspection\n")

schema <- map_dfr(
  targets$file,
  function(f) {
    x <- read_any(f)

    if (!is.data.frame(x)) {
      return(
        tibble(
          file=f,
          n_rows=NA_integer_,
          n_cols=NA_integer_,
          columns=NA_character_
        )
      )
    }

    tibble(
      file=f,
      n_rows=nrow(x),
      n_cols=ncol(x),
      columns=paste(names(x), collapse=";")
    )
  }
)

print(schema, n=Inf, width=Inf)

write_csv(
  schema,
  file.path(TAB, "21G4C4_source_schema.csv")
)

# ------------------------------------------------------------
# 3. Select most informative receiver-axis table
# ------------------------------------------------------------

cat("\nSTEP 3/8 — selecting receiver-axis patient table\n")

receiver_candidates <- targets %>%
  filter(
    str_detect(
      basename,
      "patient_receiver_axis_metrics|patient_quantitative_universe"
    )
  )

candidate_info <- map_dfr(
  receiver_candidates$file,
  function(f) {
    x <- read_any(f)

    if (!is.data.frame(x)) return(tibble())

    nm <- names(x)
    low <- str_to_lower(nm)

    tibble(
      file=f,
      n_rows=nrow(x),
      has_patient=any(str_detect(low, "^patient_id$|patient|subject")),
      has_disease=any(str_detect(low, "^disease$|condition|group")),
      has_receiver=any(str_detect(low, "^receiver$")),
      has_axis=any(str_detect(low, "^axis$")),
      n_numeric=sum(map_lgl(x, is.numeric)),
      numeric_columns=paste(
        names(x)[map_lgl(x, is.numeric)],
        collapse=";"
      )
    )
  }
)

cat("\n===== RECEIVER-AXIS CANDIDATES =====\n")
print(candidate_info, n=Inf, width=Inf)

eligible <- candidate_info %>%
  filter(
    has_patient,
    has_disease,
    has_receiver,
    has_axis
  ) %>%
  arrange(desc(n_numeric), desc(n_rows))

stopifnot(nrow(eligible) > 0)

SOURCE <- eligible$file[[1]]

cat("\nSelected source:\n", SOURCE, "\n")

x <- read_any(SOURCE)

cat("Rows:", nrow(x), "\n")
cat("Columns:", ncol(x), "\n")

# ------------------------------------------------------------
# 4. Exact column/value audit
# ------------------------------------------------------------

cat("\nSTEP 4/8 — patient / disease / receiver / axis audit\n")

cat("\n===== COLUMN NAMES =====\n")
print(names(x))

cat("\n===== PATIENTS BY DISEASE =====\n")

patient_summary <- x %>%
  distinct(patient_id, disease) %>%
  count(disease, name="n_patients") %>%
  arrange(disease)

print(patient_summary, n=Inf)

cat("\n===== PATIENT IDENTITIES =====\n")

patients <- x %>%
  distinct(patient_id, disease) %>%
  arrange(disease, patient_id)

print(patients, n=Inf)

cat("\n===== RECEIVERS =====\n")

receivers <- x %>%
  count(receiver, name="n_rows") %>%
  arrange(desc(n_rows), receiver)

print(receivers, n=Inf)

cat("\n===== AXES =====\n")

axes <- x %>%
  count(axis, name="n_rows") %>%
  arrange(desc(n_rows), axis)

print(axes, n=Inf)

# ------------------------------------------------------------
# 5. Quantitative endpoint discovery
# ------------------------------------------------------------

cat("\nSTEP 5/8 — quantitative endpoint discovery\n")

numeric_cols <- names(x)[map_lgl(x, is.numeric)]

endpoint_summary <- map_dfr(
  numeric_cols,
  function(v) {
    z <- x[[v]]

    tibble(
      endpoint=v,
      n_total=length(z),
      n_nonmissing=sum(!is.na(z)),
      n_missing=sum(is.na(z)),
      pct_missing=100 * mean(is.na(z)),
      n_unique=n_distinct(z, na.rm=TRUE),
      min=ifelse(
        all(is.na(z)),
        NA_real_,
        min(z, na.rm=TRUE)
      ),
      median=ifelse(
        all(is.na(z)),
        NA_real_,
        median(z, na.rm=TRUE)
      ),
      max=ifelse(
        all(is.na(z)),
        NA_real_,
        max(z, na.rm=TRUE)
      )
    )
  }
)

cat("\n===== NUMERIC ENDPOINTS =====\n")
print(endpoint_summary, n=Inf, width=Inf)

write_csv(
  endpoint_summary,
  file.path(TAB, "21G4C4_numeric_endpoint_summary.csv")
)

# ------------------------------------------------------------
# 6. Patient × receiver × axis coverage
# ------------------------------------------------------------

cat("\nSTEP 6/8 — patient × receiver × axis coverage\n")

key_check <- x %>%
  count(
    patient_id,
    disease,
    receiver,
    axis,
    name="n"
  )

cat(
  "Unique patient × receiver × axis keys:",
  nrow(key_check),
  "\n"
)

cat(
  "Duplicated patient × receiver × axis keys:",
  sum(key_check$n > 1),
  "\n"
)

coverage <- x %>%
  distinct(
    patient_id,
    disease,
    receiver,
    axis
  ) %>%
  count(
    disease,
    receiver,
    axis,
    name="n_patients"
  ) %>%
  arrange(
    receiver,
    axis,
    disease
  )

cat("\n===== COVERAGE DISTRIBUTION =====\n")

coverage_summary <- coverage %>%
  group_by(disease) %>%
  summarise(
    n_receiver_axis=n(),
    min_patients=min(n_patients),
    median_patients=median(n_patients),
    max_patients=max(n_patients),
    .groups="drop"
  )

print(coverage_summary, n=Inf, width=Inf)

write_csv(
  coverage,
  file.path(TAB, "21G4C4_patient_receiver_axis_coverage.csv")
)

# ------------------------------------------------------------
# 7. Matrix feasibility
# ------------------------------------------------------------

cat("\nSTEP 7/8 — matrix completeness and multivariate feasibility\n")

patient_n <- n_distinct(x$patient_id)
feature_n <- x %>%
  distinct(receiver, axis) %>%
  nrow()

full_possible <- patient_n * feature_n

observed_keys <- x %>%
  distinct(patient_id, receiver, axis) %>%
  nrow()

matrix_completeness <- observed_keys / full_possible

matrix_summary <- tibble(
  n_patients=patient_n,
  n_receiver_axis_features=feature_n,
  possible_patient_feature_cells=full_possible,
  observed_patient_feature_cells=observed_keys,
  completeness=matrix_completeness
)

cat("\n===== MATRIX FEASIBILITY =====\n")
print(matrix_summary, width=Inf)

write_csv(
  matrix_summary,
  file.path(TAB, "21G4C4_matrix_feasibility.csv")
)

# disease-balanced features
disease_counts <- x %>%
  distinct(patient_id, disease) %>%
  count(disease, name="N")

balanced <- coverage %>%
  left_join(
    disease_counts,
    by="disease"
  ) %>%
  mutate(
    coverage_fraction =
      n_patients / N
  ) %>%
  group_by(receiver, axis) %>%
  summarise(
    min_disease_coverage =
      min(coverage_fraction),
    n_diseases =
      n_distinct(disease),
    .groups="drop"
  ) %>%
  arrange(
    desc(min_disease_coverage),
    receiver,
    axis
  )

cat("\n===== RECEIVER × AXIS BALANCED COVERAGE =====\n")
print(balanced, n=Inf, width=Inf)

write_csv(
  balanced,
  file.path(TAB, "21G4C4_balanced_feature_coverage.csv")
)

# ------------------------------------------------------------
# 8. Freeze audit
# ------------------------------------------------------------

cat("\nSTEP 8/8 — formal audit\n")

formal_checks <- tribble(
  ~check, ~pass,

  "patient_source_recovered",
  nrow(x) > 0,

  "patient_id_present",
  "patient_id" %in% names(x),

  "disease_present",
  "disease" %in% names(x),

  "receiver_present",
  "receiver" %in% names(x),

  "axis_present",
  "axis" %in% names(x),

  "numeric_endpoints_present",
  length(numeric_cols) > 0,

  "no_statistical_test",
  TRUE,

  "no_threshold_refit",
  TRUE,

  "no_imputation",
  TRUE,

  "no_PCA",
  TRUE,

  "no_clustering",
  TRUE,

  "no_figure_generated",
  TRUE,

  "no_manuscript_modified",
  TRUE,

  "no_analytical_refit",
  TRUE
)

cat("\n===== FORMAL CHECKS =====\n")
print(formal_checks, n=Inf)

stopifnot(all(formal_checks$pass))

write_csv(
  formal_checks,
  file.path(TAB, "21G4C4_formal_checks.csv")
)

audit <- c(
  "21G4C4 — PATIENT-LEVEL QUANTITATIVE ARCHITECTURE AUDIT",
  "",
  paste("Source:", SOURCE),
  paste("Rows:", nrow(x)),
  paste("Patients:", patient_n),
  paste("Receivers:", n_distinct(x$receiver)),
  paste("Axes:", n_distinct(x$axis)),
  paste("Receiver-axis features:", feature_n),
  paste(
    "Patient-feature completeness:",
    signif(matrix_completeness, 4)
  ),
  "",
  "No statistical inference.",
  "No imputation.",
  "No PCA.",
  "No clustering.",
  "No figure generated.",
  "No analytical refit.",
  "",
  "PASS"
)

writeLines(
  audit,
  file.path(OUT, "21G4C4_audit.txt")
)

cat("\n============================================================\n")
cat("21G4C4 COMPLETE\n")
cat("============================================================\n")
cat("Canonical patient-level receiver-axis source recovered\n")
cat("Patient/disease structure audited\n")
cat("Receiver × axis coverage audited\n")
cat("Quantitative endpoints inventoried\n")
cat("Matrix completeness quantified\n")
cat("No PCA performed\n")
cat("No clustering performed\n")
cat("No plot generated\n")
cat("No statistical inference\n")
cat("No analytical refit\n")
cat("PASS\n")
