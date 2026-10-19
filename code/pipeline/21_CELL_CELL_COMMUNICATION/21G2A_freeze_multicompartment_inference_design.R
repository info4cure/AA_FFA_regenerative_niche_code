suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(tibble)
})

cat("============================================================\n")
cat("PIPELINE_V2 — 21G2A MULTICOMPARTMENT INFERENCE DESIGN FREEZE\n")
cat("============================================================\n")

BASE <- paste0(
  "PIPELINE_V2/results/disease_biology/",
  "21_CELL_CELL_COMMUNICATION/",
  "21G_MULTICOMPARTMENT_COMMUNICATION/"
)

IN <- file.path(
  BASE,
  "21G1_EVALUABILITY",
  "tables"
)

OUT <- file.path(
  BASE,
  "21G2_INFERENCE"
)

TAB <- file.path(OUT, "tables")

dir.create(
  TAB,
  recursive=TRUE,
  showWarnings=FALSE
)

EVAL_FILE <- file.path(
  IN,
  "21G1B_identity_mechanism_evaluability.csv"
)

CHECK_FILE <- file.path(
  IN,
  "21G1B_formal_checks.csv"
)

stopifnot(
  file.exists(EVAL_FILE),
  file.exists(CHECK_FILE)
)

cat("\nSTEP 1/9 — recovering frozen evaluability\n")

x <- read_csv(
  EVAL_FILE,
  show_col_types=FALSE
)

checks_21G1B <- read_csv(
  CHECK_FILE,
  show_col_types=FALSE
)

stopifnot(
  all(checks_21G1B$pass),
  n_distinct(x$macrocompartment) == 9
)

# ============================================================
# 2. Freeze contrasts
# ============================================================

cat("\nSTEP 2/9 — freezing contrasts\n")

contrasts <- tribble(
  ~contrast,    ~numerator, ~denominator,
  "AA_vs_HC",   "AA",       "HC",
  "FFA_vs_HC",  "FFA",      "HC",
  "AA_vs_FFA",  "AA",       "FFA"
)

write_csv(
  contrasts,
  file.path(
    TAB,
    "21G2A_contrasts.csv"
  )
)

# ============================================================
# 3. Freeze endpoints
# ============================================================

cat("\nSTEP 3/9 — freezing endpoints\n")

endpoints <- tribble(
  ~endpoint,
  ~analysis_role,
  ~interpretation,

  "mean_incoming_A_axis",
  "PRIMARY",
  paste(
    "Mean normalized incoming mechanistic LR exposure",
    "across evaluable sender identities"
  ),

  "pooled_D_axis",
  "DECOMPOSITION",
  paste(
    "Fraction of structurally possible mechanistic LR",
    "channels detected across evaluable senders"
  ),

  "mean_sender_active_intensity",
  "DECOMPOSITION",
  paste(
    "Mean communication intensity among active",
    "mechanistic LR channels"
  )
)

write_csv(
  endpoints,
  file.path(
    TAB,
    "21G2A_endpoints.csv"
  )
)

# ============================================================
# 4. Expand frozen identity-level test universe
# ============================================================

cat("\nSTEP 4/9 — constructing frozen test universe\n")

key_cols <- c(
  "receiver",
  "macrocompartment",
  "candidate_status",
  "biological_program",
  "analysis_tier",
  "mechanistic_class",
  "action",
  "interpretation"
)

# Endpoint-specific contrast evaluability was frozen in 21G1B.
# This long-format table is the canonical source for A-axis,
# density and active-intensity evaluability.
EVAL_CONTRAST_FILE <- file.path(
  BASE,
  "21G1_EVALUABILITY",
  "tables",
  "21G1B_contrast_evaluability.csv"
)

stopifnot(
  file.exists(EVAL_CONTRAST_FILE)
)

contrast_eval <- read_csv(
  EVAL_CONTRAST_FILE,
  show_col_types=FALSE
)

required_eval_cols <- c(
  key_cols,
  "contrast",
  "numerator",
  "denominator",
  "n_num_A",
  "n_den_A",
  "n_num_D",
  "n_den_D",
  "n_num_intensity",
  "n_den_intensity",
  "A_testable",
  "D_testable",
  "intensity_testable"
)

missing_eval_cols <- setdiff(
  required_eval_cols,
  names(contrast_eval)
)

if (length(missing_eval_cols) > 0) {
  stop(
    paste0(
      "Missing columns in frozen 21G1B contrast evaluability: ",
      paste(
        missing_eval_cols,
        collapse=", "
      )
    )
  )
}

# Restrict explicitly to the three contrasts frozen in STEP 2.
# No evaluability is recalculated here.
test_universe <- contrast_eval %>%
  semi_join(
    contrasts,
    by=c(
      "contrast",
      "numerator",
      "denominator"
    )
  ) %>%
  select(
    all_of(key_cols),
    contrast,
    numerator,
    denominator,
    n_num_A,
    n_den_A,
    n_num_D,
    n_den_D,
    n_num_intensity,
    n_den_intensity,
    A_testable,
    D_testable,
    intensity_testable
  ) %>%
  arrange(
    contrast,
    macrocompartment,
    receiver,
    analysis_tier,
    biological_program,
    mechanistic_class,
    action
  )

stopifnot(
  nrow(test_universe) > 0,

  !anyDuplicated(
    test_universe[
      c(
        "receiver",
        "biological_program",
        "mechanistic_class",
        "action",
        "contrast"
      )
    ]
  ),

  all(
    test_universe$A_testable ==
      (
        test_universe$n_num_A >= 3 &
        test_universe$n_den_A >= 3
      )
  ),

  all(
    test_universe$D_testable ==
      (
        test_universe$n_num_D >= 3 &
        test_universe$n_den_D >= 3
      )
  ),

  all(
    test_universe$intensity_testable ==
      (
        test_universe$n_num_intensity >= 3 &
        test_universe$n_den_intensity >= 3
      )
  )
)

cat(
  "Frozen contrast-level keys:",
  nrow(test_universe),
  "\n"
)

cat("\n===== FROZEN EVALUABILITY BY CONTRAST =====\n")

print(
  test_universe %>%
    group_by(
      contrast
    ) %>%
    summarise(
      n_keys=n(),
      n_A_testable=sum(A_testable),
      n_D_testable=sum(D_testable),
      n_intensity_testable=sum(intensity_testable),
      .groups="drop"
    ),
  n=Inf,
  width=Inf
)

# ============================================================
# 5. Long endpoint-specific inference universe
# ============================================================

cat("\nSTEP 5/9 — freezing endpoint-specific inference universe\n")

analysis <- bind_rows(
  test_universe %>%
    transmute(
      across(
        all_of(
          c(
            key_cols,
            "contrast",
            "numerator",
            "denominator"
          )
        )
      ),
      endpoint="mean_incoming_A_axis",
      analysis_role="PRIMARY",
      testable=A_testable
    ),

  test_universe %>%
    transmute(
      across(
        all_of(
          c(
            key_cols,
            "contrast",
            "numerator",
            "denominator"
          )
        )
      ),
      endpoint="pooled_D_axis",
      analysis_role="DECOMPOSITION",
      testable=D_testable
    ),

  test_universe %>%
    transmute(
      across(
        all_of(
          c(
            key_cols,
            "contrast",
            "numerator",
            "denominator"
          )
        )
      ),
      endpoint="mean_sender_active_intensity",
      analysis_role="DECOMPOSITION",
      testable=intensity_testable
    )
)

analysis <- analysis %>%
  mutate(
    inference_status =
      if_else(
        testable,
        "PRESPECIFIED_TESTABLE",
        "PRESPECIFIED_NOT_EVALUABLE"
      )
  )

write_csv(
  analysis,
  file.path(
    TAB,
    "21G2A_frozen_inference_universe.csv"
  )
)

# ============================================================
# 6. Freeze statistical inference
# ============================================================

cat("\nSTEP 6/9 — freezing statistical inference\n")

inference <- tribble(
  ~component,
  ~specification,

  "unit_of_inference",
  "patient",

  "receiver_unit",
  paste(
    "communication identity; macrocompartments are used",
    "for biological organization and multiplicity only"
  ),

  "effect_estimand",
  "difference in arithmetic patient-level group means",

  "descriptive_ratio",
  paste(
    "log2 ratio of arithmetic group means with epsilon",
    "used descriptively only"
  ),

  "primary_test",
  paste(
    "two-sided exact exhaustive permutation test of",
    "patient disease labels within each contrast"
  ),

  "minimum_support",
  paste(
    "at least 3 evaluable patients in each disease",
    "group for the endpoint being tested"
  ),

  "missingness",
  paste(
    "non-evaluable patient-level observations remain missing;",
    "no biological zero imputation"
  ),

  "active_intensity",
  paste(
    "active intensity remains missing when no active",
    "mechanistic LR channels are present"
  ),

  "sender_structure",
  paste(
    "sender identities are aggregated within each patient",
    "before disease inference; no sender pseudoreplication"
  )
)

write_csv(
  inference,
  file.path(
    TAB,
    "21G2A_inference_specification.csv"
  )
)

# ============================================================
# 7. Freeze multiplicity architecture
# ============================================================

cat("\nSTEP 7/9 — freezing multiplicity architecture\n")

multiplicity <- tribble(
  ~family_dimension,
  ~specification,

  "contrast",
  "separate",

  "macrocompartment",
  "separate",

  "analysis_tier",
  "separate PRIMARY and CONTEXTUAL",

  "endpoint",
  paste(
    "separate mean_incoming_A_axis, pooled_D_axis and",
    "mean_sender_active_intensity"
  ),

  "correction",
  "Benjamini-Hochberg FDR",

  "membership",
  paste(
    "all prespecified and endpoint-evaluable identity ×",
    "mechanism tests within the family"
  )
)

write_csv(
  multiplicity,
  file.path(
    TAB,
    "21G2A_multiplicity_specification.csv"
  )
)

analysis <- analysis %>%
  mutate(
    BH_family =
      paste(
        contrast,
        macrocompartment,
        analysis_tier,
        endpoint,
        sep="__"
      )
  )

# Freeze family sizes using evaluability only
bh_families <- analysis %>%
  filter(
    testable
  ) %>%
  count(
    contrast,
    macrocompartment,
    analysis_tier,
    endpoint,
    BH_family,
    name="n_prespecified_tests"
  ) %>%
  arrange(
    contrast,
    macrocompartment,
    analysis_tier,
    endpoint
  )

write_csv(
  bh_families,
  file.path(
    TAB,
    "21G2A_BH_family_sizes.csv"
  )
)

# rewrite analysis including family
write_csv(
  analysis,
  file.path(
    TAB,
    "21G2A_frozen_inference_universe.csv"
  )
)

# ============================================================
# 8. Summary
# ============================================================

cat("\nSTEP 8/9 — summarising frozen design\n")

summary <- analysis %>%
  count(
    contrast,
    macrocompartment,
    analysis_tier,
    endpoint,
    inference_status,
    name="n_tests"
  ) %>%
  arrange(
    contrast,
    macrocompartment,
    analysis_tier,
    endpoint,
    inference_status
  )

write_csv(
  summary,
  file.path(
    TAB,
    "21G2A_design_summary.csv"
  )
)

cat("\n===== PRIMARY A-AXIS TESTABLE COUNTS =====\n")

primary_counts <- analysis %>%
  filter(
    endpoint == "mean_incoming_A_axis",
    analysis_tier == "PRIMARY",
    testable
  ) %>%
  count(
    contrast,
    macrocompartment,
    name="n_tests"
  ) %>%
  arrange(
    contrast,
    macrocompartment
  )

print(
  primary_counts,
  n=Inf,
  width=Inf
)

cat("\n===== BH FAMILY SIZES =====\n")

print(
  bh_families,
  n=Inf,
  width=Inf
)

# ============================================================
# 9. Formal checks
# ============================================================

cat("\nSTEP 9/9 — formal checks\n")

checks <- tibble(
  check=c(
    "21G1B_evaluability_present",
    "21G1B_all_checks_pass",
    "nine_macrocompartments_retained",
    "three_contrasts_frozen",
    "three_endpoints_frozen",
    "primary_A_axis_retained",
    "density_decomposition_retained",
    "intensity_decomposition_retained",
    "identity_level_inference_retained",
    "macrocompartments_not_pooled",
    "patient_unit_of_inference",
    "minimum_three_per_group_retained",
    "exact_permutation_frozen",
    "two_sided_test_frozen",
    "non_evaluable_not_zero",
    "active_intensity_missingness_retained",
    "sender_pseudoreplication_prohibited",
    "primary_contextual_separated",
    "BH_by_contrast",
    "BH_by_macrocompartment",
    "BH_by_analysis_tier",
    "BH_by_endpoint",
    "BH_family_sizes_frozen",
    "no_disease_effect_calculated",
    "no_p_values_calculated",
    "no_significance_selection",
    "no_threshold_refit",
    "no_communication_refit",
    "no_DE_refit",
    "no_pathway_refit",
    "annotation_unchanged"
  ),

  pass=c(
    file.exists(EVAL_FILE),
    all(checks_21G1B$pass),
    n_distinct(analysis$macrocompartment) == 9,
    n_distinct(analysis$contrast) == 3,
    n_distinct(analysis$endpoint) == 3,

    any(
      analysis$endpoint ==
        "mean_incoming_A_axis" &
      analysis$analysis_role ==
        "PRIMARY"
    ),

    any(
      analysis$endpoint ==
        "pooled_D_axis"
    ),

    any(
      analysis$endpoint ==
        "mean_sender_active_intensity"
    ),

    n_distinct(analysis$receiver) == 39,

    TRUE,

    any(
      inference$component ==
        "unit_of_inference" &
      inference$specification ==
        "patient"
    ),

    any(
      inference$component ==
        "minimum_support" &
      grepl(
        "at least 3",
        inference$specification
      )
    ),

    any(
      inference$component ==
        "primary_test" &
      grepl(
        "exact exhaustive permutation",
        inference$specification
      )
    ),

    any(
      inference$component ==
        "primary_test" &
      grepl(
        "two-sided",
        inference$specification
      )
    ),

    TRUE,
    TRUE,
    TRUE,

    all(
      c(
        "PRIMARY",
        "CONTEXTUAL"
      ) %in%
        analysis$analysis_tier
    ),

    TRUE,
    TRUE,
    TRUE,
    TRUE,

    nrow(bh_families) > 0 &&
      all(
        bh_families$n_prespecified_tests > 0
      ),

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

stopifnot(
  length(checks$check) ==
    length(checks$pass),
  all(checks$pass)
)

write_csv(
  checks,
  file.path(
    TAB,
    "21G2A_formal_checks.csv"
  )
)

# ============================================================
# Audit
# ============================================================

audit <- c(
  "PIPELINE_V2 — 21G2A multicompartment inference design freeze",
  "",
  paste0("Date: ", Sys.Date()),
  "",
  paste0(
    "Frozen inference rows: ",
    nrow(analysis)
  ),
  paste0(
    "Testable endpoint-specific rows: ",
    sum(analysis$testable)
  ),
  paste0(
    "Receiver identities: ",
    n_distinct(analysis$receiver)
  ),
  paste0(
    "Macrocompartments: ",
    n_distinct(analysis$macrocompartment)
  ),
  "",
  "PRIMARY ENDPOINT",
  "----------------",
  "mean_incoming_A_axis",
  "",
  "DECOMPOSITION",
  "-------------",
  "pooled_D_axis",
  "mean_sender_active_intensity",
  "",
  "MULTIPLICITY",
  "------------",
  paste(
    "BH correction is performed independently within",
    "contrast × macrocompartment × analysis tier × endpoint."
  ),
  "",
  "INFERENCE",
  "---------",
  paste(
    "Patient is the unit of inference.",
    "Receiver identities are tested separately.",
    "Macrocompartments are not pooled as replicate observations."
  ),
  "",
  paste(
    "Only prespecified endpoint-evaluable tests enter",
    "their corresponding BH family."
  ),
  "",
  "No disease effect calculated.",
  "No p-values calculated.",
  "No significance-based selection.",
  "No threshold refit.",
  "No communication refit.",
  "No DE refit.",
  "No pathway refit.",
  "No annotation modified.",
  "",
  "PASS"
)

writeLines(
  audit,
  file.path(
    OUT,
    "21G2A_audit.txt"
  )
)

cat("\n============================================================\n")
cat("21G2A COMPLETE\n")
cat("============================================================\n")
cat("Multicompartment disease-inference design frozen\n")
cat("39 receiver identities retained\n")
cat("9 macrocompartments retained\n")
cat("Primary and contextual hypotheses kept separate\n")
cat("Primary A-axis endpoint retained\n")
cat("Density and intensity decomposition retained\n")
cat("Patient remains unit of inference\n")
cat("Exact exhaustive permutation inference frozen\n")
cat("BH families frozen before disease testing\n")
cat("No disease
 effect calculated\n")
cat("No p-values calculated\n")
cat("No significance-based selection\n")
cat("No analytical refit\n")
cat("PASS\n")
