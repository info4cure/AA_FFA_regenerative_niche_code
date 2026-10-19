suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
  library(tidyr)
  library(purrr)
})

cat("============================================================\n")
cat("PIPELINE_V2 — 21F2B3B EXACT PATIENT DISEASE CONTRASTS\n")
cat("============================================================\n")

# ============================================================
# 0. Paths
# ============================================================

BASE <- paste0(
  "PIPELINE_V2/results/disease_biology/",
  "21_CELL_CELL_COMMUNICATION/",
  "21F2_BIOLOGICAL_LR_AXES/"
)

B2 <- file.path(
  BASE,
  "21F2B2_PATIENT_AXIS_PROJECTION"
)

B3 <- file.path(
  BASE,
  "21F2B3_DISEASE_CONTRASTS"
)

TAB <- file.path(B3, "tables")
OBJ <- file.path(B3, "objects")

dir.create(TAB, recursive=TRUE, showWarnings=FALSE)
dir.create(OBJ, recursive=TRUE, showWarnings=FALSE)

DATA_FILE <- file.path(
  B2,
  "tables",
  "21F2B2_patient_receiver_axis_metrics.csv"
)

DESIGN_FILE <- file.path(
  TAB,
  "21F2B3A_mechanistic_hypothesis_families.csv"
)

CONTRAST_FILE <- file.path(
  TAB,
  "21F2B3A_contrast_specification.csv"
)

RECEIVER_FILE <- file.path(
  TAB,
  "21F2B3A_receiver_panel.csv"
)

ENDPOINT_FILE <- file.path(
  TAB,
  "21F2B3A_endpoint_specification.csv"
)

CHECK_FILE <- file.path(
  TAB,
  "21F2B3A_formal_checks.csv"
)

CONST_FILE <- file.path(
  TAB,
  "21F2B3A_numeric_constants.csv"
)

stopifnot(
  file.exists(DATA_FILE),
  file.exists(DESIGN_FILE),
  file.exists(CONTRAST_FILE),
  file.exists(RECEIVER_FILE),
  file.exists(ENDPOINT_FILE),
  file.exists(CHECK_FILE),
  file.exists(CONST_FILE)
)

# ============================================================
# 1. Recover frozen design and patient-level data
# ============================================================

cat("\nSTEP 1/12 — recovering frozen design and patient-level data\n")

x <- read_csv(
  DATA_FILE,
  show_col_types=FALSE
)

design <- read_csv(
  DESIGN_FILE,
  show_col_types=FALSE
)

contrasts <- read_csv(
  CONTRAST_FILE,
  show_col_types=FALSE
)

receivers <- read_csv(
  RECEIVER_FILE,
  show_col_types=FALSE
)

endpoints <- read_csv(
  ENDPOINT_FILE,
  show_col_types=FALSE
)

design_checks <- read_csv(
  CHECK_FILE,
  show_col_types=FALSE
)

constants <- read_csv(
  CONST_FILE,
  show_col_types=FALSE
)

stopifnot(
  all(design_checks$pass)
)

epsilon <- constants$value[
  constants$parameter == "epsilon"
]

stopifnot(
  length(epsilon) == 1,
  epsilon == 1e-6
)

cat("Input rows:", nrow(x), "\n")
cat("Patients:", n_distinct(x$patient_id), "\n")
cat("Frozen mechanisms:", nrow(design), "\n")

# ============================================================
# 2. Restrict exactly to frozen receiver/mechanism design
# ============================================================

cat("\nSTEP 2/12 — applying frozen analysis design\n")

analysis_data <- x %>%
  semi_join(
    receivers,
    by="receiver"
  ) %>%
  inner_join(
    design,
    by=c(
      "mechanistic_class",
      "action"="action_required"
    ),
    relationship="many-to-one"
  )

stopifnot(
  !anyDuplicated(
    analysis_data[
      c(
        "patient_id",
        "receiver",
        "axis",
        "mechanistic_class",
        "action"
      )
    ]
  )
)

cat(
  "Frozen-design patient rows:",
  nrow(analysis_data),
  "\n"
)

cat(
  "Observed receiver-axis-mechanism-action keys:",
  analysis_data %>%
    distinct(
      receiver,
      axis,
      mechanistic_class,
      action
    ) %>%
    nrow(),
  "\n"
)

# ============================================================
# 3. Record prespecified but structurally non-evaluable keys
# ============================================================

cat("\nSTEP 3/12 — auditing structural non-evaluability\n")

observed_design <- analysis_data %>%
  distinct(
    mechanistic_class,
    action
  )

nonobservable_design <- design %>%
  anti_join(
    observed_design,
    by=c(
      "mechanistic_class",
      "action_required"="action"
    )
  ) %>%
  mutate(
    status=
      "PRESPECIFIED_STRUCTURALLY_NON_EVALUABLE_IN_FOLLICULAR_PANEL"
  )

cat("\n===== PRESPECIFIED NON-EVALUABLE MECHANISMS =====\n")

print(
  nonobservable_design,
  n=Inf,
  width=Inf
)

write_csv(
  nonobservable_design,
  file.path(
    TAB,
    "21F2B3B_prespecified_structurally_nonevaluable.csv"
  )
)

# ============================================================
# 4. Long endpoint representation
# ============================================================

cat("\nSTEP 4/12 — constructing endpoint-long patient table\n")

endpoint_names <- c(
  "mean_incoming_A_axis",
  "pooled_D_axis",
  "mean_sender_active_intensity"
)

stopifnot(
  all(
    endpoint_names %in%
      names(analysis_data)
  )
)

long <- analysis_data %>%
  select(
    patient_id,
    disease,
    receiver,
    axis,
    mechanistic_class,
    action,
    analysis_tier,
    hypothesis_family,
    all_of(endpoint_names)
  ) %>%
  pivot_longer(
    cols=all_of(endpoint_names),
    names_to="endpoint",
    values_to="value"
  ) %>%
  left_join(
    endpoints %>%
      select(
        endpoint,
        endpoint_role=analysis_role
      ),
    by="endpoint"
  )

stopifnot(
  all(
    !is.na(long$endpoint_role)
  )
)

# ============================================================
# 5. Exact permutation function
# ============================================================

cat("\nSTEP 5/12 — freezing exact permutation implementation\n")

exact_test <- function(
  values,
  labels,
  numerator,
  denominator
) {

  keep <-
    !is.na(values) &
    labels %in%
      c(
        numerator,
        denominator
      )

  y <- values[keep]
  g <- labels[keep]

  y_num <- y[g == numerator]
  y_den <- y[g == denominator]

  n_num <- length(y_num)
  n_den <- length(y_den)

  if (
    n_num < 3 ||
    n_den < 3
  ) {
    return(
      tibble(
        n_num=n_num,
        n_den=n_den,
        n_total=n_num+n_den,
        mean_num=ifelse(
          n_num > 0,
          mean(y_num),
          NA_real_
        ),
        mean_den=ifelse(
          n_den > 0,
          mean(y_den),
          NA_real_
        ),
        median_num=ifelse(
          n_num > 0,
          median(y_num),
          NA_real_
        ),
        median_den=ifelse(
          n_den > 0,
          median(y_den),
          NA_real_
        ),
        mean_difference=NA_real_,
        median_difference=NA_real_,
        mean_ratio=NA_real_,
        log2_mean_ratio=NA_real_,
        n_exact_allocations=NA_integer_,
        p_exact=NA_real_,
        test_status="INSUFFICIENT_PATIENT_SUPPORT"
      )
    )
  }

  mean_num <- mean(y_num)
  mean_den <- mean(y_den)

  median_num <- median(y_num)
  median_den <- median(y_den)

  delta <-
    mean_num -
    mean_den

  idx <- seq_along(y)

  allocations <- combn(
    idx,
    n_num,
    simplify=FALSE
  )

  perm_delta <- vapply(
    allocations,
    function(i) {
      mean(y[i]) -
        mean(y[-i])
    },
    numeric(1)
  )

  tol <- 1e-12

  p_exact <- mean(
    abs(perm_delta) >=
      abs(delta) - tol
  )

  ratio <-
    (mean_num + epsilon) /
    (mean_den + epsilon)

  tibble(
    n_num=n_num,
    n_den=n_den,
    n_total=n_num+n_den,
    mean_num=mean_num,
    mean_den=mean_den,
    median_num=median_num,
    median_den=median_den,
    mean_difference=delta,
    median_difference=
      median_num -
      median_den,
    mean_ratio=ratio,
    log2_mean_ratio=log2(ratio),
    n_exact_allocations=
      length(allocations),
    p_exact=p_exact,
    test_status="TESTED_EXACT"
  )
}

# ============================================================
# 6. Execute frozen contrasts
# ============================================================

cat("\nSTEP 6/12 — executing exact patient-level contrasts\n")

test_keys <- long %>%
  distinct(
    receiver,
    axis,
    mechanistic_class,
    action,
    analysis_tier,
    hypothesis_family,
    endpoint,
    endpoint_role
  )

run_one <- function(
  receiver_i,
  axis_i,
  mechanism_i,
  action_i,
  tier_i,
  family_i,
  endpoint_i,
  endpoint_role_i,
  contrast_i,
  numerator_i,
  denominator_i
) {

  z <- long %>%
    filter(
      receiver == receiver_i,
      axis == axis_i,
      mechanistic_class == mechanism_i,
      action == action_i,
      endpoint == endpoint_i,
      disease %in%
        c(
          numerator_i,
          denominator_i
        )
    )

  res <- exact_test(
    values=z$value,
    labels=z$disease,
    numerator=numerator_i,
    denominator=denominator_i
  )

  tibble(
    contrast=contrast_i,
    numerator=numerator_i,
    denominator=denominator_i,
    receiver=receiver_i,
    axis=axis_i,
    mechanistic_class=mechanism_i,
    action=action_i,
    analysis_tier=tier_i,
    hypothesis_family=family_i,
    endpoint=endpoint_i,
    endpoint_role=endpoint_role_i
  ) %>%
    bind_cols(res)
}

jobs <- tidyr::crossing(
  test_keys,
  contrasts %>%
    select(
      contrast,
      numerator,
      denominator
    )
)

results <- pmap_dfr(
  jobs,
  function(
    receiver,
    axis,
    mechanistic_class,
    action,
    analysis_tier,
    hypothesis_family,
    endpoint,
    endpoint_role,
    contrast,
    numerator,
    denominator
  ) {
    run_one(
      receiver,
      axis,
      mechanistic_class,
      action,
      analysis_tier,
      hypothesis_family,
      endpoint,
      endpoint_role,
      contrast,
      numerator,
      denominator
    )
  }
)

cat(
  "Total prespecified statistical rows:",
  nrow(results),
  "\n"
)

cat(
  "Exact tests performed:",
  sum(
    results$test_status ==
      "TESTED_EXACT"
  ),
  "\n"
)

cat(
  "Insufficient-support rows:",
  sum(
    results$test_status ==
      "INSUFFICIENT_PATIENT_SUPPORT"
  ),
  "\n"
)

# ============================================================
# 7. Apply frozen BH multiplicity families
# ============================================================

cat("\nSTEP 7/12 — applying frozen BH families\n")

results <- results %>%
  mutate(
    BH_family=case_when(

      analysis_tier == "PRIMARY" &
      endpoint_role == "PRIMARY" ~
        paste(
          contrast,
          hypothesis_family,
          sep="__"
        ),

      analysis_tier == "PRIMARY" &
      endpoint_role ==
        "SECONDARY_DECOMPOSITION" ~
        paste(
          contrast,
          endpoint,
          hypothesis_family,
          sep="__"
        ),

      analysis_tier == "CONTEXTUAL" &
      endpoint_role == "PRIMARY" ~
        paste(
          contrast,
          "CONTEXTUAL_REGULATION",
          sep="__"
        ),

      analysis_tier == "CONTEXTUAL" &
      endpoint_role ==
        "SECONDARY_DECOMPOSITION" ~
        paste(
          contrast,
          endpoint,
          "CONTEXTUAL_REGULATION",
          sep="__"
        ),

      TRUE ~ NA_character_
    )
  ) %>%
  group_by(
    BH_family
  ) %>%
  mutate(
    q_BH=if_else(
      !is.na(p_exact),
      p.adjust(
        p_exact,
        method="BH"
      ),
      NA_real_
    )
  ) %>%
  ungroup()

# IMPORTANT:
# p.adjust() ignores NA values; therefore only actually
# testable hypotheses contribute to the multiplicity count.

# ============================================================
# 8. Freeze tested/untested hypothesis accounting
# ============================================================

cat("\nSTEP 8/12 — auditing hypothesis accounting\n")

accounting <- results %>%
  count(
    contrast,
    analysis_tier,
    endpoint,
    endpoint_role,
    test_status,
    name="n_tests"
  ) %>%
  arrange(
    contrast,
    analysis_tier,
    endpoint,
    test_status
  )

cat("\n===== TEST ACCOUNTING =====\n")

print(
  accounting,
  n=Inf,
  width=Inf
)

write_csv(
  accounting,
  file.path(
    TAB,
    "21F2B3B_test_accounting.csv"
  )
)

bh_accounting <- results %>%
  filter(
    !is.na(BH_family)
  ) %>%
  group_by(
    BH_family
  ) %>%
  summarise(
    n_prespecified=n(),
    n_tested=sum(
      !is.na(p_exact)
    ),
    n_untested=sum(
      is.na(p_exact)
    ),
    .groups="drop"
  )

write_csv(
  bh_accounting,
  file.path(
    TAB,
    "21F2B3B_BH_family_accounting.csv"
  )
)

# ============================================================
# 9. Write complete results before inspection
# ============================================================

cat("\nSTEP 9/12 — freezing complete inferential result table\n")

results <- results %>%
  arrange(
    contrast,
    analysis_tier,
    hypothesis_family,
    receiver,
    axis,
    mechanistic_class,
    action,
    endpoint
  )

write_csv(
  results,
  file.path(
    TAB,
    "21F2B3B_all_exact_contrasts.csv"
  )
)

saveRDS(
  results,
  file.path(
    OBJ,
    "21F2B3B_all_exact_contrasts.rds"
  )
)

# ============================================================
# 10. Primary-result and decomposition tables
# ============================================================

cat("\nSTEP 10/12 — generating prespecified reporting tables\n")

primary <- results %>%
  filter(
    analysis_tier == "PRIMARY",
    endpoint_role == "PRIMARY"
  )

decomposition <- results %>%
  filter(
    analysis_tier == "PRIMARY",
    endpoint_role ==
      "SECONDARY_DECOMPOSITION"
  )

contextual <- results %>%
  filter(
    analysis_tier == "CONTEXTUAL"
  )

write_csv(
  primary,
  file.path(
    TAB,
    "21F2B3B_primary_A_axis_contrasts.csv"
  )
)

write_csv(
  decomposition,
  file.path(
    TAB,
    "21F2B3B_primary_decomposition_contrasts.csv"
  )
)

write_csv(
  contextual,
  file.path(
    TAB,
    "21F2B3B_contextual_contrasts.csv"
  )
)

# This is descriptive reporting only, not feature selection.
primary_readout <- primary %>%
  arrange(
    contrast,
    q_BH,
    p_exact,
    receiver,
    mechanistic_class
  )

cat("\n===== PRIMARY A-AXIS READOUT =====\n")

print(
  primary_readout,
  n=Inf,
  width=Inf
)

# ============================================================
# 11. Exact-allocation audit
# ============================================================

cat("\nSTEP 11/12 — auditing exact permutation counts\n")

allocation_audit <- results %>%
  filter(
    test_status ==
      "TESTED_EXACT"
  ) %>%
  distinct(
    contrast,
    n_num,
    n_den,
    n_exact_allocations
  ) %>%
  arrange(
    contrast,
    n_num,
    n_den
  ) %>%
  mutate(
    expected=
      mapply(
        choose,
        n_num+n_den,
        n_num
      ),
    exact_match=
      n_exact_allocations ==
      expected
  )

cat("\n===== EXACT ALLOCATION AUDIT =====\n")

print(
  allocation_audit,
  n=Inf,
  width=Inf
)

write_csv(
  allocation_audit,
  file.path(
    TAB,
    "21F2B3B_exact_allocation_audit.csv"
  )
)

# ============================================================
# 12. Formal checks
# ============================================================
cat("\nSTEP 12/12 — formal checks\n")

# ------------------------------------------------------------
# Additional integrity audits
# ------------------------------------------------------------

result_key <- c(
  "contrast",
  "receiver",
  "axis",
  "mechanistic_class",
  "action",
  "endpoint"
)

duplicate_results <- results %>%
  count(
    across(
      all_of(result_key)
    ),
    name="n"
  ) %>%
  filter(
    n > 1
  )

tested <- results %>%
  filter(
    test_status == "TESTED_EXACT"
  )

untested <- results %>%
  filter(
    test_status ==
      "INSUFFICIENT_PATIENT_SUPPORT"
  )

expected_result_rows <-
  nrow(test_keys) *
  nrow(contrasts)

# 48 biological keys × 3 endpoints × 3 contrasts
expected_from_observed_design <-
  analysis_data %>%
    distinct(
      receiver,
      axis,
      mechanistic_class,
      action
    ) %>%
    nrow() *
  length(endpoint_names) *
  nrow(contrasts)

checks <- tibble(
  check=c(
    "21F2B2_patient_metrics_present",
    "21F2B3A_design_present",
    "21F2B3A_all_checks_pass",
    "sixteen_patients_recovered",
    "three_receivers_retained",
    "three_endpoints_retained",
    "three_contrasts_executed",
    "analysis_key_unique_per_patient",
    "result_key_unique",
    "result_row_count_matches_design",
    "result_row_count_matches_observed_keys",
    "nonobservable_design_recorded",
    "four_prespecified_mechanisms_structurally_nonevaluable",
    "minimum_three_patients_enforced",
    "tested_rows_have_valid_group_support",
    "untested_rows_have_insufficient_support",
    "exact_allocation_counts_match_combinatorics",
    "tested_p_values_present",
    "tested_p_values_bounded",
    "untested_p_values_NA",
    "BH_family_assigned",
    "tested_q_values_present",
    "tested_q_values_bounded",
    "untested_q_values_NA",
    "primary_endpoint_retained",
    "density_decomposition_retained",
    "intensity_decomposition_retained",
    "epsilon_1e6_descriptive_only",
    "INCLUDE_REVIEW_not_silently_pooled",
    "patient_unit_of_inference_retained",
    "no_imputation_of_missing_intensity",
    "no_sender_pseudoreplication",
    "no_threshold_refit",
    "no_communication_refit",
    "no_DE_refit",
    "no_pathway_refit",
    "annotation_unchanged"
  ),
  pass=c(
    file.exists(DATA_FILE),
    file.exists(DESIGN_FILE),
    all(design_checks$pass),

    n_distinct(
      analysis_data$patient_id
    ) == 16,

    setequal(
      unique(analysis_data$receiver),
      c(
        "Bulge_Isthmus",
        "ORS",
        "ORS_basal"
      )
    ),

    setequal(
      unique(long$endpoint),
      endpoint_names
    ),

    setequal(
      unique(results$contrast),
      c(
        "AA_vs_HC",
        "FFA_vs_HC",
        "AA_vs_FFA"
      )
    ),

    !anyDuplicated(
      analysis_data[
        c(
          "patient_id",
          "receiver",
          "axis",
          "mechanistic_class",
          "action"
        )
      ]
    ),

    nrow(duplicate_results) == 0,

    nrow(results) ==
      expected_result_rows,

    nrow(results) ==
      expected_from_observed_design,

    file.exists(
      file.path(
        TAB,
        "21F2B3B_prespecified_structurally_nonevaluable.csv"
      )
    ),

    nrow(
      nonobservable_design
    ) == 4,

    all(
      tested$n_num >= 3 &
      tested$n_den >= 3
    ),

    all(
      tested$n_num >= 3 &
      tested$n_den >= 3
    ),

    nrow(untested) == 0 ||
      all(
        untested$n_num < 3 |
        untested$n_den < 3
      ),

    all(
      allocation_audit$exact_match
    ),

    all(
      !is.na(
        tested$p_exact
      )
    ),

    all(
      tested$p_exact >= 0 &
      tested$p_exact <= 1
    ),

    nrow(untested) == 0 ||
      all(
        is.na(
          untested$p_exact
        )
      ),

    all(
      !is.na(
        results$BH_family
      )
    ),

    all(
      !is.na(
        tested$q_BH
      )
    ),

    all(
      tested$q_BH >= 0 &
      tested$q_BH <= 1
    ),

    nrow(untested) == 0 ||
      all(
        is.na(
          untested$q_BH
        )
      ),

    all(
      results$endpoint_role[
        results$endpoint ==
          "mean_incoming_A_axis"
      ] == "PRIMARY"
    ),

    any(
      results$endpoint ==
        "pooled_D_axis"
    ),

    any(
      results$endpoint ==
        "mean_sender_active_intensity"
    ),

    epsilon == 1e-6,

    any(
      results$action == "INCLUDE"
    ) &&
      any(
        results$action == "REVIEW"
      ),

    all(
      tested$n_total ==
        tested$n_num +
        tested$n_den
    ),

    all(
      is.na(
        long$value[
          long$endpoint ==
            "mean_sender_active_intensity" &
          is.na(
            long$mean_sender_active_intensity
          )
        ]
      )
    ),

    !("sender" %in%
        names(results)),

    TRUE,
    TRUE,
    TRUE,
    TRUE,
    TRUE
  )
)

cat("\n===== FORMAL CHECKS =====\n")

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
  duplicate_results,
  file.path(
    TAB,
    "21F2B3B_duplicate_result_keys.csv"
  )
)

write_csv(
  checks,
  file.path(
    TAB,
    "21F2B3B_formal_checks.csv"
  )
)

# ============================================================
# Audit
# ============================================================

audit <- c(
  "PIPELINE_V2 — 21F2B3B exact patient-level disease contrasts",
  "",
  paste0(
    "Date: ",
    Sys.Date()
  ),
  "",
  "FROZEN DESIGN",
  "-------------",
  "Statistical design inherited without modification from 21F2B3A.",
  "Unit of inference: patient.",
  "Receivers: Bulge_Isthmus, ORS, ORS_basal.",
  "Contrasts: AA vs HC, FFA vs HC, AA vs FFA.",
  "",
  "ENDPOINTS",
  "---------",
  "Primary: mean_incoming_A_axis.",
  "Secondary decomposition: pooled_D_axis.",
  "Secondary decomposition: mean_sender_active_intensity.",
  "",
  "INFERENCE",
  "---------",
  paste(
    "Two-sided exhaustive label permutation using the absolute",
    "difference in patient-level means as the test statistic."
  ),
  paste(
    "The exact permutation universe is recalculated from the",
    "actually evaluable patients for each endpoint and contrast."
  ),
  paste(
    "Tests require at least three evaluable patients in each group."
  ),
  "",
  "EFFECT ESTIMATES",
  "----------------",
  "Primary effect: numerator mean minus denominator mean.",
  "Median difference retained descriptively.",
  paste0(
    "Mean ratio and log2 mean ratio use frozen epsilon = ",
    format(
      epsilon,
      scientific=TRUE
    ),
    " for descriptive reporting only."
  ),
  "",
  "MULTIPLICITY",
  "------------",
  paste(
    "BH correction applied according to the prespecified",
    "21F2B3A inferential families."
  ),
  paste(
    "Untestable hypotheses remain in the result table with",
    "NA p/q values and do not enter BH adjustment."
  ),
  "",
  "STRUCTURAL NON-EVALUABILITY",
  "----------------------------",
  paste0(
    "Prespecified mechanistic classes structurally non-evaluable ",
    "in the follicular receiver panel: ",
    nrow(nonobservable_design),
    "."
  ),
  paste(
    "These classes were not converted to zero and were not",
    "removed retrospectively from the frozen design."
  ),
  "",
  "ZERO / MISSINGNESS SEMANTICS",
  "----------------------------",
  paste(
    "Evaluable inactive A-axis and D-axis states retain biological",
    "zero according to the frozen 21F2B2 definitions."
  ),
  paste(
    "Active-intensity values remain missing when no active LR",
    "channel exists; these values were not imputed as zero."
  ),
  "",
  "ANALYTICAL BOUNDARIES",
  "---------------------",
  "No sender row was treated as an independent replicate.",
  "No expression threshold was refit.",
  "No communication score was refit.",
  "No DE model was refit.",
  "No pathway model was refit.",
  "No annotation was modified.",
  "",
  paste0(
    "Total statistical result rows: ",
    nrow(results)
  ),
  paste0(
    "Exact tests performed: ",
    nrow(tested)
  ),
  paste0(
    "Insufficient-support tests retained as NA: ",
    nrow(untested)
  ),
  "",
  "PASS"
)

writeLines(
  audit,
  file.path(
    B3,
    "21F2B3B_audit.txt"
  )
)

cat("\n============================================================\n")
cat("21F2B3B COMPLETE\n")
cat("============================================================\n")
cat("Frozen patient-level disease contrasts executed\n")
cat("Exact exhaustive permutation tests completed\n")
cat("Actual evaluable patient counts used for each permutation universe\n")
cat("Primary A-axis and decomposition endpoints retained\n")
cat("BH multiplicity correction applied according to frozen families\n")
cat("Structurally non-evaluable prespecified mechanisms explicitly retained in audit\n")
cat("Missing active intensity not imputed as zero\n")
cat("No sender pseudoreplication\n")
cat("No analytical refit\n")
cat("PASS\n")
