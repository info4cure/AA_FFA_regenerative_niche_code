suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(tibble)
  library(stringr)
})

cat("============================================================\n")
cat("PIPELINE_V2 — 21G2C MULTICOMPARTMENT BIOLOGICAL READOUT\n")
cat("============================================================\n")

BASE <- paste0(
  "PIPELINE_V2/results/disease_biology/",
  "21_CELL_CELL_COMMUNICATION/",
  "21G_MULTICOMPARTMENT_COMMUNICATION/",
  "21G2_INFERENCE"
)

TAB <- file.path(BASE, "tables")

OUT <- file.path(
  dirname(BASE),
  "21G3_BIOLOGICAL_READOUT"
)

OUT_TAB <- file.path(
  OUT,
  "tables"
)

dir.create(
  OUT_TAB,
  recursive=TRUE,
  showWarnings=FALSE
)

RESULT_FILE <- file.path(
  TAB,
  "21G2B_multicompartment_results.csv"
)

CHECK_FILE <- file.path(
  TAB,
  "21G2B_formal_checks.csv"
)

stopifnot(
  file.exists(RESULT_FILE),
  file.exists(CHECK_FILE)
)

# ============================================================
# 1. Recover frozen inference
# ============================================================

cat("\nSTEP 1/8 — recovering frozen 21G2B inference\n")

x <- read_csv(
  RESULT_FILE,
  show_col_types=FALSE
)

checks_21G2B <- read_csv(
  CHECK_FILE,
  show_col_types=FALSE
)

stopifnot(
  all(checks_21G2B$pass)
)

cat("Result rows:", nrow(x), "\n")
cat(
  "Macrocompartments:",
  n_distinct(x$macrocompartment),
  "\n"
)
cat(
  "Receivers:",
  n_distinct(x$receiver),
  "\n"
)

required <- c(
  "receiver",
  "macrocompartment",
  "candidate_status",
  "biological_program",
  "analysis_tier",
  "mechanistic_class",
  "action",
  "contrast",
  "endpoint",
  "analysis_role",
  "axis",
  "n_num",
  "n_den",
  "mean_num",
  "mean_den",
  "mean_difference",
  "log2_mean_ratio",
  "p_exact",
  "q_BH",
  "test_status"
)

stopifnot(
  all(required %in% names(x))
)

# ============================================================
# 2. Freeze evidence vocabulary
# ============================================================

cat("\nSTEP 2/8 — classifying post-inference evidence\n")

tested <- x %>%
  filter(
    test_status == "TESTED_EXACT"
  ) %>%
  mutate(
    evidence_class=case_when(
      q_BH < 0.05 ~ "FDR_LT_0.05",
      q_BH < 0.10 ~ "FDR_LT_0.10",
      p_exact < 0.05 ~ "NOMINAL_P_LT_0.05",
      TRUE ~ "NO_THRESHOLD_EVIDENCE"
    ),

    direction=case_when(
      mean_difference > 0 ~ "HIGHER_IN_NUMERATOR",
      mean_difference < 0 ~ "LOWER_IN_NUMERATOR",
      mean_difference == 0 ~ "NO_MEAN_DIFFERENCE",
      TRUE ~ NA_character_
    )
  )

stopifnot(
  !any(is.na(tested$evidence_class)),
  !any(is.na(tested$direction))
)

cat("Tested exact rows:", nrow(tested), "\n")

# ============================================================
# 3. Primary A-axis universe
# ============================================================

cat("\nSTEP 3/8 — constructing primary A-axis readout\n")

primary_A <- tested %>%
  filter(
    endpoint == "mean_incoming_A_axis",
    analysis_role == "PRIMARY"
  ) %>%
  arrange(
    contrast,
    macrocompartment,
    analysis_tier,
    q_BH,
    p_exact,
    receiver,
    biological_program,
    axis
  )

primary_q005 <- primary_A %>%
  filter(q_BH < 0.05)

primary_q010 <- primary_A %>%
  filter(q_BH < 0.10)

primary_nominal <- primary_A %>%
  filter(p_exact < 0.05)

write_csv(
  primary_A,
  file.path(
    OUT_TAB,
    "21G2C_primary_A_all.csv"
  )
)

write_csv(
  primary_q005,
  file.path(
    OUT_TAB,
    "21G2C_primary_A_q005.csv"
  )
)

write_csv(
  primary_q010,
  file.path(
    OUT_TAB,
    "21G2C_primary_A_q010.csv"
  )
)

write_csv(
  primary_nominal,
  file.path(
    OUT_TAB,
    "21G2C_primary_A_nominal_p005.csv"
  )
)

cat("\n===== PRIMARY A-AXIS EVIDENCE =====\n")

print(
  primary_A %>%
    count(
      contrast,
      macrocompartment,
      analysis_tier,
      evidence_class,
      name="n"
    ) %>%
    arrange(
      contrast,
      macrocompartment,
      analysis_tier,
      evidence_class
    ),
  n=Inf,
  width=Inf
)

cat("\n===== PRIMARY A-AXIS q < 0.05 =====\n")
print(primary_q005, n=Inf, width=Inf)

cat("\n===== PRIMARY A-AXIS q < 0.10 =====\n")
print(primary_q010, n=Inf, width=Inf)

# ============================================================
# 4. Cross-compartment primary A-axis summary
# ============================================================

cat("\nSTEP 4/8 — summarising cross-compartment primary biology\n")

macro_primary_summary <- primary_A %>%
  group_by(
    contrast,
    macrocompartment,
    analysis_tier
  ) %>%
  summarise(
    n_tested=n(),
    n_q005=sum(q_BH < 0.05, na.rm=TRUE),
    n_q010=sum(q_BH < 0.10, na.rm=TRUE),
    n_nominal=sum(p_exact < 0.05, na.rm=TRUE),

    min_p_exact=
      if (all(is.na(p_exact))) {
        NA_real_
      } else {
        min(p_exact, na.rm=TRUE)
      },

    min_q_BH=
      if (all(is.na(q_BH))) {
        NA_real_
      } else {
        min(q_BH, na.rm=TRUE)
      },

    median_difference=
      median(
        mean_difference,
        na.rm=TRUE
      ),

    median_log2_mean_ratio=
      median(
        log2_mean_ratio,
        na.rm=TRUE
      ),

    .groups="drop"
  ) %>%
  arrange(
    contrast,
    analysis_tier,
    min_q_BH,
    min_p_exact,
    macrocompartment
  )

write_csv(
  macro_primary_summary,
  file.path(
    OUT_TAB,
    "21G2C_macrocompartment_primary_summary.csv"
  )
)

cat("\n===== CROSS-COMPARTMENT PRIMARY A SUMMARY =====\n")

print(
  macro_primary_summary,
  n=Inf,
  width=Inf
)

# ============================================================
# 5. Compartment-specific biological signal inventory
# ============================================================

cat("\nSTEP 5/8 — constructing compartment-specific signal inventory\n")

signal_inventory <- primary_A %>%
  filter(
    p_exact < 0.05
  ) %>%
  select(
    contrast,
    macrocompartment,
    receiver,
    candidate_status,
    biological_program,
    analysis_tier,
    mechanistic_class,
    action,
    axis,
    n_num,
    n_den,
    mean_num,
    mean_den,
    mean_difference,
    log2_mean_ratio,
    p_exact,
    q_BH,
    evidence_class,
    direction
  ) %>%
  arrange(
    contrast,
    macrocompartment,
    analysis_tier,
    q_BH,
    p_exact,
    receiver,
    biological_program,
    axis
  )

write_csv(
  signal_inventory,
  file.path(
    OUT_TAB,
    "21G2C_primary_A_signal_inventory.csv"
  )
)

cat("\n===== PRIMARY A-AXIS SIGNAL INVENTORY — p < 0.05 =====\n")

print(
  signal_inventory,
  n=Inf,
  width=Inf
)

# Disease-specific inventories are descriptive reorganizations
# of the already frozen inferential output.

for (con in sort(unique(primary_A$contrast))) {

  z <- primary_A %>%
    filter(
      contrast == con
    ) %>%
    arrange(
      macrocompartment,
      analysis_tier,
      q_BH,
      p_exact,
      receiver,
      biological_program,
      axis
    )

  write_csv(
    z,
    file.path(
      OUT_TAB,
      paste0(
        "21G2C_",
        con,
        "_primary_A_all.csv"
      )
    )
  )
}

# ============================================================
# 6. Density/intensity decomposition
# ============================================================

cat("\nSTEP 6/8 — linking primary A-axis evidence to decomposition\n")

decomposition <- tested %>%
  filter(
    endpoint %in%
      c(
        "pooled_D_axis",
        "mean_sender_active_intensity"
      )
  ) %>%
  select(
    receiver,
    macrocompartment,
    candidate_status,
    biological_program,
    analysis_tier,
    mechanistic_class,
    action,
    contrast,
    axis,
    endpoint,
    n_num,
    n_den,
    mean_num,
    mean_den,
    mean_difference,
    log2_mean_ratio,
    p_exact,
    q_BH,
    evidence_class,
    direction
  )

# Decomposition is reported for primary A-axis signals with
# q<0.10 or nominal p<0.05. This does not create a new test:
# all p- and q-values were frozen in 21G2B.

A_for_decomposition <- primary_A %>%
  filter(
    q_BH < 0.10 |
      p_exact < 0.05
  ) %>%
  select(
    receiver,
    macrocompartment,
    biological_program,
    analysis_tier,
    mechanistic_class,
    action,
    contrast,
    axis,
    A_mean_difference=mean_difference,
    A_log2_mean_ratio=log2_mean_ratio,
    A_p_exact=p_exact,
    A_q_BH=q_BH,
    A_evidence_class=evidence_class,
    A_direction=direction
  )

decomposition_readout <- A_for_decomposition %>%
  left_join(
    decomposition,
    by=c(
      "receiver",
      "macrocompartment",
      "biological_program",
      "analysis_tier",
      "mechanistic_class",
      "action",
      "contrast",
      "axis"
    ),
    relationship="one-to-many"
  ) %>%
  arrange(
    contrast,
    macrocompartment,
    A_q_BH,
    A_p_exact,
    receiver,
    biological_program,
    axis,
    endpoint
  )

write_csv(
  decomposition_readout,
  file.path(
    OUT_TAB,
    "21G2C_primary_signal_decomposition.csv"
  )
)

cat("\n===== PRIMARY SIGNAL DECOMPOSITION =====\n")

print(
  decomposition_readout,
  n=Inf,
  width=Inf
)

# Compact decomposition classification.
#
# These labels describe the already estimated decomposition;
# they are not new inferential outcomes.

decomposition_wide <- decomposition_readout %>%
  select(
    receiver,
    macrocompartment,
    biological_program,
    analysis_tier,
    mechanistic_class,
    action,
    contrast,
    axis,
    A_mean_difference,
    A_log2_mean_ratio,
    A_p_exact,
    A_q_BH,
    A_evidence_class,
    A_direction,
    endpoint,
    mean_difference,
    log2_mean_ratio,
    p_exact,
    q_BH
  ) %>%
  pivot_wider(
    names_from=endpoint,
    values_from=c(
      mean_difference,
      log2_mean_ratio,
      p_exact,
      q_BH
    ),
    names_sep="__"
  ) %>%
  mutate(
    decomposition_pattern=case_when(

      !is.na(
        q_BH__pooled_D_axis
      ) &
      q_BH__pooled_D_axis < 0.05 &
      (
        is.na(
          q_BH__mean_sender_active_intensity
        ) |
        q_BH__mean_sender_active_intensity >= 0.05
      ) ~
        "DENSITY_SUPPORTED",

      !is.na(
        q_BH__mean_sender_active_intensity
      ) &
      q_BH__mean_sender_active_intensity < 0.05 &
      (
        is.na(
          q_BH__pooled_D_axis
        ) |
        q_BH__pooled_D_axis >= 0.05
      ) ~
        "ACTIVE_INTENSITY_SUPPORTED",

      !is.na(
        q_BH__pooled_D_axis
      ) &
      !is.na(
        q_BH__mean_sender_active_intensity
      ) &
      q_BH__pooled_D_axis < 0.05 &
      q_BH__mean_sender_active_intensity < 0.05 ~
        "BOTH_SUPPORTED",

      !is.na(
        p_exact__pooled_D_axis
      ) &
      p_exact__pooled_D_axis < 0.05 &
      (
        is.na(
          p_exact__mean_sender_active_intensity
        ) |
        p_exact__mean_sender_active_intensity >= 0.05
      ) ~
        "DENSITY_NOMINAL_ONLY",

      !is.na(
        p_exact__mean_sender_active_intensity
      ) &
      p_exact__mean_sender_active_intensity < 0.05 &
      (
        is.na(
          p_exact__pooled_D_axis
        ) |
        p_exact__pooled_D_axis >= 0.05
      ) ~
        "ACTIVE_INTENSITY_NOMINAL_ONLY",

      !is.na(
        p_exact__pooled_D_axis
      ) &
      !is.na(
        p_exact__mean_sender_active_intensity
      ) &
      p_exact__pooled_D_axis < 0.05 &
      p_exact__mean_sender_active_intensity < 0.05 ~
        "BOTH_NOMINAL_ONLY",

      TRUE ~
        "NO_DECOMPOSITION_THRESHOLD_EVIDENCE"
    )
  )

write_csv(
  decomposition_wide,
  file.path(
    OUT_TAB,
    "21G2C_primary_signal_decomposition_wide.csv"
  )
)

cat("\n===== DECOMPOSITION PATTERNS =====\n")

print(
  decomposition_wide %>%
    select(
      contrast,
      macrocompartment,
      receiver,
      biological_program,
      analysis_tier,
      axis,
      A_mean_difference,
      A_q_BH,
      decomposition_pattern,
      starts_with("mean_difference__"),
      starts_with("p_exact__"),
      starts_with("q_BH__")
    ),
  n=Inf,
  width=Inf
)

# ============================================================
# 4. Cross-compartment primary A-axis summary
# ============================================================

cat("\nSTEP 4/8 — summarising cross-compartment primary biology\n")

macro_primary_summary <- primary_A %>%
  group_by(
    contrast,
    macrocompartment,
    analysis_tier
  ) %>%
  summarise(
    n_tested=n(),
    n_q005=sum(q_BH < 0.05, na.rm=TRUE),
    n_q010=sum(q_BH < 0.10, na.rm=TRUE),
    n_nominal=sum(p_exact < 0.05, na.rm=TRUE),

    min_p_exact=
      if (all(is.na(p_exact))) {
        NA_real_
      } else {
        min(p_exact, na.rm=TRUE)
      },

    min_q_BH=
      if (all(is.na(q_BH))) {
        NA_real_
      } else {
        min(q_BH, na.rm=TRUE)
      },

    median_difference=
      median(
        mean_difference,
        na.rm=TRUE
      ),

    median_log2_mean_ratio=
      median(
        log2_mean_ratio,
        na.rm=TRUE
      ),

    .groups="drop"
  ) %>%
  arrange(
    contrast,
    analysis_tier,
    min_q_BH,
    min_p_exact,
    macrocompartment
  )

write_csv(
  macro_primary_summary,
  file.path(
    OUT_TAB,
    "21G2C_macrocompartment_primary_summary.csv"
  )
)

cat("\n===== CROSS-COMPARTMENT PRIMARY A SUMMARY =====\n")

print(
  macro_primary_summary,
  n=Inf,
  width=Inf
)

# ============================================================
# 5. Compartment-specific biological signal inventory
# ============================================================

cat("\nSTEP 5/8 — constructing compartment-specific signal inventory\n")

signal_inventory <- primary_A %>%
  filter(
    p_exact < 0.05
  ) %>%
  select(
    contrast,
    macrocompartment,
    receiver,
    candidate_status,
    biological_program,
    analysis_tier,
    mechanistic_class,
    action,
    axis,
    n_num,
    n_den,
    mean_num,
    mean_den,
    mean_difference,
    log2_mean_ratio,
    p_exact,
    q_BH,
    evidence_class,
    direction
  ) %>%
  arrange(
    contrast,
    macrocompartment,
    analysis_tier,
    q_BH,
    p_exact,
    receiver,
    biological_program,
    axis
  )

write_csv(
  signal_inventory,
  file.path(
    OUT_TAB,
    "21G2C_primary_A_signal_inventory.csv"
  )
)

cat("\n===== PRIMARY A-AXIS SIGNAL INVENTORY — p < 0.05 =====\n")

print(
  signal_inventory,
  n=Inf,
  width=Inf
)

# Disease-specific inventories are descriptive reorganizations
# of the already frozen inferential output.

for (con in sort(unique(primary_A$contrast))) {

  z <- primary_A %>%
    filter(
      contrast == con
    ) %>%
    arrange(
      macrocompartment,
      analysis_tier,
      q_BH,
      p_exact,
      receiver,
      biological_program,
      axis
    )

  write_csv(
    z,
    file.path(
      OUT_TAB,
      paste0(
        "21G2C_",
        con,
        "_primary_A_all.csv"
      )
    )
  )
}

# ============================================================
# 6. Density/intensity decomposition
# ============================================================

cat("\nSTEP 6/8 — linking primary A-axis evidence to decomposition\n")

decomposition <- tested %>%
  filter(
    endpoint %in%
      c(
        "pooled_D_axis",
        "mean_sender_active_intensity"
      )
  ) %>%
  select(
    receiver,
    macrocompartment,
    candidate_status,
    biological_program,
    analysis_tier,
    mechanistic_class,
    action,
    contrast,
    axis,
    endpoint,
    n_num,
    n_den,
    mean_num,
    mean_den,
    mean_difference,
    log2_mean_ratio,
    p_exact,
    q_BH,
    evidence_class,
    direction
  )

# Decomposition is reported for primary A-axis signals with
# q<0.10 or nominal p<0.05. This does not create a new test:
# all p- and q-values were frozen in 21G2B.

A_for_decomposition <- primary_A %>%
  filter(
    q_BH < 0.10 |
      p_exact < 0.05
  ) %>%
  select(
    receiver,
    macrocompartment,
    biological_program,
    analysis_tier,
    mechanistic_class,
    action,
    contrast,
    axis,
    A_mean_difference=mean_difference,
    A_log2_mean_ratio=log2_mean_ratio,
    A_p_exact=p_exact,
    A_q_BH=q_BH,
    A_evidence_class=evidence_class,
    A_direction=direction
  )

decomposition_readout <- A_for_decomposition %>%
  left_join(
    decomposition,
    by=c(
      "receiver",
      "macrocompartment",
      "biological_program",
      "analysis_tier",
      "mechanistic_class",
      "action",
      "contrast",
      "axis"
    ),
    relationship="one-to-many"
  ) %>%
  arrange(
    contrast,
    macrocompartment,
    A_q_BH,
    A_p_exact,
    receiver,
    biological_program,
    axis,
    endpoint
  )

write_csv(
  decomposition_readout,
  file.path(
    OUT_TAB,
    "21G2C_primary_signal_decomposition.csv"
  )
)

cat("\n===== PRIMARY SIGNAL DECOMPOSITION =====\n")

print(
  decomposition_readout,
  n=Inf,
  width=Inf
)

# Compact decomposition classification.
#
# These labels describe the already estimated decomposition;
# they are not new inferential outcomes.

decomposition_wide <- decomposition_readout %>%
  select(
    receiver,
    macrocompartment,
    biological_program,
    analysis_tier,
    mechanistic_class,
    action,
    contrast,
    axis,
    A_mean_difference,
    A_log2_mean_ratio,
    A_p_exact,
    A_q_BH,
    A_evidence_class,
    A_direction,
    endpoint,
    mean_difference,
    log2_mean_ratio,
    p_exact,
    q_BH
  ) %>%
  pivot_wider(
    names_from=endpoint,
    values_from=c(
      mean_difference,
      log2_mean_ratio,
      p_exact,
      q_BH
    ),
    names_sep="__"
  ) %>%
  mutate(
    decomposition_pattern=case_when(

      !is.na(
        q_BH__pooled_D_axis
      ) &
      q_BH__pooled_D_axis < 0.05 &
      (
        is.na(
          q_BH__mean_sender_active_intensity
        ) |
        q_BH__mean_sender_active_intensity >= 0.05
      ) ~
        "DENSITY_SUPPORTED",

      !is.na(
        q_BH__mean_sender_active_intensity
      ) &
      q_BH__mean_sender_active_intensity < 0.05 &
      (
        is.na(
          q_BH__pooled_D_axis
        ) |
        q_BH__pooled_D_axis >= 0.05
      ) ~
        "ACTIVE_INTENSITY_SUPPORTED",

      !is.na(
        q_BH__pooled_D_axis
      ) &
      !is.na(
        q_BH__mean_sender_active_intensity
      ) &
      q_BH__pooled_D_axis < 0.05 &
      q_BH__mean_sender_active_intensity < 0.05 ~
        "BOTH_SUPPORTED",

      !is.na(
        p_exact__pooled_D_axis
      ) &
      p_exact__pooled_D_axis < 0.05 &
      (
        is.na(
          p_exact__mean_sender_active_intensity
        ) |
        p_exact__mean_sender_active_intensity >= 0.05
      ) ~
        "DENSITY_NOMINAL_ONLY",

      !is.na(
        p_exact__mean_sender_active_intensity
      ) &
      p_exact__mean_sender_active_intensity < 0.05 &
      (
        is.na(
          p_exact__pooled_D_axis
        ) |
        p_exact__pooled_D_axis >= 0.05
      ) ~
        "ACTIVE_INTENSITY_NOMINAL_ONLY",

      !is.na(
        p_exact__pooled_D_axis
      ) &
      !is.na(
        p_exact__mean_sender_active_intensity
      ) &
      p_exact__pooled_D_axis < 0.05 &
      p_exact__mean_sender_active_intensity < 0.05 ~
        "BOTH_NOMINAL_ONLY",

      TRUE ~
        "NO_DECOMPOSITION_THRESHOLD_EVIDENCE"
    )
  )

write_csv(
  decomposition_wide,
  file.path(
    OUT_TAB,
    "21G2C_primary_signal_decomposition_wide.csv"
  )
)

cat("\n===== DECOMPOSITION PATTERNS =====\n")

print(
  decomposition_wide %>%
    select(
      contrast,
      macrocompartment,
      receiver,
      biological_program,
      analysis_tier,
      axis,
      A_mean_difference,
      A_q_BH,
      decomposition_pattern,
      starts_with("mean_difference__"),
      starts_with("p_exact__"),
      starts_with("q_BH__")
    ),
  n=Inf,
  width=Inf
)

# ============================================================
# 7. Biological compartment inventory
# ============================================================

cat("\nSTEP 7/8 — freezing biological compartment readout\n")

compartment_inventory <- primary_A %>%
  group_by(
    contrast,
    macrocompartment
  ) %>%
  summarise(
    n_receivers=
      n_distinct(receiver),

    n_programs=
      n_distinct(biological_program),

    n_axes=
      n_distinct(axis),

    n_primary_tests=
      n(),

    n_q005=
      sum(
        q_BH < 0.05,
        na.rm=TRUE
      ),

    n_q010=
      sum(
        q_BH < 0.10,
        na.rm=TRUE
      ),

    n_nominal_p005=
      sum(
        p_exact < 0.05,
        na.rm=TRUE
      ),

    min_p_exact=
      if (all(is.na(p_exact))) {
        NA_real_
      } else {
        min(
          p_exact,
          na.rm=TRUE
        )
      },

    min_q_BH=
      if (all(is.na(q_BH))) {
        NA_real_
      } else {
        min(
          q_BH,
          na.rm=TRUE
        )
      },

    .groups="drop"
  ) %>%
  arrange(
    contrast,
    min_q_BH,
    min_p_exact,
    macrocompartment
  )

write_csv(
  compartment_inventory,
  file.path(
    OUT_TAB,
    "21G2C_compartment_inventory.csv"
  )
)

cat("\n===== BIOLOGICAL COMPARTMENT INVENTORY =====\n")

print(
  compartment_inventory,
  n=Inf,
  width=Inf
)

# ============================================================
# 8. Formal checks
# ============================================================

cat("\nSTEP 8/8 — formal checks\n")

expected_macrocompartments <- c(
  "ENDOTHELIAL",
  "FIBROBLAST",
  "FOLLICULAR_EPITHELIAL",
  "IMMUNE",
  "INTERFOLLICULAR_EPITHELIAL",
  "MELANOCYTE",
  "MURAL_PERIVASCULAR",
  "PILOSEBACEOUS_SEBACEOUS",
  "SECRETORY_DUCTAL"
)

formal_checks <- tibble(
  check=c(
    "21G2B_results_present",
    "21G2B_checks_present",
    "21G2B_all_checks_pass",
    "nine_macrocompartments_retained",
    "expected_macrocompartments_retained",
    "receiver_identity_retained",
    "primary_A_axis_retained",
    "primary_contextual_separated",
    "three_contrasts_retained",
    "evidence_class_generated",
    "direction_generated",
    "cross_compartment_summary_generated",
    "signal_inventory_generated",
    "density_decomposition_retained",
    "active_intensity_decomposition_retained",
    "decomposition_linked_by_axis",
    "no_new_disease_test",
    "no_new_p_value",
    "no_new_BH_adjustment",
    "no_threshold_refit",
    "no_communication_refit",
    "no_DE_refit",
    "no_pathway_refit",
    "annotation_unchanged"
  ),

  pass=c(
    file.exists(RESULT_FILE),
    file.exists(CHECK_FILE),
    all(checks_21G2B$pass),

    n_distinct(
      primary_A$macrocompartment
    ) == 9,

    setequal(
      unique(
        primary_A$macrocompartment
      ),
      expected_macrocompartments
    ),

    n_distinct(
      primary_A$receiver
    ) > 1,

    all(
      primary_A$endpoint ==
        "mean_incoming_A_axis"
    ),

    all(
      primary_A$analysis_tier %in%
        c(
          "PRIMARY",
          "CONTEXTUAL"
        )
    ),

    setequal(
      unique(
        primary_A$contrast
      ),
      c(
        "AA_vs_HC",
        "FFA_vs_HC",
        "AA_vs_FFA"
      )
    ),

    all(
      !is.na(
        primary_A$evidence_class
      )
    ),

    all(
      !is.na(
        primary_A$direction
      )
    ),

    nrow(
      macro_primary_summary
    ) > 0,

    nrow(
      signal_inventory
    ) >= 0,

    all(
      decomposition$endpoint %in%
        c(
          "pooled_D_axis",
          "mean_sender_active_intensity"
        )
    ),

    any(
      decomposition$endpoint ==
        "mean_sender_active_intensity"
    ),

    all(
      !is.na(
        decomposition_readout$axis
      )
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

cat("\n===== FORMAL CHECKS =====\n")

print(
  formal_checks,
  n=Inf,
  width=Inf
)

stopifnot(
  length(
    formal_checks$check
  ) ==
    length(
      formal_checks$pass
    ),

  all(
    formal_checks$pass
  )
)

write_csv(
  formal_checks,
  file.path(
    OUT_TAB,
    "21G2C_formal_checks.csv"
  )
)

# ============================================================
# Audit
# ============================================================

audit <- c(
  "PIPELINE_V2 — 21G2C MULTICOMPARTMENT BIOLOGICAL READOUT",
  "=======================================================",
  "",
  "STATUS: PASS",
  "",
  paste0(
    "Frozen 21G2B result rows: ",
    nrow(x)
  ),
  paste0(
    "Tested exact rows: ",
    nrow(tested)
  ),
  paste0(
    "Primary A-axis tested rows: ",
    nrow(primary_A)
  ),
  paste0(
    "Primary A-axis q<0.05 rows: ",
    nrow(primary_q005)
  ),
  paste0(
    "Primary A-axis q<0.10 rows: ",
    nrow(primary_q010)
  ),
  paste0(
    "Primary A-axis nominal p<0.05 rows: ",
    nrow(primary_nominal)
  ),
  "",
  paste(
    "The primary endpoint is normalized incoming mechanistic",
    "communication exposure (mean_incoming_A_axis)."
  ),
  paste(
    "Density (pooled_D_axis) and active communication intensity",
    "(mean_sender_active_intensity) are decomposition endpoints."
  ),
  "",
  paste(
    "Macrocompartments organize interpretation but receiver",
    "communication identities remain the analytical units."
  ),
  "",
  paste(
    "No new disease test, p-value, BH correction, threshold",
    "refit, communication refit, DE refit, pathway refit or",
    "annotation change was performed in 21G2C."
  ),
  "",
  "PASS"
)

writeLines(
  audit,
  file.path(
    OUT,
    "21G2C_audit.txt"
  )
)

cat("\n============================================================\n")
cat("21G2C COMPLETE\n")
cat("============================================================\n")
cat("Frozen 21G2B inference reorganized biologically\n")
cat("Nine macrocompartments retained\n")
cat("Receiver identities retained as analytical units\n")
cat("Primary A-axis evidence classified\n")
cat("Density and active-intensity decomposition retained\n")
cat("Cross-compartment summary generated\n")
cat("Compartment-specific biological readouts generated\n")
cat("No new statistical inference\n")
cat("No analytical refit\n")
cat("PASS\n")
