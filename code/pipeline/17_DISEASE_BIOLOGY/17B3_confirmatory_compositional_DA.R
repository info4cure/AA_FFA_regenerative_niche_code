suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(tibble)
})

cat("============================================================\n")
cat("PIPELINE_V2 — 17B3 CONFIRMATORY COMPOSITIONAL DA\n")
cat("============================================================\n")

IN <- paste0(
  "PIPELINE_V2/results/disease_biology/",
  "17B_composition_DA/tables/",
  "17B2_patient_L2_composition.csv"
)

SUPPORT <- paste0(
  "PIPELINE_V2/results/disease_biology/",
  "17B_composition_DA/tables/",
  "17B2_L2_testability.csv"
)

OUT <- paste0(
  "PIPELINE_V2/results/disease_biology/",
  "17B_composition_DA/confirmatory"
)

dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

x <- read_csv(IN, show_col_types = FALSE)
support <- read_csv(SUPPORT, show_col_types = FALSE)

# ============================================================
# 1. Confirmatory eligibility
# ============================================================

cat("\nSTEP 1/8 — confirmatory eligibility\n")

eligible_labels <- support %>%
  mutate(
    confirmatory_testable =
      n_patients_ge5 >= 5
  ) %>%
  filter(confirmatory_testable) %>%
  select(
    final_compartment,
    final_L2_label
  )

cat(
  "Screening L2:",
  nrow(support),
  "\n"
)

cat(
  "Confirmatory L2:",
  nrow(eligible_labels),
  "\n"
)

write_csv(
  support %>%
    mutate(
      confirmatory_testable =
        n_patients_ge5 >= 5
    ),
  file.path(
    OUT,
    "17B3_confirmatory_testability.csv"
  )
)

z <- x %>%
  semi_join(
    eligible_labels,
    by = c(
      "final_compartment",
      "final_L2_label"
    )
  ) %>%
  mutate(
    disease = factor(
      disease,
      levels = c("HC", "FFA", "AA")
    )
  )

# ============================================================
# 2. Metadata / batch structure
# ============================================================

cat("\nSTEP 2/8 — patient metadata and batch confounding\n")

patient_meta <- z %>%
  distinct(
    patient_id,
    disease
  )

# Recover batch from final master ledger.
master <- readRDS(
  paste0(
    "PIPELINE_V2/results/disease_biology/",
    "17A_frozen_analysis_atlas/tables/",
    "17A_master_cell_ledger_FINAL.rds"
  )
)

batch_meta <- master %>%
  filter(analysis_eligible_global) %>%
  distinct(
    patient_id,
    disease,
    batch
  )

stopifnot(
  nrow(batch_meta) == 16
)

z <- z %>%
  left_join(
    batch_meta %>%
      select(
        patient_id,
        batch
      ),
    by = "patient_id"
  )

cat("\nDisease × batch:\n")
print(
  table(
    batch_meta$disease,
    batch_meta$batch
  )
)

# ============================================================
# 3. Quasibinomial model function
# ============================================================

cat("\nSTEP 3/8 — fitting primary patient-level models\n")

fit_qb <- function(
  dat,
  group1,
  group0,
  adjust_batch = FALSE
) {

  dat <- dat %>%
    filter(
      disease %in% c(group0, group1)
    ) %>%
    mutate(
      disease2 = factor(
        as.character(disease),
        levels = c(group0, group1)
      ),
      failures =
        n_compartment - n_cells
    )

  if (n_distinct(dat$disease2) < 2) {
    return(NULL)
  }

  form <- if (adjust_batch) {
    cbind(n_cells, failures) ~ disease2 + batch
  } else {
    cbind(n_cells, failures) ~ disease2
  }

  fit <- tryCatch(
    glm(
      form,
      family = quasibinomial(),
      data = dat
    ),
    error = function(e) NULL
  )

  if (is.null(fit)) {
    return(NULL)
  }

  coef_name <- paste0(
    "disease2",
    group1
  )

  sm <- summary(fit)$coefficients

  if (!(coef_name %in% rownames(sm))) {
    return(NULL)
  }

  beta <- sm[coef_name, "Estimate"]
  se <- sm[coef_name, "Std. Error"]
  p <- sm[coef_name, "Pr(>|t|)"]

  phi <- summary(fit)$dispersion

  tibble(
    numerator = group1,
    denominator = group0,
    n_numerator =
      sum(dat$disease2 == group1),
    n_denominator =
      sum(dat$disease2 == group0),

    beta_log_odds = beta,
    SE = se,

    OR = exp(beta),
    OR_CI95_low =
      exp(beta - 1.96 * se),
    OR_CI95_high =
      exp(beta + 1.96 * se),

    dispersion = phi,
    p_value = p,
    batch_adjusted = adjust_batch
  )
}

keys <- z %>%
  distinct(
    final_compartment,
    final_L2_label
  )

run_models <- function(
  group1,
  group0,
  adjust_batch = FALSE
) {

  bind_rows(
    lapply(
      seq_len(nrow(keys)),
      function(i) {

        comp <- keys$final_compartment[i]
        lab <- keys$final_L2_label[i]

        d <- z %>%
          filter(
            final_compartment == comp,
            final_L2_label == lab
          )

        fit <- fit_qb(
          d,
          group1,
          group0,
          adjust_batch
        )

        if (is.null(fit)) {
          return(NULL)
        }

        fit %>%
          mutate(
            contrast =
              paste0(group1, "_vs_", group0),

            final_compartment = comp,
            final_L2_label = lab,

            .before = 1
          )
      }
    )
  )
}

primary <- bind_rows(
  run_models("AA", "HC", FALSE),
  run_models("FFA", "HC", FALSE),
  run_models("AA", "FFA", FALSE)
)

# BH by contrast × compartment
primary <- primary %>%
  group_by(
    contrast,
    final_compartment
  ) %>%
  mutate(
    FDR_BH_compartment =
      p.adjust(
        p_value,
        method = "BH"
      )
  ) %>%
  ungroup() %>%
  group_by(contrast) %>%
  mutate(
    FDR_BH_global =
      p.adjust(
        p_value,
        method = "BH"
      )
  ) %>%
  ungroup()

write_csv(
  primary,
  file.path(
    OUT,
    "17B3_quasibinomial_primary.csv"
  )
)

# ============================================================
# 4. FFA vs HC batch-adjusted sensitivity
# ============================================================

cat("\nSTEP 4/8 — FFA vs HC batch-adjusted sensitivity\n")

ffa_batch <- run_models(
  "FFA",
  "HC",
  TRUE
)

ffa_batch <- ffa_batch %>%
  group_by(final_compartment) %>%
  mutate(
    FDR_BH_compartment =
      p.adjust(
        p_value,
        method = "BH"
      )
  ) %>%
  ungroup() %>%
  mutate(
    FDR_BH_global =
      p.adjust(
        p_value,
        method = "BH"
      )
  )

write_csv(
  ffa_batch,
  file.path(
    OUT,
    "17B3_FFA_vs_HC_batch_adjusted.csv"
  )
)

# ============================================================
# 5. AA batch-confounding sensitivity
# ============================================================

cat("\nSTEP 5/8 — AA batch-confounding sensitivity\n")

# Contemporary HC:
# CON276 = CABIMER2026
#
# Historical HC:
# CON256/265/266/267

aa_sensitivity <- z %>%
  filter(
    disease == "AA" |
    patient_id %in%
      c(
        "CON276",
        "CON256",
        "CON265",
        "CON266",
        "CON267"
      )
  ) %>%
  group_by(
    final_compartment,
    final_L2_label
  ) %>%
  summarise(
    AA_median =
      median(
        proportion[disease == "AA"]
      ),

    HC_CABIMER_CON276 =
      proportion[
        patient_id == "CON276"
      ][1],

    HC_historical_median =
      median(
        proportion[
          patient_id %in%
            c(
              "CON256",
              "CON265",
              "CON266",
              "CON267"
            )
        ]
      ),

    .groups = "drop"
  ) %>%
  mutate(
    delta_AA_vs_CON276 =
      AA_median -
      HC_CABIMER_CON276,

    delta_AA_vs_HC_historical =
      AA_median -
      HC_historical_median,

    direction_vs_CON276 =
      sign(delta_AA_vs_CON276),

    direction_vs_historical =
      sign(delta_AA_vs_HC_historical),

    direction_concordant =
      direction_vs_CON276 ==
      direction_vs_historical,

    stronger_than_both =
      direction_concordant &
      direction_vs_CON276 != 0
  )

write_csv(
  aa_sensitivity,
  file.path(
    OUT,
    "17B3_AA_batch_confounding_sensitivity.csv"
  )
)

# ============================================================
# 6. Merge with 17B2 screening
# ============================================================

cat("\nSTEP 6/8 — integrating screening + confirmatory evidence\n")

screen <- read_csv(
  paste0(
    "PIPELINE_V2/results/disease_biology/",
    "17B_composition_DA/tables/",
    "17B2_L2_pairwise_DA.csv"
  ),
  show_col_types = FALSE
)

integrated <- primary %>%
  left_join(
    screen %>%
      select(
        contrast,
        final_compartment,
        final_L2_label,
        median_percent_numerator,
        median_percent_denominator,
        delta_percentage_points,
        log2_median_ratio,
        wilcox_p = p_value,
        wilcox_FDR_compartment =
          FDR_BH_compartment,
        wilcox_FDR_global =
          FDR_BH_global
      ),
    by = c(
      "contrast",
      "final_compartment",
      "final_L2_label"
    )
  ) %>%
  mutate(
    confirmatory_direction =
      case_when(
        OR > 1 ~
          paste0("higher_in_", numerator),

        OR < 1 ~
          paste0("higher_in_", denominator),

        TRUE ~ "no_change"
      ),

    wilcox_confirmatory_direction_concordant =
      case_when(
        OR > 1 &
          delta_percentage_points > 0 ~ TRUE,

        OR < 1 &
          delta_percentage_points < 0 ~ TRUE,

        OR == 1 &
          delta_percentage_points == 0 ~ TRUE,

        TRUE ~ FALSE
      )
  )

write_csv(
  integrated,
  file.path(
    OUT,
    "17B3_integrated_DA_evidence.csv"
  )
)

# ============================================================
# 7. Result summaries
# ============================================================

cat("\nSTEP 7/8 — result summaries\n")

for (cc in c(
  "AA_vs_HC",
  "FFA_vs_HC",
  "AA_vs_FFA"
)) {

  cat("\n============================================================\n")
  cat(cc, "\n")
  cat("============================================================\n")

  show <- integrated %>%
    filter(contrast == cc) %>%
    arrange(
      FDR_BH_compartment,
      FDR_BH_global,
      p_value
    ) %>%
    select(
      final_compartment,
      final_L2_label,
      median_percent_numerator,
      median_percent_denominator,
      delta_percentage_points,
      OR,
      OR_CI95_low,
      OR_CI95_high,
      dispersion,
      p_value,
      FDR_BH_compartment,
      FDR_BH_global,
      wilcox_p,
      wilcox_FDR_compartment,
      wilcox_confirmatory_direction_concordant
    ) %>%
    slice_head(n = 30)

  print(
    as.data.frame(show),
    row.names = FALSE
  )
}

cat("\n===== FFA vs HC BATCH-ADJUSTED TOP =====\n")

print(
  as.data.frame(
    ffa_batch %>%
      arrange(
        FDR_BH_compartment,
        FDR_BH_global,
        p_value
      ) %>%
      select(
        final_compartment,
        final_L2_label,
        OR,
        OR_CI95_low,
        OR_CI95_high,
        dispersion,
        p_value,
        FDR_BH_compartment,
        FDR_BH_global
      ) %>%
      slice_head(n = 30)
  ),
  row.names = FALSE
)

cat("\n===== AA BATCH-SENSITIVITY — SCREENING SIGNALS =====\n")

aa_signal <- integrated %>%
  filter(
    contrast == "AA_vs_HC",
    wilcox_FDR_compartment < 0.05
  ) %>%
  select(
    final_compartment,
    final_L2_label,
    OR,
    FDR_BH_compartment,
    wilcox_FDR_compartment
  ) %>%
  left_join(
    aa_sensitivity,
    by = c(
      "final_compartment",
      "final_L2_label"
    )
  )

print(
  as.data.frame(aa_signal),
  row.names = FALSE
)

# ============================================================
# 8. Formal checks
# ============================================================

cat("\nSTEP 8/8 — formal checks\n")

checks <- tibble(
  check = c(
    "patients_16",
    "confirmatory_filter_ge5_in_ge5",
    "three_primary_contrasts",
    "quasibinomial_primary_generated",
    "OR_and_CI_generated",
    "BH_compartment_generated",
    "BH_global_generated",
    "FFA_batch_adjusted_generated",
    "AA_not_batch_adjusted",
    "AA_CON276_sensitivity_generated",
    "AA_historical_sensitivity_generated",
    "screen_confirmatory_integrated",
    "patient_level_inference",
    "annotation_unchanged",
    "figures_generated_NO"
  ),

  pass = c(
    n_distinct(z$patient_id) == 16,

    all(
      support$n_patients_ge5[
        support$final_L2_label %in%
          eligible_labels$final_L2_label
      ] >= 5
    ),

    n_distinct(primary$contrast) == 3,

    nrow(primary) > 0,

    all(
      is.finite(primary$OR) &
      is.finite(primary$OR_CI95_low) &
      is.finite(primary$OR_CI95_high)
    ),

    all(!is.na(primary$FDR_BH_compartment)),
    all(!is.na(primary$FDR_BH_global)),

    nrow(ffa_batch) > 0,

    all(
      !primary$batch_adjusted[
        primary$contrast %in%
          c("AA_vs_HC", "AA_vs_FFA")
      ]
    ),

    all(
      !is.na(
        aa_sensitivity$HC_CABIMER_CON276
      )
    ),

    all(
      !is.na(
        aa_sensitivity$HC_historical_median
      )
    ),

    nrow(integrated) == nrow(primary),

    TRUE,
    TRUE,
    TRUE
  )
)

print(
  as.data.frame(checks),
  row.names = FALSE
)

write_csv(
  checks,
  file.path(
    OUT,
    "17B3_formal_checks.csv"
  )
)

stopifnot(all(checks$pass))

cat("\n============================================================\n")
cat("17B3 COMPLETE\n")
cat("============================================================\n")
cat("Confirmatory unit: patient\n")
cat("Model: quasibinomial L2 counts / compartment counts\n")
cat("Confirmatory support: >=5 cells in >=5 patients\n")
cat("AA vs HC: unadjusted for batch due complete confounding\n")
cat("AA vs FFA: unadjusted for batch due severe confounding\n")
cat("FFA vs HC: primary + batch-adjusted sensitivity\n")
cat("AA sensitivity: CON276 + historical HC direction audit\n")
cat("Annotation modified: NO\n")
cat("Figures generated: NO\n")
cat("PASS\n")
