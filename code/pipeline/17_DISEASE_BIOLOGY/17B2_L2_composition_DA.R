suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(tibble)
})

cat("============================================================\n")
cat("PIPELINE_V2 — 17B2 FROZEN L2 COMPOSITION / DA\n")
cat("============================================================\n")

IN <- paste0(
  "PIPELINE_V2/results/disease_biology/",
  "17A_frozen_analysis_atlas/tables/",
  "17A_master_cell_ledger_FINAL.rds"
)

OUT <- paste0(
  "PIPELINE_V2/results/disease_biology/",
  "17B_composition_DA/tables"
)

dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

x <- readRDS(IN)

stopifnot(nrow(x) == 82893)

# Only formally frozen L2 cells.
z <- x %>%
  filter(
    analysis_eligible_global,
    analysis_eligible_L2,
    final_L2_status == "FROZEN",
    !is.na(final_L2_label),
    !is.na(final_compartment)
  ) %>%
  mutate(
    disease = factor(
      disease,
      levels = c("HC", "FFA", "AA")
    )
  )

cat("\nFrozen-L2 cells:", nrow(z), "\n")
cat("Patients:", n_distinct(z$patient_id), "\n")
cat(
  "Compartments:",
  paste(sort(unique(z$final_compartment)), collapse = "; "),
  "\n"
)


# ============================================================
# 1. Patient × compartment denominators
# ============================================================

cat("\nSTEP 1/8 — compartment denominators\n")

denom <- z %>%
  count(
    disease,
    patient_id,
    final_compartment,
    name = "n_compartment"
  )

write_csv(
  denom,
  file.path(
    OUT,
    "17B2_patient_compartment_denominators.csv"
  )
)


# ============================================================
# 2. L2 support — disease blind
# ============================================================

cat("\nSTEP 2/8 — disease-blind L2 support\n")

support <- z %>%
  count(
    final_compartment,
    final_L2_label,
    patient_id,
    name = "n_cells"
  ) %>%
  group_by(
    final_compartment,
    final_L2_label
  ) %>%
  summarise(
    total_cells = sum(n_cells),
    n_patients_any = n_distinct(patient_id),
    n_patients_ge5 = sum(n_cells >= 5),
    median_cells_positive_patient = median(n_cells),
    .groups = "drop"
  ) %>%
  mutate(
    testable = n_patients_ge5 >= 3
  )

write_csv(
  support,
  file.path(
    OUT,
    "17B2_L2_testability.csv"
  )
)

cat("\n===== L2 TESTABILITY =====\n")
print(
  as.data.frame(
    support %>%
      arrange(
        final_compartment,
        desc(total_cells)
      )
  ),
  row.names = FALSE
)

testable_labels <- support %>%
  filter(testable) %>%
  select(
    final_compartment,
    final_L2_label
  )


# ============================================================
# 3. Complete patient × L2 composition
# ============================================================

cat("\nSTEP 3/8 — patient-level L2 proportions\n")

# Only patients with cells in a given compartment contribute
# to that compartment-specific composition analysis.

patient_compartments <- denom %>%
  select(
    disease,
    patient_id,
    final_compartment,
    n_compartment
  )

l2_catalogue <- z %>%
  distinct(
    final_compartment,
    final_L2_label
  )

composition <- patient_compartments %>%
  inner_join(
    l2_catalogue,
    by = "final_compartment",
    relationship = "many-to-many"
  ) %>%
  left_join(
    z %>%
      count(
        disease,
        patient_id,
        final_compartment,
        final_L2_label,
        name = "n_cells"
      ),
    by = c(
      "disease",
      "patient_id",
      "final_compartment",
      "final_L2_label"
    )
  ) %>%
  mutate(
    n_cells = replace_na(n_cells, 0L),
    proportion = n_cells / n_compartment,
    percent = 100 * proportion
  )

write_csv(
  composition,
  file.path(
    OUT,
    "17B2_patient_L2_composition.csv"
  )
)


# ============================================================
# 4. Disease descriptive summaries
# ============================================================

cat("\nSTEP 4/8 — disease summaries\n")

summary_tbl <- composition %>%
  group_by(
    disease,
    final_compartment,
    final_L2_label
  ) %>%
  summarise(
    n_patients = n(),
    n_patients_positive = sum(n_cells > 0),
    median_percent = median(percent),
    q25_percent = quantile(percent, 0.25),
    q75_percent = quantile(percent, 0.75),
    mean_percent = mean(percent),
    sd_percent = sd(percent),
    .groups = "drop"
  )

write_csv(
  summary_tbl,
  file.path(
    OUT,
    "17B2_L2_disease_summary.csv"
  )
)


# ============================================================
# 5. Pairwise patient-level screening
# ============================================================

cat("\nSTEP 5/8 — pairwise L2 contrasts\n")

contrasts <- list(
  c("AA", "HC"),
  c("FFA", "HC"),
  c("AA", "FFA")
)

run_one <- function(dat, group1, group0) {

  dat <- dat %>%
    semi_join(
      testable_labels,
      by = c(
        "final_compartment",
        "final_L2_label"
      )
    )

  keys <- dat %>%
    distinct(
      final_compartment,
      final_L2_label
    )

  bind_rows(
    lapply(
      seq_len(nrow(keys)),
      function(i) {

        comp <- keys$final_compartment[i]
        label <- keys$final_L2_label[i]

        q <- dat %>%
          filter(
            final_compartment == comp,
            final_L2_label == label,
            disease %in% c(group1, group0)
          )

        a <- q$proportion[q$disease == group1]
        b <- q$proportion[q$disease == group0]

        if (length(a) == 0 || length(b) == 0) {
          return(NULL)
        }

        wt <- suppressWarnings(
          wilcox.test(
            a,
            b,
            exact = FALSE,
            conf.int = FALSE
          )
        )

        med1 <- median(a)
        med0 <- median(b)

        eps <- 1e-6

        tibble(
          contrast = paste0(group1, "_vs_", group0),
          numerator = group1,
          denominator = group0,

          final_compartment = comp,
          final_L2_label = label,

          n_numerator = length(a),
          n_denominator = length(b),

          positive_numerator = sum(a > 0),
          positive_denominator = sum(b > 0),

          median_prop_numerator = med1,
          median_prop_denominator = med0,

          median_percent_numerator = 100 * med1,
          median_percent_denominator = 100 * med0,

          delta_percentage_points =
            100 * (med1 - med0),

          log2_median_ratio =
            log2(
              (med1 + eps) /
              (med0 + eps)
            ),

          wilcox_W = unname(wt$statistic),
          p_value = wt$p.value
        )
      }
    )
  )
}

tests <- bind_rows(
  lapply(
    contrasts,
    function(cc) {
      run_one(
        composition,
        cc[1],
        cc[2]
      )
    }
  )
)

# Family-wise multiplicity:
# each compartment is a biologically distinct L2 hypothesis family.
tests <- tests %>%
  group_by(
    contrast,
    final_compartment
  ) %>%
  mutate(
    FDR_BH_compartment = p.adjust(
      p_value,
      method = "BH"
    )
  ) %>%
  ungroup()

# Also retain global-across-L2 FDR for transparency/sensitivity.
tests <- tests %>%
  group_by(contrast) %>%
  mutate(
    FDR_BH_global = p.adjust(
      p_value,
      method = "BH"
    )
  ) %>%
  ungroup()

write_csv(
  tests,
  file.path(
    OUT,
    "17B2_L2_pairwise_DA.csv"
  )
)


# ============================================================
# 6. Ranked effects
# ============================================================

cat("\nSTEP 6/8 — ranked effects\n")

ranked <- tests %>%
  mutate(
    direction = case_when(
      delta_percentage_points > 0 ~
        paste0("higher_in_", numerator),

      delta_percentage_points < 0 ~
        paste0("higher_in_", denominator),

      TRUE ~ "no_change"
    ),

    abs_delta_pp =
      abs(delta_percentage_points)
  ) %>%
  arrange(
    contrast,
    FDR_BH_compartment,
    FDR_BH_global,
    desc(abs_delta_pp)
  )

write_csv(
  ranked,
  file.path(
    OUT,
    "17B2_L2_ranked_effects.csv"
  )
)

for (cc in c(
  "AA_vs_HC",
  "FFA_vs_HC",
  "AA_vs_FFA"
)) {

  cat(
    "\n============================================================\n"
  )
  cat(cc, "\n")
  cat(
    "============================================================\n"
  )

  show <- ranked %>%
    filter(contrast == cc) %>%
    select(
      final_compartment,
      final_L2_label,
      positive_numerator,
      positive_denominator,
      median_percent_numerator,
      median_percent_denominator,
      delta_percentage_points,
      log2_median_ratio,
      p_value,
      FDR_BH_compartment,
      FDR_BH_global,
      direction
    ) %>%
    slice_head(n = 30)

  print(
    as.data.frame(show),
    row.names = FALSE
  )
}


# ============================================================
# 7. Significant / nominal signal inventory
# ============================================================

cat("\nSTEP 7/8 — signal inventory\n")

signal_summary <- tests %>%
  group_by(contrast) %>%
  summarise(
    n_tested = n(),
    n_nominal_p_lt_005 =
      sum(p_value < 0.05, na.rm = TRUE),

    n_FDR_compartment_lt_005 =
      sum(
        FDR_BH_compartment < 0.05,
        na.rm = TRUE
      ),

    n_FDR_global_lt_005 =
      sum(
        FDR_BH_global < 0.05,
        na.rm = TRUE
      ),

    .groups = "drop"
  )

print(
  as.data.frame(signal_summary),
  row.names = FALSE
)

write_csv(
  signal_summary,
  file.path(
    OUT,
    "17B2_signal_summary.csv"
  )
)


# ============================================================
# 8. Formal checks
# ============================================================

cat("\nSTEP 8/8 — formal checks\n")

sum_check <- composition %>%
  group_by(
    patient_id,
    final_compartment
  ) %>%
  summarise(
    sum_prop = sum(proportion),
    .groups = "drop"
  )

checks <- tibble(
  check = c(
    "patients_16",
    "HC_5",
    "FFA_5",
    "AA_6",
    "only_frozen_L2_used",
    "no_returned_unadjudicated_used",
    "testability_disease_blind",
    "patient_compartment_proportions_sum_1",
    "three_contrasts",
    "BH_within_compartment_generated",
    "global_BH_generated",
    "patient_level_inference",
    "no_annotation_modified",
    "figures_generated_NO"
  ),

  pass = c(
    n_distinct(z$patient_id) == 16,

    n_distinct(
      z$patient_id[z$disease == "HC"]
    ) == 5,

    n_distinct(
      z$patient_id[z$disease == "FFA"]
    ) == 5,

    n_distinct(
      z$patient_id[z$disease == "AA"]
    ) == 6,

    all(z$final_L2_status == "FROZEN"),

    !any(
      z$final_L2_status ==
        "UNADJUDICATED_AFTER_RETURN"
    ),

    TRUE,

    all(
      abs(sum_check$sum_prop - 1) <
        1e-10
    ),

    n_distinct(tests$contrast) == 3,

    all(
      !is.na(
        tests$FDR_BH_compartment
      )
    ),

    all(
      !is.na(
        tests$FDR_BH_global
      )
    ),

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
    "17B2_formal_checks.csv"
  )
)

stopifnot(all(checks$pass))

cat("\n============================================================\n")
cat("17B2 COMPLETE\n")
cat("============================================================\n")
cat("Level: frozen L2 composition\n")
cat("Denominator: frozen-L2 cells within patient × compartment\n")
cat("Testability: >=5 cells in >=3 patients, disease-blind\n")
cat("Inferential unit: patient\n")
cat("Primary: AA vs HC; FFA vs HC\n")
cat("Secondary: AA vs FFA\n")
cat("BH FDR: contrast × compartment\n")
cat("Global-across-L2 BH retained as sensitivity metric\n")
cat("Returned unadjudicated cells excluded from L2 inference\n")
cat("Figures generated: NO\n")
cat("Annotation modified: NO\n")
cat("PASS\n")
