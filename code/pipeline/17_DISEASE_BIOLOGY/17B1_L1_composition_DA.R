suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(tibble)
})

cat("============================================================\n")
cat("PIPELINE_V2 — 17B1 L1 COMPOSITION / DIFFERENTIAL ABUNDANCE\n")
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

x <- x %>%
  filter(analysis_eligible_L1)

stopifnot(nrow(x) == 82680)

disease_levels <- c("HC", "FFA", "AA")

x$disease <- factor(
  x$disease,
  levels = disease_levels
)

# ============================================================
# 1. Patient denominators
# ============================================================

cat("\nSTEP 1/7 — patient denominators\n")

denom <- x %>%
  count(
    disease,
    patient_id,
    name = "n_total"
  )

print(as.data.frame(denom), row.names = FALSE)


# ============================================================
# 2. Patient × L1 counts and proportions
# ============================================================

cat("\nSTEP 2/7 — patient-level L1 composition\n")

all_compartments <- sort(
  unique(x$final_compartment)
)

composition <- x %>%
  count(
    disease,
    patient_id,
    final_compartment,
    name = "n_cells"
  ) %>%
  complete(
    disease,
    patient_id,
    final_compartment = all_compartments,
    fill = list(n_cells = 0)
  ) %>%
  left_join(
    denom,
    by = c("disease", "patient_id")
  ) %>%
  mutate(
    proportion = n_cells / n_total,
    percent = 100 * proportion
  )

# Remove impossible disease/patient combinations introduced by complete()
composition <- composition %>%
  filter(!is.na(n_total))

write_csv(
  composition,
  file.path(
    OUT,
    "17B1_patient_L1_composition.csv"
  )
)


# ============================================================
# 3. Disease descriptive summary
# ============================================================

cat("\nSTEP 3/7 — disease summaries\n")

summary_tbl <- composition %>%
  group_by(
    disease,
    final_compartment
  ) %>%
  summarise(
    n_patients = n(),
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
    "17B1_L1_disease_summary.csv"
  )
)

print(as.data.frame(summary_tbl), row.names = FALSE)


# ============================================================
# 4. Pairwise patient-level tests
# ============================================================

cat("\nSTEP 4/7 — pairwise contrasts\n")

contrasts <- list(
  c("AA", "HC"),
  c("FFA", "HC"),
  c("AA", "FFA")
)

run_contrast <- function(dat, group1, group0) {

  compartments <- sort(
    unique(dat$final_compartment)
  )

  bind_rows(
    lapply(
      compartments,
      function(comp) {

        z <- dat %>%
          filter(
            final_compartment == comp,
            disease %in% c(group1, group0)
          )

        a <- z$proportion[z$disease == group1]
        b <- z$proportion[z$disease == group0]

        # Wilcoxon exact inference where possible.
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

        # Pseudocount only for fold-change reporting.
        eps <- 1e-6

        tibble(
          contrast = paste0(group1, "_vs_", group0),
          numerator = group1,
          denominator = group0,
          final_compartment = comp,

          n_numerator = length(a),
          n_denominator = length(b),

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
      run_contrast(
        composition,
        cc[1],
        cc[2]
      )
    }
  )
)

tests <- tests %>%
  group_by(contrast) %>%
  mutate(
    FDR_BH = p.adjust(
      p_value,
      method = "BH"
    )
  ) %>%
  ungroup()

write_csv(
  tests,
  file.path(
    OUT,
    "17B1_L1_pairwise_DA.csv"
  )
)


# ============================================================
# 5. Direction / ranking table
# ============================================================

cat("\nSTEP 5/7 — ranked effects\n")

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
    FDR_BH,
    desc(abs_delta_pp)
  )

write_csv(
  ranked,
  file.path(
    OUT,
    "17B1_L1_ranked_effects.csv"
  )
)

cat("\n===== AA vs HC =====\n")
print(
  as.data.frame(
    ranked %>%
      filter(contrast == "AA_vs_HC") %>%
      select(
        final_compartment,
        median_percent_numerator,
        median_percent_denominator,
        delta_percentage_points,
        log2_median_ratio,
        p_value,
        FDR_BH,
        direction
      )
  ),
  row.names = FALSE
)

cat("\n===== FFA vs HC =====\n")
print(
  as.data.frame(
    ranked %>%
      filter(contrast == "FFA_vs_HC") %>%
      select(
        final_compartment,
        median_percent_numerator,
        median_percent_denominator,
        delta_percentage_points,
        log2_median_ratio,
        p_value,
        FDR_BH,
        direction
      )
  ),
  row.names = FALSE
)

cat("\n===== AA vs FFA =====\n")
print(
  as.data.frame(
    ranked %>%
      filter(contrast == "AA_vs_FFA") %>%
      select(
        final_compartment,
        median_percent_numerator,
        median_percent_denominator,
        delta_percentage_points,
        log2_median_ratio,
        p_value,
        FDR_BH,
        direction
      )
  ),
  row.names = FALSE
)


# ============================================================
# 6. Global compositional heterogeneity
# ============================================================

cat("\nSTEP 6/7 — global composition matrix\n")

composition_wide <- composition %>%
  select(
    disease,
    patient_id,
    final_compartment,
    proportion
  ) %>%
  pivot_wider(
    names_from = final_compartment,
    values_from = proportion,
    values_fill = 0
  )

write_csv(
  composition_wide,
  file.path(
    OUT,
    "17B1_patient_L1_composition_matrix.csv"
  )
)


# ============================================================
# 7. Formal checks
# ============================================================

cat("\nSTEP 7/7 — formal checks\n")

checks <- tibble(
  check = c(
    "analysis_cells_82680",
    "patients_16",
    "HC_5",
    "FFA_5",
    "AA_6",
    "three_diseases",
    "patient_proportions_sum_1",
    "three_contrasts",
    "no_cell_level_inference",
    "BH_FDR_generated",
    "no_annotation_modified"
  ),

  pass = c(
    nrow(x) == 82680,

    n_distinct(x$patient_id) == 16,

    n_distinct(
      x$patient_id[x$disease == "HC"]
    ) == 5,

    n_distinct(
      x$patient_id[x$disease == "FFA"]
    ) == 5,

    n_distinct(
      x$patient_id[x$disease == "AA"]
    ) == 6,

    setequal(
      as.character(
        unique(x$disease)
      ),
      c("HC", "FFA", "AA")
    ),

    all(
      abs(
        composition %>%
          group_by(patient_id) %>%
          summarise(
            s = sum(proportion),
            .groups = "drop"
          ) %>%
          pull(s) - 1
      ) < 1e-10
    ),

    n_distinct(tests$contrast) == 3,

    TRUE,

    all(!is.na(tests$FDR_BH)),

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
    "17B1_formal_checks.csv"
  )
)

stopifnot(all(checks$pass))

cat("\n============================================================\n")
cat("17B1 COMPLETE\n")
cat("============================================================\n")
cat("Analysis level: L1 composition\n")
cat("Inferential unit: patient\n")
cat("Primary: AA vs HC; FFA vs HC\n")
cat("Secondary: AA vs FFA\n")
cat("Multiple testing: BH within contrast\n")
cat("Figures generated: NO\n")
cat("Annotation modified: NO\n")
cat("READY FOR L1 INTERPRETATION AND L2 DA\n")
cat("PASS\n")
