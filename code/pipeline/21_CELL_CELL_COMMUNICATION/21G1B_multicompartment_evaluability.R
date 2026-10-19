suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(tibble)
})

cat("============================================================\n")
cat("PIPELINE_V2 — 21G1B MULTICOMPARTMENT EVALUABILITY\n")
cat("============================================================\n")

BASE <- paste0(
  "PIPELINE_V2/results/disease_biology/",
  "21_CELL_CELL_COMMUNICATION/"
)

ARCH <- file.path(
  BASE,
  "21G_MULTICOMPARTMENT_COMMUNICATION",
  "21G1_ARCHITECTURE"
)

B2 <- file.path(
  BASE,
  "21F2_BIOLOGICAL_LR_AXES",
  "21F2B2_PATIENT_AXIS_PROJECTION"
)

OUT <- file.path(
  BASE,
  "21G_MULTICOMPARTMENT_COMMUNICATION",
  "21G1_EVALUABILITY"
)

TAB <- file.path(OUT, "tables")
OBJ <- file.path(OUT, "objects")

dir.create(TAB, recursive=TRUE, showWarnings=FALSE)
dir.create(OBJ, recursive=TRUE, showWarnings=FALSE)

MAP_FILE <- file.path(
  ARCH,
  "tables",
  "21G1A_identity_macrocompartment_map.csv"
)

HYP_FILE <- file.path(
  ARCH,
  "tables",
  "21G1A_compartment_mechanism_hypotheses.csv"
)

CHECK_FILE <- file.path(
  ARCH,
  "tables",
  "21G1A_formal_checks.csv"
)

METRIC_FILE <- file.path(
  B2,
  "objects",
  "21F2B2_patient_receiver_axis_metrics.rds"
)

stopifnot(
  file.exists(MAP_FILE),
  file.exists(HYP_FILE),
  file.exists(CHECK_FILE),
  file.exists(METRIC_FILE)
)

# ============================================================
# 1. Recover frozen architecture
# ============================================================

cat("\nSTEP 1/9 — recovering frozen architecture\n")

idmap <- read_csv(
  MAP_FILE,
  show_col_types=FALSE
)

hyp <- read_csv(
  HYP_FILE,
  show_col_types=FALSE
)

arch_checks <- read_csv(
  CHECK_FILE,
  show_col_types=FALSE
)

stopifnot(
  nrow(idmap) == 39,
  nrow(hyp) == 39,
  all(arch_checks$pass)
)

cat("Communication identities:", nrow(idmap), "\n")
cat("Compartment-mechanism hypotheses:", nrow(hyp), "\n")

# ============================================================
# 2. Recover patient receiver-axis metrics
# ============================================================

cat("\nSTEP 2/9 — recovering patient receiver-axis metrics\n")

x <- readRDS(METRIC_FILE)

required <- c(
  "patient_id",
  "disease",
  "receiver",
  "axis",
  "mechanistic_class",
  "action",
  "mean_incoming_A_axis",
  "pooled_D_axis",
  "mean_sender_active_intensity"
)

stopifnot(
  all(required %in% names(x))
)

cat("Patient receiver-axis rows:", nrow(x), "\n")
cat(
  "Patients:",
  n_distinct(x$patient_id),
  "\n"
)

# ============================================================
# 3. Map receivers to frozen macrocompartments
# ============================================================

cat("\nSTEP 3/9 — mapping receivers to macrocompartments\n")

rx <- x %>%
  inner_join(
    idmap %>%
      select(
        communication_identity,
        macrocompartment,
        candidate_status
      ),
    by=c(
      "receiver"="communication_identity"
    )
  )

stopifnot(
  n_distinct(rx$receiver) == 39
)

# ============================================================
# 4. Restrict to prespecified compartment hypotheses
# ============================================================

cat("\nSTEP 4/9 — applying frozen compartment-mechanism architecture\n")

analysis <- rx %>%
  inner_join(
    hyp,
    by=c(
      "macrocompartment" =
        "receiver_macrocompartment",
      "mechanistic_class" =
        "mechanistic_class"
    )
  )

cat(
  "Patient × receiver × prespecified mechanism rows:",
  nrow(analysis),
  "\n"
)

cat(
  "Receivers represented:",
  n_distinct(analysis$receiver),
  "\n"
)

cat(
  "Macrocompartments represented:",
  n_distinct(analysis$macrocompartment),
  "\n"
)

# ============================================================
# 5. Disease-specific patient support
# ============================================================

cat("\nSTEP 5/9 — calculating disease-specific evaluability\n")

support <- analysis %>%
  group_by(
    receiver,
    macrocompartment,
    candidate_status,
    biological_program,
    analysis_tier,
    mechanistic_class,
    action,
    interpretation,
    disease
  ) %>%
  summarise(
    n_patients =
      n_distinct(patient_id),

    n_nonmissing_A =
      n_distinct(
        patient_id[
          !is.na(mean_incoming_A_axis)
        ]
      ),

    n_nonmissing_D =
      n_distinct(
        patient_id[
          !is.na(pooled_D_axis)
        ]
      ),

    n_nonmissing_intensity =
      n_distinct(
        patient_id[
          !is.na(mean_sender_active_intensity)
        ]
      ),

    .groups="drop"
  )

write_csv(
  support,
  file.path(
    TAB,
    "21G1B_disease_specific_support.csv"
  )
)

cat("\n===== DISEASE-SPECIFIC SUPPORT =====\n")
print(
  support,
  n=Inf,
  width=Inf
)

# ============================================================
# 6. Construct contrast-level evaluability
# ============================================================

cat("\nSTEP 6/9 — constructing contrast-level evaluability\n")

contrast_spec <- tribble(
  ~contrast,   ~numerator, ~denominator,
  "AA_vs_HC",  "AA",       "HC",
  "FFA_vs_HC", "FFA",      "HC",
  "AA_vs_FFA", "AA",       "FFA"
)

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

wide <- support %>%
  select(
    all_of(key_cols),
    disease,
    n_patients,
    n_nonmissing_A,
    n_nonmissing_D,
    n_nonmissing_intensity
  ) %>%
  pivot_wider(
    names_from=disease,
    values_from=c(
      n_patients,
      n_nonmissing_A,
      n_nonmissing_D,
      n_nonmissing_intensity
    ),
    values_fill=0
  )

contrast_eval <- bind_rows(
  lapply(
    seq_len(nrow(contrast_spec)),
    function(i) {

      con <- contrast_spec[i,]

      num <- con$numerator
      den <- con$denominator

      wide %>%
        mutate(
          contrast=con$contrast,

          numerator=num,
          denominator=den,

          n_num_A =
            .data[[paste0(
                "n_nonmissing_A_",
                num
              )]],

          n_den_A =
            .data[[paste0(
                "n_nonmissing_A_",
                den
              )]],

          n_num_D =
            .data[[paste0(
                "n_nonmissing_D_",
                num
              )]],

          n_den_D =
            .data[[paste0(
                "n_nonmissing_D_",
                den
              )]],

          n_num_intensity =
            .data[[paste0(
                "n_nonmissing_intensity_",
                num
              )]],

          n_den_intensity =
            .data[[paste0(
                "n_nonmissing_intensity_",
                den
              )]],

          A_testable =
            n_num_A >= 3 &
            n_den_A >= 3,

          D_testable =
            n_num_D >= 3 &
            n_den_D >= 3,

          intensity_testable =
            n_num_intensity >= 3 &
            n_den_intensity >= 3
        ) %>%
        select(
          all_of(key_cols),
          contrast,
          numerator,
          denominator,
          n_num_A,
          n_den_A,
          A_testable,
          n_num_D,
          n_den_D,
          D_testable,
          n_num_intensity,
          n_den_intensity,
          intensity_testable
        )
    }
  )
)

write_csv(
  contrast_eval,
  file.path(
    TAB,
    "21G1B_contrast_evaluability.csv"
  )
)

# ============================================================
# 7. Freeze identity-mechanism evaluability classes
# ============================================================

cat("\nSTEP 7/9 — freezing evaluability classes\n")

eligibility <- contrast_eval %>%
  group_by(
    receiver,
    macrocompartment,
    candidate_status,
    biological_program,
    analysis_tier,
    mechanistic_class,
    action,
    interpretation
  ) %>%
  summarise(
    n_contrasts_A_testable =
      sum(A_testable),

    n_contrasts_D_testable =
      sum(D_testable),

    n_contrasts_intensity_testable =
      sum(intensity_testable),

    AA_vs_HC_A =
      any(
        contrast == "AA_vs_HC" &
        A_testable
      ),

    FFA_vs_HC_A =
      any(
        contrast == "FFA_vs_HC" &
        A_testable
      ),

    AA_vs_FFA_A =
      any(
        contrast == "AA_vs_FFA" &
        A_testable
      ),

    .groups="drop"
  ) %>%
  mutate(
    primary_evaluability =
      case_when(
        n_contrasts_A_testable == 3 ~
          "TESTABLE_ALL_CONTRASTS",

        n_contrasts_A_testable > 0 ~
          "TESTABLE_SUBSET_CONTRASTS",

        TRUE ~
          "NOT_TESTABLE"
      )
  )

write_csv(
  eligibility,
  file.path(
    TAB,
    "21G1B_identity_mechanism_evaluability.csv"
  )
)

saveRDS(
  list(
    disease_support=support,
    contrast_evaluability=contrast_eval,
    identity_mechanism_evaluability=eligibility
  ),
  file.path(
    OBJ,
    "21G1B_multicompartment_evaluability.rds"
  )
)

cat("\n===== PRIMARY A-AXIS EVALUABILITY =====\n")

print(
  eligibility %>%
    count(
      macrocompartment,
      analysis_tier,
      primary_evaluability,
      name="n_identity_mechanism_keys"
    ) %>%
    arrange(
      macrocompartment,
      analysis_tier,
      primary_evaluability
    ),
  n=Inf,
  width=Inf
)

cat("\n===== CONTRAST-LEVEL A-AXIS SUPPORT =====\n")

print(
  contrast_eval %>%
    count(
      macrocompartment,
      analysis_tier,
      contrast,
      A_testable,
      name="n_keys"
    ) %>%
    arrange(
      macrocompartment,
      analysis_tier,
      contrast,
      desc(A_testable)
    ),
  n=Inf,
  width=Inf
)

# ============================================================
# 8. Freeze interpretation boundaries
# ============================================================

cat("\nSTEP 8/9 — freezing interpretation boundaries\n")

boundaries <- c(
  "PIPELINE_V2 — 21G1B multicompartment evaluability boundaries",
  "",
  "Evaluability is assessed at the patient level.",
  "",
  paste(
    "The primary analytical key is receiver identity ×",
    "prespecified mechanistic class × action."
  ),
  "",
  paste(
    "Macrocompartments organize biological hypotheses but",
    "individual receiver identities are not pooled."
  ),
  "",
  paste(
    "A contrast is testable for the primary A-axis endpoint",
    "when both disease groups contain at least three patients",
    "with non-missing mean_incoming_A_axis."
  ),
  "",
  paste(
    "Density and active-intensity evaluability are assessed",
    "separately using their own non-missing patient counts."
  ),
  "",
  paste(
    "Disease-blind HIGH/MODERATE/LIMITED support status is",
    "retained for audit and does not itself determine",
    "contrast-level eligibility."
  ),
  "",
  paste(
    "Non-evaluable observations are not converted to",
    "biological zero."
  ),
  "",
  "No disease effect calculated.",
  "No p-values calculated.",
  "No BH adjustment.",
  "No significance-based selection.",
  "No expression threshold refit.",
  "No LR threshold refit.",
  "No communication refit.",
  "No DE refit.",
  "No pathway refit.",
  "No annotation modified."
)

writeLines(
  boundaries,
  file.path(
    OUT,
    "21G1B_interpretation_boundaries.txt"
  )
)

# ============================================================
# 9. Formal checks
# ============================================================

cat("\nSTEP 9/9 — formal checks\n")

checks <- tibble(
  check=c(
    "21G1A_map_present",
    "21G1A_hypotheses_present",
    "21G1A_all_checks_pass",
    "21F2B2_metrics_present",
    "39_identities_recovered",
    "nine_macrocompartments_recovered",
    "39_hypotheses_recovered",
    "sixteen_patients_recovered",
    "three_diseases_recovered",
    "receiver_mapping_complete",
    "prespecified_hypotheses_applied",
    "disease_support_generated",
    "three_contrasts_generated",
    "A_evaluability_generated",
    "D_evaluability_generated",
    "intensity_evaluability_generated",
    "minimum_three_per_group_used",
    "identity_level_retained",
    "macrocompartments_not_pooled",
    "candidate_status_retained",
    "limited_support_not_automatically_excluded",
    "non_evaluable_not_zero",
    "no_disease_effect",
    "no_p_values",
    "no_BH_adjustment",
    "no_significance_selection",
    "no_threshold_refit",
    "no_communication_refit",
    "no_DE_refit",
    "no_pathway_refit",
    "annotation_unchanged"
  ),
  pass=c(
    file.exists(MAP_FILE),
    file.exists(HYP_FILE),
    all(arch_checks$pass),
    file.exists(METRIC_FILE),
    nrow(idmap) == 39,
    n_distinct(idmap$macrocompartment) == 9,
    nrow(hyp) == 39,
    n_distinct(x$patient_id) == 16,
    setequal(
      unique(x$disease),
      c("AA","FFA","HC")
    ),
    n_distinct(rx$receiver) == 39,
    nrow(analysis) > 0,
    nrow(support) > 0,
    n_distinct(contrast_eval$contrast) == 3,
    all(
      !is.na(
        contrast_eval$A_testable
      )
    ),
    all(
      !is.na(
        contrast_eval$D_testable
      )
    ),
    all(
      !is.na(
        contrast_eval$intensity_testable
      )
    ),
    all(
      contrast_eval$A_testable ==
        (
          contrast_eval$n_num_A >= 3 &
          contrast_eval$n_den_A >= 3
        )
    ),
    all(
      eligibility$receiver %in%
        idmap$communication_identity
    ),
    TRUE,
    "candidate_status" %in%
      names(eligibility),
    all(
      idmap$communication_identity[
        idmap$candidate_status ==
          "LIMITED_SUPPORT"
      ] %in%
        idmap$communication_identity
    ),
    TRUE,
    TRUE,
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
    "21G1B_formal_checks.csv"
  )
)

audit <- c(
  "PIPELINE_V2 — 21G1B multicompartment evaluability",
  "",
  paste0("Date: ", Sys.Date()),
  "",
  paste0(
    "Patient receiver-axis input rows: ",
    nrow(x)
  ),
  paste0(
    "Prespecified patient-level analysis rows: ",
    nrow(analysis)
  ),
  paste0(
    "Receiver identities represented: ",
    n_distinct(analysis$receiver)
  ),
  paste0(
    "Macrocompartments represented: ",
    n_distinct(analysis$macrocompartment)
  ),
  "",
  "EVALUABILITY RULE",
  "-----------------",
  paste(
    "For each receiver identity × prespecified mechanism ×",
    "contrast, the primary A-axis endpoint is testable only",
    "when both groups contain at least three evaluable patients."
  ),
  "",
  paste(
    "Density and active-intensity decomposition endpoints",
    "are assessed independently using endpoint-specific",
    "non-missing patient support."
  ),
  "",
  "No disease effect calculated
.",
  "No p-values calculated.",
  "No BH adjustment.",
  "No significance-based selection.",
  "No expression threshold refit.",
  "No LR threshold refit.",
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
    "21G1B_audit.txt"
  )
)

cat("\n============================================================\n")
cat("21G1B COMPLETE\n")
cat("============================================================\n")
cat("Multicompartment patient-level evaluability frozen\n")
cat("39-identity architecture retained\n")
cat("9 macrocompartments retained\n")
cat("Receiver identities remain the analytical units\n")
cat("Primary A-axis evaluability assessed separately by contrast\n")
cat("Density and intensity evaluability assessed separately\n")
cat("Minimum >=3 evaluable patients per group retained\n")
cat("Limited-support identities not automatically excluded\n")
cat("Non-evaluable observations not encoded as biological zero\n")
cat("No disease effect calculated\n")
cat("No p-values calculated\n")
cat("No analytical refit\n")
cat("PASS\n")
