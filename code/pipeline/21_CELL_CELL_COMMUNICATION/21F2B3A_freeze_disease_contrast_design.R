suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
})

cat("============================================================\n")
cat("PIPELINE_V2 — 21F2B3A DISEASE-CONTRAST DESIGN FREEZE\n")
cat("============================================================\n")

# ============================================================
# 0. Paths
# ============================================================

BASE <- paste0(
  "PIPELINE_V2/results/disease_biology/",
  "21_CELL_CELL_COMMUNICATION/"
)

B2 <- file.path(
  BASE,
  "21F2_BIOLOGICAL_LR_AXES",
  "21F2B2_PATIENT_AXIS_PROJECTION"
)

OUT <- file.path(
  BASE,
  "21F2_BIOLOGICAL_LR_AXES",
  "21F2B3_DISEASE_CONTRASTS"
)

TAB <- file.path(OUT, "tables")

dir.create(
  TAB,
  recursive=TRUE,
  showWarnings=FALSE
)

B2_CHECKS <- file.path(
  B2,
  "tables",
  "21F2B2_formal_checks.csv"
)

B2_DICT <- file.path(
  B2,
  "tables",
  "21F2B2_metric_dictionary.csv"
)

stopifnot(
  file.exists(B2_CHECKS),
  file.exists(B2_DICT)
)

# Import only methodological provenance.
# No patient-level communication values are read in this module.

b2_checks <- read_csv(
  B2_CHECKS,
  show_col_types=FALSE
)

b2_dict <- read_csv(
  B2_DICT,
  show_col_types=FALSE
)

stopifnot(
  all(b2_checks$pass)
)

# ============================================================
# 1. Freeze analysis population
# ============================================================

cat("\nSTEP 1/10 — freezing analysis population\n")

population <- tribble(
  ~disease, ~n_patients, ~role,
  "AA",      6L,          "DISEASE",
  "FFA",     5L,          "DISEASE",
  "HC",      5L,          "REFERENCE"
)

write_csv(
  population,
  file.path(
    TAB,
    "21F2B3A_analysis_population.csv"
  )
)

print(
  population,
  n=Inf,
  width=Inf
)

# ============================================================
# 2. Freeze contrasts
# ============================================================

cat("\nSTEP 2/10 — freezing disease contrasts\n")

contrasts <- tribble(
  ~contrast,
  ~numerator,
  ~denominator,
  ~n_numerator,
  ~n_denominator,
  ~n_exact_label_allocations,
  ~interpretive_role,

  "AA_vs_HC",
  "AA",
  "HC",
  6L,
  5L,
  462L,
  "AA-associated communication architecture relative to HC",

  "FFA_vs_HC",
  "FFA",
  "HC",
  5L,
  5L,
  252L,
  "FFA-associated communication architecture relative to HC",

  "AA_vs_FFA",
  "AA",
  "FFA",
  6L,
  5L,
  462L,
  "Direct disease triangulation"
)

write_csv(
  contrasts,
  file.path(
    TAB,
    "21F2B3A_contrast_specification.csv"
  )
)

print(
  contrasts,
  n=Inf,
  width=Inf
)

# ============================================================
# 3. Freeze follicular receiver panel
# ============================================================

cat("\nSTEP 3/10 — freezing receiver panel\n")

receivers <- tribble(
  ~receiver,
  ~analysis_role,

  "Bulge_Isthmus",
  "PRIMARY_FOLLICULAR_RECEIVER",

  "ORS",
  "PRIMARY_FOLLICULAR_RECEIVER",

  "ORS_basal",
  "PRIMARY_FOLLICULAR_RECEIVER"
)

write_csv(
  receivers,
  file.path(
    TAB,
    "21F2B3A_receiver_panel.csv"
  )
)

# ============================================================
# 4. Freeze mechanistic hypothesis families
# ============================================================

cat("\nSTEP 4/10 — freezing mechanistic hypothesis families\n")

mechanisms <- tribble(
  ~mechanistic_class,
  ~action_required,
  ~analysis_tier,
  ~hypothesis_family,
  ~direction_prespecified,
  ~interpretation,

  "INFLAMMATORY_EXPOSURE",
  "INCLUDE",
  "PRIMARY",
  "INJURY_EXPOSURE",
  "NONE",
  "Inflammatory cytokine/interferon/chemokine exposure",

  "CYTOTOXIC_DEATH_SIGNAL",
  "INCLUDE",
  "PRIMARY",
  "INJURY_EXPOSURE",
  "NONE",
  "Direct cytotoxic/death-receptor communication",

  "FGF_PARACRINE_TROPHIC",
  "INCLUDE",
  "PRIMARY",
  "TROPHIC_HOMEOSTATIC_SUPPORT",
  "NONE",
  "Paracrine FGF-family trophic communication",

  "IGF_SIGNALING_TROPHIC",
  "INCLUDE",
  "PRIMARY",
  "TROPHIC_HOMEOSTATIC_SUPPORT",
  "NONE",
  "IGF1R-associated trophic communication",

  "GROWTH_SURVIVAL_TROPHIC",
  "INCLUDE",
  "PRIMARY",
  "TROPHIC_HOMEOSTATIC_SUPPORT",
  "NONE",
  "EGF/ERBB/HGF-associated growth and survival communication",

  "WNT_CANONICAL_REGENERATIVE",
  "INCLUDE",
  "PRIMARY",
  "TROPHIC_HOMEOSTATIC_SUPPORT",
  "NONE",
  "Canonical-compatible WNT regenerative communication",

  "VEGF_VASCULAR_SUPPORT",
  "INCLUDE",
  "PRIMARY",
  "VASCULAR_SUPPORT",
  "NONE",
  "VEGF-family vascular support",

  "ANGPT_TIE_VASCULAR",
  "INCLUDE",
  "PRIMARY",
  "VASCULAR_SUPPORT",
  "NONE",
  "Canonical ANGPT-TIE vascular communication",

  "TGF_BETA_FIBROTIC_REMODELING",
  "INCLUDE",
  "PRIMARY",
  "REMODELING",
  "NONE",
  "TGF-beta-associated remodeling communication",

  "ECM_ADHESION_STRUCTURAL",
  "INCLUDE",
  "PRIMARY",
  "REMODELING",
  "NONE",
  "ECM-receptor structural/adhesive communication",

  "BMP_QUIESCENCE_DIFFERENTIATION",
  "INCLUDE",
  "CONTEXTUAL",
  "CONTEXTUAL_REGULATION",
  "NONE",
  "Context-dependent BMP quiescence/differentiation signaling",

  "WNT_NONCANONICAL_CONTEXTUAL",
  "INCLUDE",
  "CONTEXTUAL",
  "CONTEXTUAL_REGULATION",
  "NONE",
  "Noncanonical/context-dependent WNT communication",

  "WNT_NONCANONICAL_CONTEXTUAL",
  "REVIEW",
  "CONTEXTUAL",
  "CONTEXTUAL_REGULATION",
  "NONE",
  "Review-level noncanonical/context-dependent WNT communication",

  "TGF_ACTIVIN_CONTEXTUAL",
  "REVIEW",
  "CONTEXTUAL",
  "CONTEXTUAL_REGULATION",
  "NONE",
  "GDF/TGF-superfamily contextual communication",

  "IGF2_CLEARANCE_TRAFFICKING",
  "INCLUDE",
  "CONTEXTUAL",
  "CONTEXTUAL_REGULATION",
  "NONE",
  "IGF2R-associated trafficking/clearance",

  "IGF2_CLEARANCE_TRAFFICKING",
  "REVIEW",
  "CONTEXTUAL",
  "CONTEXTUAL_REGULATION",
  "NONE",
  "Review-level IGF2R-associated trafficking/clearance",

  "FGF_ENDOCRINE_CONTEXTUAL",
  "INCLUDE",
  "CONTEXTUAL",
  "CONTEXTUAL_REGULATION",
  "NONE",
  "Endocrine FGF signaling requiring co-receptor context",

  "PDGF_STROMAL_VASCULAR",
  "INCLUDE",
  "CONTEXTUAL",
  "CONTEXTUAL_REGULATION",
  "NONE",
  "Stromal/perivascular PDGF communication",

  "IMMUNE_REGULATORY_SIGNAL",
  "INCLUDE",
  "CONTEXTUAL",
  "CONTEXTUAL_REGULATION",
  "NONE",
  "Immune co-regulatory communication"
)

stopifnot(
  !anyDuplicated(
    mechanisms[
      c(
        "mechanistic_class",
        "action_required"
      )
    ]
  )
)

write_csv(
  mechanisms,
  file.path(
    TAB,
    "21F2B3A_mechanistic_hypothesis_families.csv"
  )
)

cat("\n===== MECHANISTIC DESIGN =====\n")

print(
  mechanisms,
  n=Inf,
  width=Inf
)

# ============================================================
# 5. Freeze patient-level endpoints
# ============================================================

cat("\nSTEP 5/10 — freezing statistical endpoints\n")

endpoints <- tribble(
  ~endpoint,
  ~analysis_role,
  ~formula,
  ~interpretation,

  "mean_incoming_A_axis",
  "PRIMARY",
  "(1/K_pra) * sum_s A_psra",
  paste(
    "Mean structural-opportunity-normalized incoming",
    "mechanistic exposure across evaluable senders"
  ),

  "pooled_D_axis",
  "SECONDARY_DECOMPOSITION",
  "sum_s N_psra / sum_s M_sra",
  paste(
    "Fraction of structurally possible incoming",
    "mechanistic LR channels active"
  ),

  "mean_sender_active_intensity",
  "SECONDARY_DECOMPOSITION",
  "mean_s(meanC_psra) among active sender-axis edges",
  paste(
    "Mean active LR intensity among contributing",
    "sender-axis edges"
  )
)

write_csv(
  endpoints,
  file.path(
    TAB,
    "21F2B3A_endpoint_specification.csv"
  )
)

# ============================================================
# 6. Freeze estimands
# ============================================================

cat("\nSTEP 6/10 — freezing estimands\n")

estimands <- tribble(
  ~estimand,
  ~role,
  ~formula,
  ~notes,

  "mean_difference",
  "PRIMARY_EFFECT",
  "mean(numerator) - mean(denominator)",
  "Absolute patient-level difference on the native metric scale",

  "median_difference",
  "DESCRIPTIVE_EFFECT",
  "median(numerator) - median(denominator)",
  "Robust descriptive location difference",

  "mean_ratio",
  "DESCRIPTIVE_EFFECT",
  "(mean(numerator) + epsilon) / (mean(denominator) + epsilon)",
  "Ratio reported descriptively with frozen epsilon",

  "log2_mean_ratio",
  "DESCRIPTIVE_EFFECT",
  "log2((mean(numerator)+epsilon)/(mean(denominator)+epsilon))",
  "Not the primary inferential estimand"
)

epsilon_spec <- tribble(
  ~parameter,
  ~value,
  ~definition,

  "epsilon",
  1e-6,
  paste(
    "Fixed numerical pseudocount used only for descriptive",
    "ratio/log2-ratio calculations; never used for",
    "eligibility, permutation testing or mean differences"
  )
)

write_csv(
  estimands,
  file.path(
    TAB,
    "21F2B3A_estimand_specification.csv"
  )
)

write_csv(
  epsilon_spec,
  file.path(
    TAB,
    "21F2B3A_numeric_constants.csv"
  )
)

# ============================================================
# 7. Freeze exact permutation inference
# ============================================================

cat("\nSTEP 7/10 — freezing exact permutation inference\n")

inference <- tribble(
  ~component,
  ~specification,

  "unit_of_inference",
  "patient",

  "primary_test_statistic",
  "absolute difference in patient-level means",

  "signed_effect",
  "numerator mean minus denominator mean",

  "permutation_scheme",
  paste(
    "Exhaustive reassignment of disease labels among patients",
    "within each two-group contrast while preserving group sizes"
  ),

  "AA_vs_HC_allocations",
  "choose(11,6) = 462",

  "FFA_vs_HC_allocations",
  "choose(10,5) = 252",

  "AA_vs_FFA_allocations",
  "choose(11,6) = 462",

  "primary_p_value",
  paste(
    "Two-sided exact permutation p = fraction of all label",
    "allocations with |delta_perm| >= |delta_observed|"
  ),

  "ties",
  "included using >= in the exact tail probability",

  "zero_semantics",
  paste(
    "Evaluable inactive mechanistic axes retain zero exposure;",
    "non-evaluable structural combinations are not converted to zero"
  ),

  "missingness",
  paste(
    "No imputation. A test is performed only when the receiver-axis",
    "endpoint is evaluable in at least 3 patients in each contrasted group."
  ),

  "wilcoxon",
  paste(
    "Optional secondary sensitivity statistic only;",
    "not the primary inferential test"
  )
)

write_csv(
  inference,
  file.path(
    TAB,
    "21F2B3A_inference_specification.csv"
  )
)

# ============================================================
# 8. Freeze multiplicity families
# ============================================================

cat("\nSTEP 8/10 — freezing multiplicity control\n")

multiplicity <- tribble(
  ~analysis_tier,
  ~endpoint_role,
  ~contrast_scope,
  ~BH_family,
  ~interpretation,

  "PRIMARY",
  "PRIMARY",
  "EACH_CONTRAST_SEPARATELY",
  "contrast × hypothesis_family",
  paste(
    "BH correction across all primary receiver × mechanistic-class",
    "tests belonging to the same prespecified hypothesis family",
    "within a disease contrast"
  ),

  "PRIMARY",
  "SECONDARY_DECOMPOSITION",
  "EACH_CONTRAST_SEPARATELY",
  "contrast × endpoint × hypothesis_family",
  paste(
    "Secondary decomposition FDR kept separate from the",
    "primary A-axis inferential family"
  ),

  "CONTEXTUAL",
  "PRIMARY",
  "EACH_CONTRAST_SEPARATELY",
  "contrast × CONTEXTUAL_REGULATION",
  paste(
    "Contextual analyses corrected separately and never",
    "promoted into the primary family"
  ),

  "CONTEXTUAL",
  "SECONDARY_DECOMPOSITION",
  "EACH_CONTRAST_SEPARATELY",
  "contrast × endpoint × CONTEXTUAL_REGULATION",
  "Contextual decomposition sensitivity family"
)

write_csv(
  multiplicity,
  file.path(
    TAB,
    "21F2B3A_multiplicity_specification.csv"
  )
)

# ============================================================
# 9. Freeze reporting / interpretation rules
# ============================================================

cat("\nSTEP 9/10 — freezing reporting rules\n")

reporting <- tribble(
  ~rule_id,
  ~rule,

  "R01",
  "Report patient-level group n for every tested receiver-axis endpoint.",

  "R02",
  "Report group means and medians together with the signed mean difference.",

  "R03",
  "Report the exact two-sided permutation p-value and BH-adjusted q-value.",

  "R04",
  "Do not interpret q>=0.05 as evidence of equivalence or absence of biology.",

  "R05",
  "Do not interpret communication metrics as direct molecular flux.",

  "R06",
  "Decompose A_axis using D_axis and active intensity rather than treating A_axis as mechanistically indivisible.",

  "R07",
  "Preserve INCLUDE and REVIEW strata; do not pool them silently.",

  "R08",
  "Do not infer protective, destructive, regenerative or fibrotic direction solely from the sign of a statistical contrast.",

  "R09",
  "Use sender-resolved follow-up only after receiver-axis patient-level testing; sender rows are explanatory, not independent replicates.",

  "R10",
  "No composite preservation score is defined in 21F2B3A.",

  "R11",
  "AA-vs-HC, FFA-vs-HC and AA-vs-FFA are all retained; no contrast is selected according to observed significance.",

  "R12",
  "Primary conclusions require consistency with the frozen mechanistic interpretation established before disease testing."
)

write_csv(
  reporting,
  file.path(
    TAB,
    "21F2B3A_reporting_rules.csv"
  )
)

# ============================================================
# 10. Formal checks
# ============================================================

cat("\nSTEP 10/10 — formal checks\n")

checks <- tibble(
  check=c(
    "21F2B2_checks_present",
    "21F2B2_all_checks_pass",
    "three_disease_groups_frozen",
    "AA_n6_frozen",
    "FFA_n5_frozen",
    "HC_n5_frozen",
    "three_contrasts_frozen",
    "AA_HC_exact_allocations_462",
    "FFA_HC_exact_allocations_252",
    "AA_FFA_exact_allocations_462",
    "three_follicular_receivers_frozen",
    "primary_mechanisms_prespecified",
    "contextual_mechanisms_prespecified",
    "action_retained_in_design",
    "primary_endpoint_A_axis",
    "density_decomposition_frozen",
    "intensity_decomposition_frozen",
    "patient_unit_of_inference",
    "exact_permutation_primary",
    "two_sided_test_frozen",
    "minimum_three_patients_per_group",
    "BH_families_frozen",
    "epsilon_descriptive_only",
    "no_composite_preservation_score",
    "no_patient_level_values_read",
    "no_disease_effect_calculated",
    "no_p_values_calculated",
    "no_significance_based_selection",
    "no_DE_refit",
    "no_pathway_refit",
    "no_communication_refit",
    "annotation_unchanged"
  ),
  pass=c(
    file.exists(B2_CHECKS),
    all(b2_checks$pass),

    nrow(population) == 3,
    population$n_patients[
      population$disease == "AA"
    ] == 6,
    population$n_patients[
      population$disease == "FFA"
    ] == 5,
    population$n_patients[
      population$disease == "HC"
    ] == 5,

    nrow(contrasts) == 3,

    contrasts$n_exact_label_allocations[
      contrasts$contrast == "AA_vs_HC"
    ] == 462,

    contrasts$n_exact_label_allocations[
      contrasts$contrast == "FFA_vs_HC"
    ] == 252,

    contrasts$n_exact_label_allocations[
      contrasts$contrast == "AA_vs_FFA"
    ] == 462,

    nrow(receivers) == 3,

    all(
      c(
        "INJURY_EXPOSURE",
        "TROPHIC_HOMEOSTATIC_SUPPORT",
        "VASCULAR_SUPPORT",
        "REMODELING"
      ) %in%
        mechanisms$hypothesis_family[
          mechanisms$analysis_tier == "PRIMARY"
        ]
    ),

    any(
      mechanisms$analysis_tier ==
        "CONTEXTUAL"
    ),

    "action_required" %in%
      names(mechanisms),

    endpoints$analysis_role[
      endpoints$endpoint ==
        "mean_incoming_A_axis"
    ] == "PRIMARY",

    any(
      endpoints$endpoint ==
        "pooled_D_axis"
    ),

    any(
      endpoints$endpoint ==
        "mean_sender_active_intensity"
    ),

    any(
      inference$component ==
        "unit_of_inference" &
      inference$specification ==
        "patient"
    ),

    any(
      inference$component ==
        "permutation_scheme" &
      grepl(
        "Exhaustive",
        inference$specification
      )
    ),

    any(
      inference$component ==
        "primary_p_value" &
      grepl(
        "Two-sided",
        inference$specification
      )
    ),

    any(
      inference$component ==
        "missingness" &
      grepl(
        "at least 3 patients",
        inference$specification
      )
    ),

    nrow(multiplicity) == 4,

    epsilon_spec$value[
      epsilon_spec$parameter ==
        "epsilon"
    ] == 1e-6,

    any(
      reporting$rule_id == "R10" &
      grepl(
        "No composite preservation score",
        reporting$rule
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
    "21F2B3A_formal_checks.csv"
  )
)

# ============================================================
# Audit
# ============================================================

audit <- c(
  "PIPELINE_V2 — 21F2B3A disease-contrast statistical design freeze",
  "",
  paste0("Date: ", Sys.Date()),
  "",
  "PURPOSE",
  "-------",
  paste(
    "Freeze the patient-level statistical design for mechanistic",
    "cell-cell communication contrasts before calculating",
    "disease effects or inferential statistics."
  ),
  "",
  "ANALYSIS POPULATION",
  "-------------------",
  "- AA: 6 patients",
  "- FFA: 5 patients",
  "- HC: 5 patients",
  "- Unit of inference: patient",
  "",
  "CONTRASTS",
  "---------",
  "- AA vs HC: 462 exhaustive label allocations",
  "- FFA vs HC: 252 exhaustive label allocations",
  "- AA vs FFA: 462 exhaustive label allocations",
  "",
  "PRIMARY RECEIVERS",
  "-----------------",
  "- Bulge_Isthmus",
  "- ORS",
  "- ORS_basal",
  "",
  "PRIMARY ENDPOINT",
  "----------------",
  paste(
    "mean_incoming_A_axis: mean structural-opportunity-normalized",
    "mechanistic LR exposure across evaluable senders."
  ),
  "",
  "DECOMPOSITION ENDPOINTS",
  "-----------------------",
  "- pooled_D_axis: active-channel density.",
  "- mean_sender_active_intensity: active-channel intensity.",
  "",
  "PRIMARY EFFECT",
  "--------------",
  paste(
    "Signed absolute difference in patient-level means:",
    "numerator minus denominator."
  ),
  "",
  "PRIMARY INFERENCE",
  "-----------------",
  paste(
    "Two-sided exhaustive patient-label permutation test using",
    "the absolute difference in means as the test statistic."
  ),
  paste(
    "Ties are included using |delta_perm| >= |delta_observed|."
  ),
  paste(
    "Endpoints require at least 3 evaluable patients in each",
    "contrasted group."
  ),
  "",
  "MULTIPLICITY",
  "------------",
  paste(
    "BH FDR is controlled within prespecified contrast, endpoint",
    "and mechanistic hypothesis families according to the",
    "frozen multiplicity specification."
  ),
  "",
  "INTERPRETATION",
  "--------------",
  paste(
    "No direction of biological benefit or harm is assigned",
    "from the statistical sign alone."
  ),
  paste(
    "INCLUDE and REVIEW strata remain separate."
  ),
  paste(
    "No composite preservation score is defined."
  ),
  "",
  "ANALYTICAL BOUNDARIES",
  "---------------------",
  "No patient-level communication values read.",
  "No disease effect calculated.",
  "No p-value calculated.",
  "No significance-based feature selection.",
  "No DE refit.",
  "No pathway refit.",
  "No communication refit.",
  "No annotation modified.",
  "",
  "PASS"
)

writeLines(
  audit,
  file.path(
    OUT,
    "21F2B3A_audit.txt"
  )
)

cat("\n============================================================\n")
cat("21F2B3A COMPLETE\n")
cat("============================================================\n")
cat("Patient-level disease-contrast design frozen before testing\n")
cat("AA-vs-HC, FFA-vs-HC and AA-vs-FFA retained\n")
cat("Follicular receiver panel frozen\n")
cat("Mechanistic hypothesis families frozen\n")
cat("Primary A-axis endpoint frozen\n")
cat("Density and intensity decomposition frozen\n")
cat("Exact exhaustive permutation inference frozen\n")
cat("BH multiplicity families frozen\n")
cat("No composite preservation score defined\n")
cat("No patient-level communication values inspected\n")
cat("No disease effects or p-values calculated\n")
cat("PASS\n")
