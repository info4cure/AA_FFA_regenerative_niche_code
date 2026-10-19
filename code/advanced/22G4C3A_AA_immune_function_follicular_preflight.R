
# ============================================================
# 22G4C3A — AA IMMUNE FUNCTION × FOLLICULAR TOPOLOGY PREFLIGHT
#
# PURPOSE
# Build the exact AA-only patient-level universe linking
# pre-specified immune functional programmes to frozen
# follicular / mesenchymal / upper comparator structures.
#
# IMPORTANT
# - AA ONLY
# - patient is the statistical unit
# - no correlations
# - no regression
# - no significance testing
# - no feature selection from observed associations
# - no re-estimation of immune programmes
# - no re-estimation of structural identities
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(stringr)
})

cat("\n============================================================\n")
cat("22G4C3A — AA IMMUNE FUNCTION × FOLLICULAR TOPOLOGY PREFLIGHT\n")
cat("============================================================\n")

ROOT <- "PIPELINE_V2"

BASE <- file.path(
  ROOT,
  "results/disease_biology",
  "22_MULTICOMPARTMENT_DISEASE_AXES",
  "22G_BULGE_CENTERED_ECOSYSTEM",
  "22G4_IMMUNE_ECOSYSTEM",
  "22G4C_AA_HC_LOWER_FOLLICLE"
)

OUT <- file.path(
  BASE,
  "22G4C3_IMMUNE_FUNCTION_FOLLICULAR_TOPOLOGY",
  "22G4C3A_PREFLIGHT"
)

TAB <- file.path(
  OUT,
  "tables"
)

OBJ <- file.path(
  OUT,
  "objects"
)

AUD <- file.path(
  OUT,
  "audit"
)

dir.create(
  TAB,
  recursive=TRUE,
  showWarnings=FALSE
)

dir.create(
  OBJ,
  recursive=TRUE,
  showWarnings=FALSE
)

dir.create(
  AUD,
  recursive=TRUE,
  showWarnings=FALSE
)

# ============================================================
# 1. IMMUNE FUNCTION SOURCE
# ============================================================

IMMUNE_FILE <- file.path(
  BASE,
  "22G4C2_IMMUNE_FUNCTION",
  "tables",
  "22G4C2B_targeted_immune_functional_scores.csv"
)

stopifnot(
  file.exists(
    IMMUNE_FILE
  )
)

immune_all <- read.csv(
  IMMUNE_FILE,
  stringsAsFactors=FALSE,
  check.names=FALSE
) |>
  as_tibble()

cat("\n===== IMMUNE FUNCTION SOURCE =====\n")

cat(
  "Rows:",
  nrow(immune_all),
  "\n"
)

cat(
  "Patients:",
  n_distinct(
    immune_all$patient_id
  ),
  "\n"
)

# ============================================================
# 2. PRE-SPECIFIED IMMUNE PROGRAMMES
# ============================================================

immune_hypotheses <- tribble(
  ~immune_L2,
  ~programme,
  ~immune_axis,
  ~priority,

  "CD8 T",
  "CYTOTOXICITY",
  "CD8_CYTOTOXICITY",
  "PRIMARY",

  "CD8 T",
  "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "CD8_IFNG",
  "PRIMARY",

  "CD8 T",
  "HALLMARK_INTERFERON_ALPHA_RESPONSE",
  "CD8_IFNA",
  "PRIMARY",

  "NK",
  "CYTOTOXICITY",
  "NK_CYTOTOXICITY",
  "PRIMARY",

  "NK",
  "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "NK_IFNG",
  "PRIMARY",

  "Macrophage/monocyte",
  "HALLMARK_INTERFERON_ALPHA_RESPONSE",
  "MAC_IFNA",
  "SECONDARY",

  "Macrophage/monocyte",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "MAC_TNFA_NFKB",
  "SECONDARY"
)

cat("\n===== PRE-SPECIFIED IMMUNE AXES =====\n")

print(
  immune_hypotheses,
  n=Inf,
  width=Inf
)

# ============================================================
# 3. AA IMMUNE FUNCTIONAL MASTER
# ============================================================

immune_AA <- immune_all |>
  filter(
    disease=="AA"
  ) |>
  inner_join(
    immune_hypotheses,
    by=c(
      "final_L2_label"="immune_L2",
      "programme"="programme"
    )
  ) |>
  dplyr::select(
    patient_id,
    disease,
    batch,
    immune_L2=final_L2_label,
    immune_axis,
    priority,
    programme,
    immune_score=score,
    immune_n_cells=n_cells,
    immune_support_tier=support_tier,
    analysis_population
  ) |>
  arrange(
    immune_axis,
    patient_id
  )

cat("\n===== AA IMMUNE FUNCTIONAL SUPPORT =====\n")

immune_support <- immune_AA |>
  group_by(
    immune_L2,
    immune_axis,
    priority
  ) |>
  summarise(
    n_patients=
      n_distinct(
        patient_id
      ),

    patients=
      paste(
        sort(
          unique(
            patient_id
          )
        ),
        collapse=","
      ),

    total_cells=
      sum(
        immune_n_cells
      ),

    min_cells=
      min(
        immune_n_cells
      ),

    median_cells=
      median(
        immune_n_cells
      ),

    max_cells=
      max(
        immune_n_cells
      ),

    n_ge5=
      sum(
        immune_n_cells >= 5
      ),

    n_ge10=
      sum(
        immune_n_cells >= 10
      ),

    n_ge20=
      sum(
        immune_n_cells >= 20
      ),

    .groups="drop"
  )

print(
  immune_support,
  n=Inf,
  width=Inf
)

# ============================================================
# 4. LOCATE 22G4C1 STRUCTURAL OUTPUTS
# ============================================================

cat("\n===== SEARCH 22G4C1 STRUCTURAL SOURCES =====\n")

candidate_files <- list.files(
  BASE,
  pattern="22G4C1.*\\.csv$",
  recursive=TRUE,
  full.names=TRUE
)

cat(
  "Candidate 22G4C1 CSV files:",
  length(candidate_files),
  "\n"
)

if(
  length(candidate_files)>0
) {
  cat(
    paste(
      candidate_files,
      collapse="\n"
    ),
    "\n"
  )
}

# ============================================================
# 5. USE ORIGINAL FROZEN L2 COMPOSITION AS STRUCTURAL SOURCE
# ============================================================

# This is the same frozen patient-level compositional source
# used in the previous topology modules and avoids deriving
# structural values from model output.

STRUCT_FILE <- file.path(
  ROOT,
  "results/disease_biology",
  "20_CELLULAR_COMPOSITION",
  "20C_PATIENT_LEVEL_COMPOSITION",
  "tables",
  "20C_L2_patient_frequencies_within_L1.csv"
)

stopifnot(
  file.exists(
    STRUCT_FILE
  )
)

struct_all <- read.csv(
  STRUCT_FILE,
  stringsAsFactors=FALSE,
  check.names=FALSE
) |>
  as_tibble()

# ============================================================
# 6. PRE-SPECIFIED STRUCTURAL TARGETS
# ============================================================

structural_def <- tribble(
  ~anatomical_axis,
  ~target_L1,
  ~target_L2,
  ~target_variable,
  ~role,

  "LOWER_FOLLICLE",
  "Epithelial",
  "ORS",
  "ORS",
  "PRIMARY",

  "LOWER_FOLLICLE",
  "Epithelial",
  "ORS_basal",
  "ORS_basal",
  "PRIMARY",

  "LOWER_FOLLICLE",
  "Epithelial",
  "IRS_companion",
  "IRS_companion",
  "PRIMARY",

  "LOWER_FOLLICLE",
  "Epithelial",
  "Hair_shaft",
  "Hair_shaft",
  "PRIMARY",

  "LOWER_MESENCHYME",
  "Fibroblast",
  "COCH/AQP4 dermal papilla",
  "COCH_AQP4_dermal_papilla",
  "PRIMARY",

  "LOWER_MESENCHYME",
  "Fibroblast",
  "COL11A1/DPEP1 dermal-sheath",
  "COL11A1_DPEP1_dermal_sheath",
  "PRIMARY",

  "UPPER_COMPARATOR",
  "Epithelial",
  "IFE_Infundibulum_differentiating",
  "IFE_differentiating",
  "COMPARATOR",

  "UPPER_COMPARATOR",
  "Epithelial",
  "IFE_Infundibulum_suprabasal",
  "IFE_suprabasal",
  "COMPARATOR",

  "UPPER_COMPARATOR",
  "Epithelial",
  "IFE_Infundibulum_terminal",
  "IFE_terminal",
  "COMPARATOR"
)

cat("\n===== PRE-SPECIFIED STRUCTURAL TARGETS =====\n")

print(
  structural_def,
  n=Inf,
  width=Inf
)

# ============================================================
# 7. AA STRUCTURAL MASTER
# ============================================================

struct_AA <- struct_all |>
  filter(
    disease=="AA"
  ) |>
  inner_join(
    structural_def,
    by=c(
      "L1"="target_L1",
      "L2"="target_L2"
    )
  ) |>
  transmute(
    patient_id,
    disease,
    anatomical_axis,
    target_L1=L1,
    target_L2=L2,
    target_variable,
    role,
    structural_frequency=frequency,
    structural_percent=percent,
    structural_n_cells=n,
    structural_denominator=
      denominator_L1
  ) |>
  arrange(
    anatomical_axis,
    target_variable,
    patient_id
  )

cat("\n===== AA STRUCTURAL SUPPORT =====\n")

struct_support <- struct_AA |>
  group_by(
    anatomical_axis,
    target_L1,
    target_L2,
    target_variable,
    role
  ) |>
  summarise(
    n_patients=
      n_distinct(
        patient_id
      ),

    patients=
      paste(
        sort(
          unique(
            patient_id
          )
        ),
        collapse=","
      ),

    total_cells=
      sum(
        structural_n_cells,
        na.rm=TRUE
      ),

    zero_patients=
      sum(
        structural_n_cells==0,
        na.rm=TRUE
      ),

    median_frequency=
      median(
        structural_frequency,
        na.rm=TRUE
      ),

    .groups="drop"
  )

print(
  struct_support,
  n=Inf,
  width=Inf
)

# ============================================================
# 8. BUILD AA IMMUNE × STRUCTURAL PAIR UNIVERSE
# ============================================================

cat("\n===== BUILD AA IMMUNE × STRUCTURAL UNIVERSE =====\n")

pair_master <- immune_AA |>
  inner_join(
    struct_AA,
    by=c(
      "patient_id",
      "disease"
    ),
    relationship="many-to-many"
  ) |>
  arrange(
    immune_axis,
    anatomical_axis,
    target_variable,
    patient_id
  )

cat(
  "Patient-level pair rows:",
  nrow(pair_master),
  "\n"
)

# ============================================================
# 9. CONNECTION-LEVEL COVERAGE
# ============================================================

connection_support <- pair_master |>
  group_by(
    immune_L2,
    immune_axis,
    priority,
    anatomical_axis,
    target_L1,
    target_L2,
    target_variable,
    role
  ) |>
  summarise(
    n_common=
      n_distinct(
        patient_id
      ),

    common_patients=
      paste(
        sort(
          unique(
            patient_id
          )
        ),
        collapse=","
      ),

    immune_total_cells=
      sum(
        immune_n_cells
      ),

    immune_min_cells=
      min(
        immune_n_cells
      ),

    immune_median_cells=
      median(
        immune_n_cells
      ),

    immune_n_ge5=
      sum(
        immune_n_cells >= 5
      ),

    immune_n_ge10=
      sum(
        immune_n_cells >= 10
      ),

    structural_total_cells=
      sum(
        structural_n_cells
      ),

    structural_zero_patients=
      sum(
        structural_n_cells==0
      ),

    immune_unique_scores=
      n_distinct(
        immune_score
      ),

    structural_unique_values=
      n_distinct(
        structural_frequency
      ),

    .groups="drop"
  ) |>
  mutate(
    evaluable_ALL=
      n_common >= 4 &
      immune_unique_scores >= 3 &
      structural_unique_values >= 3,

    evaluable_GE5=
      immune_n_ge5 >= 4,

    evaluable_GE10=
      immune_n_ge10 >= 4,

    analysis_tier=case_when(

      immune_axis %in%c(
        "CD8_CYTOTOXICITY",
        "CD8_IFNG",
        "CD8_IFNA"
      ) &
        evaluable_ALL
        ~ "PRIMARY_CD8",

      immune_axis %in%c(
        "NK_CYTOTOXICITY",
        "NK_IFNG"
      ) &
        evaluable_ALL
        ~ "PRIMARY_NK_LOW_SUPPORT",

      priority=="SECONDARY" &
        evaluable_ALL
        ~ "SECONDARY",

      TRUE
        ~ "NOT_EVALUABLE"
    )
  )

cat("\n===== CONNECTION SUPPORT =====\n")

print(
  connection_support,
  n=Inf,
  width=Inf
)

# ============================================================
# 10. EXPECTED CONNECTION COUNTS
# ============================================================

cat("\n===== CONNECTION COUNTS =====\n")

connection_counts <- connection_support |>
  count(
    immune_axis,
    priority,
    analysis_tier,
    name="n_connections"
  ) |>
  arrange(
    priority,
    immune_axis,
    analysis_tier
  )

print(
  connection_counts,
  n=Inf,
  width=Inf
)

cat("\n===== CONNECTIONS BY ANATOMICAL AXIS =====\n")

anatomical_counts <- connection_support |>
  count(
    anatomical_axis,
    analysis_tier,
    name="n_connections"
  ) |>
  arrange(
    anatomical_axis,
    analysis_tier
  )

print(
  anatomical_counts,
  n=Inf,
  width=Inf
)

# ============================================================
# 11. AA PATIENT COVERAGE MATRIX
# ============================================================

cat("\n===== AA PATIENT COVERAGE MATRIX =====\n")

immune_patient_matrix <- immune_AA |>
  distinct(
    patient_id,
    immune_axis,
    immune_n_cells
  ) |>
  mutate(
    value=
      paste0(
        "n=",
        immune_n_cells
      )
  ) |>
  dplyr::select(
    patient_id,
    immune_axis,
    value
  ) |>
  pivot_wider(
    names_from=immune_axis,
    values_from=value,
    values_fill="ABSENT"
  ) |>
  arrange(
    patient_id
  )

print(
  immune_patient_matrix,
  n=Inf,
  width=Inf
)

# ============================================================
# 12. SAVE TABLES
# ============================================================

write.csv(
  immune_hypotheses,
  file.path(
    TAB,
    "22G4C3A_prespecified_immune_axes.csv"
  ),
  row.names=FALSE
)

write.csv(
  structural_def,
  file.path(
    TAB,
    "22G4C3A_prespecified_structural_targets.csv"
  ),
  row.names=FALSE
)

write.csv(
  immune_AA,
  file.path(
    TAB,
    "22G4C3A_AA_immune_functional_master.csv"
  ),
  row.names=FALSE
)

write.csv(
  struct_AA,
  file.path(
    TAB,
    "22G4C3A_AA_structural_master.csv"
  ),
  row.names=FALSE
)

write.csv(
  pair_master,
  file.path(
    TAB,
    "22G4C3A_AA_patient_pair_master.csv"
  ),
  row.names=FALSE
)

write.csv(
  connection_support,
  file.path(
    TAB,
    "22G4C3A_connection_support.csv"
  ),
  row.names=FALSE
)

write.csv(
  immune_support,
  file.path(
    TAB,
    "22G4C3A_immune_support.csv"
  ),
  row.names=FALSE
)

write.csv(
  struct_support,
  file.path(
    TAB,
    "22G4C3A_structural_support.csv"
  ),
  row.names=FALSE
)

write.csv(
  immune_patient_matrix,
  file.path(
    TAB,
    "22G4C3A_AA_patient_immune_coverage_matrix.csv"
  ),
  row.names=FALSE
)

# ============================================================
# 13. FREEZE
# ============================================================

freeze <- list(
  module="22G4C3A",

  purpose=paste(
    "AA-only preflight for
 immune functional programmes and frozen follicular",
    "structural states prior to patient-level coupling analysis."
  ),

  disease_universe="AA only",

  statistical_unit="patient",

  immune_source=
    IMMUNE_FILE,

  structural_source=
    STRUCT_FILE,

  immune_hypotheses=
    immune_hypotheses,

  structural_targets=
    structural_def,

  immune_support=
    immune_support,

  structural_support=
    struct_support,

  connection_support=
    connection_support,

  patient_pair_master=
    pair_master,

  interpretation_constraints=c(
    "No correlation testing performed.",
    "No regression testing performed.",
    "No significance-based feature selection performed.",
    "Immune functional programmes were pre-specified before coupling analysis.",
    "Structural targets were pre-specified before coupling analysis.",
    "Patient is the statistical unit.",
    "AA-only analysis avoids the AA-versus-HC disease/batch confounding for within-AA coupling.",
    "NK analyses retain low-cell-support warnings.",
    "Association in subsequent modules will not establish causal direction."
  )
)

saveRDS(
  freeze,
  file.path(
    OBJ,
    "22G4C3A_AA_immune_function_follicular_preflight_freeze.rds"
  )
)

# ============================================================
# 14. FORMAL CHECKS
# ============================================================

cat("\n===== FORMAL CHECKS =====\n")

expected_immune_axes <- c(
  "CD8_CYTOTOXICITY",
  "CD8_IFNG",
  "CD8_IFNA",
  "NK_CYTOTOXICITY",
  "NK_IFNG",
  "MAC_IFNA",
  "MAC_TNFA_NFKB"
)

expected_structural_targets <- c(
  "ORS",
  "ORS_basal",
  "IRS_companion",
  "Hair_shaft",
  "COCH_AQP4_dermal_papilla",
  "COL11A1_DPEP1_dermal_sheath",
  "IFE_differentiating",
  "IFE_suprabasal",
  "IFE_terminal"
)

expected_connections <-
  length(
    expected_immune_axes
  ) *
  length(
    expected_structural_targets
  )

checks <- tibble(
  check=c(
    "AA-only immune master",
    "AA-only structural master",
    "Seven pre-specified immune axes",
    "Three CD8 functional axes",
    "Two NK functional axes",
    "Two macrophage functional axes",
    "Nine structural targets",
    "Four lower-follicle targets",
    "Two lower-mesenchymal targets",
    "Three upper comparators",
    "63 pre-specified immune-structure connections",
    "CD8 axes evaluable",
    "NK axes retained",
    "Macrophage axes evaluable",
    "No inferential testing performed",
    "Freeze exists"
  ),

  pass=c(
    nrow(immune_AA)>0 &&
      all(
        immune_AA$disease=="AA"
      ),

    nrow(struct_AA)>0 &&
      all(
        struct_AA$disease=="AA"
      ),

    setequal(
      unique(
        immune_hypotheses$immune_axis
      ),
      expected_immune_axes
    ),

    sum(
      grepl(
        "^CD8_",
        expected_immune_axes
      )
    )==3,

    sum(
      grepl(
        "^NK_",
        expected_immune_axes
      )
    )==2,

    sum(
      grepl(
        "^MAC_",
        expected_immune_axes
      )
    )==2,

    setequal(
      unique(
        structural_def$target_variable
      ),
      expected_structural_targets
    ),

    sum(
      structural_def$anatomical_axis==
        "LOWER_FOLLICLE"
    )==4,

    sum(
      structural_def$anatomical_axis==
        "LOWER_MESENCHYME"
    )==2,

    sum(
      structural_def$anatomical_axis==
        "UPPER_COMPARATOR"
    )==3,

    nrow(
      connection_support
    )==
      expected_connections,

    all(
      connection_support$evaluable_ALL[
        grepl(
          "^CD8_",
          connection_support$immune_axis
        )
      ]
    ),

    all(
      c(
        "NK_CYTOTOXICITY",
        "NK_IFNG"
      ) %in%
        connection_support$immune_axis
    ),

    all(
      connection_support$evaluable_ALL[
        grepl(
          "^MAC_",
          connection_support$immune_axis
        )
      ]
    ),

    TRUE,

    file.exists(
      file.path(
        OBJ,
        "22G4C3A_AA_immune_function_follicular_preflight_freeze.rds"
      )
    )
  )
)

print(
  checks,
  n=Inf,
  width=Inf
)

write.csv(
  checks,
  file.path(
    AUD,
    "22G4C3A_formal_checks.csv"
  ),
  row.names=FALSE
)

# ============================================================
# 15. ANALYSIS-TIER SUMMARY
# ============================================================

cat("\n===== ANALYSIS-TIER SUMMARY =====\n")

tier_summary <- connection_support |>
  count(
    analysis_tier,
    name="n_connections"
  ) |>
  arrange(
    desc(
      n_connections
    )
  )

print(
  tier_summary,
  n=Inf,
  width=Inf
)

write.csv(
  tier_summary,
  file.path(
    TAB,
    "22G4C3A_analysis_tier_summary.csv"
  ),
  row.names=FALSE
)

# ============================================================
# 16. AUDIT
# ============================================================

audit_lines <- c(
  "============================================================",
  "22G4C3A — AA IMMUNE FUNCTION × FOLLICULAR TOPOLOGY PREFLIGHT",
  "============================================================",
  "",
  "Disease universe: AA only",
  "Statistical unit: patient",
  "",
  paste0(
    "Pre-specified immune functional axes: ",
    n_distinct(
      immune_hypotheses$immune_axis
    )
  ),
  paste0(
    "Pre-specified structural targets: ",
    n_distinct(
      structural_def$target_variable
    )
  ),
  paste0(
    "Pre-specified immune-structure connections: ",
    nrow(
      connection_support
    )
  ),
  "",
  "PRIMARY CD8 axes:",
  "- CD8 cytotoxicity",
  "- CD8 IFN-gamma response",
  "- CD8 IFN-alpha response",
  "",
  "PRIMARY NK axes:",
  "- NK cytotoxicity",
  "- NK IFN-gamma response",
  "",
  "SECONDARY macrophage axes:",
  "- macrophage IFN-alpha response",
  "- macrophage TNFA/NFKB signaling",
  "",
  "Primary lower-follicular structures:",
  "- ORS",
  "- ORS basal",
  "- IRS/companion",
  "- Hair shaft",
  "",
  "Primary lower-mesenchymal structures:",
  "- COCH/AQP4 dermal papilla",
  "- COL11A1/DPEP1 dermal sheath",
  "",
  "Upper anatomical comparators:",
  "- IFE/Infundibulum differentiating",
  "- IFE/Infundibulum suprabasal",
  "- IFE/Infundibulum terminal",
  "",
  paste(
    "Within-AA coupling is evaluated separately from the",
    "AA-versus-HC functional-state comparison."
  ),
  paste(
    "Because all AA samples belong to CABIMER2026, within-AA",
    "patient coupling is not confounded by the cross-disease",
    "historical-versus-CABIMER batch structure."
  ),
  "",
  paste(
    "NK is retained as a biologically prioritized low-support",
    "population; downstream NK results require explicit",
    "cell-support sensitivity analysis."
  ),
  "",
  "No correlations were calculated.",
  "No regression models were fitted.",
  "No P values or FDR values were generated.",
  "No features were selected from observed associations.",
  "",
  paste0(
    "FORMAL CHECKS: ",
    sum(
      checks$pass
    ),
    "/",
    nrow(
      checks
    ),
    " PASS"
  ),
  "",
  if(
    all(
      checks$pass
    )
  ) {
    "STATUS=PASS"
  } else {
    "STATUS=FAIL"
  }
)

writeLines(
  audit_lines,
  file.path(
    AUD,
    "22G4C3A_audit.txt"
  )
)

cat(
  paste(
    audit_lines,
    collapse="\n"
  ),
  "\n"
)

# ============================================================
# 17. FINAL STATUS
# ============================================================

if(
  !all(
    checks$pass
  )
) {
  stop(
    "22G4C3A formal checks failed."
  )
}

cat("\n============================================================\n")
cat("22G4C3A COMPLETE\n")
cat("STATUS=PASS\n")
cat("============================================================\n")

