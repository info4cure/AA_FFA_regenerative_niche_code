
# ============================================================
# 22G4C — AA vs HC IMMUNE–LOWER-FOLLICLE PREFLIGHT
#
# Purpose:
# Test an AA-specific lower-follicular inflammatory topology
# without FFA.
#
# Patients:
#   AA = 6
#   HC = 5
#
# Candidate architecture:
#
# Immune
#   |
#   +-- Bulge (comparator)
#   |
#   +-- ORS / ORS_basal
#   +-- IRS_companion
#   +-- Hair_shaft
#   |
#   +-- COCH/AQP4 dermal papilla
#   +-- COL11A1/DPEP1 dermal sheath
#
# Upper IFE states are retained as anatomical comparators.
#
# NO inferential modelling is performed here.
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
})

cat("\n============================================================\n")
cat("22G4C — AA vs HC IMMUNE–LOWER-FOLLICLE PREFLIGHT\n")
cat("============================================================\n")

ROOT <- "PIPELINE_V2"

# ============================================================
# 1. INPUTS
# ============================================================

L2_FILE <- file.path(
  ROOT,
  "results/disease_biology",
  "20_CELLULAR_COMPOSITION",
  "20C_PATIENT_LEVEL_COMPOSITION",
  "tables",
  "20C_L2_patient_frequencies_within_L1.csv"
)

BULGE_FILE <- file.path(
  ROOT,
  "results/disease_biology",
  "22_MULTICOMPARTMENT_DISEASE_AXES",
  "22F_BULGE_STEM_NICHE_INTEGRITY",
  "22F2_FUNCTIONAL_STATE_DENSITY",
  "22F2A_SCORE_FREEZE",
  "tables",
  "22F2A_Bulge_patient_functional_scores.csv"
)

CFILE <- file.path(
  ROOT,
  "results/disease_biology",
  "22_MULTICOMPARTMENT_DISEASE_AXES",
  "22F_BULGE_STEM_NICHE_INTEGRITY",
  "22F4_BULGE_MONOCLE3_TRAJECTORIES",
  "22F4C_TRUNK_PSEUDOTIME",
  "objects",
  "22F4C_trunk_definition_freeze.rds"
)

OUT <- file.path(
  ROOT,
  "results/disease_biology",
  "22_MULTICOMPARTMENT_DISEASE_AXES",
  "22G_BULGE_CENTERED_ECOSYSTEM",
  "22G4_IMMUNE_ECOSYSTEM",
  "22G4C_AA_HC_LOWER_FOLLICLE"
)

TAB <- file.path(OUT,"tables")
AUD <- file.path(OUT,"audit")

stopifnot(
  file.exists(L2_FILE),
  file.exists(BULGE_FILE),
  file.exists(CFILE)
)

# ============================================================
# 2. AA + HC ONLY
# ============================================================

l2 <- read.csv(
  L2_FILE,
  stringsAsFactors=FALSE
) |>
  as_tibble() |>
  filter(
    disease %in% c(
      "HC",
      "AA"
    )
  )

cat("\n===== PATIENTS =====\n")

patient_audit <- l2 |>
  distinct(
    patient_id,
    disease
  ) |>
  count(
    disease,
    name="n_patients"
  )

print(
  patient_audit,
  n=Inf
)

stopifnot(
  n_distinct(
    l2$patient_id
  ) == 11
)

# ============================================================
# 3. IMMUNE UNIVERSE — REAUDIT IN AA + HC
# ============================================================

IMMUNE_L1 <- c(
  "T_NK",
  "Myeloid_DC"
)

immune <- l2 |>
  filter(
    L1 %in%
      IMMUNE_L1
  )

immune_support <- immune |>
  group_by(
    L1,
    L2
  ) |>
  summarise(
    n_patients=
      n_distinct(
        patient_id
      ),

    total_cells=
      sum(
        n,
        na.rm=TRUE
      ),

    positive_patients=
      sum(
        n > 0,
        na.rm=TRUE
      ),

    positive_HC=
      sum(
        disease=="HC" &
          n > 0,
        na.rm=TRUE
      ),

    positive_AA=
      sum(
        disease=="AA" &
          n > 0,
        na.rm=TRUE
      ),

    median_n=
      median(
        n,
        na.rm=TRUE
      ),

    .groups="drop"
  ) |>
  mutate(
    feature=paste(
      L1,
      L2,
      sep=" :: "
    ),

    # AA/HC-specific eligibility.
    # Slightly adapted from pan-disease rule because n=11.
    eligible_AA_HC=
      total_cells >= 30 &
      positive_patients >= 8 &
      positive_HC >= 3 &
      positive_AA >= 3,

    exclusion_reason=case_when(
      eligible_AA_HC
        ~ "INCLUDE",

      total_cells < 30
        ~ "<30 AA+HC cells",

      positive_patients < 8
        ~ "<8 positive AA+HC patients",

      positive_HC < 3
        ~ "<3 positive HC",

      positive_AA < 3
        ~ "<3 positive AA",

      TRUE
        ~ "OTHER"
    )
  ) |>
  arrange(
    L1,
    L2
  )

cat("\n===== AA+HC IMMUNE SUPPORT =====\n")

print(
  immune_support,
  n=Inf
)

eligible_immune <- immune_support |>
  filter(
    eligible_AA_HC
  )

cat("\n===== ELIGIBLE AA+HC IMMUNE STATES =====\n")

print(
  eligible_immune,
  n=Inf
)

# ============================================================
# 4. LOWER-FOLLICLE TARGET PANEL
# ============================================================

target_panel <- tribble(
  ~anatomical_axis, ~L1, ~L2, ~variable, ~role,

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

# ============================================================
# 5. TARGET SUPPORT
# ============================================================

target_source <- target_panel |>
  inner_join(
    l2,
    by=c(
      "L1",
      "L2"
    ),
    relationship="many-to-many"
  )

target_support <- target_source |>
  group_by(
    anatomical_axis,
    L1,
    L2,
    variable,
    role
  ) |>
  summarise(
    n_patients=
      n_distinct(
        patient_id
      ),

    total_cells=
      sum(
        n,
        na.rm=TRUE
      ),

    positive_patients=
      sum(
        n > 0,
        na.rm=TRUE
      ),

    positive_HC=
      sum(
        disease=="HC" &
          n > 0,
        na.rm=TRUE
      ),

    positive_AA=
      sum(
        disease=="AA" &
          n > 0,
        na.rm=TRUE
      ),

    median_n=
      median(
        n,
        na.rm=TRUE
      ),

    .groups="drop"
  ) |>
  mutate(
    evaluable=
      n_patients >= 10 &
      positive_HC >= 2 &
      positive_AA >= 3
  )

cat("\n===== LOWER-FOLLICLE / COMPARATOR SUPPORT =====\n")

print(
  target_support,
  n=Inf
)

# ============================================================
# 6. BULGE COMPARATOR
# ============================================================

C <- readRDS(
  CFILE
)

bulge_position <- as_tibble(
  C$patient_position
) |>
  filter(
    disease %in%
      c(
        "HC",
        "AA"
      )
  ) |>
  dplyr::select(
    patient_id,
    disease,
    median_position
  )

bulge_scores <- read.csv(
  BULGE_FILE,
  stringsAsFactors=FALSE
) |>
  as_tibble() |>
  filter(
    disease %in%
      c(
        "HC",
        "AA"
      )
  ) |>
  dplyr::select(
    patient_id,
    disease,
    stem_niche_mean,
    MHCII_mean,
    cycling_mean,
    CD200_z_mean
  )

bulge <- bulge_position |>
  left_join(
    bulge_scores,
    by=c(
      "patient_id",
      "disease"
    )
  )

cat("\n===== BULGE AA+HC SUPPORT =====\n")

print(
  bulge |>
    count(
      disease,
      name="n_patients"
    ),
  n=Inf
)

# ============================================================
# 7. DISEASE DISTRIBUTION — DESCRIPTIVE PREFLIGHT
# ============================================================

target_disease <- target_source |>
  group_by(
    anatomical_axis,
    variable,
    disease
  ) |>
  summarise(
    median_frequency=
      median(
        frequency,
        na.rm=TRUE
      ),

    mean_frequency=
      mean(
        frequency,
        na.rm=TRUE
      ),

    total_cells=
      sum(
        n,
        na.rm=TRUE
      ),

    .groups="drop"
  )

cat("\n===== TARGET DISEASE DISTRIBUTION =====\n")

print(
  target_disease,
  n=Inf
)

# ============================================================
# 8. SAVE
# ============================================================

write.csv(
  immune_support,
  file.path(
    TAB,
    "22G4C_AA_HC_immune_support.csv"
  ),
  row.names=FALSE
)

write.csv(
  eligible_immune,
  file.path(
    TAB,
    "22G4C_AA_HC_eligible_immune_states.csv"
  ),
  row.names=FALSE
)

write.csv(
  target_panel,
  file.path(
    TAB,
    "22G4C_lower_follicle_target_panel.csv"
  ),
  row.names=FALSE
)

write.csv(
  target_support,
  file.path(
    TAB,
    "22G4C_lower_follicle_target_support.csv"
  ),
  row.names=FALSE
)

write.csv(
  target_source,
  file.path(
    TAB,
    "22G4C_lower_follicle_patient_source.csv"
  ),
  row.names=FALSE
)

write.csv(
  bulge,
  file.path(
    TAB,
    "22G4C_AA_HC_Bulge_comparator.csv"
  ),
  row.names=FALSE
)

write.csv(
  target_disease,
  file.path(
    TAB,
    "22G4C_target_disease_distribution.csv"
  ),
  row.names=FALSE
)

# ============================================================
# 9. FORMAL CHECKS
# ============================================================

checks <- tibble(
  check=c(
    "11 AA+HC patients",
    "6 AA patients",
    "5 HC patients",
    "Immune universe non-empty",
    "Eligible AA+HC immune states > 0",
    "ORS present",
    "ORS_basal present",
    "IRS_companion present",
    "Hair_shaft present",
    "COCH/AQP4 dermal papilla present",
    "COL11A1/DPEP1 dermal sheath present",
    "Three upper IFE comparators present",
    "Bulge comparator has 11 patients"
  ),

  pass=c(
    n_distinct(
      l2$patient_id
    )==11,

    n_distinct(
      l2$patient_id[
        l2$disease=="AA"
      ]
    )==6,

    n_distinct(
      l2$patient_id[
        l2$disease=="HC"
      ]
    )==5,

    nrow(
      immune_support
    )>0,

    nrow(
      eligible_immune
    )>0,

    "ORS" %in%
      target_support$variable,

    "ORS_basal" %in%
      target_support$variable,

    "IRS_companion" %in%
      target_support$variable,

    "Hair_shaft" %in%
      target_support$variable,

    "COCH_AQP4_dermal_papilla" %in%
      target_support$variable,

    "COL11A1_DPEP1_dermal_sheath" %in%
      target_support$variable,

    sum(
      target_support$anatomical_axis==
        "UPPER_COMPARATOR"
    )==3,

    nrow(
      bulge
    )==11
  )
)

cat("\n===== FORMAL CHECKS =====\n")

print(
  checks,
  n=Inf
)

write.csv(
  checks,
  file.path(
    AUD,
    "22G4C_preflight_checks.csv"
  ),
  row.names=FALSE
)

# ============================================================
# 10. AUDIT
# ============================================================

audit_lines <- c(
  "============================================================",
  "22G4C — AA vs HC IMMUNE–LOWER-FOLLICLE PREFLIGHT",
  "============================================================",
  "",
  "Disease universe: AA + HC only; FFA excluded.",
  "AA patients: 6",
  "HC patients: 5",
  "",
  paste0(
    "Eligible immune states: ",
    nrow(
      eligible_immune
    )
  ),
  "",
  "Primary lower-follicular targets:",
  "- ORS",
  "- ORS_basal",
  "- IRS_companion",
  "- Hair_shaft",
  "- COCH/AQP4 dermal papilla",
  "- COL11A1/DPEP1 dermal sheath",
  "",
  "Upper-follicular comparators:",
  "- IFE/Infundibulum differentiating",
  "- IFE/Infundibulum suprabasal",
  "- IFE/Infundibulum terminal",
  "",
  "Bulge retained as comparator rather than assumed AA epicentre.",
  "",
  paste0(
    "FORMAL CHECKS: ",
    sum(checks$pass),
    "/",
    nrow(checks),
    " PASS"
  ),
  "",
  if(all(checks$pass)) {
    "STATUS=PASS"
  } else {
    "STATUS=FAIL"
  }
)

writeLines(
  audit_lines,
  file.path(
    AUD,
    "22G4C_preflight_audit.txt"
  )
)

cat(
  paste(
    audit_lines,
    collapse="\n"
  ),
  "\n"
)

if(!all(checks$pass)) {
  stop(
    "22G4C preflight failed."
  )
}

cat("\n============================================================\n")
cat("22G4C PREFLIGHT COMPLETE\n")
cat("STATUS=PASS\n")
cat("============================================================\n")

