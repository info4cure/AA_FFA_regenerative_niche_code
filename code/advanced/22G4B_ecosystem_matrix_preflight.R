
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
})

cat("\n============================================================\n")
cat("22G4B — IMMUNE × ECOSYSTEM MATRIX PREFLIGHT\n")
cat("============================================================\n")

ROOT <- "PIPELINE_V2"

# ============================================================
# INPUTS
# ============================================================

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

SCORE_FILE <- file.path(
  ROOT,
  "results/disease_biology",
  "22_MULTICOMPARTMENT_DISEASE_AXES",
  "22F_BULGE_STEM_NICHE_INTEGRITY",
  "22F2_FUNCTIONAL_STATE_DENSITY",
  "22F2A_SCORE_FREEZE",
  "tables",
  "22F2A_Bulge_patient_functional_scores.csv"
)

L2_FILE <- file.path(
  ROOT,
  "results/disease_biology",
  "20_CELLULAR_COMPOSITION",
  "20C_PATIENT_LEVEL_COMPOSITION",
  "tables",
  "20C_L2_patient_frequencies_within_L1.csv"
)

IMMUNE_FILE <- file.path(
  ROOT,
  "results/disease_biology",
  "22_MULTICOMPARTMENT_DISEASE_AXES",
  "22G_BULGE_CENTERED_ECOSYSTEM",
  "22G4_IMMUNE_ECOSYSTEM",
  "22G4A_PREFLIGHT",
  "tables",
  "22G4A_eligible_immune_states.csv"
)

OUT <- file.path(
  ROOT,
  "results/disease_biology",
  "22_MULTICOMPARTMENT_DISEASE_AXES",
  "22G_BULGE_CENTERED_ECOSYSTEM",
  "22G4_IMMUNE_ECOSYSTEM",
  "22G4B_ECOSYSTEM_COUPLING"
)

TAB <- file.path(OUT,"tables")
AUD <- file.path(OUT,"audit")

stopifnot(
  file.exists(CFILE),
  file.exists(SCORE_FILE),
  file.exists(L2_FILE),
  file.exists(IMMUNE_FILE)
)

# ============================================================
# BULGE
# ============================================================

C <- readRDS(CFILE)

pos <- as_tibble(
  C$patient_position
) |>
  dplyr::select(
    patient_id,
    disease,
    median_position
  )

scores <- read.csv(
  SCORE_FILE,
  stringsAsFactors=FALSE
) |>
  as_tibble() |>
  dplyr::select(
    patient_id,
    disease,
    stem_niche_mean,
    MHCII_mean,
    cycling_mean,
    CD200_z_mean
  )

bulge <- pos |>
  left_join(
    scores,
    by=c(
      "patient_id",
      "disease"
    )
  )

# ============================================================
# FROZEN L2
# ============================================================

l2 <- read.csv(
  L2_FILE,
  stringsAsFactors=FALSE
) |>
  as_tibble()

# ============================================================
# ECOSYSTEM TARGET DEFINITION
# ============================================================

targets <- tribble(
  ~layer, ~L1, ~L2, ~variable,

  "FOLLICULAR_EPITHELIUM",
  "Epithelial",
  "ORS",
  "ORS",

  "FOLLICULAR_EPITHELIUM",
  "Epithelial",
  "IRS_companion",
  "IRS_companion",

  "FOLLICULAR_EPITHELIUM",
  "Epithelial",
  "Hair_shaft",
  "Hair_shaft",

  "SUPERFICIAL_EPITHELIUM",
  "Epithelial",
  "IFE_Infundibulum_differentiating",
  "IFE_differentiating",

  "SUPERFICIAL_EPITHELIUM",
  "Epithelial",
  "IFE_Infundibulum_suprabasal",
  "IFE_suprabasal",

  "SUPERFICIAL_EPITHELIUM",
  "Epithelial",
  "IFE_Infundibulum_terminal",
  "IFE_terminal",

  "MESENCHYMAL_NICHE",
  "Fibroblast",
  "PTGDS/IGF1 fibroblast",
  "PTGDS_IGF1",

  "MESENCHYMAL_NICHE",
  "Fibroblast",
  "COL11A1/DPEP1 dermal-sheath",
  "COL11A1_DPEP1",

  "MESENCHYMAL_NICHE",
  "Fibroblast",
  "APOE/CCL19 immune-perivascular",
  "APOE_CCL19",

  "MESENCHYMAL_NICHE",
  "Vascular_EC",
  "Epithelial-like EC",
  "Vascular_epithelial_like"
)

# ============================================================
# BUILD ECOSYSTEM L2 MATRIX
# ============================================================

eco_long <- targets |>
  inner_join(
    l2,
    by=c(
      "L1",
      "L2"
    ),
    relationship="many-to-many"
  ) |>
  dplyr::select(
    patient_id,
    disease,
    layer,
    variable,
    L1,
    L2,
    n,
    denominator_L1,
    frequency
  )

eco_wide <- eco_long |>
  dplyr::select(
    patient_id,
    disease,
    variable,
    frequency
  ) |>
  pivot_wider(
    names_from=variable,
    values_from=frequency
  )

# ============================================================
# ADD BULGE FEATURES
# ============================================================

master <- bulge |>
  left_join(
    eco_wide,
    by=c(
      "patient_id",
      "disease"
    )
  )

# ============================================================
# IMMUNE UNIVERSE
# ============================================================

immune_states <- read.csv(
  IMMUNE_FILE,
  stringsAsFactors=FALSE
) |>
  as_tibble() |>
  dplyr::select(
    L1,
    L2,
    feature
  )

immune_long <- l2 |>
  inner_join(
    immune_states,
    by=c(
      "L1",
      "L2"
    ),
    relationship="many-to-many"
  ) |>
  dplyr::select(
    patient_id,
    disease,
    L1,
    L2,
    feature,
    n,
    denominator_L1,
    frequency
  )

# ============================================================
# AUDIT ECOSYSTEM FEATURES
# ============================================================

bulge_features <- tibble(
  layer="BULGE",
  variable=c(
    "median_position",
    "stem_niche_mean",
    "MHCII_mean",
    "cycling_mean",
    "CD200_z_mean"
  )
)

other_features <- targets |>
  dplyr::select(
    layer,
    variable
  )

feature_manifest <- bind_rows(
  bulge_features,
  other_features
)

feature_support <- lapply(
  seq_len(
    nrow(feature_manifest)
  ),
  function(i) {

    v <- feature_manifest$variable[i]

    tibble(
      layer=
        feature_manifest$layer[i],

      variable=v,

      n_valid=
        sum(
          is.finite(
            master[[v]]
          )
        ),

      n_unique=
        n_distinct(
          master[[v]][
            is.finite(
              master[[v]]
            )
          ]
        )
    )
  }
) |>
  bind_rows()

cat("\n===== ECOSYSTEM FEATURE SUPPORT =====\n")

print(
  feature_support,
  n=Inf
)

cat("\n===== IMMUNE STATE SUPPORT =====\n")

immune_support <- immune_long |>
  group_by(
    L1,
    L2,
    feature
  ) |>
  summarise(
    n_valid=n_distinct(patient_id),
    total_cells=sum(n),
    zero_patients=sum(n==0),
    .groups="drop"
  )

print(
  immune_support,
  n=Inf
)

# ============================================================
# SAVE
# ============================================================

write.csv(
  master,
  file.path(
    TAB,
    "22G4B_ecosystem_patient_master.csv"
  ),
  row.names=FALSE
)

write.csv(
  feature_manifest,
  file.path(
    TAB,
    "22G4B_ecosystem_feature_manifest.csv"
  ),
  row.names=FALSE
)

write.csv(
  feature_support,
  file.path(
    TAB,
    "22G4B_ecosystem_feature_support.csv"
  ),
  row.names=FALSE
)

write.csv(
  immune_long,
  file.path(
    TAB,
    "22G4B_immune_patient_source.csv"
  ),
  row.names=FALSE
)

write.csv(
  immune_support,
  file.path(
    TAB,
    "22G4B_immune_state_support.csv"
  ),
  row.names=FALSE
)

# ============================================================
# CHECKS
# ============================================================

checks <- tibble(
  check=c(
    "16 patients in ecosystem master",
    "15 ecosystem features",
    "9 eligible immune states",
    "All Bulge features have 16 observations",
    "All epithelial/niche features have >=15 observations",
    "All immune states represented"
  ),

  pass=c(
    nrow(master)==16,

    nrow(feature_manifest)==15,

    n_distinct(
      immune_states$feature
    )==9,

    all(
      feature_support$n_valid[
        feature_support$layer=="BULGE"
      ]==16
    ),

    all(
      feature_support$n_valid[
        feature_support$layer!="BULGE"
      ]>=15
    ),

    nrow(
      immune_support
    )==9
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
    "22G4B_preflight_checks.csv"
  ),
  row.names=FALSE
)

if(!all(checks$pass)) {
  stop(
    "22G4B ecosystem matrix preflight failed."
  )
}

cat("\n============================================================\n")
cat("22G4B MATRIX PREFLIGHT COMPLETE\n")
cat("STATUS=PASS\n")
cat("============================================================\n")

