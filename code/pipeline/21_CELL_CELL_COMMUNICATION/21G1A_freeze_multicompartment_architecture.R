suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
})

cat("============================================================\n")
cat("PIPELINE_V2 — 21G1A MULTICOMPARTMENT ARCHITECTURE FREEZE\n")
cat("============================================================\n")

BASE <- paste0(
  "PIPELINE_V2/results/disease_biology/",
  "21_CELL_CELL_COMMUNICATION/"
)

G0 <- file.path(
  BASE,
  "21G_MULTICOMPARTMENT_COMMUNICATION",
  "21G0_PREFLIGHT"
)

OUT <- file.path(
  BASE,
  "21G_MULTICOMPARTMENT_COMMUNICATION",
  "21G1_ARCHITECTURE"
)

TAB <- file.path(OUT, "tables")
dir.create(TAB, recursive=TRUE, showWarnings=FALSE)

UNIVERSE_FILE <- file.path(
  G0,
  "tables",
  "21G0_candidate_multicompartment_universe.csv"
)

stopifnot(file.exists(UNIVERSE_FILE))

u <- read_csv(
  UNIVERSE_FILE,
  show_col_types=FALSE
)

stopifnot(nrow(u) == 39)

cat("\nSTEP 1/7 — recovering disease-blind identity universe\n")
cat("Identities:", nrow(u), "\n")

# ============================================================
# 2. Frozen identity -> macrocompartment map
# ============================================================

cat("\nSTEP 2/7 — freezing macrocompartment map\n")

map <- tribble(
  ~communication_identity, ~macrocompartment,

  "Bulge_Isthmus",
  "FOLLICULAR_EPITHELIAL",

  "ORS",
  "FOLLICULAR_EPITHELIAL",

  "ORS_basal",
  "FOLLICULAR_EPITHELIAL",

  "IRS_companion",
  "FOLLICULAR_EPITHELIAL",

  "Hair_shaft",
  "FOLLICULAR_EPITHELIAL",

  "Basal_epithelial",
  "INTERFOLLICULAR_EPITHELIAL",

  "IFE_Infundibulum_differentiating",
  "INTERFOLLICULAR_EPITHELIAL",

  "IFE_Infundibulum_suprabasal",
  "INTERFOLLICULAR_EPITHELIAL",

  "IFE_Infundibulum_terminal",
  "INTERFOLLICULAR_EPITHELIAL",

  "Basal pilosebaceous progenitor-like",
  "PILOSEBACEOUS_SEBACEOUS",

  "Junctional pilosebaceous progenitor-like",
  "PILOSEBACEOUS_SEBACEOUS",

  "Early sebocyte commitment",
  "PILOSEBACEOUS_SEBACEOUS",

  "Terminal lipid-rich inflammatory sebocyte",
  "PILOSEBACEOUS_SEBACEOUS",

  "Clear secretory cell",
  "SECRETORY_DUCTAL",

  "Dark/apocrine-like secretory cell",
  "SECRETORY_DUCTAL",

  "Ductal/clear transitional cell",
  "SECRETORY_DUCTAL",

  "Inflammatory secretory/ductal cell",
  "SECRETORY_DUCTAL",

  "APOE/CCL19 immune-perivascular",
  "FIBROBLAST",

  "COL11A1/DPEP1 dermal-sheath",
  "FIBROBLAST",

  "COL6A5 papillary",
  "FIBROBLAST",

  "PTGDS/IGF1 fibroblast",
  "FIBROBLAST",

  "SFRP2/DPP4 reticular",
  "FIBROBLAST",

  "ZNF804B/ADGRB3 fibroblast",
  "FIBROBLAST",

  "CD4 T",
  "IMMUNE",

  "CD8 T",
  "IMMUNE",

  "Treg",
  "IMMUNE",

  "Macrophage/monocyte",
  "IMMUNE",

  "Langerhans",
  "IMMUNE",

  "ACKR1+ venular EC",
  "ENDOTHELIAL",

  "CA4/RGCC capillary EC",
  "ENDOTHELIAL",

  "GJA5/FBLN5 arterial EC",
  "ENDOTHELIAL",

  "SELE/CXCL10 activated venular EC",
  "ENDOTHELIAL",

  "SLC10A1/CD83 endothelial state",
  "ENDOTHELIAL",

  "Low-complexity EC",
  "ENDOTHELIAL",

  "Contractile mural cell",
  "MURAL_PERIVASCULAR",

  "RGS5/PDGFRB ECM-associated pericyte",
  "MURAL_PERIVASCULAR",

  "Antigen-presentation-high melanocyte",
  "MELANOCYTE",

  "Differentiated/homeostatic melanocyte",
  "MELANOCYTE",

  "Stress/neural-crest-like melanocyte",
  "MELANOCYTE"
)

stopifnot(
  nrow(map) == 39,
  !anyDuplicated(map$communication_identity),
  setequal(
    map$communication_identity,
    u$communication_identity
  )
)

identity_architecture <- u %>%
  left_join(
    map,
    by="communication_identity"
  ) %>%
  select(
    communication_identity,
    macrocompartment,
    candidate_status,
    n_patients_evaluable,
    n_patients_sender,
    n_patients_receiver,
    total_cells,
    n_axis_keys,
    n_keys_ge5_patients
  ) %>%
  arrange(
    macrocompartment,
    communication_identity
  )

write_csv(
  identity_architecture,
  file.path(
    TAB,
    "21G1A_identity_macrocompartment_map.csv"
  )
)

cat("\n===== MACROCOMPARTMENT INVENTORY =====\n")

macro_summary <- identity_architecture %>%
  group_by(
    macrocompartment
  ) %>%
  summarise(
    n_identities=n(),
    n_high_support=
      sum(candidate_status == "HIGH_SUPPORT"),
    n_moderate_support=
      sum(candidate_status == "MODERATE_SUPPORT"),
    n_limited_support=
      sum(candidate_status == "LIMITED_SUPPORT"),
    total_cells=sum(total_cells),
    .groups="drop"
  ) %>%
  arrange(desc(total_cells))

print(
  macro_summary,
  n=Inf,
  width=Inf
)

write_csv(
  macro_summary,
  file.path(
    TAB,
    "21G1A_macrocompartment_inventory.csv"
  )
)

# ============================================================
# 3. Freeze biological mechanism vocabulary
# ============================================================

cat("\nSTEP 3/7 — freezing biological question architecture\n")

mechanisms <- tribble(
  ~biological_program,
  ~mechanistic_class,
  ~interpretation,

  "INFLAMMATORY_INJURY",
  "INFLAMMATORY_EXPOSURE",
  "IFN/cytokine/chemokine-associated exposure",

  "CYTOTOXIC_INJURY",
  "CYTOTOXIC_DEATH_SIGNAL",
  "Death-receptor/cytotoxic communication",

  "IMMUNE_REGULATION",
  "IMMUNE_REGULATORY_SIGNAL",
  "Immune co-regulatory communication",

  "TGF_REMODELING",
  "TGF_BETA_FIBROTIC_REMODELING",
  "TGF-beta-associated remodeling",

  "ECM_ADHESION",
  "ECM_ADHESION_STRUCTURAL",
  "ECM-receptor structural/adhesive communication",

  "PDGF_STROMAL",
  "PDGF_STROMAL_VASCULAR",
  "PDGF-associated stromal/perivascular communication",

  "VEGF_VASCULAR",
  "VEGF_VASCULAR_SUPPORT",
  "VEGF-associated vascular communication",

  "ANGPT_VASCULAR",
  "ANGPT_TIE_VASCULAR",
  "ANGPT-TIE vascular communication",

  "FGF_TROPHIC",
  "FGF_PARACRINE_TROPHIC",
  "Paracrine FGF trophic communication",

  "IGF_TROPHIC",
  "IGF_SIGNALING_TROPHIC",
  "IGF1R-associated trophic communication",

  "GROWTH_SURVIVAL",
  "GROWTH_SURVIVAL_TROPHIC",
  "EGF/ERBB/HGF-associated growth-survival communication",

  "WNT_CANONICAL",
  "WNT_CANONICAL_REGENERATIVE",
  "Canonical-compatible WNT communication",

  "WNT_CONTEXTUAL",
  "WNT_NONCANONICAL_CONTEXTUAL",
  "Noncanonical/context-dependent WNT communication",

  "BMP_CONTEXTUAL",
  "BMP_QUIESCENCE_DIFFERENTIATION",
  "BMP-associated quiescence/differentiation communication",

  "IGF2R_CONTEXTUAL",
  "IGF2_CLEARANCE_TRAFFICKING",
  "IGF2R-associated trafficking/clearance",

  "TGF_ACTIVIN_CONTEXTUAL",
  "TGF_ACTIVIN_CONTEXTUAL",
  "GDF/TGF-superfamily contextual communication"
)

write_csv(
  mechanisms,
  file.path(
    TAB,
    "21G1A_mechanistic_programs.csv"
  )
)

# ============================================================
# 4. Freeze compartment × program hypotheses
# ============================================================

cat("\nSTEP 4/7 — freezing compartment × mechanism hypotheses\n")

hypotheses <- tribble(
  ~receiver_macrocompartment,
  ~biological_program,
  ~analysis_tier,

  "FOLLICULAR_EPITHELIAL",
  "INFLAMMATORY_INJURY",
  "PRIMARY",

  "FOLLICULAR_EPITHELIAL",
  "CYTOTOXIC_INJURY",
  "PRIMARY",

  "FOLLICULAR_EPITHELIAL",
  "FGF_TROPHIC",
  "PRIMARY",

  "FOLLICULAR_EPITHELIAL",
  "IGF_TROPHIC",
  "PRIMARY",

  "FOLLICULAR_EPITHELIAL",
  "GROWTH_SURVIVAL",
  "PRIMARY",

  "FOLLICULAR_EPITHELIAL",
  "WNT_CANONICAL",
  "PRIMARY",

  "FOLLICULAR_EPITHELIAL",
  "BMP_CONTEXTUAL",
  "CONTEXTUAL",

  "FOLLICULAR_EPITHELIAL",
  "IGF2R_CONTEXTUAL",
  "CONTEXTUAL",

  "FIBROBLAST",
  "TGF_REMODELING",
  "PRIMARY",

  "FIBROBLAST",
  "ECM_ADHESION",
  "PRIMARY",

  "FIBROBLAST",
  "PDGF_STROMAL",
  "PRIMARY",

  "FIBROBLAST",
  "INFLAMMATORY_INJURY",
  "PRIMARY",

  "FIBROBLAST",
  "WNT_CONTEXTUAL",
  "CONTEXTUAL",

  "FIBROBLAST",
  "BMP_CONTEXTUAL",
  "CONTEXTUAL",

  "IMMUNE",
  "INFLAMMATORY_INJURY",
  "PRIMARY",

  "IMMUNE",
  "CYTOTOXIC_INJURY",
  "PRIMARY",

  "IMMUNE",
  "IMMUNE_REGULATION",
  "PRIMARY",

  "ENDOTHELIAL",
  "INFLAMMATORY_INJURY",
  "PRIMARY",

  "ENDOTHELIAL",
  "VEGF_VASCULAR",
  "PRIMARY",

  "ENDOTHELIAL",
  "ANGPT_VASCULAR",
  "PRIMARY",

  "ENDOTHELIAL",
  "ECM_ADHESION",
  "CONTEXTUAL",

  "MURAL_PERIVASCULAR",
  "PDGF_STROMAL",
  "PRIMARY",

  "MURAL_PERIVASCULAR",
  "TGF_REMODELING",
  "PRIMARY",

  "MURAL_PERIVASCULAR",
  "ECM_ADHESION",
  "PRIMARY",

  "MURAL_PERIVASCULAR",
  "VEGF_VASCULAR",
  "CONTEXTUAL",

  "PILOSEBACEOUS_SEBACEOUS",
  "INFLAMMATORY_INJURY",
  "PRIMARY",

  "PILOSEBACEOUS_SEBACEOUS",
  "GROWTH_SURVIVAL",
  "PRIMARY",

  "PILOSEBACEOUS_SEBACEOUS",
  "IGF_TROPHIC",
  "PRIMARY",

  "PILOSEBACEOUS_SEBACEOUS",
  "WNT_CANONICAL",
  "PRIMARY",

  "PILOSEBACEOUS_SEBACEOUS",
  "BMP_CONTEXTUAL",
  "CONTEXTUAL",

  "SECRETORY_DUCTAL",
  "INFLAMMATORY_INJURY",
  "PRIMARY",

  "SECRETORY_DUCTAL",
  "GROWTH_SURVIVAL",
  "PRIMARY",

  "SECRETORY_DUCTAL",
  "ECM_ADHESION",
  "CONTEXTUAL",

  "INTERFOLLICULAR_EPITHELIAL",
  "INFLAMMATORY_INJURY",
  "PRIMARY",

  "INTERFOLLICULAR_EPITHELIAL",
  "GROWTH_SURVIVAL",
  "PRIMARY",

  "INTERFOLLICULAR_EPITHELIAL",
  "WNT_CANONICAL",
  "CONTEXTUAL",

  "MELANOCYTE",
  "INFLAMMATORY_INJURY",
  "PRIMARY",

  "MELANOCYTE",
  "IMMUNE_REGULATION",
  "CONTEXTUAL",

  "MELANOCYTE",
  "GROWTH_SURVIVAL",
  "CONTEXTUAL"
)

hypotheses <- hypotheses %>%
  left_join(
    mechanisms,
    by="biological_program"
  )

stopifnot(
  !any(is.na(hypotheses$mechanistic_class))
)

write_csv(
  hypotheses,
  file.path(
    TAB,
    "21G1A_compartment_mechanism_hypotheses.csv"
  )
)

cat("\n===== FROZEN MULTICOMPARTMENT HYPOTHESES =====\n")

print(
  hypotheses,
  n=Inf,
  width=Inf
)

# ============================================================
# 5. Freeze inference hierarchy
# ============================================================

cat("\nSTEP 5/7 — freezing inference hierarchy\n")

inference <- tribble(
  ~component,
  ~specification,

  "unit_of_inference",
  "patient",

  "primary_receiver_level",
  "individual communication identity",

  "macrocompartment_role",
  "biological organization and multiplicity family; identities are not pooled before inference",

  "primary_endpoint",
  "mean_incoming_A_axis",

  "decomposition_endpoint_1",
  "pooled_D_axis",

  "decomposition_endpoint_2",
  "mean_sender_active_intensity",

  "primary_contrasts",
  "AA_vs_HC;FFA_vs_HC;AA_vs_FFA",

  "minimum_support",
  "at least 3 evaluable patients in each contrast group",

  "test",
  "two-sided exhaustive exact permutation test",

  "missingness",
  "non-evaluable identity-axis values remain missing and are never converted to biological zero",

  "sender_followup",
  "sender decomposition only after a receiver-level signal is retained; descriptive unless separately prespecified"
)

write_csv(
  inference,
  file.path(
    TAB,
    "21G1A_inference_hierarchy.csv"
  )
)

# ============================================================
# 6. Freeze multiplicity architecture
# ============================================================

cat("\nSTEP 6/7 — freezing multiplicity architecture\n")

multiplicity <- tribble(
  ~family,
  ~scope,

  "PRIMARY",
  paste(
    "BH correction separately within contrast ×",
    "receiver_macrocompartment × PRIMARY endpoint"
  ),

  "CONTEXTUAL",
  paste(
    "BH correction separately within contrast ×",
    "receiver_macrocompartment × CONTEXTUAL endpoint"
  ),

  "DECOMPOSITION",
  paste(
    "Density and active-intensity endpoints interpreted",
    "only as decomposition of retained A-axis findings"
  )
)

write_csv(
  multiplicity,
  file.path(
    TAB,
    "21G1A_multiplicity_plan.csv"
  )
)

# ============================================================
# 7. Formal checks
# ============================================================

cat("\nSTEP 7/7 — formal checks\n")

checks <- tibble(
  check=c(
    "21G0_universe_present",
    "39_identities_recovered",
    "39_identities_mapped_once",
    "nine_macrocompartments_defined",
    "support_status_retained",
    "limited_support_not_deleted",
    "mechanistic_programs_defined",
    "compartment_hypotheses_defined",
    "primary_contextual_separated",
    "patient_unit_of_inference",
    "identity_level_inference_retained",
    "macrocompartments_not_pooled_for_testing",
   
    "primary_endpoint_A_axis",
    "density_decomposition_retained",
    "intensity_decomposition_retained",
    "three_contrasts_prespecified",
    "minimum_support_three_per_group",
    "exact_permutation_prespecified",
    "missing_not_zero",
    "sender_followup_prespecified",
    "multiplicity_families_defined",
    "no_disease_values_used_for_architecture",
    "no_disease_test",
    "no_AA_FFA_test",
    "no_p_values",
    "no_significance_based_selection",
    "no_threshold_refit",
    "no_communication_refit",
    "no_DE_refit",
    "no_pathway_refit",
    "annotation_unchanged"
  ),

  pass=c(
    file.exists(UNIVERSE_FILE),

    nrow(u) == 39,

    nrow(map) == 39 &&
      !anyDuplicated(
        map$communication_identity
      ) &&
      setequal(
        map$communication_identity,
        u$communication_identity
      ),

    n_distinct(
      map$macrocompartment
    ) == 9,

    all(
      !is.na(
        identity_architecture$candidate_status
      )
    ),

    all(
      u$communication_identity[
        u$candidate_status ==
          "LIMITED_SUPPORT"
      ] %in%
        identity_architecture$communication_identity
    ),

    nrow(mechanisms) > 0 &&
      !anyDuplicated(
        mechanisms$mechanistic_class
      ),

    nrow(hypotheses) > 0,

    all(
      c(
        "PRIMARY",
        "CONTEXTUAL"
      ) %in%
        hypotheses$analysis_tier
    ),

    any(
      inference$component ==
        "unit_of_inference" &
      inference$specification ==
        "patient"
    ),

    any(
      inference$component ==
        "primary_receiver_level" &
      grepl(
        "individual communication identity",
        inference$specification,
        fixed=TRUE
      )
    ),

    any(
      inference$component ==
        "macrocompartment_role" &
      grepl(
        "not pooled",
        inference$specification,
        fixed=TRUE
      )
    ),

    any(
      inference$component ==
        "primary_endpoint" &
      inference$specification ==
        "mean_incoming_A_axis"
    ),

    any(
      inference$component ==
        "decomposition_endpoint_1" &
      inference$specification ==
        "pooled_D_axis"
    ),

    any(
      inference$component ==
        "decomposition_endpoint_2" &
      inference$specification ==
        "mean_sender_active_intensity"
    ),

    any(
      inference$component ==
        "primary_contrasts" &
      inference$specification ==
        "AA_vs_HC;FFA_vs_HC;AA_vs_FFA"
    ),

    any(
      inference$component ==
        "minimum_support" &
      grepl(
        "at least 3 evaluable patients",
        inference$specification,
        fixed=TRUE
      )
    ),

    any(
      inference$component ==
        "test" &
      grepl(
        "exhaustive exact permutation",
        inference$specification,
        fixed=TRUE
      )
    ),

    any(
      inference$component ==
        "missingness" &
      grepl(
        "never converted to biological zero",
        inference$specification,
        fixed=TRUE
      )
    ),

    any(
      inference$component ==
        "sender_followup" &
      grepl(
        "sender decomposition",
        inference$specification,
        fixed=TRUE
      )
    ),

    nrow(multiplicity) == 3 &&
      all(
        c(
          "PRIMARY",
          "CONTEXTUAL",
          "DECOMPOSITION"
        ) %in%
          multiplicity$family
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
    "21G1A_formal_checks.csv"
  )
)

# ============================================================
# Audit
# ============================================================

audit <- c(
  "PIPELINE_V2 — 21G1A multicompartment architecture freeze",
  "",
  paste0(
    "Date: ",
    Sys.Date()
  ),
  "",
  paste0(
    "Frozen communication identities: ",
    nrow(identity_architecture)
  ),
  paste0(
    "Frozen macrocompartments: ",
    n_distinct(
      identity_architecture$macrocompartment
    )
  ),
  paste0(
    "Frozen compartment-mechanism hypotheses: ",
    nrow(hypotheses)
  ),
  "",
  "ARCHITECTURE",
  "------------",
  paste(
    "All 39 frozen communication identities are retained",
    "and mapped exactly once to a biological macrocompartment."
  ),
  paste(
    "Original communication identities remain the primary",
    "receiver-level analytical units."
  ),
  paste(
    "Macrocompartments organize biological hypotheses and",
    "multiplicity families and are not pooled before inference."
  ),
  "",
  "SUPPORT",
  "-------",
  paste(
    "Disease-blind support classifications from 21G0 are retained."
  ),
  paste(
    "Limited-support identities are not deleted; non-evaluable",
    "patient-level observations remain missing rather than",
    "being encoded as biological zero."
  ),
  "",
  "MECHANISTIC SCOPE",
  "-----------------",
  paste(
    "The frozen architecture includes follicular and",
    "interfollicular epithelium, pilosebaceous/sebaceous",
    "states, secretory/ductal states, fibroblasts, immune",
    "cells, endothelium, mural/perivascular cells and melanocytes."
  ),
  paste(
    "Mechanistic programs include inflammatory and cytotoxic",
    "injury, immune regulation, TGF-beta remodeling, ECM",
    "adhesion, PDGF, VEGF/ANGPT vascular communication,",
    "FGF/IGF/growth-survival signaling, WNT and BMP."
  ),
  "",
  "INFERENCE",
  "---------",
  paste(
    "The patient remains the unit of inference."
  ),
  paste(
    "The primary endpoint remains mean incoming normalized",
    "axis exposure (mean_incoming_A_axis)."
  ),
  paste(
    "Pooled axis density and mean sender active intensity",
    "are retained as decomposition endpoints."
  ),
  paste(
    "AA-vs-HC, FFA-vs-HC and AA-vs-FFA are prespecified."
  ),
  paste(
    "Exact exhaustive two-sided permutation testing and",
    "minimum support of three evaluable patients per group",
    "are prespecified for subsequent inference."
  ),
  "",
  "ANALYTICAL BOUNDARIES",
  "---------------------",
  "No disease effect inspected to define this architecture.",
  "No disease comparison performed.",
  "No p-values calculated.",
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
    "21G1A_audit.txt"
  )
)

cat("\n============================================================\n")
cat("21G1A COMPLETE\n")
cat("============================================================\n")
cat("39 communication identities retained\n")
cat("9 biological macrocompartments frozen\n")
cat("Identity-level inference retained\n")
cat("Macrocompartments used for biological organization\n")
cat("Multicompartment mechanistic hypotheses frozen disease-blind\n")
cat("Patient remains unit of inference\n")
cat("Primary A-axis endpoint retained\n")
cat("Density and intensity decomposition retained\n")
cat("No disease comparison\n")
cat("No analytical refit\n")
cat("PASS\n")
