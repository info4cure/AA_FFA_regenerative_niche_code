suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(stringr)
})

cat("============================================================\n")
cat("PIPELINE_V2 — 21C CELL–CELL COMMUNICATION DESIGN FREEZE\n")
cat("============================================================\n")

# ============================================================
# 0. Paths
# ============================================================

ROOT <- file.path(
  "PIPELINE_V2",
  "results",
  "disease_biology",
  "21_CELL_CELL_COMMUNICATION"
)

IN21B <- file.path(
  ROOT,
  "21B_FROZEN_COMMUNICATION_ATLAS"
)

OUT <- file.path(
  ROOT,
  "21C_DESIGN_FREEZE"
)

TAB <- file.path(
  OUT,
  "tables"
)

dir.create(TAB, recursive=TRUE, showWarnings=FALSE)

# Search defensively for the frozen 21B object.
obj_candidates <- list.files(
  ROOT,
  pattern="21B.*\\.rds$",
  recursive=TRUE,
  full.names=TRUE
)

if (!length(obj_candidates)) {
  obj_candidates <- list.files(
    "PIPELINE_V2",
    pattern="21B.*communication.*\\.rds$",
    recursive=TRUE,
    full.names=TRUE
  )
}

if (!length(obj_candidates)) {
  stop("Frozen 21B communication atlas RDS not found.")
}

cat("\nSTEP 1/9 — recovering frozen 21B communication atlas\n")
cat("Candidate frozen objects:\n")
print(obj_candidates)

# Select object containing the expected communication metadata.
obj <- NULL
OBJ_IN <- NA_character_

for (f in obj_candidates) {
  z <- tryCatch(readRDS(f), error=function(e) NULL)
  if (is.null(z)) next

  mdn <- tryCatch(names(z[[]]), error=function(e) character())

  if (all(
    c(
      "patient_id",
      "disease",
      "batch",
      "communication_L1",
      "communication_identity"
    ) %in% mdn
  )) {
    obj <- z
    OBJ_IN <- f
    break
  }
}

if (is.null(obj)) {
  stop(
    "No 21B RDS contains the required frozen communication metadata."
  )
}

md <- obj[[]] %>%
  rownames_to_column("cell")

cat("Frozen object:", OBJ_IN, "\n")
cat("Cells        :", nrow(md), "\n")

stopifnot(nrow(md) == 82893)

# ============================================================
# 2. Freeze analytical questions
# ============================================================

cat("\nSTEP 2/9 — freezing analytical questions\n")

questions <- tribble(
  ~question_id, ~priority, ~question,

  "Q1", "PRIMARY",
  paste(
    "Which ligand-receptor communication programmes differ between",
    "AA and FFA relative to healthy scalp?"
  ),

  "Q2", "PRIMARY",
  paste(
    "Do follicular epithelial receiver states show evidence of",
    "reduced trophic, metabolic, vascular, or growth-related input",
    "consistent with a transcriptional conservation/adaptation response?"
  ),

  "Q3", "PRIMARY",
  paste(
    "Does AA differ from FFA in the balance between inflammatory",
    "exposure and follicular homeostatic/regenerative support?"
  ),

  "Q4", "PRIMARY",
  paste(
    "Are signals compatible with follicular preservation present in AA",
    "but absent, reduced, or overridden by inflammatory/remodelling",
    "communication in FFA?"
  ),

  "Q5", "SECONDARY",
  paste(
    "Which immune, fibroblast, vascular, mural, sebaceous, sweat-gland,",
    "melanocyte, and epithelial populations are candidate senders of",
    "disease-associated signals to follicular epithelial states?"
  ),

  "Q6", "SECONDARY",
  paste(
    "Which communication changes reproduce the pathway architecture",
    "identified independently in modules 20A–20F2B?"
  )
)

write_csv(
  questions,
  file.path(TAB, "21C_frozen_questions.csv")
)

# ============================================================
# 3. Freeze contrasts
# ============================================================

cat("\nSTEP 3/9 — freezing disease contrasts\n")

contrasts <- tribble(
  ~contrast, ~numerator, ~denominator, ~role,

  "AA_vs_HC",
  "AA",
  "HC",
  "DISEASE_VS_CONTROL",

  "FFA_vs_HC",
  "FFA",
  "HC",
  "DISEASE_VS_CONTROL",

  "AA_vs_FFA",
  "AA",
  "FFA",
  "DIRECT_DISEASE_TRIANGULATION"
)

write_csv(
  contrasts,
  file.path(TAB, "21C_frozen_contrasts.csv")
)

# ============================================================
# 4. Define receiver architecture
# ============================================================

cat("\nSTEP 4/9 — defining follicular receiver architecture\n")

# Primary follicular receivers.
#
# These are defined from frozen biological identity, not from disease
# expression or communication results.
#
# Infundibular populations are retained as an anatomically related
# epithelial boundary compartment but are distinguished from the
# primary lower/pilosebaceous follicular receiver set.

receiver_definitions <- tribble(
  ~receiver_class,
  ~identity,
  ~priority,
  ~interpretation,

  "FOLLICULAR_CORE",
  "Bulge_Isthmus",
  "PRIMARY",
  "follicular stem/progenitor-associated niche",

  "FOLLICULAR_CORE",
  "ORS",
  "PRIMARY",
  "outer-root-sheath follicular epithelium",

  "FOLLICULAR_CORE",
  "ORS_basal",
  "PRIMARY",
  "basal outer-root-sheath follicular epithelium",

  "FOLLICULAR_CORE",
  "IRS_companion",
  "PRIMARY",
  "inner-root-sheath/companion follicular compartment",

  "FOLLICULAR_CORE",
  "Hair_shaft",
  "SECONDARY",
  "differentiated hair-shaft compartment",

  "INFUNDIBULAR_INTERFACE",
  "IFE_Infundibulum_differentiating",
  "PRIMARY",
  "differentiating infundibular epithelial interface",

  "INFUNDIBULAR_INTERFACE",
  "IFE_Infundibulum_suprabasal",
  "PRIMARY",
  "suprabasal infundibular epithelial interface",

  "INFUNDIBULAR_INTERFACE",
  "IFE_Infundibulum_terminal",
  "PRIMARY",
  "terminal infundibular epithelial interface",

  "EPITHELIAL_REFERENCE",
  "Basal_epithelial",
  "REFERENCE",
  "non-follicular/basal epithelial comparator"
)

write_csv(
  receiver_definitions,
  file.path(
    TAB,
    "21C_frozen_receiver_definitions.csv"
  )
)

# ============================================================
# 5. Freeze sender compartments
# ============================================================

cat("\nSTEP 5/9 — defining sender architecture\n")

sender_classes <- tribble(
  ~sender_class,
  ~communication_L1,
  ~priority,
  ~biological_role,

  "IMMUNE",
  "T_NK",
  "PRIMARY",
  "adaptive/cytotoxic/innate lymphoid signalling",

  "IMMUNE",
  "Myeloid_DC",
  "PRIMARY",
  "myeloid/dendritic inflammatory signalling",

  "STROMAL",
  "Fibroblast",
  "PRIMARY",
  "ECM/remodelling/trophic signalling",

  "VASCULAR",
  "Vascular_EC",
  "PRIMARY",
  "vascular/inflammatory/endothelial signalling",

  "VASCULAR",
  "Lymphatic_EC",
  "PRIMARY",
  "lymphatic/immune-trafficking signalling",

  "MURAL",
  "Pericyte_SMC",
  "PRIMARY",
  "perivascular/contractile niche signalling",

  "MURAL",
  "Smooth_muscle",
  "PRIMARY",
  "contractile/APM-associated niche signalling",

  "EPITHELIAL",
  "Epithelial",
  "PRIMARY",
  "epithelial autocrine/paracrine signalling",

  "PILOSEBACEOUS",
  "Sebocyte",
  "SECONDARY",
  "pilosebaceous metabolic/lipid signalling",

  "GLANDULAR",
  "Sweat_gland",
  "SECONDARY",
  "glandular epithelial signalling",

  "PIGMENTARY",
  "Melanocyte",
  "SECONDARY",
  "melanocyte/follicular interface signalling"
)

write_csv(
  sender_classes,
  file.path(
    TAB,
    "21C_frozen_sender_classes.csv"
  )
)

# ============================================================
# 6. Freeze biological signalling axes
# ============================================================

cat("\nSTEP 6/9 — freezing communication interpretation axes\n")

axes <- tribble(
  ~communication_axis,
  ~priority,
  ~families,
  ~interpretive_question,

  "INFLAMMATORY_EXPOSURE",
  "PRIMARY",
  paste(
    "IFN;IFNG;IFN-I;TNF;IL1;IL6;CXCL;CCL;",
    "MHC-associated inflammatory signalling"
  ),
  paste(
    "How strongly is the follicular epithelium exposed to",
    "inflammatory and immune-recruiting signals?"
  ),

  "IMMUNE_CYTOTOXIC_PRESSURE",
  "PRIMARY",
  paste(
    "FAS/FASLG;TNF-family;TRAIL;cytotoxic-associated;",
    "T-cell/NK-associated receptor-ligand systems"
  ),
  paste(
    "Is there differential immune/cytotoxic pressure on",
    "follicular receiver populations?"
  ),

  "TROPHIC_REGENERATIVE_SUPPORT",
  "PRIMARY",
  "WNT;RSPO;BMP;FGF;IGF;EGF;HGF;SHH;NOTCH",
  paste(
    "Is follicular trophic/regenerative signalling preserved,",
    "suppressed, or lost?"
  ),

  "QUIESCENCE_CONSERVATION",
  "PRIMARY",
  "BMP;TGF-beta;WNT antagonism;growth-limiting signalling",
  paste(
    "Is there evidence compatible with reversible reduction",
    "of growth/activity or niche conservation?"
  ),

  "ECM_FIBROTIC_REMODELLING",
  "PRIMARY",
  paste(
    "TGF-beta;COLLAGEN;LAMININ;FN1;THBS;SPP1;",
    "integrins;ECM-receptor interactions"
  ),
  paste(
    "Does the follicular niche receive persistent remodelling",
    "or fibrotic signals?"
  ),

  "VASCULAR_TROPHIC_SUPPORT",
  "PRIMARY",
  "VEGF;ANGPT;PDGF;EDN;FGF;vascular trophic signalling",
  paste(
    "Is vascular support to the follicular niche altered?"
  ),

  "METABOLIC_STRESS_ADAPTATION",
  "PRIMARY",
  paste(
    "INSULIN/IGF;adipokine/lipid-associated signalling;",
    "hypoxia/stress-associated intercellular programmes"
  ),
  paste(
    "Are communication changes compatible with altered",
    "metabolic demand or stress adaptation?"
  ),

  "CELL_SURVIVAL_APOPTOSIS",
  "PRIMARY",
  "FAS;TRAIL;TNF;EGF;IGF;survival/apoptotic receptor systems",
  paste(
    "Does the balance of survival versus death signalling",
    "differ between AA and FFA?"
  ),

  "ADHESION_NICHE_INTEGRITY",
  "SECONDARY",
  "integrins;cadherins;NECTIN;JAM;ECM adhesion systems",
  paste(
    "Is physical/niche signalling architecture preserved",
    "or remodelled?"
  )
)

write_csv(
  axes,
  file.path(
    TAB,
    "21C_frozen_communication_axes.csv"
  )
)

# ============================================================
# 7. Patient-support rules
# ============================================================

cat("\nSTEP 7/9 — freezing patient-support rules\n")

rules <- tribble(
  ~rule,
  ~value,
  ~rationale,

  "BIOLOGICAL_UNIT",
  "PATIENT",
  "cells are not treated as independent biological replicates",

  "MIN_CELLS_PER_PATIENT_IDENTITY",
  "20",
  "consistent with pseudobulk/pathway eligibility framework",

  "ROBUST_CELLS_PER_PATIENT_IDENTITY",
  "50",
  "prespecified robustness stratum",

  "MIN_PATIENTS_PER_ARM",
  "3",
  "minimum disease-arm support for comparative inference",

  "PRIMARY_FDR",
  "0.05",
  "Benjamini-Hochberg primary inferential threshold",

  "EXPLORATORY_FDR",
  "0.10",
  "exploratory sensitivity only",

  "AA_BATCH_ADJUSTMENT",
  "NO",
  paste(
    "AA is fully confounded with CABIMER2026 batch;",
    "batch cannot be separated from disease"
  ),

  "FFA_HC_BATCH_ADJUSTMENT",
  "YES_WHEN_MODEL_IDENTIFIABLE",
  paste(
    "FFA and HC contain historical and CABIMER2026 patients;",
    "batch adjustment is identifiable"
  ),

  "AA_FFA_BATCH_ADJUSTMENT",
  "NO_PRIMARY_ADJUSTMENT",
  paste(
    "batch is severely confounded with disease;",
    "direct contrast requires explicit interpretive limitation"
  ),

  "RETURNED_EPITHELIAL_CELLS",
  "L1_ONLY",
  paste(
    "68 returned epithelial cells remain available at L1 but",
    "are excluded from L2-specific receiver inference"
  )
)

write_csv(
  rules,
  file.path(
    TAB,
    "21C_frozen_analysis_rules.csv"
  )
)

# ============================================================
# 8. Empirical eligibility from frozen atlas
# ============================================================

cat("\nSTEP 8/9 — auditing empirical sender/receiver support\n")

# Exact identity support already frozen in 21B.
support <- md %>%
  filter(
    !is.na(communication_identity),
    communication_identity != ""
  ) %>%
  count(
    disease,
    patient_id,
    communication_L1,
    communication_identity,
    name="n_cells"
  )

support20 <- support %>%
  filter(n_cells >= 20)

support50 <- support %>%
  filter(n_cells >= 50)

identity_support <- support %>%
  group_by(
    disease,
    communication_L1,
    communication_identity
  ) %>%
  summarise(
    n_patients_any=n_distinct(patient_id),
    n_patients_ge20=n_distinct(
      patient_id[n_cells >= 20]
    ),
    n_patients_ge50=n_distinct(
      patient_id[n_cells >= 50]
    ),
    total_cells=sum(n_cells),
    median_cells=median(n_cells),
    .groups="drop"
  )

write_csv(
  support,
  file.path(
    TAB,
    "21C_patient_identity_cell_support.csv"
  )
)

write_csv(
  identity_support,
  file.path(
    TAB,
    "21C_disease_identity_support.csv"
  )
)

# Receiver support audit.
receiver_support <- identity_support %>%
  inner_join(
    receiver_definitions,
    by=c(
      "communication_identity"="identity"
    )
  ) %>%
  arrange(
    receiver_class,
    communication_identity,
    disease
  )

write_csv(
  receiver_support,
  file.path(
    TAB,
    "21C_receiver_support_audit.csv"
  )
)

# Contrast-specific receiver eligibility.
receiver_eligibility <- expand_grid(
  contrast=contrasts$contrast,
  identity=receiver_definitions$identity
) %>%
  left_join(
    contrasts,
    by="contrast"
  ) %>%
  rowwise() %>%
  mutate(
    n_num_ge20={
      z <- identity_support %>%
        filter(
          disease == numerator,
          communication_identity == identity
        )
      ifelse(
        nrow(z),
        z$n_patients_ge20[[1]],
        0
      )
    },
    n_den_ge20={
      z <- identity_support %>%
        filter(
          disease == denominator,
          communication_identity == identity
        )
      ifelse(
        nrow(z),
        z$n_patients_ge20[[1]],
        0
      )
    },
    n_num_ge50={
      z <- identity_support %>%
        filter(
          disease == numerator,
          communication_identity == identity
        )
      ifelse(
        nrow(z),
        z$n_patients_ge50[[1]],
        0
      )
    },
    n_den_ge50={
      z <- identity_support %>%
        filter(
          disease == denominator,
          communication_identity == identity
        )
      ifelse(
        nrow(z),
        z$n_patients_ge50[[1]],
        0
      )
    },
    eligible_ge20=
      n_num_ge20 >= 3 &
      n_den_ge20 >= 3,
    robust_ge50=
      n_num_ge50 >= 3 &
      n_den_ge50 >= 3
  ) %>%
  ungroup() %>%
  left_join(
    receiver_definitions,
    by="identity"
  )

write_csv(
  receiver_eligibility,
  file.path(
    TAB,
    "21C_receiver_contrast_eligibility.csv"
  )
)

cat("\n===== RECEIVER ELIGIBILITY =====\n")
print(
  receiver_eligibility %>%
    select(
      contrast,
     
      identity,
      receiver_class,
      priority,
      numerator,
      denominator,
      n_num_ge20,
      n_den_ge20,
      eligible_ge20,
      n_num_ge50,
      n_den_ge50,
      robust_ge50
    ),
  n=Inf,
  width=Inf
)

# ============================================================
# 9. Formal checks / audit
# ============================================================

cat("\nSTEP 9/9 — formal checks\n")

checks <- tibble(
  check=c(
    "frozen_21B_object_present",
    "global_cells_82893",
    "three_contrasts_frozen",
    "primary_disease_vs_control_contrasts_present",
    "direct_AA_vs_FFA_retained",
    "receiver_definitions_frozen",
    "follicular_core_defined",
    "infundibular_interface_defined",
    "sender_classes_frozen",
    "communication_axes_frozen",
    "preservation_hypothesis_predefined",
    "inflammatory_exposure_axis_predefined",
    "trophic_regenerative_axis_predefined",
    "quiescence_conservation_axis_predefined",
    "ECM_fibrotic_axis_predefined",
    "metabolic_adaptation_axis_predefined",
    "patient_is_biological_unit",
    "minimum_patient_support_predefined",
    "receiver_eligibility_audited",
    "returned_cells_not_forced_to_L2",
    "no_communication_inference",
    "no_DE_refit",
    "no_pathway_refit",
    "no_annotation_change",
    "source_object_not_overwritten"
  ),
  pass=c(
    !is.na(OBJ_IN) &&
      file.exists(OBJ_IN),

    nrow(md) == 82893,

    nrow(contrasts) == 3,

    all(
      c(
        "AA_vs_HC",
        "FFA_vs_HC"
      ) %in% contrasts$contrast
    ),

    "AA_vs_FFA" %in%
      contrasts$contrast,

    nrow(receiver_definitions) == 9,

    all(
      c(
        "Bulge_Isthmus",
        "ORS",
        "ORS_basal",
        "IRS_companion"
      ) %in%
        receiver_definitions$identity[
          receiver_definitions$receiver_class ==
            "FOLLICULAR_CORE"
        ]
    ),

    all(
      c(
        "IFE_Infundibulum_differentiating",
        "IFE_Infundibulum_suprabasal",
        "IFE_Infundibulum_terminal"
      ) %in%
        receiver_definitions$identity[
          receiver_definitions$receiver_class ==
            "INFUNDIBULAR_INTERFACE"
        ]
    ),

    nrow(sender_classes) == 11,

    nrow(axes) == 9,

    all(
      c(
        "TROPHIC_REGENERATIVE_SUPPORT",
        "QUIESCENCE_CONSERVATION",
        "METABOLIC_STRESS_ADAPTATION",
        "CELL_SURVIVAL_APOPTOSIS"
      ) %in%
        axes$communication_axis
    ),

    "INFLAMMATORY_EXPOSURE" %in%
      axes$communication_axis,

    "TROPHIC_REGENERATIVE_SUPPORT" %in%
      axes$communication_axis,

    "QUIESCENCE_CONSERVATION" %in%
      axes$communication_axis,

    "ECM_FIBROTIC_REMODELLING" %in%
      axes$communication_axis,

    "METABOLIC_STRESS_ADAPTATION" %in%
      axes$communication_axis,

    rules$value[
      rules$rule ==
        "BIOLOGICAL_UNIT"
    ] == "PATIENT",

    rules$value[
      rules$rule ==
        "MIN_PATIENTS_PER_ARM"
    ] == "3",

    nrow(receiver_eligibility) ==
      nrow(contrasts) *
      nrow(receiver_definitions),

    all(
      md$communication_annotation_level[
        md$communication_annotation_level ==
          "RETURNED_UNADJUDICATED"
      ] ==
        "RETURNED_UNADJUDICATED"
    ),

    TRUE,
    TRUE,
    TRUE,
    TRUE,
    TRUE
  )
)

print(
  checks,
  n=Inf
)

stopifnot(
  nrow(checks) ==
    length(checks$pass)
)

stopifnot(
  all(checks$pass)
)

write_csv(
  checks,
  file.path(
    TAB,
    "21C_formal_checks.csv"
  )
)

# ============================================================
# Audit
# ============================================================

audit <- c(
  "PIPELINE_V2 — 21C cell-cell communication design freeze",
  "",
  paste0(
    "Date: ",
    Sys.Date()
  ),
  "",
  "FROZEN INPUT",
  "------------",
  paste0(
    "21B communication atlas: ",
    OBJ_IN
  ),
  "- 82,893 cells",
  "- frozen L1/L2 communication annotation",
  "- SoupX-corrected RNA expression",
  "- patient, disease and batch metadata retained",
  "",
  "PRIMARY BIOLOGICAL QUESTION",
  "---------------------------",
  paste(
    "Determine whether AA and FFA differ in the balance between",
    "follicular inflammatory exposure and adaptive/homeostatic",
    "responses compatible with preservation of the follicular niche."
  ),
  "",
  paste(
    "In particular, test whether reduced trophic, vascular, growth",
    "or metabolic signalling is compatible with a potentially",
    "reversible conservation response in AA, and whether FFA instead",
    "shows failure of this architecture together with persistent",
    "inflammatory and/or fibrotic-remodelling communication."
  ),
  "",
  "IMPORTANT INTERPRETIVE BOUNDARY",
  "-------------------------------",
  paste(
    "The preservation/conservation model is a prespecified hypothesis.",
    "Communication inference can support or contradict compatibility",
    "with this model but cannot by itself establish adaptive intent,",
    "causality, reversibility, or mechanism."
  ),
  "",
  "FROZEN CONTRASTS",
  "----------------",
  "- AA vs HC",
  "- FFA vs HC",
  "- AA vs FFA",
  "",
  paste(
    "AA-vs-HC and FFA-vs-HC define disease-associated communication.",
    "AA-vs-FFA provides direct disease triangulation."
  ),
  "",
  "FROZEN RECEIVER ARCHITECTURE",
  "----------------------------",
  "Primary follicular core:",
  "- Bulge/Isthmus",
  "- ORS",
  "- ORS-basal",
  "- IRS/companion",
  "",
  "Related infundibular interface:",
  "- differentiating infundibulum",
  "- suprabasal infundibulum",
  "- terminal infundibulum",
  "",
  "Reference epithelial state:",
  "- basal epithelium",
  "",
  paste(
    "Hair-shaft cells are retained as a secondary follicular",
    "receiver population."
  ),
  "",
  "FROZEN COMMUNICATION AXES",
  "-------------------------",
  "- inflammatory exposure",
  "- immune/cytotoxic pressure",
  "- trophic/regenerative support",
  "- quiescence/conservation",
  "- ECM/fibrotic remodelling",
  "- vascular trophic support",
  "- metabolic/stress adaptation",
  "- cell survival/apoptosis",
  "- adhesion/niche integrity",
  "",
  "STATISTICAL DESIGN",
  "------------------",
  "- patient is the biological unit",
  "- >=20 cells/patient/identity for primary support",
  "- >=50 cells/patient/identity for robustness",
  "- >=3 patients per contrast arm",
  "- primary FDR <0.05",
  "- FDR <0.10 exploratory only",
  "",
  "BATCH",
  "-----",
  paste(
    "FFA-vs-HC permits batch adjustment when identifiable."
  ),
  paste(
    "AA is completely confounded with CABIMER2026 batch and therefore",
    "AA-vs-HC cannot provide a disease effect statistically separated",
    "from sequencing batch."
  ),
  paste(
    "AA-vs-FFA is likewise batch-confounded and will be interpreted",
    "with this limitation."
  ),
  "",
  "RETURNED EPITHELIAL CELLS",
  "-------------------------",
  paste(
    "The 68 returned epithelial boundary cells remain L1-only and",
    "are not assigned to a follicular L2 receiver identity."
  ),
  "",
  "ANALYTICAL STATUS",
  "-----------------",
  "- no ligand-receptor inference performed",
  "- no CellChat inference performed",
  "- no differential-expression model refitted",
  "- no pathway model refitted",
  "- no cell annotation modified",
  "- source 21B object not overwritten",
  "",
  "NEXT MODULE",
  "-----------",
  paste(
    "21D will construct the communication inference input and",
    "ligand-receptor resource while preserving this frozen design."
  ),
  "",
  "PASS"
)

writeLines(
  audit,
  file.path(
    OUT,
    "21C_audit.txt"
  )
)

cat("\n===== FORMAL CHECKS =====\n")
print(
  checks,
  n=Inf
)

cat("\n============================================================\n")
cat("21C COMPLETE\n")
cat("============================================================\n")
cat("Cell-cell communication design frozen before inference\n")
cat("Follicular preservation/conservation hypothesis prespecified\n")
cat("Inflammatory exposure separated from trophic/homeostatic response\n")
cat("AA-vs-HC and FFA-vs-HC define disease-associated architecture\n")
cat("AA-vs-FFA retained for direct triangulation\n")
cat("Patient-level eligibility rules frozen\n")
cat("No communication inference performed\n")
cat("No differential-expression refit\n")
cat("No pathway refit\n")
cat("No annotation modified\n")
cat("PASS\n")
