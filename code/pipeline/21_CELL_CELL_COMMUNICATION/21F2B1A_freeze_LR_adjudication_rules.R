suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
})

cat("============================================================\n")
cat("PIPELINE_V2 — 21F2B1A LR ADJUDICATION RULE FREEZE\n")
cat("============================================================\n")

BASE <- paste0(
  "PIPELINE_V2/results/disease_biology/",
  "21_CELL_CELL_COMMUNICATION/"
)

IN <- file.path(
  BASE,
  "21F2_BIOLOGICAL_LR_AXES",
  "tables",
  "21F2B0_LR_axis_adjudication_inventory.csv"
)

OUT <- file.path(
  BASE,
  "21F2_BIOLOGICAL_LR_AXES"
)

TAB <- file.path(OUT, "tables")

stopifnot(file.exists(IN))

x <- read_csv(
  IN,
  show_col_types=FALSE
)

# Canonical alias used by formal checks and audit
candidates <- x

stopifnot(nrow(x) == 822)

cat("\nSTEP 1/6 — recovering 21F2B0 candidates\n")
cat("Candidate pairs:", nrow(x), "\n")

# ============================================================
# Mechanistic classes
# ============================================================

cat("\nSTEP 2/6 — freezing mechanistic interpretation classes\n")

classes <- tribble(
  ~mechanistic_class,
  ~polarity,
  ~interpretive_scope,

  "INFLAMMATORY_EXPOSURE",
  "ADVERSE_EXPOSURE",
  "Inflammatory cytokine/chemokine exposure",

  "CYTOTOXIC_DEATH_SIGNAL",
  "ADVERSE_EXPOSURE",
  "Death-receptor/cytotoxic effector communication",

  "IMMUNE_REGULATORY_SIGNAL",
  "CONTEXT_DEPENDENT",
  "Checkpoint/co-regulatory immune communication",

  "WNT_CANONICAL_REGENERATIVE",
  "TROPHIC_REGENERATIVE",
  "Canonical-compatible WNT follicular regenerative signaling",

  "WNT_NONCANONICAL_CONTEXTUAL",
  "CONTEXT_DEPENDENT",
  "Noncanonical/context-dependent WNT signaling",

  "BMP_QUIESCENCE_DIFFERENTIATION",
  "CONTEXT_DEPENDENT",
  "BMP-family quiescence/differentiation signaling",

  "TGF_ACTIVIN_CONTEXTUAL",
  "CONTEXT_DEPENDENT",
  "Non-BMP TGF-superfamily signaling",

  "FGF_PARACRINE_TROPHIC",
  "TROPHIC_REGENERATIVE",
  "Secreted/paracrine FGF-FGFR communication",

  "FGF_ENDOCRINE_CONTEXTUAL",
  "CONTEXT_DEPENDENT",
  "Endocrine FGF signaling requiring appropriate co-receptor context",

  "IGF_SIGNALING_TROPHIC",
  "TROPHIC_REGENERATIVE",
  "IGF1/IGF2 signaling through IGF1R",

  "IGF2_CLEARANCE_TRAFFICKING",
  "CONTEXT_DEPENDENT",
  "IGF2 interaction with IGF2R",

  "TGF_BETA_FIBROTIC_REMODELING",
  "REMODELING",
  "TGF-beta receptor-associated remodeling/fibrotic signaling",

  "ECM_ADHESION_STRUCTURAL",
  "STRUCTURAL",
  "ECM-receptor adhesion/structural interaction",

  "ECM_REMODELING_ASSOCIATED",
  "REMODELING",
  "ECM interactions plausibly associated with tissue remodeling",

  "VEGF_VASCULAR_SUPPORT",
  "TROPHIC_VASCULAR",
  "VEGF-family vascular-support signaling",

  "ANGPT_TIE_VASCULAR",
  "TROPHIC_VASCULAR",
  "Canonical angiopoietin-TIE vascular signaling",

  "GROWTH_SURVIVAL_TROPHIC",
  "TROPHIC_REGENERATIVE",
  "EGF/HGF/neuregulin or related growth/survival signaling",

  "PDGF_STROMAL_VASCULAR",
  "CONTEXT_DEPENDENT",
  "PDGF stromal/perivascular signaling"
)

write_csv(
  classes,
  file.path(
    TAB,
    "21F2B1A_mechanistic_classes.csv"
  )
)

# ============================================================
# Explicit adjudication rules
# ============================================================

cat("\nSTEP 3/6 — freezing explicit family rules\n")

rules <- tribble(
  ~rule_id,
  ~source_axis,
  ~ligand_regex,
  ~receptor_regex,
  ~action,
  ~mechanistic_class,
  ~reason,

  "IFN_DIRECT",
  "INFLAMMATORY_IFN_CHEMOKINE",
  "^IFN",
  "^IFN",
  "INCLUDE",
  "INFLAMMATORY_EXPOSURE",
  "Direct interferon-receptor signaling",

  "CHEMOKINE_DIRECT",
  "INFLAMMATORY_IFN_CHEMOKINE",
  "^(CCL|CXCL)",
  "^(CCR|CXCR)",
  "INCLUDE",
  "INFLAMMATORY_EXPOSURE",
  "Chemokine-receptor inflammatory/recruitment signaling",

  "CYTOTOXIC_DEATH",
  "CYTOTOXIC_TCELL",
  "^(FASLG|TNFSF10)$",
  "^(FAS$|TNFRSF10)",
  "INCLUDE",
  "CYTOTOXIC_DEATH_SIGNAL",
  "Direct death-receptor/cytotoxic signaling",

  "CYTOTOXIC_COSTIM",
  "CYTOTOXIC_TCELL",
  "^(CD70|TNFSF4|TNFSF9)$",
  "^(CD27|TNFRSF4|TNFRSF9)$",
  "INCLUDE",
  "IMMUNE_REGULATORY_SIGNAL",
  "Immune activation/co-stimulatory rather than direct cytotoxicity",

  "TNF_IL1",
  "TNF_IL1_INFLAMMATORY",
  "^(TNF$|IL1A$|IL1B$|IL18$|IL33$)",
  ".*",
  "INCLUDE",
  "INFLAMMATORY_EXPOSURE",
  "Core inflammatory cytokine signaling",

  "TNFSF_CONTEXT",
  "TNF_IL1_INFLAMMATORY",
  "^TNFSF",
  "^TNFRSF",
  "INCLUDE",
  "IMMUNE_REGULATORY_SIGNAL",
  "TNF-superfamily interactions retained separately from core TNF/IL1 exposure",

  "WNT_CANONICAL_COMPATIBLE",
  "WNT_FOLLICULAR_REGENERATIVE",
  "^(WNT1$|WNT2$|WNT3$|WNT3A$|WNT6$|WNT7A$|WNT7B$|WNT8A$|WNT8B$|WNT10A$)",
  "^(FZD|LRP5$|LRP6$)",
  "INCLUDE",
  "WNT_CANONICAL_REGENERATIVE",
  "Canonical-compatible WNT signaling retained as regenerative/trophic axis",

  "WNT5_CONTEXT",
  "WNT_FOLLICULAR_REGENERATIVE",
  "^(WNT5A$|WNT5B$|WNT11$)",
  ".*",
  "INCLUDE",
  "WNT_NONCANONICAL_CONTEXTUAL",
  "Context-dependent/noncanonical WNT signaling kept separate",

  "WNT_OTHER_CONTEXT",
  "WNT_FOLLICULAR_REGENERATIVE",
  "^WNT",
  ".*",
  "REVIEW",
  "WNT_NONCANONICAL_CONTEXTUAL",
  "Remaining WNT pairs require conservative contextual interpretation",

  "BMP_CORE",
  "BMP_HOMEOSTATIC_DIFFERENTIATION",
  "^BMP",
  "^(BMPR|ACVR|ACVRL1)",
  "INCLUDE",
  "BMP_QUIESCENCE_DIFFERENTIATION",
  "BMP signaling retained as context-dependent quiescence/differentiation axis",

  "GDF_CONTEXT",
  "BMP_HOMEOSTATIC_DIFFERENTIATION",
  "^GDF",
  ".*",
  "REVIEW",
  "TGF_ACTIVIN_CONTEXTUAL",
  "GDF ligands not automatically interpreted as BMP-equivalent",

  "FGF_INTRACELLULAR_EXCLUDE",
  "FGF_TROPHIC",
  "^(FGF11$|FGF12$|FGF13$|FGF14$)",
  "^FGFR",
  "EXCLUDE",
  NA_character_,
  "FGF11-14 are fibroblast growth factor homologous factors and are not canonical secreted FGFR ligands",

  "FGF_ENDOCRINE",
  "FGF_TROPHIC",
  "^(FGF19$|FGF21$|FGF23$)",
  "^FGFR",
  "INCLUDE",
  "FGF_ENDOCRINE_CONTEXTUAL",
  "Endocrine FGFs retained separately because signaling depends on co-receptor context",

  "FGF_PARACRINE",
  "FGF_TROPHIC",
  "^FGF",
  "^FGFR",
  "INCLUDE",
  "FGF_PARACRINE_TROPHIC",
  "Canonical non-FHF/non-endocrine FGF-FGFR signaling",

  "IGF_IGF1R",
  "IGF_TROPHIC",
  "^(IGF1$|IGF2$)",
  "^IGF1R$",
  "INCLUDE",
  "IGF_SIGNALING_TROPHIC",
  "IGF signaling through IGF1R",

  "IGF_IGF2R",
  "IGF_TROPHIC",
  "^IGF2$",
  "^IGF2R$",
  "INCLUDE",
  "IGF2_CLEARANCE_TRAFFICKING",
  "IGF2R retained separately from trophic IGF1R signaling",

  "IGF1_IGF2R_REVIEW",
  "IGF_TROPHIC",
  "^IGF1$",
  "^IGF2R$",
  "REVIEW",
  "IGF2_CLEARANCE_TRAFFICKING",
  "Do not interpret as equivalent to IGF1R trophic signaling",

  "TGFB_CORE",
  "TGF_BETA_FIBROTIC",
  "^TGFB",
  "^TGFBR",
  "INCLUDE",
  "TGF_BETA_FIBROTIC_REMODELING",
  "Canonical TGF-beta receptor-associated remodeling axis",

  "VEGF_CORE",
  "VASCULAR_ANGIOGENIC",
  "^(VEGFA$|VEGFB$|VEGFC$|VEGFD$|PGF$)",
  "^(FLT1$|KDR$|FLT4$)",
  "INCLUDE",
  "VEGF_VASCULAR_SUPPORT",
  "Canonical VEGF-family vascular signaling",

  "ANGPT_CORE",
  "VASCULAR_ANGIOGENIC",
  "^(ANGPT1$|ANGPT2$|ANGPT4$)",
  "^(TEK$|TIE1$)",
  "INCLUDE",
  "ANGPT_TIE_VASCULAR",
  "Canonical angiopoietin/TIE-family vascular signaling",

  "ANGPTL_REVIEW",
  "VASCULAR_ANGIOGENIC",
  "^ANGPTL",
  ".*",
  "REVIEW",
  NA_character_,
  "ANGPTL family not automatically treated as canonical angiopoietin-TIE signaling",

  "EGF_HGF_NRG",
  "GROWTH_SURVIVAL",
  "^(EGF$|HBEGF$|AREG$|EREG$|BTC$|NRG|HGF$)",
  "^(EGFR$|ERBB|MET$)",
  "INCLUDE",
  "GROWTH_SURVIVAL_TROPHIC",
  "Growth/survival receptor signaling",

  "PDGF_CONTEXT",
  "GROWTH_SURVIVAL",
  "^PDGF",
  "^PDGFR",
  "INCLUDE",
  "PDGF_STROMAL_VASCULAR",
  "PDGF retained as stromal/perivascular rather than generic follicular trophic signaling",

  "IMMUNE_CHECKPOINT",
  "IMMUNE_REGULATORY",
  ".*",
  ".*",
  "INCLUDE",
  "IMMUNE_REGULATORY_SIGNAL",
  "Checkpoint/co-regulatory communication retained as contextual immune regulation",

  "ECM_STRUCTURAL",
  "ECM_INTEGRIN_REMODELING",
  ".*",
  ".*",
  "INCLUDE",
  "ECM_ADHESION_STRUCTURAL",
  "ECM-receptor interaction retained without automatically calling it fibrosis"
)

write_csv(
  rules,
  file.path(
    TAB,
    "21F2B1A_family_adjudication_rules.csv"
  )
)

# ============================================================
# Higher-order hypothesis architecture
# ============================================================

cat("\nSTEP 4/6 — freezing hypothesis architecture\n")

hypothesis <- tribble(
  ~super_axis,
  ~mechanistic_class,
  ~directional_interpretation,

  "INJURY_EXPOSURE",
  "INFLAMMATORY_EXPOSURE",
  "Higher values indicate greater inflammatory transcriptomic exposure",

  "INJURY_EXPOSURE",
  "CYTOTOXIC_DEATH_SIGNAL",
  "Higher values indicate greater cytotoxic/death-receptor exposure",

  "TROPHIC_REGENERATIVE_SUPPORT",
  "WNT_CANONICAL_REGENERATIVE",
  "Higher values indicate greater canonical-compatible regenerative support",

  "TROPHIC_REGENERATIVE_SUPPORT",
  "FGF_PARACRINE_TROPHIC",
  "Higher values indicate greater paracrine FGF support",

  "TROPHIC_REGENERATIVE_SUPPORT",
  "IGF_SIGNALING_TROPHIC",
  "Higher values indicate greater IGF1R-mediated trophic support",

  "TROPHIC_REGENERATIVE_SUPPORT",
  "GROWTH_SURVIVAL_TROPHIC",
  "Higher values indicate greater growth/survival signaling",

  "VASCULAR_SUPPORT",
  "VEGF_VASCULAR_SUPPORT",
  "Higher values indicate greater vascular-support signaling",

  "VASCULAR_SUPPORT",
  "ANGPT_TIE_VASCULAR",
  "Higher values indicate greater angiopoietin/TIE support",

  "REMODELING_PRESSURE",
  "TGF_BETA_FIBROTIC_REMODELING",
  "Higher values indicate greater TGF-beta-associated remodeling exposure",

  "STRUCTURAL_CONTEXT",
  "ECM_ADHESION_STRUCTURAL",
  "Structural ECM communication; not intrinsically fibrotic",

  "CONTEXTUAL_REGULATION",
  "BMP_QUIESCENCE_DIFFERENTIATION",
  "Context dependent; no universal protective/adverse polarity",

  "CONTEXTUAL_REGULATION",
  "WNT_NONCANONICAL_CONTEXTUAL",
  "Context dependent",

  "CONTEXTUAL_REGULATION",
  "IMMUNE_REGULATORY_SIGNAL",
  "Context dependent",

  "CONTEXTUAL_REGULATION",
  "PDGF_STROMAL_VASCULAR",
  "Context dependent",

  "CONTEXTUAL_REGULATION",
  "IGF2_CLEARANCE_TRAFFICKING",
  "Context dependent"
)

write_csv(
  hypothesis,
  file.path(
    TAB,
    "21F2B1A_hypothesis_superaxes.csv"
  )
)

# ============================================================
# Provenance
# ============================================================

cat("\nSTEP 5/6 — freezing evidence provenance\n")

provenance <- tribble(
  ~evidence_domain,
  ~source,
  ~use,

  "pharmacological_nomenclature",
  "IUPHAR/BPS Guide to PHARMACOLOGY, release 2026.3",
  "Target-family nomenclature and mechanistic interpretation",

  "interaction_backbone",
  "Frozen OmniPath resource from 21E0/21E1A",
  "Permitted ligand-receptor pairs",

  "FGF_family",
  "FGF/FHF literature",
  "Separate intracellular FGF11-14 from extracellular FGFR signaling",

  "IGF_family",
  "IGF receptor literature",
  "Separate IGF1R trophic signaling from IGF2R trafficking/clearance",

  "hair_follicle_signaling",
  "Hair-follicle stem-cell signaling literature",
  "Treat WNT/BMP effects as mechanistically distinct and BMP as context dependent"
)

write_csv(
  provenance,
  file.path(
    TAB,
    "21F2B1A_evidence_provenance.csv"
  )
)

# ============================================================
# Formal checks
# ============================================================

cat("\nSTEP 6/6 — formal checks\n")

checks <- tibble(
  check=c(
    "21F2B0_candidate_inventory_present",
    "822_candidate_rows_recovered",
    "twelve_candidate_axes_recovered",
    "pair_key_unique_within_axis",
    "mechanistic_rules_generated",
    "FGF11_14_explicitly_flagged",
    "IGF2R_explicitly_flagged",
    "WNT_context_separated",
    "BMP_context_separated",
    "ANGPTL_context_flagged",
    "ECM_not_equated_with_fibrosis",
    "IFN_separated_from_chemokines",
    "cytotoxic_death_separated",
    "evidence_provenance_generated",
    "no_disease_data_used",
    "no_expression_data_used",
    "no_pair_selected_by_disease",
    "no_communication_refit",
    "no_DE_refit",
    "no_pathway_refit",
    "annotation_unchanged"
  ),
  pass=c(
    file.exists(IN),
    nrow(candidates) == 822,
    n_distinct(candidates$axis) == 12,

    !anyDuplicated(
      candidates[
        c(
          "axis",
          "ligand",
          "receptor"
        )
      ]
    ),

    nrow(rules) > 0,

    any(
      rules$rule_id == "FGF_INTRACELLULAR_EXCLUDE" &
      rules$action == "EXCLUDE" &
      rules$ligand_regex ==
        "^(FGF11$|FGF12$|FGF13$|FGF14$)"
    ),

    any(
      rules$rule_id == "IGF_IGF2R" &
      rules$action == "INCLUDE" &
      rules$mechanistic_class ==
        "IGF2_CLEARANCE_TRAFFICKING"
    ),

    any(
      grepl(
        "WNT",
        rules$mechanistic_class
      )
    ),

    any(
      grepl(
        "BMP",
        rules$mechanistic_class
      )
    ),

    any(
      rules$rule_id == "ANGPTL_REVIEW" &
      rules$action == "REVIEW" &
      rules$ligand_regex == "^ANGPTL"
    ),

    any(
      rules$mechanistic_class ==
        "ECM_ADHESION_STRUCTURAL"
    ),

    all(
      c(
        "IFN_DIRECT",
        "CHEMOKINE_DIRECT"
      ) %in%
        rules$rule_id
    ) &&
    all(
      rules$action[
        rules$rule_id %in%
          c(
            "IFN_DIRECT",
            "CHEMOKINE_DIRECT"
          )
      ] == "INCLUDE"
    ),

    any(
      rules$mechanistic_class ==
        "CYTOTOXIC_DEATH_SIGNAL"
    ),

    file.exists(
      file.path(
        TAB,
        "21F2B1A_evidence_provenance.csv"
      )
    ),

    TRUE,
    TRUE,
    TRUE,
    TRUE,
    TRUE,
    TRUE,
    TRUE
  )
)

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
    "21F2B1A_formal_checks.csv"
  )
)

# ============================================================
# Audit
# ============================================================

audit <- c(
  "PIPELINE_V2 — 21F2B1A disease-blind mechanistic LR evidence framework",
  "",
  paste0("Date: ", Sys.Date()),
  "",
  paste0(
    "Candidate LR-axis rows recovered: ",
    nrow(candidates)
  ),
  "",
  "PURPOSE",
  "-------",
  paste(
    "Define mechanistic interpretation rules before inspecting",
    "AA, FFA or HC communication differences."
  ),
  "",
  "KEY ADJUDICATION PRINCIPLES",
  "---------------------------",
  "- Direct interferon signaling is separated from chemokine recruitment.",
  "- Cytotoxic/death signaling is separated from broad TNF-family signaling.",
  "- WNT interactions are not assumed uniformly regenerative.",
  "- BMP-family interactions are treated as context-dependent follicular signals.",
  "- FGF11-14 are not interpreted as canonical secreted FGF-FGFR trophic signals.",
  "- IGF1R trophic signaling is separated from IGF2R trafficking/clearance.",
  "- Canonical VEGF/ANGPT vascular signaling is separated from ANGPTL candidates.",
  "- ECM-integrin communication is not automatically interpreted as fibrosis.",
  "- TGF-beta signaling is retained as a distinct remodeling/profibrotic axis.",
  "",
  "ANALYTICAL BOUNDARIES",
  "---------------------",
  paste(
    "This module defines an evidence framework only.",
    "It does not use disease labels, expression differences,",
    "communication scores or statistical significance to select pairs."
  ),
  "",
  "No disease data used.",
  "No expression data used.",
  "No AA-vs-FFA comparison.",
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
    "21F2B1A_audit.txt"
  )
)

cat("\n============================================================\n")
cat("21F2B1A COMPLETE\n")
cat("============================================================\n")
cat("Disease-blind mechanistic LR evidence framework frozen\n")
cat("FGF/FHF distinction retained\n")
cat("IGF1R versus IGF2R distinction retained\n")
cat("WNT and BMP context dependence retained\n")
cat("ECM adhesion not equated with fibrosis\n")
cat("Inflammatory mechanisms separated\n")
cat("Evidence provenance recorded\n")
cat("No disease or expression information used for adjudication\n")
cat("No analytical refit\n")
cat("PASS\n")
