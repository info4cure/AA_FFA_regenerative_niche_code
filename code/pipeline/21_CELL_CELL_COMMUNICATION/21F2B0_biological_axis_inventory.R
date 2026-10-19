suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(tidyr)
  library(tibble)
})

cat("============================================================\n")
cat("PIPELINE_V2 — 21F2B0 BIOLOGICAL LR-AXIS INVENTORY\n")
cat("============================================================\n")

BASE <- paste0(
  "PIPELINE_V2/results/disease_biology/",
  "21_CELL_CELL_COMMUNICATION/"
)

E1 <- file.path(
  BASE,
  "21E1_OBSERVABLE_UNIVERSE"
)

F1 <- file.path(
  BASE,
  "21F_PATIENT_LEVEL_COMMUNICATION"
)

OUT <- file.path(
  BASE,
  "21F2_BIOLOGICAL_LR_AXES"
)

TAB <- file.path(OUT, "tables")

dir.create(
  TAB,
  recursive=TRUE,
  showWarnings=FALSE
)

LR_FILE <- file.path(
  E1,
  "tables",
  "21E1A_observable_LR_universe.csv"
)

DEN_FILE <- file.path(
  F1,
  "tables",
  "21F1B_sender_receiver_structural_denominator.csv"
)

stopifnot(
  file.exists(LR_FILE),
  file.exists(DEN_FILE)
)

# ============================================================
# 1. Recover unique frozen LR universe
# ============================================================

cat("\nSTEP 1/7 — recovering frozen symbol-level LR universe\n")

lr0 <- read_csv(
  LR_FILE,
  show_col_types=FALSE
)

lr <- lr0 %>%
  distinct(
    source_genesymbol,
    target_genesymbol
  ) %>%
  rename(
    ligand=source_genesymbol,
    receptor=target_genesymbol
  )

cat("Observable resource rows:", nrow(lr0), "\n")
cat("Unique LR pairs:", nrow(lr), "\n")

stopifnot(
  nrow(lr) == 4547
)

# ============================================================
# 2. Freeze candidate biological-axis vocabulary
# ============================================================

cat("\nSTEP 2/7 — freezing candidate biological-axis vocabulary\n")

axes <- tribble(
  ~axis,
  ~purpose,

  "INFLAMMATORY_IFN_CHEMOKINE",
  "Interferon-associated and chemokine-mediated inflammatory exposure",

  "CYTOTOXIC_TCELL",
  "T/NK recruitment, activation and cytotoxic immune communication",

  "TNF_IL1_INFLAMMATORY",
  "TNF/IL1-family inflammatory signaling",

  "WNT_FOLLICULAR_REGENERATIVE",
  "WNT-associated follicular regenerative/homeostatic signaling",

  "BMP_HOMEOSTATIC_DIFFERENTIATION",
  "BMP-associated epithelial/follicular homeostasis and differentiation",

  "FGF_TROPHIC",
  "FGF-family trophic/growth signaling",

  "IGF_TROPHIC",
  "IGF-family trophic/survival signaling",

  "TGF_BETA_FIBROTIC",
  "TGF-beta-family profibrotic/remodeling signaling",

  "ECM_INTEGRIN_REMODELING",
  "ECM-receptor/integrin-associated structural and remodeling communication",

  "VASCULAR_ANGIOGENIC",
  "VEGF/angiogenic and vascular-support signaling",

  "GROWTH_SURVIVAL",
  "Additional growth-factor and survival signaling",

  "IMMUNE_REGULATORY",
  "Immune-regulatory/checkpoint or suppressive communication"
)

write_csv(
  axes,
  file.path(
    TAB,
    "21F2B0_candidate_axis_vocabulary.csv"
  )
)

# ============================================================
# 3. Freeze candidate gene-family rules
# ============================================================

cat("\nSTEP 3/7 — defining transparent candidate-family rules\n")

rules <- tribble(
  ~axis,
  ~endpoint,
  ~regex,

  "INFLAMMATORY_IFN_CHEMOKINE",
  "LIGAND",
  "^(IFNA|IFNB|IFNG|CXCL|CCL)",

  "INFLAMMATORY_IFN_CHEMOKINE",
  "RECEPTOR",
  "^(IFNAR|IFNGR|CXCR|CCR)",

  "CYTOTOXIC_TCELL",
  "LIGAND",
  "^(FASLG|TNFSF10|TNFSF9|TNFSF4|CD70)$",

  "CYTOTOXIC_TCELL",
  "RECEPTOR",
  "^(FAS|TNFRSF10|TNFRSF9|TNFRSF4|CD27)",

  "TNF_IL1_INFLAMMATORY",
  "LIGAND",
  "^(TNF$|TNFSF|IL1A$|IL1B$|IL18$|IL33$)",

  "TNF_IL1_INFLAMMATORY",
  "RECEPTOR",
  "^(TNFRSF|IL1R|IL18R|IL1RL)",

  "WNT_FOLLICULAR_REGENERATIVE",
  "LIGAND",
  "^WNT",

  "WNT_FOLLICULAR_REGENERATIVE",
  "RECEPTOR",
  "^(FZD|LRP5$|LRP6$|ROR1$|ROR2$|RYK$)",

  "BMP_HOMEOSTATIC_DIFFERENTIATION",
  "LIGAND",
  "^(BMP|GDF)",

  "BMP_HOMEOSTATIC_DIFFERENTIATION",
  "RECEPTOR",
  "^(BMPR|ACVR)",

  "FGF_TROPHIC",
  "LIGAND",
  "^FGF",

  "FGF_TROPHIC",
  "RECEPTOR",
  "^FGFR",

  "IGF_TROPHIC",
  "LIGAND",
  "^(IGF1$|IGF2$)",

  "IGF_TROPHIC",
  "RECEPTOR",
  "^(IGF1R$|IGF2R$)",

  "TGF_BETA_FIBROTIC",
  "LIGAND",
  "^(TGFB1$|TGFB2$|TGFB3$)",

  "TGF_BETA_FIBROTIC",
  "RECEPTOR",
  "^(TGFBR1$|TGFBR2$|TGFBR3$)",

  "ECM_INTEGRIN_REMODELING",
  "LIGAND",
  "^(COL[0-9]|FN1$|LAMA|LAMB|LAMC|THBS|VCAN$|POSTN$|TNC$)",

  "ECM_INTEGRIN_REMODELING",
  "RECEPTOR",
  "^(ITGA|ITGB|CD44$|SDC)",

  "VASCULAR_ANGIOGENIC",
  "LIGAND",
  "^(VEGFA$|VEGFB$|VEGFC$|VEGFD$|PGF$|ANGPT)",

  "VASCULAR_ANGIOGENIC",
  "RECEPTOR",
  "^(FLT1$|KDR$|FLT4$|TEK$|TIE1$)",

  "GROWTH_SURVIVAL",
  "LIGAND",
  "^(EGF$|HBEGF$|AREG$|EREG$|BTC$|NRG|HGF$|PDGFA$|PDGFB$|PDGFC$|PDGFD$)",

  "GROWTH_SURVIVAL",
  "RECEPTOR",
  "^(EGFR$|ERBB|MET$|PDGFRA$|PDGFRB$)",

  "IMMUNE_REGULATORY",
  "LIGAND",
  "^(CD274$|PDCD1LG2$|LGALS9$|CD80$|CD86$|ICOSLG$)",

  "IMMUNE_REGULATORY",
  "RECEPTOR",
  "^(PDCD1$|CTLA4$|HAVCR2$|ICOS$|CD28$)"
)

write_csv(
  rules,
  file.path(
    TAB,
    "21F2B0_candidate_gene_family_rules.csv"
  )
)

# ============================================================
# 4. Candidate endpoint inventory
# ============================================================

cat("\nSTEP 4/7 — inventorying candidate endpoints\n")

genes <- sort(
  unique(
    c(
      lr$ligand,
      lr$receptor
    )
  )
)

candidate_endpoints <- bind_rows(
  lapply(
    seq_len(nrow(rules)),
    function(i) {

      rr <- rules[i, ]

      universe <- if (
        rr$endpoint == "LIGAND"
      ) {
        unique(lr$ligand)
      } else {
        unique(lr$receptor)
      }

      hit <- universe[
        str_detect(
          universe,
          regex(
            rr$regex,
            ignore_case=FALSE
          )
        )
      ]

      tibble(
        axis=rr$axis,
        endpoint=rr$endpoint,
        regex=rr$regex,
        gene=sort(unique(hit))
      )
    }
  )
) %>%
  distinct()

endpoint_summary <- candidate_endpoints %>%
  count(
    axis,
    endpoint,
    name="n_candidate_genes"
  )

cat("\n===== CANDIDATE ENDPOINTS =====\n")
print(
  endpoint_summary,
  n=Inf,
  width=Inf
)

write_csv(
  candidate_endpoints,
  file.path(
    TAB,
    "21F2B0_candidate_endpoints.csv"
  )
)

write_csv(
  endpoint_summary,
  file.path(
    TAB,
    "21F2B0_candidate_endpoint_summary.csv"
  )
)

# ============================================================
# 5. Candidate LR-pair inventory
# ============================================================

cat("\nSTEP 5/7 — inventorying candidate LR pairs\n")

candidate_pairs <- bind_rows(
  lapply(
    axes$axis,
    function(ax) {

      lig <- candidate_endpoints %>%
        filter(
          axis == ax,
          endpoint == "LIGAND"
        ) %>%
        pull(gene) %>%
        unique()

      rec <- candidate_endpoints %>%
        filter(
          axis == ax,
          endpoint == "RECEPTOR"
        ) %>%
        pull(gene) %>%
        unique()

      lr %>%
        filter(
          ligand %in% lig,
          receptor %in% rec
        ) %>%
        mutate(
          axis=ax
        )
    }
  )
) %>%
  distinct(
    axis,
    ligand,
    receptor
  )

pair_summary <- candidate_pairs %>%
  group_by(axis) %>%
  summarise(
    n_pairs=n(),
    n_ligands=n_distinct(ligand),
    n_receptors=n_distinct(receptor),
    .groups="drop"
  ) %>%
  right_join(
    axes %>%
      select(axis),
    by="axis"
  ) %>%
  mutate(
    across(
      c(
        n_pairs,
        n_ligands,
        n_receptors
      ),
      ~coalesce(.x, 0L)
    )
  )

cat("\n===== CANDIDATE LR PAIRS BY AXIS =====\n")
print(
  pair_summary,
  n=Inf,
  width=Inf
)

write_csv(
  candidate_pairs,
  file.path(
    TAB,
    "21F2B0_candidate_LR_pairs.csv"
  )
)

write_csv(
  pair_summary,
  file.path(
    TAB,
    "21F2B0_candidate_LR_pair_summary.csv"
  )
)

# ============================================================
# 6. Explicit pair listing for adjudication
# ============================================================

cat("\nSTEP 6/7 — writing adjudication inventory\n")

adjudication_inventory <- candidate_pairs %>%
  left_join(
    axes,
    by="axis"
  ) %>%
  arrange(
    axis,
    ligand,
    receptor
  ) %>%
  mutate(
    adjudication_status=
      "UNADJUDICATED_CANDIDATE",
    include_final=NA,
    adjudication_note=NA_character_
  )

write_csv(
  adjudication_inventory,
  file.path(
    TAB,
    "21F2B0_LR_axis_adjudication_inventory.csv"
  )
)

cat("\n===== FULL CANDIDATE LR INVENTORY =====\n")
print(
  adjudication_inventory,
  n=Inf,
  width=Inf
)

# ============================================================
# 7. Formal checks
# ============================================================

cat("\nSTEP 7/7 — formal checks\n")

checks <- tibble(
  check=c(
    "21E1A_LR_input_present",
    "4547_unique_symbol_pairs",
    "twelve_candidate_axes_defined",
    "rules_explicit",
    "candidate_endpoints_generated",
    "candidate_pairs_generated",
    "adjudication_inventory_generated",
    "all_candidate_pairs_in_frozen_OmniPath",
    "no_axis_finalized",
    "no_disease_filtering",
    "no_expression_filtering",
    "no_disease_test",
    "no_AA_FFA_test",
    "no_communication_refit",
    "no_DE_refit",
    "no_pathway_refit",
    "annotation_unchanged"
  ),
  pass=c(
    file.exists(LR_FILE),
    nrow(lr) == 4547,
    nrow(axes) == 12,
    nrow(rules) > 0,
    nrow(candidate_endpoints) > 0,
    nrow(candidate_pairs) > 0,
    nrow(adjudication_inventory) ==
      nrow(candidate_pairs),
    all(
      paste(
        candidate_pairs$ligand,
        candidate_pairs$receptor
      ) %in%
        paste(
          lr$ligand,
          lr$receptor
        )
    ),
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

print(
  checks,
  n=Inf,
  width=Inf
)

stopifnot(
  all(checks$pass)
)

write_csv(
  checks,
  file.path(
    TAB,
    "21F2B0_formal_checks.csv"
  )
)

audit <- c(
  "PIPELINE_V2 — 21F2B0 biological LR-axis inventory",
  "",
  paste0("Date: ", Sys.Date()),
  "",
  "Purpose:",
  paste(
    "Generate a disease-blind candidate inventory of",
    "biologically interpretable LR axes before disease testing."
  ),
  "",
  "Twelve candidate biological axes were prespecified.",
  "",
  paste(
    "Regex rules are discovery aids only and do not constitute",
    "final biological adjudication."
  ),
  "",
  paste(
    "Every candidate LR pair must already exist in the frozen",
    "OmniPath symbol-level interaction universe."
  ),
  "",
  "No axis finalized.",
  "No disease filtering.",
  "No disease comparison.",
  "No communication refit.",
  "No DE refit.",
  "No pathway refit.",
  "",
  "PASS"
)

writeLines(
  audit,
  file.path(
    OUT,
    "21F2B0_audit.txt"
  )
)

cat("\n============================================================\n")
cat("21F2B0 COMPLETE\n")
cat("============================================================\n")
cat("Disease-blind biological LR candidate inventory generated\n")
cat("All candidate pairs anchored to frozen OmniPath resource\n")
cat("Regex used for discovery only\n")
cat("No biological axis finalized\n")
cat("No disease comparison\n")
cat("PASS\n")
