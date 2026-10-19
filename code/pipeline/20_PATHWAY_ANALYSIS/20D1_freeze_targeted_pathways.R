suppressPackageStartupMessages({
  library(msigdbr)
  library(dplyr)
  library(readr)
  library(tibble)
  library(stringr)
  library(purrr)
})

cat("============================================================\n")
cat("PIPELINE_V2 — 20D1 FREEZE REACTOME / GO-BP PATHWAYS\n")
cat("============================================================\n")

BASE <- file.path(
  "PIPELINE_V2","results","disease_biology",
  "20_PATHWAY_ANALYSIS"
)

OUT <- file.path(
  BASE,"20D_TARGETED_PATHWAYS"
)

GS <- file.path(OUT,"gene_sets")
TAB <- file.path(OUT,"tables")

dir.create(GS,recursive=TRUE,showWarnings=FALSE)
dir.create(TAB,recursive=TRUE,showWarnings=FALSE)

# ============================================================
# 1. Retrieve exact MSigDB collections
# ============================================================

cat("\nSTEP 1/6 — retrieving MSigDB 2026.1.Hs collections\n")

reactome <- msigdbr(
  species="Homo sapiens",
  collection="C2",
  subcollection="CP:REACTOME"
)

gobp <- msigdbr(
  species="Homo sapiens",
  collection="C5",
  subcollection="GO:BP"
)

stopifnot(
  length(unique(reactome$gs_name))==1839,
  length(unique(gobp$gs_name))==7538
)

cat("Reactome:",n_distinct(reactome$gs_name),"\n")
cat("GO:BP   :",n_distinct(gobp$gs_name),"\n")

# ============================================================
# 2. Freeze complete collections
# ============================================================

cat("\nSTEP 2/6 — freezing complete collections\n")

R <- reactome %>%
  transmute(
    database="REACTOME",
    pathway=gs_name,
    gene=gene_symbol
  ) %>%
  distinct()

G <- gobp %>%
  transmute(
    database="GO_BP",
    pathway=gs_name,
    gene=gene_symbol
  ) %>%
  distinct()

write_csv(
  R,
  file.path(GS,"20D1_Reactome_2026.1.Hs.csv")
)

write_csv(
  G,
  file.path(GS,"20D1_GO_BP_2026.1.Hs.csv")
)

saveRDS(
  split(R$gene,R$pathway),
  file.path(GS,"20D1_Reactome_2026.1.Hs.rds")
)

saveRDS(
  split(G$gene,G$pathway),
  file.path(GS,"20D1_GO_BP_2026.1.Hs.rds")
)

# ============================================================
# 3. Predefine biological axes
# ============================================================

cat("\nSTEP 3/6 — defining biological axes\n")

axes <- tribble(
  ~axis, ~pattern,

  "IFN_JAK_STAT_ANTIGEN",
  "INTERFERON|JAK_STAT|JAK.STAT|STAT_SIGNAL|ANTIGEN_PROCESS|ANTIGEN_PRESENT|MHC_CLASS|MAJOR_HISTOCOMPAT",

  "T_CELL_CYTOTOXIC",
  "T_CELL|T.CELL|CYTOTOX|NATURAL_KILLER|NK_CELL|TCR|T_CELL_RECEPTOR|LYMPHOCYTE_MEDIATED",

  "ECM_COLLAGEN_FIBROSIS",
  "EXTRACELLULAR_MATRIX|ECM_|COLLAGEN|FIBRO|MATRIX_ORGAN|MATRIX_ASSEMB|MATRIX_DEGRAD|INTEGRIN|TGF_BETA|TGFB",

  "WNT_HEDGEHOG_FOLLICLE",
  "WNT|BETA_CATENIN|CTNNB|HEDGEHOG|HAIR_FOLLIC|HAIR_CYCLE|PILOSEBAC|EPIDERMAL_APPENDAGE",

  "EPITHELIAL_DIFFERENTIATION",
  "KERATIN|CORNIF|EPIDERM|EPITHELIAL_DIFFERENT|CELL_JUNCTION|TIGHT_JUNCTION|ADHERENS|DESMOSOM|EPITHELIAL_CELL",

  "VASCULAR_ANGIOGENESIS",
  "ANGIOGEN|VASCULAR|ENDOTHEL|BLOOD_VESSEL|VEGF|VASCULOGEN",

  "METABOLISM_OXPHOS_ROS",
  "OXIDATIVE_PHOSPHORYL|RESPIRATORY_CHAIN|MITOCHONDR|GLYCOLY|REACTIVE_OXYGEN|OXIDATIVE_STRESS|FATTY_ACID|LIPID_METAB",

  "CELL_CYCLE_DNA_REPAIR",
  "CELL_CYCLE|MITOTIC|MITOSIS|DNA_REPAIR|DNA_DAMAGE|DNA_REPLICATION|G1_S|G2_M|CHECKPOINT"
)

write_csv(
  axes,
  file.path(TAB,"20D1_axis_definitions.csv")
)

# ============================================================
# 4. Map pathways to axes by pathway name
# ============================================================

cat("\nSTEP 4/6 — mapping pathways to axes\n")

inventory <- bind_rows(
  R %>% distinct(database,pathway),
  G %>% distinct(database,pathway)
)

axis_map <- map_dfr(
  seq_len(nrow(axes)),
  function(i) {

    inventory %>%
      filter(
        str_detect(
          pathway,
          regex(
            axes$pattern[i],
            ignore_case=TRUE
          )
        )
      ) %>%
      mutate(
        axis=axes$axis[i]
      )
  }
) %>%
  distinct(
    database,
    pathway,
    axis
  )

write_csv(
  axis_map,
  file.path(TAB,"20D1_targeted_pathway_inventory.csv")
)

# ============================================================
# 5. Coverage summaries
# ============================================================

cat("\nSTEP 5/6 — coverage summaries\n")

coverage <- axis_map %>%
  count(
    database,
    axis,
    name="n_pathways"
  ) %>%
  arrange(database,axis)

write_csv(
  coverage,
  file.path(TAB,"20D1_axis_coverage.csv")
)

multi_axis <- axis_map %>%
  count(
    database,
    pathway,
    name="n_axes"
  ) %>%
  filter(n_axes>1) %>%
  arrange(desc(n_axes),pathway)

write_csv(
  multi_axis,
  file.path(TAB,"20D1_multi_axis_pathways.csv")
)

cat("\n===== AXIS COVERAGE =====\n")
print(as.data.frame(coverage),row.names=FALSE)

cat("\nPathways assigned to >1 axis:",
    nrow(multi_axis),"\n")

# ============================================================
# 6. Formal checks / audit
# ============================================================

cat("\nSTEP 6/6 — formal checks\n")

checks <- tibble(
  check=c(
    "Reactome_1839",
    "GO_BP_7538",
    "eight_axes",
    "all_axes_have_pathways",
    "collections_frozen",
    "inventory_generated",
    "no_GSEA_performed",
    "no_DE_refit",
    "annotation_unchanged"
  ),
  pass=c(
    n_distinct(R$pathway)==1839,
    n_distinct(G$pathway)==7538,
    nrow(axes)==8,
    all(axes$axis %in% coverage$axis),
    file.exists(
      file.path(GS,"20D1_Reactome_2026.1.Hs.rds")
    ) &&
      file.exists(
        file.path(GS,"20D1_GO_BP_2026.1.Hs.rds")
      ),
    nrow(axis_map)>0,
    TRUE,
    TRUE,
    TRUE
  )
)

print(as.data.frame(checks),row.names=FALSE)
stopifnot(all(checks$pass))

audit <- c(
  "PIPELINE_V2 — 20D1 Reactome / GO-BP pathway freeze",
  "",
  paste0("Date: ",Sys.Date()),
  paste0(
    "msigdbr version: ",
    as.character(packageVersion("msigdbr"))
  ),
  "MSigDB version: 2026.1.Hs",
  "",
  "Complete collections frozen before enrichment:",
  "- Reactome: 1,839 pathways",
  "- GO Biological Process: 7,538 pathways",
  "",
  "Eight interpretation axes were predefined by pathway-name patterns:",
  "- IFN/JAK-STAT/antigen presentation",
  "- T-cell/cytotoxic immunity",
  "- ECM/collagen/fibrosis",
  "- WNT/Hedgehog/hair follicle",
  "- epithelial differentiation/junction",
  "- vascular/angiogenesis",
  "- metabolism/OXPHOS/ROS",
  "- cell cycle/DNA repair",
  "",
  "The complete Reactome and GO-BP collections are retained.",
  "Axis assignment is for targeted interpretation and does not determine",
  "which pathways will be tested in full enrichment.",
  "",
  "No GSEA performed in 20D1.",
  "No differential expression refit.",
  "No annotation modified."
)

writeLines(
  audit,
  file.path(OUT,"20D1_audit.txt")
)

cat("\n============================================================\n")
cat("20D1 COMPLETE\n")
cat("============================================================\n")
cat("Reactome / GO-BP collections frozen\n")
cat("Targeted interpretation inventory frozen\n")
cat("No enrichment performed\n")
cat("PASS\n")
