
# ============================================================
# 22F4F6A — FUNCTIONAL ANNOTATION OF FROZEN K5 MODULES
# ORA Hallmark + GO Biological Process
#
# Query universe:
#   2774 frozen dynamic genes
#
# IMPORTANT:
#   No new trajectory inference
#   No new DE testing
#   No module reassignment
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(msigdbr)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
})

cat("\n============================================================\n")
cat("22F4F6A — FROZEN K5 FUNCTIONAL ANNOTATION\n")
cat("============================================================\n")

# ============================================================
# 1. PATHS
# ============================================================

BASE <- file.path(
  "PIPELINE_V2",
  "results",
  "disease_biology",
  "22_MULTICOMPARTMENT_DISEASE_AXES",
  "22F_BULGE_STEM_NICHE_INTEGRITY",
  "22F4_BULGE_MONOCLE3_TRAJECTORIES"
)

D3_FILE <- file.path(
  BASE,
  "22F4D_DYNAMIC_EXPRESSION",
  "22F4D3B_K5_FREEZE",
  "objects",
  "22F4D3B_K5_dynamic_modules_freeze.rds"
)

OUT <- file.path(
  BASE,
  "22F4F_PUBLICATION_FIGURES",
  "22F4F6_INTEGRATED_TRAJECTORY"
)

TAB <- file.path(
  OUT,
  "tables"
)

OBJ <- file.path(
  OUT,
  "objects"
)

AUDIT <- file.path(
  OUT,
  "audit"
)

for (d in c(
  TAB,
  OBJ,
  AUDIT
)) {
  dir.create(
    d,
    recursive=TRUE,
    showWarnings=FALSE
  )
}

stopifnot(
  file.exists(D3_FILE)
)

# ============================================================
# 2. LOAD FROZEN MODULES
# ============================================================

cat("\n===== LOAD FROZEN K5 =====\n")

D3 <- readRDS(
  D3_FILE
)

assignments <- D3$gene_module_assignments |>
  mutate(
    gene=as.character(gene),
    module=as.character(module)
  )

stopifnot(
  nrow(assignments) == 2774,
  D3$K_final == 5,
  all(
    paste0("M",1:5) %in%
      unique(assignments$module)
  )
)

universe_symbol <- unique(
  assignments$gene
)

cat(
  "Dynamic-gene universe:",
  length(universe_symbol),
  "\n"
)

cat("\nGenes/module:\n")

print(
  assignments |>
    count(
      module,
      name="n_genes"
    ) |>
    arrange(module),
  n=Inf
)

# ============================================================
# 3. HALLMARK GENE SETS
# ============================================================

cat("\n===== HALLMARK GENE SETS =====\n")

hallmark_raw <- msigdbr(
  species="Homo sapiens",
  collection="H"
)

# msigdbr versions differ slightly in column names.
# Resolve robustly.

gene_col <- intersect(
  c(
    "gene_symbol",
    "human_gene_symbol"
  ),
  colnames(hallmark_raw)
)[1]

set_col <- intersect(
  c(
    "gs_name"
  ),
  colnames(hallmark_raw)
)[1]

if (
  is.na(gene_col) ||
  is.na(set_col)
) {
  stop(
    "Could not identify msigdbr gene/set columns."
  )
}

hallmark <- hallmark_raw |>
  transmute(
    term=.data[[set_col]],
    gene=.data[[gene_col]]
  ) |>
  filter(
    gene %in% universe_symbol
  ) |>
  distinct()

cat(
  "Hallmark sets represented:",
  n_distinct(hallmark$term),
  "\n"
)

cat(
  "Hallmark universe genes:",
  n_distinct(hallmark$gene),
  "\n"
)

# ============================================================
# 4. GENERIC ORA FUNCTION
# ============================================================

run_ora <- function(
  query_genes,
  universe_genes,
  term2gene
) {

  enr <- clusterProfiler::enricher(
    gene=query_genes,
    universe=universe_genes,
    TERM2GENE=term2gene,
    pAdjustMethod="BH",
    pvalueCutoff=1,
    qvalueCutoff=1,
    minGSSize=5,
    maxGSSize=5000
  )

  if (
    is.null(enr) ||
    nrow(as.data.frame(enr)) == 0
  ) {
    return(
      tibble()
    )
  }

  as.data.frame(enr) |>
    as_tibble()
}

# ============================================================
# 5. HALLMARK ORA
# ============================================================

cat("\n===== HALLMARK ORA =====\n")

hallmark_results <- lapply(
  paste0("M",1:5),
  function(m) {

    genes <- assignments |>
      filter(
        module == m
      ) |>
      pull(gene) |>
      unique()

    x <- run_ora(
      query_genes=genes,
      universe_genes=universe_symbol,
      term2gene=hallmark
    )

    if (nrow(x) == 0) {
      return(
        tibble(
          module=character()
        )
      )
    }

    x |>
      mutate(
        module=m,
        collection="HALLMARK",
        .before=1
      )
  }
) |>
  bind_rows()

if (nrow(hallmark_results) > 0) {

  hallmark_results <- hallmark_results |>
    mutate(
      pathway=gsub(
        "^HALLMARK_",
        "",
        ID
      ),
      pathway=gsub(
        "_",
        " ",
        pathway
      )
    ) |>
    arrange(
      module,
      p.adjust,
      desc(Count)
    )
}

write.csv(
  hallmark_results,
  file.path(
    TAB,
    "22F4F6A_K5_Hallmark_ORA.csv"
  ),
  row.names=FALSE
)

cat("\nTop Hallmark results:\n")

if (nrow(hallmark_results) > 0) {
  print(
    hallmark_results |>
      group_by(module) |>
      slice_head(n=10) |>
      ungroup() |>
      dplyr::select(
        module,
        pathway,
        Count,
        GeneRatio,
        p.adjust
      ),
    n=Inf
  )
}

# ============================================================
# 6. SYMBOL -> ENTREZ
# ============================================================

cat("\n===== ENTREZ MAPPING =====\n")

mapping <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys=universe_symbol,
  keytype="SYMBOL",
  columns=c(
    "SYMBOL",
    "ENTREZID"
  )
) |>
  as_tibble() |>
  filter(
    !is.na(ENTREZID)
  ) |>
  distinct(
    SYMBOL,
    ENTREZID
  )

cat(
  "Mapped symbols:",
  n_distinct(mapping$SYMBOL),
  "/",
  length(universe_symbol),
  "\n"
)

universe_entrez <- unique(
  mapping$ENTREZID
)

# ============================================================
# 7. GO BIOLOGICAL PROCESS ORA
# ============================================================

cat("\n===== GO BIOLOGICAL PROCESS ORA =====\n")

go_results <- lapply(
  paste0("M",1:5),
  function(m) {

    symbols <- assignments |>
      filter(
        module == m
      ) |>
      pull(gene) |>
      unique()

    entrez <- mapping |>
      filter(
        SYMBOL %in% symbols
      ) |>
      pull(ENTREZID) |>
      unique()

    if (length(entrez) < 5) {
      return(
        tibble()
      )
    }

    ego <- clusterProfiler::enrichGO(
      gene=entrez,
      universe=universe_entrez,
      OrgDb=org.Hs.eg.db,
      keyType="ENTREZID",
      ont="BP",
      pAdjustMethod="BH",
      pvalueCutoff=1,
      qvalueCutoff=1,
      minGSSize=10,
      maxGSSize=500,
      readable=TRUE
    )

    if (
      is.null(ego) ||
      nrow(as.data.frame(ego)) == 0
    ) {
      return(
        tibble()
      )
    }

    as.data.frame(ego) |>
      as_tibble() |>
      mutate(
        module=m,
        collection="GO_BP",
        .before=1
      )
  }
) |>
  bind_rows()

if (nrow(go_results) > 0) {

  go_results <- go_results |>
    arrange(
      module,
      p.adjust,
      desc(Count)
    )
}

write.csv(
  go_results,
  file.path(
    TAB,
    "22F4F6A_K5_GO_BP_ORA.csv"
  ),
  row.names=FALSE
)

cat("\nTop GO:BP results:\n")

if (nrow(go_results) > 0) {

  print(
    go_results |>
      group_by(module) |>
      slice_head(n=10) |>
      ungroup() |>
      dplyr::select(
        module,
        Description,
        Count,
        GeneRatio,
        p.adjust
      ),
    n=Inf
  )
}

# ============================================================
# 8. DISPLAY ANNOTATION CANDIDATES
# ============================================================

cat("\n===== BUILD DISPLAY ANNOTATION CANDIDATES =====\n")

hallmark_display <- hallmark_results |>
  filter(
    p.adjust < 0.10
  ) |>
  group_by(module) |>
  arrange(
    p.adjust,
    desc(Count),
    .by_group=TRUE
  ) |>
  slice_head(
    n=3
  ) |>
  ungroup() |>
  transmute(
    module,
    collection,
    term=pathway,
    FDR=p.adjust,
    Count
  )

go_display <- go_results |>
  filter(
    p.adjust < 0.05
  ) |>
  group_by(module) |>
  arrange(
    p.adjust,
    desc(Count),
    .by_group=TRUE
  ) |>
  slice_head(
    n=5
  ) |>
  ungroup() |>
  transmute(
    module,
    collection,
    term=Description,
    FDR=p.adjust,
    Count
  )

display_candidates <- bind_rows(
  hallmark_display,
  go_display
) |>
  arrange(
    module,
    collection,
    FDR
  )

write.csv(
  display_candidates,
  file.path(
    TAB,
    "22F4F6A_K5_functional_annotation_candidates.csv"
  ),
  row.names=FALSE
)

cat("\nDISPLAY CANDIDATES:\n")

print(
  display_candidates,
  n=Inf
)

# ============================================================
# 9. MODULE SHAPE SUMMARY
# ============================================================

shape_summary <- D3$module_summary |>
  dplyr::select(
    module,
    n_genes,
    centroid_peak_position,
    centroid_trough_position,
    centroid_amplitude
  ) |>
  arrange(module)

write.csv(
  shape_summary,
  file.path(
    TAB,
    "22F4F6A_K5_shape_summary.csv"
  ),
  row.names=FALSE
)

cat("\n===== MODULE SHAPES =====\n")

print(
  shape_summary,
  n=Inf
)

# ============================================================
# 10. FREEZE
# ============================================================

freeze <- list(
  module="22F4F6A",
  method="ORA",
  universe_definition=
    "2774 frozen 22F4D3B dynamic genes",
  K_final=D3$K_final,
  assignments=assignments,
  module_summary=D3$module_summary,
  hallmark_results=hallmark_results,
  go_bp_results=go_results,
  display_candidates=display_candidates
)

saveRDS(
  freeze,
  file.path(
    OBJ,
    "22F4F6A_K5_functional_annotation_freeze.rds"
  )
)

# ============================================================
# 11. FORMAL CHECKS
# ============================================================

checks <- tibble(
  check=c(
    "dynamic universe = 2774",
    "K = 5",
    "all M1-M5 represented",
    "Hallmark result table generated",
    "GO BP result table generated",
    "freeze exists"
  ),
  pass=c(
    length(universe_symbol) == 2774,
    D3$K_final == 5,
    all(
      paste0("M",1:5) %in%
        unique(assignments$module)
    ),
    file.exists(
      file.path(
        TAB,
        "22F4F6A_K5_Hallmark_ORA.csv"
      )
    ),
    file.exists(
      file.path(
        TAB,
        "22F4F6A_K5_GO_BP_ORA.csv"
      )
    ),
    file.exists(
      file.path(
        OBJ,
        "22F4F6A_K5_functional_annotation_freeze.rds"
      )
    )
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
    AUDIT,
    "22F4F6A_formal_checks.csv"
  ),
  row.names=FALSE
)

if (!all(checks$pass)) {
  stop(
    "22F4F6A formal checks failed."
  )
}

cat("\n============================================================\n")
cat("22F4F6A COMPLETE\n")
cat("STATUS=PASS\n")
cat("============================================================\n")

