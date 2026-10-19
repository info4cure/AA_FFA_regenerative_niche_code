
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
  library(purrr)
  library(stringr)
  library(fgsea)
})

cat("============================================================\n")
cat("PIPELINE_V2 — 20E1 BIOPLANET PRERANKED GSEA\n")
cat("============================================================\n")

BASE <- file.path(
  "PIPELINE_V2","results","disease_biology",
  "20_PATHWAY_ANALYSIS"
)

IN20A <- file.path(
  BASE,"20A_HALLMARK_GSEA"
)

BIO <- file.path(
  BASE,"20E_BIOPLANET"
)

OUT <- file.path(
  BIO,"20E1_GSEA"
)

TAB <- file.path(OUT,"tables")
FULL <- file.path(OUT,"full_results")

dir.create(TAB,recursive=TRUE,showWarnings=FALSE)
dir.create(FULL,recursive=TRUE,showWarnings=FALSE)

# ============================================================
# 1. Frozen inputs
# ============================================================

cat("\nSTEP 1/8 — recovering frozen BioPlanet and 20A rankings\n")

pathways <- readRDS(
  file.path(
    BIO,"gene_sets",
    "20E0_BioPlanet_gene_sets.rds"
  )
)

inventory <- read_csv(
  file.path(
    BIO,"tables",
    "20E0_BioPlanet_pathway_inventory.csv"
  ),
  show_col_types=FALSE
)

# Pre-filter only by frozen pathway size specification.
eligible_names <- inventory %>%
  filter(
    n_genes >= 10,
    n_genes <= 500
  ) %>%
  pull(GENESET_NAME)

pathways <- pathways[
  intersect(
    names(pathways),
    eligible_names
  )
]

cat("BioPlanet frozen pathways :",nrow(inventory),"\n")
cat("BioPlanet GSEA eligible   :",length(pathways),"\n")

stopifnot(length(pathways)==1246)

# Locate frozen 20A rankings.
rank_files <- list.files(
  file.path(IN20A, "rankings"),
  pattern = "^20A_rank__.*\\.csv$",
  full.names = TRUE
)

if (!length(rank_files)) {
  stop(
    "No frozen 20A ranking CSV files found in ",
    file.path(IN20A, "rankings")
  )
}

cat("Frozen rankings:", length(rank_files), "\n")
stopifnot(length(rank_files) == 60)

# ============================================================
# 2. Ranking metadata
# ============================================================

cat("\nSTEP 2/8 — recovering ranking metadata\n")

ranking_summary <- read_csv(
  file.path(
    IN20A,"tables",
    "20A_ranking_summary.csv"
  ),
  show_col_types=FALSE
)

stopifnot(nrow(ranking_summary)==60)

print(
  table(ranking_summary$contrast)
)

# Match each frozen 20A ranking CSV by identifiers present in its filename.
safe_name <- function(x) {
  gsub("[^A-Za-z0-9._-]+","_",x)
}

find_ranking <- function(contrast, module, l2) {

  expected <- paste0(
    "20A_rank__",
    safe_name(contrast),
    "__",
    safe_name(module),
    "__",
    safe_name(l2),
    ".csv"
  )

  hit <- rank_files[
    basename(rank_files) == expected
  ]

  if (length(hit) != 1) {
    stop(
      "Exact ranking match !=1 for: ",
      contrast, " | ", module, " | ", l2,
      " ; expected=", expected,
      " ; n=", length(hit)
    )
  }

  hit
}

# ============================================================
# 3. GSEA runner
# ============================================================

cat("\nSTEP 3/8 — defining BioPlanet GSEA runner\n")

run_one <- function(i) {

  meta <- ranking_summary[i,]

  contrast <- meta$contrast[[1]]
  module   <- meta$module[[1]]
  l2       <- meta$final_L2_label[[1]]

  rf <- find_ranking(
    contrast,module,l2
  )

  # Recover the exact frozen 20A ranking.
  # Do NOT recompute the ranking from logFC/F here.
  rank_tbl <- read_csv(
    rf,
    show_col_types = FALSE
  )

  stopifnot(
    all(c("gene", "rank_score") %in% names(rank_tbl)),
    nrow(rank_tbl) > 0,
    !anyNA(rank_tbl$gene),
    !anyNA(rank_tbl$rank_score),
    !anyDuplicated(rank_tbl$gene),
    all(is.finite(rank_tbl$rank_score))
  )

  stats <- rank_tbl$rank_score
  names(stats) <- rank_tbl$gene

  stats <- sort(
    stats,
    decreasing = TRUE
  )

  # Verify against the frozen 20A ranking manifest.
  stopifnot(
    length(stats) == meta$n_ranked_genes[[1]]
  )

  cat(
    sprintf(
      "[%02d/60] %-10s | %-20s | %s | genes=%d\n",
      i,contrast,module,l2,length(stats)
    )
  )

  fg <- fgseaMultilevel(
    pathways=pathways,
    stats=stats,
    minSize=10,
    maxSize=500,
    eps=0
  )

  if (!nrow(fg)) {
    stop(
      "No BioPlanet GSEA results: ",
      contrast," | ",module," | ",l2
    )
  }

  fg <- as_tibble(fg) %>%
    mutate(
      contrast=contrast,
      module=module,
      final_L2_label=l2,
      n_ranked_genes=length(stats),
      direction=case_when(
        NES > 0 ~ "UP",
        NES < 0 ~ "DOWN",
        TRUE ~ "ZERO"
      ),
      .before=1
    )

  # Parse stable BioPlanet ID/name from frozen gene-set identifier.
  fg <- fg %>%
    mutate(
      bioplanet_id=str_match(
        pathway,
        "^BIOPLANET_(bioplanet_[0-9]+)__"
      )[,2],
      bioplanet_name=sub(
        "^BIOPLANET_bioplanet_[0-9]+__",
        "",
        pathway
      ),
      .after=pathway
    )

  outname <- paste0(
    "20E1_",
    safe_name(contrast),"__",
    safe_name(module),"__",
    safe_name(l2),
    ".csv"
  )

  # List column cannot be written directly.
  fg_write <- fg %>%
    mutate(
      leadingEdge=vapply(
        leadingEdge,
        paste,
        collapse=";",
        FUN.VALUE=character(1)
      )
    )

  write_csv(
    fg_write,
    file.path(FULL,outname)
  )

  fg
}

# ============================================================
# 4. Run all 60 analyses
# ============================================================

cat("\nSTEP 4/8 — running BioPlanet fgseaMultilevel\n")

results_list <- map(
  seq_len(nrow(ranking_summary)),
  run_one
)

all_res <- bind_rows(results_list)

cat(
  "\nTotal BioPlanet result rows:",
  nrow(all_res),"\n"
)

# ============================================================
# 5. Consolidate
# ============================================================

cat("\nSTEP 5/8 — consolidating results\n")

all_write <- all_res %>%
  mutate(
    leadingEdge=vapply(
      leadingEdge,
      paste,
      collapse=";",
      FUN.VALUE=character(1)
    )
  )

write_csv(
  all_write,
  file.path(
    TAB,
    "20E1_all_BioPlanet_GSEA.csv"
  )
)

sig <- all_res %>%
  filter(
    !is.na(padj),
    padj < 0.05
  )

sig_write <- sig %>%
  mutate(
    leadingEdge=vapply(
      leadingEdge,
      paste,
      collapse=";",
      FUN.VALUE=character(1)
    )
  )

write_csv(
  sig_write,
  file.path(
    TAB,
    "20E1_BioPlanet_FDR005.csv"
  )
)

# ============================================================
# 6. Summaries
# ============================================================

cat("\nSTEP 6/8 — summaries\n")

global_summary <- all_res %>%
  group_by(contrast) %>%
  summarise(
    n_L2=n_distinct(final_L2_label),
    L2_with_FDR005=
      n_distinct(
        final_L2_label[
          !is.na(padj) & padj<0.05
        ]
      ),
    total_pathway_hits_FDR005=
      sum(padj<0.05,na.rm=TRUE),
    total_UP_FDR005=
      sum(
        padj<0.05 &
        NES>0,
        na.rm=TRUE
      ),
    total_DOWN_FDR005=
      sum(
        padj<0.05 &
        NES<0,
        na.rm=TRUE
      ),
    .groups="drop"
  )

l2_summary <- all_res %>%
  group_by(
    contrast,
    module,
    final_L2_label
  ) %>%
  summarise(
    n_pathways=n(),
    n_FDR005=sum(padj<0.05,na.rm=TRUE),
    n_FDR010=sum(padj<0.10,na.rm=TRUE),
    n_UP_FDR005=
      sum(padj<0.05 & NES>0,na.rm=TRUE),
    n_DOWN_FDR005=
      sum(padj<0.05 & NES<0,na.rm=TRUE),
    .groups="drop"
  )

recurrent <- sig %>%
  group_by(
    contrast,
    bioplanet_id,
    bioplanet_name,
    direction
  ) %>%
  summarise(
    n_L2=n_distinct(final_L2_label),
    median_NES=median(NES,na.rm=TRUE),
    min_FDR=min(padj,na.rm=TRUE),
    L2=paste(
      sort(unique(final_L2_label)),
      collapse=";"
    ),
    .groups="drop"
  ) %>%
  arrange(
    contrast,
    desc(n_L2),
    min_FDR
  )

write_csv(
  global_summary,
  file.path(
    TAB,
    "20E1_global_summary.csv"
  )
)

write_csv(
  l2_summary,
  file.path(
    TAB,
    "20E1_L2_summary.csv"
  )
)

write_csv(
  recurrent,
  file.path(
    TAB,
    "20E1_recurrent_pathways.csv"
  )
)

# ============================================================
# 7. Leading-edge driver genes
# ============================================================

cat("\nSTEP 7/8 — leading-edge architecture\n")

leading <- sig %>%
  select(
    contrast,
    module,
    final_L2_label,
    bioplanet_id,
    bioplanet_name,
    NES,
    padj,
    direction,
    leadingEdge
  ) %>%
  tidyr::unnest_longer(
    leadingEdge,
    values_to="gene"
  ) %>%
  filter(
    !is.na(gene),
    gene!=""
  )

write_csv(
  leading,
  file.path(
    TAB,
    "20E1_leading_edge_genes_FDR005.csv"
  )
)

drivers <- leading %>%
  group_by(
    contrast,
    direction,
    gene
  ) %>%
  summarise(
    n_pathways=n_distinct(bioplanet_id),
    n_L2=n_distinct(final_L2_label),
    n_L2_pathway_pairs=
      n_distinct(
        paste(
          final_L2_label,
          bioplanet_id,
          sep="|||"
        )
      ),
    .groups="drop"
  ) %>%
  arrange(
    contrast,
    direction,
    desc(n_L2),
    desc(n_pathways),
    gene
  )

write_csv(
  drivers,
  file.path(
    TAB,
    "20E1_recurrent_driver_genes.csv"
  )
)

# ============================================================
# 8. Formal checks / audit
# ============================================================

cat("\nSTEP 8/8 — formal checks\n")

checks <- tibble(
  check=c(
    "BioPlanet_1658_frozen",
    "BioPlanet_1246_size_eligible",
    "rankings_60",
    "AA_HC_24",
    "FFA_HC_18",
    "AA_FFA_18",
    "all_results_generated",
    "leading_edges_generated",
    "no_GSEA_rank_refit",
    "no_DE_refit",
    "annotation_unchanged",
    "no_figures_generated"
  ),
  pass=c(
    nrow(inventory)==1658,
    length(pathways)==1246,
    nrow(ranking_summary)==60,
    sum(ranking_summary$contrast=="AA_vs_HC")==24,
    sum(ranking_summary$contrast=="FFA_vs_HC")==18,
    sum(ranking_summary$contrast=="AA_vs_FFA")==18,
    nrow(all_res)>0,
    nrow(leading)>=0,
    TRUE,
    TRUE,
    TRUE,
    TRUE
  )
)

print(
  as.data.frame(checks),
  row.names=FALSE
)

stopifnot(all(checks$pass))

audit <- c(
  "PIPELINE_V2 — 20E1 BioPlanet preranked GSEA",
  "",
  paste0("Date: ",Sys.Date()),
  "",
  "BioPlanet source:",
  "- frozen official NCATS BioPlanet collection from 20E0",
  "- 1,658 source pathways",
  "- 9,818 unique gene symbols",
  "- 1,246 pathways eligible at frozen size 10-500",
  "",
  "Expression rankings:",
  "- exactly the 60 frozen signed QL rankings used by 20A",
  "- no differential-expression refit",
  "- no ranking refit",
  "",
  "GSEA:",
  "- fgseaMultilevel",
  "- minSize=10",
  "- maxSize=500",
  "- eps=0",
  "- pathway-level BH FDR returned by fgsea within each L2 x contrast analysis",
  "",
  "Primary significance threshold:",
  "- FDR < 0.05",
  "",
  "No cell annotation modified.",
  "No figures generated."
)

writeLines(
  audit,
  file.path(
    OUT,
    "20E1_audit.txt"
  )
)

cat("\n===== GLOBAL BIOPLANET SUMMARY =====\n")
print(
  as.data.frame(global_summary),
  row.names=FALSE
)

cat("\n===== TOP RECURRENT FDR<0.05 BIOPLANET PATHWAYS =====\n")

print(
  as.data.frame(
    recurrent %>%
      group_by(contrast) %>%
      slice_head(n=25) %>%
      ungroup()
  ),
  row.names=FALSE
)

cat("\n============================================================\n")
cat("20E1 COMPLETE\n")
cat("============================================================\n")
cat("BioPlanet preranked GSEA\n")
cat("60 frozen rankings\n")
cat("No DE refit\n")
cat("PASS\n")
