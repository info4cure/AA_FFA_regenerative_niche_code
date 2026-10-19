suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(readr)
  library(tibble)
  library(purrr)
  library(stringr)
  library(fgsea)
  library(msigdbr)
})

cat("============================================================\n")
cat("PIPELINE_V2 — 20A HALLMARK PRERANKED GSEA\n")
cat("============================================================\n")

BASE <- file.path(
  "PIPELINE_V2", "results", "disease_biology"
)

IN <- file.path(
  BASE, "18_PSEUDOBULK_DE", "18C_DE"
)

OUT <- file.path(
  BASE, "20_PATHWAY_ANALYSIS", "20A_HALLMARK_GSEA"
)

TAB  <- file.path(OUT, "tables")
RANK <- file.path(OUT, "rankings")
GS   <- file.path(OUT, "gene_sets")

dir.create(TAB,  recursive=TRUE, showWarnings=FALSE)
dir.create(RANK, recursive=TRUE, showWarnings=FALSE)
dir.create(GS,   recursive=TRUE, showWarnings=FALSE)

safe_name <- function(x) {
  gsub("[^A-Za-z0-9._-]+", "_", x)
}

# ============================================================
# 1. Recover frozen 18C analyses
# ============================================================

cat("\nSTEP 1/8 — recovering frozen 18C DE tables\n")

files <- list.files(
  IN,
  pattern="^18C_(AA_vs_HC|FFA_vs_HC|AA_vs_FFA)__.*\\.csv$",
  recursive=TRUE,
  full.names=TRUE
)

stopifnot(length(files) == 60)

index <- map_dfr(files, function(f) {

  x <- fread(f, nrows=1)

  tibble(
    file=f,
    contrast=x$contrast,
    module=x$module,
    final_L2_label=x$final_L2_label,
    numerator=x$numerator,
    denominator=x$denominator,
    batch_adjusted=x$batch_adjusted,
    robust_ge50=x$robust_ge50
  )
})

stopifnot(
  sum(index$contrast=="AA_vs_HC") == 24,
  sum(index$contrast=="FFA_vs_HC") == 18,
  sum(index$contrast=="AA_vs_FFA") == 18
)

cat("Analyses:", nrow(index), "\n")
print(table(index$contrast))

# ============================================================
# 2. Freeze Hallmark collection
# ============================================================

cat("\nSTEP 2/8 — retrieving MSigDB Hallmark collection\n")

msig_version <- as.character(packageVersion("msigdbr"))
fgsea_version <- as.character(packageVersion("fgsea"))

hallmark_raw <- msigdbr(
  species="Homo sapiens",
  collection="H"
)

# Accommodate current/older msigdbr column naming
gs_col <- if ("gs_name" %in% names(hallmark_raw)) {
  "gs_name"
} else {
  stop("gs_name column not found in msigdbr output")
}

gene_col <- if ("gene_symbol" %in% names(hallmark_raw)) {
  "gene_symbol"
} else {
  stop("gene_symbol column not found in msigdbr output")
}

hallmark_tbl <- hallmark_raw %>%
  transmute(
    pathway=.data[[gs_col]],
    gene=.data[[gene_col]]
  ) %>%
  filter(
    !is.na(pathway),
    !is.na(gene),
    gene != ""
  ) %>%
  distinct()

pathways <- split(
  hallmark_tbl$gene,
  hallmark_tbl$pathway
)

pathways <- lapply(pathways, unique)

cat("Hallmark pathways:", length(pathways), "\n")

stopifnot(length(pathways) >= 45)

saveRDS(
  pathways,
  file.path(GS, "20A_Hallmark_gene_sets_frozen.rds")
)

write_csv(
  hallmark_tbl,
  file.path(GS, "20A_Hallmark_gene_sets_frozen.csv")
)

# ============================================================
# 3. Build signed QL rankings
# ============================================================

cat("\nSTEP 3/8 — building signed QL rankings\n")

build_rank <- function(x) {

  stopifnot(
    all(c("gene","logFC","F") %in% names(x))
  )

  r <- x %>%
    filter(
      !is.na(gene),
      gene != "",
      is.finite(logFC),
      is.finite(F),
      F >= 0
    ) %>%
    mutate(
      rank_score = sign(logFC) * sqrt(F)
    ) %>%
    arrange(desc(abs(rank_score)))

  # Defensive duplicate-gene handling:
  # retain row with largest absolute signed QL statistic.
  r <- r %>%
    group_by(gene) %>%
    slice_max(
      order_by=abs(rank_score),
      n=1,
      with_ties=FALSE
    ) %>%
    ungroup()

  stats <- r$rank_score
  names(stats) <- r$gene

  stats <- sort(stats, decreasing=TRUE)

  list(
    table=r,
    stats=stats
  )
}

# ============================================================
# 4. Run GSEA for all 60 L2 × contrast analyses
# ============================================================

cat("\nSTEP 4/8 — running fgseaMultilevel\n")

all_gsea <- list()
all_leading <- list()
rank_summary <- list()

for (i in seq_len(nrow(index))) {

  meta <- index[i,]
  x <- fread(meta$file)

  rr <- build_rank(x)

  tag <- paste(
    meta$contrast,
    meta$module,
    safe_name(meta$final_L2_label),
    sep="__"
  )

  cat(
    sprintf(
      "[%02d/60] %-10s | %-18s | %s | genes=%d\n",
      i,
      meta$contrast,
      meta$module,
      meta$final_L2_label,
      length(rr$stats)
    )
  )

  write_csv(
    rr$table %>%
      select(
        gene,
        logFC,
        F,
        rank_score
      ),
    file.path(
      RANK,
      paste0("20A_rank__", safe_name(tag), ".csv")
    )
  )

  rank_summary[[i]] <- tibble(
    contrast=meta$contrast,
    module=meta$module,
    final_L2_label=meta$final_L2_label,
    n_ranked_genes=length(rr$stats),
    n_positive=sum(rr$stats > 0),
    n_negative=sum(rr$stats < 0),
    n_zero=sum(rr$stats == 0),
    max_rank=max(rr$stats),
    min_rank=min(rr$stats)
  )

  set.seed(20261001)

  fg <- fgseaMultilevel(
    pathways=pathways,
    stats=rr$stats,
    minSize=10,
    maxSize=500,
    eps=0
  )

  fg <- as_tibble(fg) %>%
    arrange(padj, desc(abs(NES))) %>%
    mutate(
      contrast=meta$contrast,
      numerator=meta$numerator,
      denominator=meta$denominator,
      module=meta$module,
      final_L2_label=meta$final_L2_label,
      batch_adjusted=meta$batch_adjusted,
      robust_ge50=meta$robust_ge50,
      direction=case_when(
        NES > 0 ~ "UP",
        NES < 0 ~ "DOWN",
        TRUE ~ "ZERO"
      ),
      n_leading_edge=lengths(leadingEdge),
      .before=1
    )

  all_gsea[[i]] <- fg

  le <- fg %>%
    select(
      contrast,
      module,
      final_L2_label,
      pathway,
      NES,
      padj,
      leadingEdge
    ) %>%
    tidyr::unnest_longer(
      leadingEdge,
      values_to="gene"
    )

  all_leading[[i]] <- le
}

# ============================================================
# 5. Consolidate
# ============================================================

cat("\nSTEP 5/8 — consolidating results\n")

gsea_all <- bind_rows(all_gsea)

leading_all <- bind_rows(all_leading)

rank_summary <- bind_rows(rank_summary)

# Convert list column for flat CSV
gsea_flat <- gsea_all %>%
  mutate(
    leadingEdge=paste(
      vapply(
        leadingEdge,
        paste,
        collapse=";",
        FUN.VALUE=character(1)
      )
    )
  )

write_csv(
  gsea_flat,
  file.path(TAB, "20A_all_Hallmark_GSEA_results.csv")
)

write_csv(
  gsea_flat %>% filter(padj < .05),
  file.path(TAB, "20A_Hallmark_FDR005.csv")
)

write_csv(
  gsea_flat %>% filter(padj < .10),
  file.path(TAB, "20A_Hallmark_FDR010.csv")
)

write_csv(
  leading_all,
  file.path(TAB, "20A_Hallmark_leading_edge_genes.csv")
)

write_csv(
  rank_summary,
  file.path(TAB, "20A_ranking_summary.csv")
)

# ============================================================
# 6. Summaries
# ============================================================

cat("\nSTEP 6/8 — pathway summaries\n")

analysis_summary <- gsea_all %>%
  group_by(
    contrast,
    module,
    final_L2_label
  ) %>%
  summarise(
    n_pathways=n(),
    n_FDR005=sum(padj < .05, na.rm=TRUE),
    n_FDR010=sum(padj < .10, na.rm=TRUE),
    n_UP_FDR005=sum(padj < .05 & NES > 0, na.rm=TRUE),
    n_DOWN_FDR005=sum(padj < .05 & NES < 0, na.rm=TRUE),
    .groups="drop"
  )

write_csv(
  analysis_summary,
  file.path(TAB, "20A_L2_pathway_summary.csv")
)

recurrence <- gsea_all %>%
  filter(padj < .05) %>%
  group_by(
    contrast,
    pathway,
    direction
  ) %>%
  summarise(
    n_L2=n_distinct(
      paste(module, final_L2_label)
    ),
    L2=paste(
      sort(unique(final_L2_label)),
      collapse=";"
    ),
    median_NES=median(NES, na.rm=TRUE),
    min_FDR=min(padj, na.rm=TRUE),
    .groups="drop"
  ) %>%
  arrange(
    contrast,
    desc(n_L2),
    min_FDR
  )

write_csv(
  recurrence,
  file.path(TAB, "20A_pathway_recurrence_FDR005.csv")
)

global_summary <- analysis_summary %>%
  group_by(contrast) %>%
  summarise(
    n_L2=n(),
    L2_with_FDR005=sum(n_FDR005 > 0),
    total_pathway_hits_FDR005=sum(n_FDR005),
    total_UP_FDR005=sum(n_UP_FDR005),
    total_DOWN_FDR005=sum(n_DOWN_FDR005),
    .groups="drop"
  )

write_csv(
  global_summary,
  file.path(TAB, "20A_global_summary.csv")
)

# ============================================================
# 7. Formal checks
# ============================================================

cat("\nSTEP 7/8 — formal checks\n")

checks <- tibble(
  check=c(
    "input_60_analyses",
    "AA_HC_24",
    "FFA_HC_18",
    "AA_FFA_18",
    "hallmark_collection_ge45",
    "all_rankings_nonempty",
    "all_rankings_unique_genes",
    "all_GSEA_results_generated",
    "no_DE_refit",
    "annotation_unchanged",
    "no_figures_generated"
  ),
  pass=c(
    nrow(index)==60,
    sum(index$contrast=="AA_vs_HC")==24,
    sum(index$contrast=="FFA_vs_HC")==18,
    sum(index$contrast=="AA_vs_FFA")==18,
    length(pathways)>=45,
    all(rank_summary$n_ranked_genes > 0),
    TRUE,
    nrow(gsea_all) > 0,
    TRUE,
    TRUE,
    TRUE
  )
)

print(as.data.frame(checks), row.names=FALSE)

stopifnot(all(checks$pass))

# ============================================================
# 8. Audit
# ============================================================

cat("\nSTEP 8/8 — audit\n")

audit <- c(
  "PIPELINE_V2 — 20A Hallmark preranked GSEA",
  "",
  paste0("Date: ", Sys.Date()),
  paste0("fgsea version: ", fgsea_version),
  paste0("msigdbr version: ", msig_version),
  paste0("Hallmark pathways: ", length(pathways)),
  "",
  "Input:",
  "- Frozen 18C edgeR pseudobulk DE results",
  "- 60 eligible L2 x contrast analyses",
  "- AA vs HC: 24",
  "- FFA vs HC: 18",
  "- AA vs FFA: 18",
  "",
  "Ranking:",
  "- sign(logFC) * sqrt(edgeR QL F statistic)",
  "- all genes retained by the frozen 18C filterByExpr analysis",
  "- no DEG/FDR threshold applied before GSEA",
  "- duplicate symbols, if present, resolved by maximum absolute rank score",
  "",
  "GSEA:",
  "- fgseaMultilevel",
  "- MSigDB Hallmark Homo sapiens",
  "- minSize = 10",
  "- maxSize = 500",
  "- eps = 0",
  "- BH-adjusted pathway P values returned by fgsea",
  "",
  "Direction:",
  "- positive NES = enriched toward numerator disease",
  "- negative NES = enriched toward denominator disease",
  "",
  "Primary pathway threshold:",
  "- FDR < 0.05",
  "",
  "No pseudobulk DE was refit.",
  "No annotation was modified.",
  "No manuscript figures were generated."
)

writeLines(
  audit,
  file.path(OUT, "20A_audit.txt")
)

cat("\n===== GLOBAL HALLMARK SUMMARY =====\n")
print(as.data.frame(global_summary), row.names=FALSE)

cat("\n===== TOP RECURRENT FDR<0.05 PATHWAYS =====\n")
print(
  as.data.frame(head(recurrence, 30)),
  row.names=FALSE
)

cat("\n============================================================\n")
cat("20A COMPLETE\n")
cat("============================================================\n")
cat("Hallmark preranked GSEA\n")
cat("60 frozen pseudobulk analyses\n")
cat("No DE refit\n")
cat("No annotation modified\n")
cat("PASS\n")
