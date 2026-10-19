suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(readr)
  library(tibble)
  library(purrr)
  library(stringr)
  library(tidyr)
  library(fgsea)
})

cat("============================================================\n")
cat("PIPELINE_V2 — 20B2 FFA HALLMARK BATCH SENSITIVITY\n")
cat("============================================================\n")

BASE <- file.path(
  "PIPELINE_V2","results","disease_biology"
)

IN18 <- file.path(
  BASE,"18_PSEUDOBULK_DE","18D_interpretation",
  "18D3_FFA_batch_sensitivity","full_results"
)

GSFILE <- file.path(
  BASE,"20_PATHWAY_ANALYSIS","20A_HALLMARK_GSEA",
  "gene_sets","20A_Hallmark_gene_sets_frozen.rds"
)

OUT <- file.path(
  BASE,"20_PATHWAY_ANALYSIS",
  "20B_HALLMARK_ARCHITECTURE",
  "20B2_FFA_batch_sensitivity"
)

TAB <- file.path(OUT,"tables")
dir.create(TAB,recursive=TRUE,showWarnings=FALSE)

stopifnot(
  dir.exists(IN18),
  file.exists(GSFILE)
)

pathways <- readRDS(GSFILE)

# ============================================================
# 1. Recover 18D3 full results
# ============================================================

cat("\nSTEP 1/6 — recovering 18D3 models\n")

files <- list.files(
  IN18,
  pattern="^18D3_(adjusted|unadjusted|historical)__.*\\.csv$",
  full.names=TRUE
)

cat("Files:",length(files),"\n")

index <- map_dfr(files,function(f){

  x <- fread(f,nrows=1)

  model <- sub(
    "^18D3_(adjusted|unadjusted|historical)__.*$",
    "\\1",
    basename(f)
  )

  tibble(
    file=f,
    model=model,
    module=x$module,
    final_L2_label=x$final_L2_label
  )
})

print(table(index$model))

stopifnot(
  sum(index$model=="adjusted")==18,
  sum(index$model=="unadjusted")==18
)

# historical can be 17 because one L2 lacked eligibility

# ============================================================
# 2. Build rankings and GSEA
# ============================================================

cat("\nSTEP 2/6 — Hallmark GSEA for sensitivity models\n")

run_one <- function(f,model,module,l2){

  x <- fread(f)

  r <- x %>%
    filter(
      !is.na(gene),
      gene!="",
      is.finite(logFC),
      is.finite(F),
      F>=0
    ) %>%
    mutate(
      rank_score=sign(logFC)*sqrt(F)
    ) %>%
    arrange(desc(abs(rank_score))) %>%
    group_by(gene) %>%
    slice_max(
      order_by=abs(rank_score),
      n=1,
      with_ties=FALSE
    ) %>%
    ungroup()

  stats <- r$rank_score
  names(stats) <- r$gene
  stats <- sort(stats,decreasing=TRUE)

  set.seed(20261001)

  fg <- fgseaMultilevel(
    pathways=pathways,
    stats=stats,
    minSize=10,
    maxSize=500,
    eps=0
  )

  as_tibble(fg) %>%
    mutate(
      model=model,
      module=module,
      final_L2_label=l2,
      direction=if_else(NES>0,"UP","DOWN"),
      .before=1
    )
}

res <- pmap(
  index,
  function(file,model,module,final_L2_label){
    cat(
      model,"|",module,"|",final_L2_label,"\n"
    )
    run_one(
      file,model,module,final_L2_label
    )
  }
)

gsea <- bind_rows(res)

# ============================================================
# 3. Flat output
# ============================================================

cat("\nSTEP 3/6 — saving results\n")

flat <- gsea %>%
  mutate(
    leadingEdge=vapply(
      leadingEdge,
      paste,
      collapse=";",
      FUN.VALUE=character(1)
    )
  )

write_csv(
  flat,
  file.path(TAB,"20B2_all_FFA_batch_sensitivity_GSEA.csv")
)

# ============================================================
# 4. Pair model NES
# ============================================================

cat("\nSTEP 4/6 — comparing NES\n")

wide <- gsea %>%
  select(
    model,module,final_L2_label,
    pathway,NES,padj
  ) %>%
  pivot_wider(
    names_from=model,
    values_from=c(NES,padj)
  )

write_csv(
  wide,
  file.path(TAB,"20B2_FFA_Hallmark_model_comparison.csv")
)

pair_stats <- function(df,a,b){

  ok <- !is.na(df[[a]]) & !is.na(df[[b]])

  tibble(
    comparison=paste(a,b,sep="_vs_"),
    n_pairs=sum(ok),
    spearman_NES=cor(
      df[[a]][ok],
      df[[b]][ok],
      method="spearman"
    ),
    pearson_NES=cor(
      df[[a]][ok],
      df[[b]][ok],
      method="pearson"
    ),
    direction_concordance=mean(
      sign(df[[a]][ok])==
        sign(df[[b]][ok])
    )
  )
}

comparison_summary <- bind_rows(
  pair_stats(
    wide,
    "NES_adjusted",
    "NES_unadjusted"
  ),
  pair_stats(
    wide,
    "NES_adjusted",
    "NES_historical"
  ),
  pair_stats(
    wide,
    "NES_unadjusted",
    "NES_historical"
  )
)

write_csv(
  comparison_summary,
  file.path(TAB,"20B2_NES_concordance.csv")
)

# ============================================================
# 5. Significant pathway retention
# ============================================================

cat("\nSTEP 5/6 — pathway retention\n")

retention <- wide %>%
  mutate(
    sig_adjusted=
      !is.na(padj_adjusted) &
      padj_adjusted<.05,

    sig_unadjusted=
      !is.na(padj_unadjusted) &
      padj_unadjusted<.05,

    sig_historical=
      !is.na(padj_historical) &
      padj_historical<.05,

    same_direction_AU=
      sign(NES_adjusted)==
      sign(NES_unadjusted),

    same_direction_AH=
      if_else(
        !is.na(NES_historical),
        sign(NES_adjusted)==
          sign(NES_historical),
        NA
      )
  )

retention_summary <- tibble(
  adjusted_hits=
    sum(retention$sig_adjusted,na.rm=TRUE),

  unadjusted_hits=
    sum(retention$sig_unadjusted,na.rm=TRUE),

  historical_hits=
    sum(retention$sig_historical,na.rm=TRUE),

  adjusted_hits_same_direction_unadjusted=
    sum(
      retention$sig_adjusted &
      retention$same_direction_AU,
      na.rm=TRUE
    ),

  adjusted_hits_same_direction_historical=
    sum(
      retention$sig_adjusted &
      retention$same_direction_AH,
      na.rm=TRUE
    ),

  adjusted_hits_FDR005_historical=
    sum(
      retention$sig_adjusted &
      retention$sig_historical &
      retention$same_direction_AH,
      na.rm=TRUE
    )
)

write_csv(
  retention,
  file.path(TAB,"20B2_pathway_retention.csv")
)

write_csv(
  retention_summary,
  file.path(TAB,"20B2_pathway_retention_summary.csv")
)

# ============================================================
# 6. Checks + audit
# ============================================================

cat("\nSTEP 6/6 — formal checks\n")

checks <- tibble(
  check=c(
    "adjusted_18",
    "unadjusted_18",
    "historical_ge16",
    "hallmark_frozen_used",
    "GSEA_generated",
    "comparison_generated",
    "no_DE_refit",
    "annotation_unchanged",
    "no_figures_generated"
  ),
  pass=c(
    sum(index$model=="adjusted")==18,
    sum(index$model=="unadjusted")==18,
    sum(index$model=="historical")>=16,
    file.exists(GSFILE),
    nrow(gsea)>0,
    nrow(comparison_summary)==3,
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
  "PIPELINE_V2 — 20B2 FFA Hallmark batch sensitivity",
  "",
  "Purpose:",
  "- Evaluate stability of FFA-vs-HC Hallmark enrichment across batch specifications.",
  "",
  "Models:",
  "- adjusted: all eligible FFA/HC, ~ batch + disease",
  "- unadjusted: same patients, ~ disease",
  "- historical: historical-only patients, ~ disease",
  "",
  "Input:",
  "- Frozen 18D3 edgeR sensitivity results",
  "- Frozen 20A Hallmark gene sets",
  "",
  "Ranking:",
  "- sign(logFC) * sqrt(edgeR QL F)",
  "",
  "GSEA:",
  "- fgseaMultilevel",
  "- minSize=10",
  "- maxSize=500",
  "- eps=0",
  "",
  "This is a sensitivity analysis.",
  "Primary FFA-vs-HC inference remains the batch-adjusted 18C/20A analysis.",
  "Historical-only is not an independent validation cohort.",
  "No DE was refit.",
  "No annotation modified.",
  "No figures generated."
)

writeLines(
  audit,
  file.path(OUT,"20B2_audit.txt")
)

cat("\n===== NES CONCORDANCE =====\n")
print(
  as.data.frame(comparison_summary),
  row.names=FALSE
)

cat("\n===== PATHWAY RETENTION =====\n")
print(
  as.data.frame(retention_summary),
  row.names=FALSE
)

cat("\n============================================================\n")
cat("20B2 COMPLETE\n")
cat("============================================================\n")
cat("FFA Hallmark batch sensitivity\n")
cat("Primary 20A unchanged\n")
cat("PASS\n")
