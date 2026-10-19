suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(readr)
  library(tibble)
  library(purrr)
  library(stringr)
  library(fgsea)
})

cat("============================================================\n")
cat("PIPELINE_V2 — 20D2 FULL REACTOME + GO-BP GSEA\n")
cat("============================================================\n")

BASE <- file.path(
  "PIPELINE_V2","results","disease_biology",
  "20_PATHWAY_ANALYSIS"
)

RANKDIR <- file.path(
  BASE,"20A_HALLMARK_GSEA","rankings"
)

D20 <- file.path(
  BASE,"20D_TARGETED_PATHWAYS"
)

GS <- file.path(D20,"gene_sets")
TAB <- file.path(D20,"tables")
FULL <- file.path(D20,"full_results")

dir.create(FULL,recursive=TRUE,showWarnings=FALSE)

RFILE <- file.path(
  GS,"20D1_Reactome_2026.1.Hs.rds"
)

GFILE <- file.path(
  GS,"20D1_GO_BP_2026.1.Hs.rds"
)

AXISFILE <- file.path(
  TAB,"20D1_targeted_pathway_inventory.csv"
)

stopifnot(
  file.exists(RFILE),
  file.exists(GFILE),
  file.exists(AXISFILE)
)

collections <- list(
  REACTOME=readRDS(RFILE),
  GO_BP=readRDS(GFILE)
)

axis_map <- read_csv(
  AXISFILE,
  show_col_types=FALSE
)

# ============================================================
# 1. Recover the 60 frozen rankings
# ============================================================

cat("\nSTEP 1/7 — recovering frozen 20A rankings\n")

rank_files <- list.files(
  RANKDIR,
  pattern="^20A_rank__.*\\.csv$",
  full.names=TRUE
)

stopifnot(length(rank_files)==60)

# Metadata are recovered from the filename and ranking source.
# Use frozen 20A ranking summary/result table to avoid ambiguity.

GSEA20A <- file.path(
  BASE,"20A_HALLMARK_GSEA","tables",
  "20A_all_Hallmark_GSEA_results.csv"
)

meta20 <- read_csv(
  GSEA20A,
  show_col_types=FALSE
) %>%
  distinct(
    contrast,
    module,
    final_L2_label
  )

stopifnot(nrow(meta20)==60)

safe_name <- function(x) {
  gsub("[^A-Za-z0-9._-]+","_",x)
}

expected_file <- function(ct,m,l2) {
  file.path(
    RANKDIR,
    paste0(
      "20A_rank__",
      safe_name(
        paste(ct,m,safe_name(l2),sep="__")
      ),
      ".csv"
    )
  )
}

index <- meta20 %>%
  mutate(
    rank_file=pmap_chr(
      list(
        contrast,
        module,
        final_L2_label
      ),
      expected_file
    )
  )

missing <- index$rank_file[
  !file.exists(index$rank_file)
]

if (length(missing)) {
  stop(
    paste(
      "Missing frozen ranking files:",
      paste(missing,collapse="\n")
    )
  )
}

cat("Rankings:",nrow(index),"\n")
print(table(index$contrast))

# ============================================================
# 2. GSEA runner
# ============================================================

cat("\nSTEP 2/7 — defining full GSEA runner\n")

run_one <- function(rank_file,pathways,database,ct,m,l2) {

  r <- fread(rank_file)

  stopifnot(
    all(c("gene","rank_score") %in% names(r))
  )

  r <- r[
    !is.na(gene) &
    gene!="" &
    is.finite(rank_score)
  ]

  if (anyDuplicated(r$gene)) {
    stop(
      paste("Duplicated ranking genes:",rank_file)
    )
  }

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
      database=database,
      contrast=ct,
      module=m,
      final_L2_label=l2,
      direction=if_else(NES>0,"UP","DOWN"),
      n_leading_edge=lengths(leadingEdge),
      .before=1
    )
}

# ============================================================
# 3. Run complete collections
# ============================================================

cat("\nSTEP 3/7 — running complete Reactome and GO-BP collections\n")

all_results <- list()
k <- 0L

for (db in names(collections)) {

  pathways <- collections[[db]]

  cat(
    "\nDATABASE:",db,
    "| pathways:",length(pathways),"\n"
  )

  for (i in seq_len(nrow(index))) {

    ct <- index$contrast[i]
    m  <- index$module[i]
    l2 <- index$final_L2_label[i]

    cat(
      sprintf(
        "[%s %02d/60] %-10s | %-18s | %s\n",
        db,i,ct,m,l2
      )
    )

    fg <- run_one(
      index$rank_file[i],
      pathways,
      db,
      ct,m,l2
    )

    k <- k+1L
    all_results[[k]] <- fg

    flat <- fg %>%
      mutate(
        leadingEdge=vapply(
          leadingEdge,
          paste,
          collapse=";",
          FUN.VALUE=character(1)
        )
      )

    outfile <- file.path(
      FULL,
      paste0(
        "20D2_",db,"__",
        safe_name(ct),"__",
        safe_name(m),"__",
        safe_name(l2),
        ".csv"
      )
    )

    write_csv(flat,outfile)
  }
}

res <- bind_rows(all_results)

# ============================================================
# 4. Consolidate and targeted-axis annotation
# ============================================================

cat("\nSTEP 4/7 — consolidating and annotating targeted axes\n")

flat_all <- res %>%
  mutate(
    leadingEdge=vapply(
      leadingEdge,
      paste,
      collapse=";",
      FUN.VALUE=character(1)
    )
  )

write_csv(
  flat_all,
  file.path(
    TAB,
    "20D2_all_Reactome_GOBP_GSEA.csv"
  )
)

sig005 <- flat_all %>%
  filter(padj<.05)

write_csv(
  sig005,
  file.path(
    TAB,
    "20D2_FDR005_all.csv"
  )
)

sig010 <- flat_all %>%
  filter(padj<.10)

write_csv(
  sig010,
  file.path(
    TAB,
    "20D2_FDR010_all.csv"
  )
)

targeted <- flat_all %>%
  inner_join(
    axis_map,
    by=c(
      "database",
      "pathway"
    )
  )

write_csv(
  targeted,
  file.path(
    TAB,
    "20D2_targeted_axes_all_results.csv"
  )
)

write_csv(
  targeted %>% filter(padj<.05),
  file.path(
    TAB,
    "20D2_targeted_axes_FDR005.csv"
  )
)

# ============================================================
# 5. Summaries
# ============================================================

cat("\nSTEP 5/7 — summaries\n")

global_summary <- flat_all %>%
  group_by(
    database,
    contrast
  ) %>%
  summarise(
    n_L2=n_distinct(
      paste(module,final_L2_label)
    ),
    L2_with_FDR005=n_distinct(
      paste(
        module[padj<.05],
        final_L2_label[padj<.05]
      )
    ),
    total_FDR005=sum(
      padj<.05,
      na.rm=TRUE
    ),
    total_UP_FDR005=sum(
      padj<.05 & NES>0,
      na.rm=TRUE
    ),
    total_DOWN_FDR005=sum(
      padj<.05 & NES<0,
      na.rm=TRUE
    ),
    .groups="drop"
  )

write_csv(
  global_summary,
  file.path(
    TAB,
    "20D2_global_summary.csv"
  )
)

axis_summary <- targeted %>%
  filter(padj<.05) %>%
  group_by(
    database,
    contrast,
    axis,
    direction
  ) %>%
  summarise(
    n_hits=n(),
    n_pathways=n_distinct(pathway),
    n_L2=n_distinct(final_L2_label),
    median_NES=median(NES,na.rm=TRUE),
    min_FDR=min(padj,na.rm=TRUE),
    .groups="drop"
  ) %>%
  arrange(
    database,
    contrast,
    axis,
    direction
  )

write_csv(
  axis_summary,
  file.path(
    TAB,
    "20D2_targeted_axis_summary_FDR005.csv"
  )
)

recurrence <- targeted %>%
  filter(padj<.05) %>%
  group_by(
    database,
    contrast,
    axis,
    pathway,
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
    database,
    contrast,
    axis,
    desc(n_L2),
    min_FDR
  )

write_csv(
  recurrence,
  file.path(
    TAB,
    "20D2_targeted_pathway_recurrence_FDR005.csv"
  )
)

# ============================================================
# 6. Leading-edge table for significant targeted pathways
# ============================================================

cat("\nSTEP 6/7 — targeted leading edges\n")

leading <- res %>%
  filter(padj<.05) %>%
  inner_join(
    axis_map,
    by=c(
      "database",
      "pathway"
    )
  ) %>%
  select(
    database,
    contrast,
    module,
    final_L2_label,
    axis,
    pathway,
    NES,
    padj,
    leadingEdge
  ) %>%
  tidyr::unnest_longer(
    leadingEdge,
    values_to="gene"
  )

write_csv(
  leading,
  file.path(
    TAB,
    "20D2_targeted_leading_edge_genes_FDR005.csv"
  )
)

# ============================================================
# 7. Checks + audit
# ============================================================

cat("\nSTEP 7/7 — formal checks\n")

checks <- tibble(
  check=c(
    "rankings_60",
    "Reactome_1839_frozen",
    "GO_BP_7538_frozen",
    "two_databases_run",
    "all_results_generated",
    "targeted_axes_annotated",
    "leading_edges_generated",
    "no_DE_refit",
    "no_annotation_change",
    "no_figures_generated"
  ),
  pass=c(
    nrow(index)==60,
    length(collections$REACTOME)==1839,
    length(collections$GO_BP)==7538,
    all(
      c("REACTOME","GO_BP")
      %in% unique(flat_all$database)
    ),
    nrow(flat_all)>0,
    nrow(targeted)>0,
    nrow(leading)>0,
    TRUE,
    TRUE,
    TRUE
  )
)

print(as.data.frame(checks),row.names=FALSE)
stopifnot(all(checks$pass))

audit <- c(
  "PIPELINE_V2 — 20D2 full Reactome + GO-BP preranked GSEA",
  "",
  paste0("Date: ",Sys.Date()),
  paste0(
    "fgsea version: ",
    as.character(packageVersion("fgsea"))
  ),
  "MSigDB version: 2026.1.Hs",
  "",
  "Inputs:",
  "- 60 frozen rankings generated in 20A",
  "- Reactome 2026.1.Hs: 1,839 pathways",
  "- GO Biological Process 2026.1.Hs: 7,538 pathways",
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
  "Multiple testing:",
  "- BH FDR returned by fgsea within each database x L2 x contrast analysis",
  "- primary pathway threshold FDR < 0.05",
  "",
  "Interpretation:",
  "- full collections tested without targeted prefiltering",
  "- eight axes frozen in 20D1 used only for targeted interpretation",
  "- leading-edge genes retained for significant targeted pathways",
  "",
  "No pseudobulk DE refit.",
  "No cell annotation modified.",
  "No manuscript figures generated."
)

writeLines(
  audit,
  file.path(
    D20,
    "20D2_audit.txt"
  )
)

cat("\n===== GLOBAL SUMMARY =====\n")
print(
  as.data.frame(global_summary),
  row.names=FALSE
)

cat("\n===== TARGETED AXIS SUMMARY =====\n")
print(
  as.data.frame(axis_summary),
  row.names=FALSE
)

cat("\n============================================================\n")
cat("20D2 COMPLETE\n")
cat("============================================================\n")
cat("Full Reactome + GO-BP GSEA\n")
cat("60 frozen rankings\n")
cat("Targeted axes annotated post hoc from frozen 20D1 inventory\n")
cat("PASS\n")
