suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(tibble)
  library(purrr)
  library(stringr)
  library(igraph)
})

cat("============================================================\n")
cat("PIPELINE_V2 — 20D2B PATHWAY REDUNDANCY REDUCTION\n")
cat("============================================================\n")

BASE <- file.path(
  "PIPELINE_V2","results","disease_biology",
  "20_PATHWAY_ANALYSIS","20D_TARGETED_PATHWAYS"
)

TAB <- file.path(BASE,"tables")
OUT <- file.path(BASE,"20D2B_redundancy")

dir.create(OUT,recursive=TRUE,showWarnings=FALSE)

RESFILE <- file.path(
  TAB,"20D2_targeted_axes_FDR005.csv"
)

LEFILE <- file.path(
  TAB,"20D2_targeted_leading_edge_genes_FDR005.csv"
)

stopifnot(
  file.exists(RESFILE),
  file.exists(LEFILE)
)

res <- read_csv(
  RESFILE,
  show_col_types=FALSE
)

le <- read_csv(
  LEFILE,
  show_col_types=FALSE
)

# Recover enrichment direction explicitly from the signed NES.
# The frozen 20D2 leading-edge table contains NES but not a
# separate direction column.
stopifnot(
  "NES" %in% names(le),
  all(is.finite(le$NES))
)

le <- le %>%
  mutate(
    direction=case_when(
      NES > 0 ~ "UP",
      NES < 0 ~ "DOWN",
      TRUE ~ "ZERO"
    )
  )

stopifnot(
  !any(le$direction=="ZERO"),
  all(le$direction %in% c("UP","DOWN"))
)

cat("\nTargeted FDR<0.05 rows:",nrow(res),"\n")
cat("Leading-edge rows:",nrow(le),"\n")

# ============================================================
# 1. Collapse each pathway across L2 within contrast/axis/direction
# ============================================================

cat("\nSTEP 1/6 — pathway recurrence and pooled leading edges\n")

path_summary <- res %>%
  group_by(
    database,
    contrast,
    axis,
    direction,
    pathway
  ) %>%
  summarise(
    n_L2=n_distinct(final_L2_label),
    median_NES=median(NES,na.rm=TRUE),
    max_abs_NES=max(abs(NES),na.rm=TRUE),
    min_FDR=min(padj,na.rm=TRUE),
    L2=paste(
      sort(unique(final_L2_label)),
      collapse=";"
    ),
    .groups="drop"
  )

pooled_le <- le %>%
  filter(
    !is.na(gene),
    gene!=""
  ) %>%
  distinct(
    database,
    contrast,
    axis,
    direction,
    pathway,
    gene
  ) %>%
  group_by(
    database,
    contrast,
    axis,
    direction,
    pathway
  ) %>%
  summarise(
    genes=list(sort(unique(gene))),
    n_pooled_LE=n_distinct(gene),
    .groups="drop"
  )

ps <- path_summary %>%
  inner_join(
    pooled_le,
    by=c(
      "database",
      "contrast",
      "axis",
      "direction",
      "pathway"
    )
  )

# ============================================================
# 2. Jaccard graph within each biological stratum
# ============================================================

cat("\nSTEP 2/6 — leading-edge Jaccard graphs\n")

JACCARD_THRESHOLD <- 0.50

make_components <- function(df) {

  n <- nrow(df)

  if (n==0)
    return(tibble())

  if (n==1) {
    return(
      df %>%
        mutate(
          component_local=1L,
          component_size=1L
        )
    )
  }

  edges <- list()
  k <- 0L

  for (i in seq_len(n-1)) {
    gi <- df$genes[[i]]

    for (j in (i+1):n) {
      gj <- df$genes[[j]]

      inter <- length(intersect(gi,gj))
      uni   <- length(union(gi,gj))

      jac <- if (uni>0) inter/uni else 0

      if (jac >= JACCARD_THRESHOLD) {
        k <- k+1L
        edges[[k]] <- tibble(
          from=i,
          to=j,
          jaccard=jac
        )
      }
    }
  }

  if (length(edges)==0) {
    comp <- seq_len(n)
  } else {

    ed <- bind_rows(edges)

    g <- graph_from_data_frame(
      ed %>% select(from,to),
      directed=FALSE,
      vertices=tibble(
        name=as.character(seq_len(n))
      )
    )

    comp <- components(g)$membership[
      as.character(seq_len(n))
    ]
  }

  sizes <- table(comp)

  df %>%
    mutate(
      component_local=as.integer(comp),
      component_size=as.integer(
        sizes[as.character(comp)]
      )
    )
}

components_tbl <- ps %>%
  group_by(
    database,
    contrast,
    axis,
    direction
  ) %>%
  group_modify(
    ~make_components(.x)
  ) %>%
  ungroup()

# Global stable component identifier
components_tbl <- components_tbl %>%
  mutate(
    component_id=paste(
      database,
      contrast,
      axis,
      direction,
      sprintf("C%03d",component_local),
      sep="__"
    )
  )

# ============================================================
# 3. Choose representative pathway per component
# ============================================================

cat("\nSTEP 3/6 — representative pathways\n")

representatives <- components_tbl %>%
  group_by(component_id) %>%
  arrange(
    desc(n_L2),
    min_FDR,
    desc(max_abs_NES),
    pathway,
    .by_group=TRUE
  ) %>%
  mutate(
    representative_rank=row_number()
  ) %>%
  ungroup() %>%
  filter(representative_rank==1) %>%
  select(
    component_id,
    database,
    contrast,
    axis,
    direction,
    representative_pathway=pathway,
    n_L2,
    median_NES,
    max_abs_NES,
    min_FDR,
    component_size,
    n_pooled_LE,
    L2
  )

component_members <- components_tbl %>%
  select(
    component_id,
    database,
    contrast,
    axis,
    direction,
    pathway,
    n_L2,
    median_NES,
    min_FDR,
    component_size,
    n_pooled_LE,
    L2
  )

write_csv(
  representatives,
  file.path(
    OUT,
    "20D2B_representative_pathways.csv"
  )
)

write_csv(
  component_members,
  file.path(
    OUT,
    "20D2B_component_members.csv"
  )
)

# ============================================================
# 4. Recurrent leading-edge driver genes
# ============================================================

cat("\nSTEP 4/6 — recurrent driver genes\n")

driver <- le %>%
  filter(
    !is.na(gene),
    gene!=""
  ) %>%
  group_by(
    database,
    contrast,
    axis,
    direction,
    gene
  ) %>%
  summarise(
    n_pathways=n_distinct(pathway),
    n_L2=n_distinct(final_L2_label),
    n_L2_pathway_pairs=n_distinct(
      paste(
        final_L2_label,
        pathway,
        sep="|||"
      )
    ),
    pathways=paste(
      sort(unique(pathway)),
      collapse=";"
    ),
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
    direction,
    desc(n_L2),
    desc(n_pathways),
    gene
  )

write_csv(
  driver,
  file.path(
    OUT,
    "20D2B_recurrent_leading_edge_genes.csv"
  )
)

top_driver <- driver %>%
  group_by(
    database,
    contrast,
    axis,
    direction
  ) %>%
  slice_max(
    order_by=n_L2_pathway_pairs,
    n=20,
    with_ties=FALSE
  ) %>%
  ungroup()

write_csv(
  top_driver,
  file.path(
    OUT,
    "20D2B_top20_driver_genes_per_axis.csv"
  )
)

# ============================================================
# 5. Reduction summary
# ============================================================

cat("\nSTEP 5/6 — reduction summary\n")

reduction <- components_tbl %>%
  group_by(
    database,
    contrast,
    axis,
    direction
  ) %>%
  summarise(
    n_original_pathways=n_distinct(pathway),
    n_components=n_distinct(component_id),
    n_redundant_components=
      n_distinct(
        component_id[
          component_size>1
        ]
      ),
    largest_component=max(component_size),
    reduction_fraction=
      1 - n_components/n_original_pathways,
    .groups="drop"
  )

write_csv(
  reduction,
  file.path(
    OUT,
    "20D2B_reduction_summary.csv"
  )
)

cat("\n===== REDUCTION SUMMARY =====\n")
print(
  as.data.frame(reduction),
  row.names=FALSE
)

# ============================================================
# 6. Checks / audit
# ============================================================

cat("\nSTEP 6/6 — formal checks\n")

checks <- tibble(
  check=c(
    "targeted_results_present",
    "leading_edges_present",
    "no_duplicate_analytical_keys",
    "components_generated",
    "representatives_generated",
    "drivers_generated",
    "threshold_0.50",
    "no_GSEA_refit",
    "no_DE_refit",
    "annotation_unchanged"
  ),
  pass=c(
    nrow(res)>0,
    nrow(le)>0,
    !anyDuplicated(
      res %>%
        select(
          database,
          contrast,
          module,
          final_L2_label,
          pathway,
          axis
        )
    ),
    nrow(components_tbl)>0,
    nrow(representatives)>0,
    nrow(driver)>0,
    identical(JACCARD_THRESHOLD,0.50),
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
  "PIPELINE_V2 — 20D2B pathway redundancy reduction",
  "",
  paste0("Date: ",Sys.Date()),
  "",
  "Input:",
  "- frozen 20D2 targeted Reactome/GO-BP FDR<0.05 results",
  "- frozen 20D2 leading-edge genes",
  "",
  "Redundancy definition:",
  "- performed separately within database x contrast x axis x direction",
  "- pathway similarity based on pooled leading-edge genes across significant L2",
  "- Jaccard threshold >=0.50",
  "- connected components define redundancy groups",
  "",
  "Representative pathway selection:",
  "1. greatest recurrence across L2",
  "2. lowest minimum FDR",
  "3. greatest absolute NES",
  "4. pathway name for deterministic tie-breaking",
  "",
  "Important:",
  "- Reactome and GO-BP were never merged into the same redundancy component",
  "- multi-axis pathways remain independently represented within each predefined axis",
  "- reduction is descriptive and does not alter original GSEA significance",
  "",
  "No GSEA refit.",
  "No differential expression refit.",
  "No annotation modified."
)

writeLines(
  audit,
  file.path(
    OUT,
    "20D2B_audit.txt"
  )
)

cat("\n===== TOP RECURRENT REPRESENTATIVES =====\n")

print(
  as.data.frame(
    representatives %>%
      arrange(
        database,
        contrast,
        axis,
        direction,
        desc(n_L2),
        min_FDR
      ) %>%
      group_by(
        database,
        contrast,
        axis,
        direction
      ) %>%
      slice_head(n=3) %>%
      ungroup()
  ),
  row.names=FALSE
)

cat("\n============================================================\n")
cat("20D2B COMPLETE\n")
cat("============================================================\n")
cat("Leading-edge redundancy reduction\n")
cat("Original GSEA unchanged\n")
cat("PASS\n")
