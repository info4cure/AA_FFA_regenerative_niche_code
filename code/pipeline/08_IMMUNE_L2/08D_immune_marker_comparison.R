suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(readr)
})

cat("\n============================================================\n")
cat("PIPELINE_V2 — 08D IMMUNE MARKER COMPARISON\n")
cat("============================================================\n\n")

IN <- "PIPELINE_V2/objects/annotated/immune_v2_08B_reclustered.rds"

OUT <- "PIPELINE_V2/results/annotation/08_IMMUNE_L2/08D_markers"
dir.create(OUT,recursive=TRUE,showWarnings=FALSE)

x <- readRDS(IN)
DefaultAssay(x) <- "RNA"

stopifnot(
  ncol(x)==4762,
  nrow(x)==23412
)

cat("Cells:",ncol(x),"\n")
cat("Genes:",nrow(x),"\n")

# ============================================================
# helper
# ============================================================

run_markers <- function(obj,col,res){

  cat("\nMarkers resolution",res,"\n")

  Idents(obj) <- obj@meta.data[[col]]

  mk <- FindAllMarkers(
    obj,
    assay="RNA",
    only.pos=TRUE,
    min.pct=0.10,
    logfc.threshold=0.25,
    test.use="wilcox",
    verbose=FALSE
  )

  mk$cluster <- as.character(mk$cluster)

  mk <- mk %>%
    arrange(
      as.numeric(cluster),
      desc(avg_log2FC)
    )

  top20 <- mk %>%
    group_by(cluster) %>%
    slice_max(
      avg_log2FC,
      n=20,
      with_ties=FALSE
    ) %>%
    ungroup()

  write_csv(
    mk,
    file.path(
      OUT,
      paste0(
        "08D_markers_res",
        res,
        "_all.csv"
      )
    )
  )

  write_csv(
    top20,
    file.path(
      OUT,
      paste0(
        "08D_markers_res",
        res,
        "_top20.csv"
      )
    )
  )

  list(
    all=mk,
    top20=top20
  )
}

# ============================================================
# 1. RESOLUTION 0.4
# ============================================================

cat("\nSTEP 1/6 — markers resolution 0.4\n")

m04 <- run_markers(
  x,
  "immune_snn_res.0.4",
  "0.4"
)

# ============================================================
# 2. RESOLUTION 0.5
# ============================================================

cat("\nSTEP 2/6 — markers resolution 0.5\n")

m05 <- run_markers(
  x,
  "immune_snn_res.0.5",
  "0.5"
)

# ============================================================
# 3. RESOLUTION 0.6
# ============================================================

cat("\nSTEP 3/6 — markers resolution 0.6\n")

m06 <- run_markers(
  x,
  "immune_snn_res.0.6",
  "0.6"
)

# ============================================================
# 4. 0.5 -> 0.6 TRANSITIONS
# ============================================================

cat("\nSTEP 4/6 — 0.5 -> 0.6 transitions\n")

md <- x@meta.data %>%
  mutate(
    parent=as.character(
      .data[["immune_snn_res.0.5"]]
    ),
    child=as.character(
      .data[["immune_snn_res.0.6"]]
    )
  )

tr <- md %>%
  count(parent,child,name="n_cells") %>%
  group_by(parent) %>%
  mutate(
    parent_total=sum(n_cells),
    fraction_parent=n_cells/parent_total
  ) %>%
  ungroup()

write_csv(
  tr,
  file.path(
    OUT,
    "08D_transition_res0.5_to_0.6.csv"
  )
)

cat("\n=== 0.5 -> 0.6 TRANSITIONS >=5% ===\n")

print(
  as.data.frame(
    tr %>%
      filter(fraction_parent>=0.05) %>%
      arrange(
        as.numeric(parent),
        desc(fraction_parent)
      )
  ),
  row.names=FALSE
)

# ============================================================
# 5. DIRECT MARKERS FOR TRUE SPLITS
# ============================================================

cat("\nSTEP 5/6 — direct markers for true splits\n")

splits <- tr %>%
  filter(fraction_parent>=0.05) %>%
  count(parent,name="n_children") %>%
  filter(n_children>=2)

cat("\nParents with >=2 children contributing >=5%:\n")
print(as.data.frame(splits),row.names=FALSE)

split_summary <- list()

for(parent_id in splits$parent){

  children <- tr %>%
    filter(
      parent==parent_id,
      fraction_parent>=0.05
    ) %>%
    arrange(desc(fraction_parent)) %>%
    pull(child)

  cells_parent <- rownames(md)[
    md$parent==parent_id
  ]

  sub <- subset(
    x,
    cells=cells_parent
  )

  Idents(sub) <-
    as.character(
      sub@meta.data[["immune_snn_res.0.6"]]
    )

  # pairwise comparisons between all major children
  if(length(children)>=2){

    cmb <- combn(
      children,
      2,
      simplify=FALSE
    )

    for(pair in cmb){

      a <- pair[1]
      b <- pair[2]

      cat(
        "Parent",parent_id,
        ": child",a,
        "vs",b,"\n"
      )

      de <- FindMarkers(
        sub,
        ident.1=a,
        ident.2=b,
        assay="RNA",
        min.pct=0.05,
        logfc.threshold=0,
        test.use="wilcox",
        verbose=FALSE
      )

      de$gene <- rownames(de)

      de <- de %>%
        arrange(desc(abs(avg_log2FC)))

      fn <- paste0(
        "08D_split_parent",
        parent_id,
        "_child",
        a,
        "_vs_",
        b,
        ".csv"
      )

      write_csv(
        de,
        file.path(OUT,fn)
      )

      top <- de %>%
        slice_head(n=30) %>%
        mutate(
          parent=parent_id,
          child_A=a,
          child_B=b
        )

      split_summary[[length(split_summary)+1]] <-
        top
    }
  }
}

if(length(split_summary)>0){

  split_top <- bind_rows(split_summary)

  write_csv(
    split_top,
    file.path(
      OUT,
      "08D_split_comparison_top30.csv"
    )
  )
}

# ============================================================
# 6. COMPACT MARKER REPORT
# ============================================================

cat("\nSTEP 6/6 — compact marker report\n")

print_compact <- function(top,res){

  cat(
    "\n================ RESOLUTION ",
    res,
    " ================\n",
    sep=""
  )

  clusters <- sort(
    unique(as.character(top$cluster)),
    na.last=TRUE
  )

  # numeric ordering
  clusters <- clusters[
    order(as.numeric(clusters))
  ]

  for(cl in clusters){

    z <- top %>%
      filter(
        as.character(cluster)==cl
      ) %>%
      arrange(desc(avg_log2FC))

    cat(
      "C",cl," | ",
      paste(z$gene,collapse=", "),
      "\n",
      sep=""
    )
  }
}

print_compact(m04$top20,"0.4")
print_compact(m05$top20,"0.5")
print_compact(m06$top20,"0.6")

# ============================================================
# AUDIT
# ============================================================

audit <- c(
  "PIPELINE_V2 — 08D immune marker comparison",
  paste0("Generated: ",Sys.time()),
  "Input: immune_v2_08B_reclustered.rds",
  "Cells: 4762",
  "Markers calculated independently at resolutions 0.4, 0.5 and 0.6",
  "FindAllMarkers: Wilcoxon; only.pos=TRUE; min.pct=0.10; logfc.threshold=0.25",
  "0.5->0.6 cluster transitions explicitly audited",
  "Major splits defined descriptively as children >=5% of parent",
  "Direct pairwise marker tests generated for major split children",
  "No resolution selected",
  "No immune identities assigned",
  "No object modified"
)

writeLines(
  audit,
  file.path(
    OUT,
    "08D_audit.txt"
  )
)

cat("\n============================================================\n")
cat("08D COMPLETE\n")
cat("============================================================\n")
cat("No resolution selected.\n")
cat("No immune identities assigned.\n")
cat("PASS — de novo immune marker evidence generated.\n")
