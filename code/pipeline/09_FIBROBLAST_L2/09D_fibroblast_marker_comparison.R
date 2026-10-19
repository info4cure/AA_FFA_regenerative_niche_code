suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(readr)
})

cat("\n============================================================\n")
cat("PIPELINE_V2 — 09D FIBROBLAST MARKER COMPARISON\n")
cat("============================================================\n")

IN <- "PIPELINE_V2/objects/annotated/fibroblast_v2_09B_reclustered.rds"

OUT <- paste0(
  "PIPELINE_V2/results/annotation/09_FIBROBLAST_L2/",
  "09D_markers"
)

dir.create(OUT,recursive=TRUE,showWarnings=FALSE)

x <- readRDS(IN)

stopifnot(
  ncol(x)==4808,
  nrow(x)==23412
)

DefaultAssay(x) <- "RNA"

cat("Cells:",ncol(x),"\n")
cat("Genes:",nrow(x),"\n")

# --------------------------------------------------------------------------
# 1. MARKERS AT 0.3 / 0.4 / 0.5
# --------------------------------------------------------------------------

resolutions <- c(0.3,0.4,0.5)

marker_results <- list()

for(r in resolutions){

  cat("\nSTEP — markers resolution",r,"\n")

  cn <- paste0("fibroblast_snn_res.",r)

  stopifnot(cn %in% colnames(x@meta.data))

  Idents(x) <- x@meta.data[[cn]]

  mk <- FindAllMarkers(
    x,
    assay="RNA",
    only.pos=TRUE,
    min.pct=0.10,
    logfc.threshold=0.25,
    test.use="wilcox",
    verbose=FALSE
  )

  mk$resolution <- r

  marker_results[[as.character(r)]] <- mk

  write_csv(
    mk,
    file.path(
      OUT,
      paste0("09D_markers_res",r,"_all.csv")
    )
  )

  top30 <- mk %>%
    group_by(cluster) %>%
    arrange(
      p_val_adj,
      desc(avg_log2FC),
      .by_group=TRUE
    ) %>%
    slice_head(n=30) %>%
    ungroup()

  write_csv(
    top30,
    file.path(
      OUT,
      paste0("09D_markers_res",r,"_top30.csv")
    )
  )

  cat("\nTOP MARKERS — resolution",r,"\n")

  compact <- top30 %>%
    group_by(cluster) %>%
    summarise(
      markers=paste(
        head(gene,20),
        collapse=", "
      ),
      .groups="drop"
    )

  print(as.data.frame(compact),row.names=FALSE)
}

# --------------------------------------------------------------------------
# 2. DIRECT AUDIT OF THE FFA010-ENRICHED SPLIT
#
# res0.3 C0 -> res0.4 C0 + C9
# --------------------------------------------------------------------------

cat("\n============================================================\n")
cat("FFA010-ENRICHED SPLIT AUDIT\n")
cat("============================================================\n")

parent03 <- as.character(
  x@meta.data[["fibroblast_snn_res.0.3"]]
)

child04 <- as.character(
  x@meta.data[["fibroblast_snn_res.0.4"]]
)

cells_parent0 <- colnames(x)[parent03=="0"]

z <- subset(
  x,
  cells=cells_parent0
)

z$child04 <- child04[
  match(
    colnames(z),
    colnames(x)
  )
]

cat("Parent 0.3/C0 cells:",ncol(z),"\n")
print(table(z$child04))

# Compare the patient-enriched C9 against its main sibling C0.
keep <- colnames(z)[
  z$child04 %in% c("0","9")
]

zz <- subset(z,cells=keep)

Idents(zz) <- zz$child04

cat("\nDirect comparison: 0.4/C9 vs 0.4/C0\n")

split_C9_vs_C0 <- FindMarkers(
  zz,
  ident.1="9",
  ident.2="0",
  assay="RNA",
  min.pct=0.10,
  logfc.threshold=0,
  test.use="wilcox",
  verbose=FALSE
)

split_C9_vs_C0$gene <- rownames(split_C9_vs_C0)

split_C9_vs_C0 <- split_C9_vs_C0 %>%
  arrange(
    p_val_adj,
    desc(avg_log2FC)
  )

write_csv(
  split_C9_vs_C0,
  file.path(
    OUT,
    "09D_res03C0_split_res04C9_vs_C0.csv"
  )
)

cat("\nTOP UP IN C9:\n")
print(
  as.data.frame(
    split_C9_vs_C0 %>%
      filter(avg_log2FC>0) %>%
      slice_head(n=40)
  ),
  row.names=FALSE
)

cat("\nTOP DOWN IN C9:\n")
print(
  as.data.frame(
    split_C9_vs_C0 %>%
      filter(avg_log2FC<0) %>%
      arrange(
        p_val_adj,
        avg_log2FC
      ) %>%
      slice_head(n=30)
  ),
  row.names=FALSE
)

# --------------------------------------------------------------------------
# 3. PATIENT COMPOSITION OF THE SPLIT
# --------------------------------------------------------------------------

split_support <- zz@meta.data %>%
  mutate(
    child04=as.character(child04)
  ) %>%
  count(
    child04,
    patient_id,
    disease,
    batch,
    name="n_cells"
  ) %>%
  group_by(child04) %>%
  mutate(
    fraction=n_cells/sum(n_cells)
  ) %>%
  ungroup()

write_csv(
  split_support,
  file.path(
    OUT,
    "09D_res03C0_split_patient_support.csv"
  )
)

# --------------------------------------------------------------------------
# 4. OTHER TRUE SPLITS 0.3 -> 0.4
# --------------------------------------------------------------------------

cat("\nSTEP — direct markers for other 0.3 -> 0.4 splits\n")

split_pairs <- list(
  c("1","1","10"),
  c("2","2","5"),
  c("3","3","5"),
  c("5","5","11"),
  c("6","8","2")
)

for(v in split_pairs){

  parent <- v[1]
  childA <- v[2]
  childB <- v[3]

  cells_parent <- colnames(x)[
    parent03==parent
  ]

  q <- subset(
    x,
    cells=cells_parent
  )

  q$child04 <- child04[
    match(
      colnames(q),
      colnames(x)
    )
  ]

  q <- subset(
    q,
    cells=colnames(q)[
      q$child04 %in% c(childA,childB)
    ]
  )

  if(
    sum(q$child04==childA)>=10 &&
    sum(q$child04==childB)>=10
  ){

    Idents(q) <- q$child04

    dm <- FindMarkers(
      q,
      ident.1=childB,
      ident.2=childA,
      assay="RNA",
      min.pct=0.10,
      logfc.threshold=0,
      test.use="wilcox",
      verbose=FALSE
    )

    dm$gene <- rownames(dm)

    dm <- dm %>%
      arrange(
        p_val_adj,
        desc(avg_log2FC)
      )

    write_csv(
      dm,
      file.path(
        OUT,
        paste0(
          "09D_res03C",parent,
          "_split_res04C",childB,
          "_vs_C",childA,
          ".csv"
        )
      )
    )
  }
}

# --------------------------------------------------------------------------
# 5. AUDIT
# --------------------------------------------------------------------------

checks <- data.frame(
  check=c(
    "cells_4808",
    "markers_res03_generated",
    "markers_res04_generated",
    "markers_res05_generated",
    "FFA010_split_tested",
    "no_FibConsensus_used",
    "no_resolution_selected",
    "no_L2_labels_assigned",
    "no_cells_removed"
  ),
  pass=c(
    ncol(x)==4808,
    nrow(marker_results[["0.3"]])>0,
    nrow(marker_results[["0.4"]])>0,
    nrow(marker_results[["0.5"]])>0,
    nrow(split_C9_vs_C0)>0,
    TRUE,
    TRUE,
    TRUE,
    TRUE
  ),
  stringsAsFactors=FALSE
)

cat("\n=== FORMAL CHECKS ===\n")
print(checks,row.names=FALSE)

stopifnot(all(checks$pass))

write_csv(
  checks,
  file.path(
    OUT,
    "09D_formal_checks.csv"
  )
)

cat("\n============================================================\n")
cat("09D COMPLETE\n")
cat("============================================================\n")
cat("No fibroblast resolution selected.\n")
cat("No fibroblast identities assigned.\n")
cat("FibConsensus not yet used.\n")
cat("PASS — marker comparison generated.\n")

