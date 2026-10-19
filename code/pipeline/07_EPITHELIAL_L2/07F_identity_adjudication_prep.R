suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(readr)
})

cat("\n============================================================\n")
cat("PIPELINE_V2 — 07F L2 IDENTITY ADJUDICATION PREP\n")
cat("============================================================\n\n")

IN <- "PIPELINE_V2/objects/annotated/epithelial_v2_07B_reclustered.rds"
OUT <- "PIPELINE_V2/results/annotation/07_EPITHELIAL_L2/07F_adjudication"
dir.create(OUT,recursive=TRUE,showWarnings=FALSE)

x <- readRDS(IN)
DefaultAssay(x) <- "RNA"
CC <- "epi_snn_res.0.5"

stopifnot(
  ncol(x)==63893,
  CC %in% colnames(x@meta.data)
)

# ------------------------------------------------------------
# Discriminating panel
# ------------------------------------------------------------

panels <- list(

  IFE_suprabasal=c(
    "KRT1","KRT10","KRT2","IVL",
    "KRTDAP","DSC1","CALML5"
  ),

  IFE_terminal=c(
    "FLG","FLG2","CDSN","LOR",
    "SLURP1","SLURP2","ASPRV1"
  ),

  Basal_ORS=c(
    "KRT5","KRT14","KRT15","KRT19",
    "COL17A1","ITGA6","ITGB4","TP63"
  ),

  Bulge_stem=c(
    "KRT15","SOX9","LHX2",
    "LGR5","LGR6","NFATC1"
  ),

  Isthmus_progenitor=c(
    "LGR5","LGR6","FGF18",
    "KRT15","KRT19","CD34"
  ),

  Hair_shaft=c(
    "KRT31","KRT32","KRT33A","KRT33B",
    "KRT35","KRT36","KRT81","KRT83","KRT85"
  ),

  IRS_companion=c(
    "KRT25","KRT27","KRT28",
    "KRT71","KRT72","KRT73","KRT74","KRT75",
    "TCHH"
  ),

  Sebaceous=c(
    "PPARG","FASN","SCD",
    "PLIN2","MGST1"
  ),

  IFN_MHCII=c(
    "CXCL9","CXCL10","CXCL11",
    "ISG15","IFIT2","IFIT3","RSAD2",
    "HLA-DRA","HLA-DRB1","HLA-DQA1",
    "HLA-DPA1","CD74"
  ),

  Stress=c(
    "KRT6A","KRT6B","KRT16","KRT17",
    "S100A7","S100A8","S100A9",
    "SERPINB3","SERPINB4"
  ),

  Cycling=c(
    "MKI67","TOP2A","UBE2C",
    "RRM2","CENPF"
  )
)

genes <- intersect(
  unique(unlist(panels)),
  rownames(x)
)

mat <- LayerData(
  x,
  assay="RNA",
  layer="data"
)[genes,,drop=FALSE]

md <- x@meta.data
md$cell <- rownames(md)
md$cluster <- as.character(md[[CC]])

# ------------------------------------------------------------
# 1. Gene-level cluster evidence
# ------------------------------------------------------------

gene_rows <- list()

for(cl in as.character(0:18)){

  cells <- md$cell[md$cluster==cl]

  for(g in genes){

    v <- mat[g,cells]

    gene_rows[[length(gene_rows)+1]] <- data.frame(
      cluster=cl,
      gene=g,
      mean_expression=mean(v),
      pct_detected=100*mean(v>0)
    )
  }
}

gene_ev <- bind_rows(gene_rows)

write_csv(
  gene_ev,
  file.path(OUT,"07F_gene_cluster_evidence.csv")
)

# ------------------------------------------------------------
# 2. Programme evidence
# ------------------------------------------------------------

programme_rows <- list()

for(cl in as.character(0:18)){

  cells <- md$cell[md$cluster==cl]

  for(p in names(panels)){

    gg <- intersect(panels[[p]],rownames(mat))

    v <- Matrix::colMeans(
      mat[gg,cells,drop=FALSE]
    )

    programme_rows[[length(programme_rows)+1]] <-
      data.frame(
        cluster=cl,
        programme=p,
        n_genes=length(gg),
        mean_programme=mean(v),
        pct_positive=100*mean(v>0)
      )
  }
}

prog <- bind_rows(programme_rows)

prog <- prog %>%
  group_by(programme) %>%
  mutate(
    z={
      s <- sd(mean_programme)
      if(is.finite(s) && s>0)
        as.numeric(scale(mean_programme))
      else rep(0,n())
    }
  ) %>%
  ungroup()

write_csv(
  prog,
  file.path(OUT,"07F_programme_cluster_evidence.csv")
)

# ------------------------------------------------------------
# 3. Patient reproducibility of programmes
# ------------------------------------------------------------

patient_rows <- list()

for(cl in as.character(0:18)){

  z <- md[md$cluster==cl,,drop=FALSE]

  for(pt in unique(z$patient_id)){

    cells <- z$cell[z$patient_id==pt]

    for(p in names(panels)){

      gg <- intersect(panels[[p]],rownames(mat))

      v <- Matrix::colMeans(
        mat[gg,cells,drop=FALSE]
      )

      patient_rows[[length(patient_rows)+1]] <-
        data.frame(
          cluster=cl,
          patient_id=pt,
          disease=unique(
            z$disease[z$patient_id==pt]
          )[1],
          programme=p,
          n_cells=length(cells),
          mean_programme=mean(v)
        )
    }
  }
}

patient_prog <- bind_rows(patient_rows)

write_csv(
  patient_prog,
  file.path(OUT,"07F_patient_programme_evidence.csv")
)

# ------------------------------------------------------------
# 4. Focused report
# ------------------------------------------------------------

focus <- c(
  "0","1","2","3","4","5","6","7","8",
  "9","10","11","12","13","14","15","16","17","18"
)

cat("============================================================\n")
cat("FOCUSED IDENTITY EVIDENCE\n")
cat("============================================================\n")

for(cl in focus){

  cat("\n------------------------------------------------------------\n")
  cat("CLUSTER",cl,"\n")

  pp <- prog %>%
    filter(cluster==cl) %>%
    arrange(desc(z))

  print(
    as.data.frame(
      pp %>%
        select(
          programme,
          z,
          mean_programme,
          pct_positive
        )
    ),
    row.names=FALSE
  )

  topgenes <- gene_ev %>%
    filter(cluster==cl) %>%
    arrange(desc(mean_expression)) %>%
    slice_head(n=20)

  cat("TOP PANEL GENES:\n")

  print(
    as.data.frame(topgenes),
    row.names=FALSE
  )
}

# ------------------------------------------------------------
# 5. Special patient-enriched C14
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("CLUSTER 14 PATIENT AUDIT\n")
cat("============================================================\n")

c14 <- md %>%
  filter(cluster=="14") %>%
  count(
    patient_id,
    disease,
    name="n_cells"
  ) %>%
  mutate(
    fraction=n_cells/sum(n_cells)
  ) %>%
  arrange(desc(n_cells))

print(
  as.data.frame(c14),
  row.names=FALSE
)

write_csv(
  c14,
  file.path(OUT,"07F_cluster14_patient_audit.csv")
)

audit <- c(
  "PIPELINE_V2 — 07F L2 identity adjudication prep",
  paste0("Generated: ",Sys.time()),
  "Working resolution: 0.5",
  "Cells: 63893",
  "Discriminating programmes combine V2 evidence with historical A09 marker framework",
  "Historical labels were not transferred",
  "Disease was not used to define identity",
  "No L2 labels assigned",
  "No resolution frozen"
)

writeLines(
  audit,
  file.path(OUT,"07F_audit.txt")
)

cat("\n============================================================\n")
cat("07F COMPLETE\n")
cat("============================================================\n")
cat("No labels assigned.\n")
cat("PASS — adjudication evidence generated.\n")
