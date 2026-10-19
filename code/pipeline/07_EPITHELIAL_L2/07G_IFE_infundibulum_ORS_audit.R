suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(tidyr)
  library(readr)
})

cat("\n============================================================\n")
cat("PIPELINE_V2 — 07G IFE / INFUNDIBULUM / ORS BOUNDARY AUDIT\n")
cat("============================================================\n\n")

IN <- "PIPELINE_V2/objects/annotated/epithelial_v2_07B_reclustered.rds"
OUT <- "PIPELINE_V2/results/annotation/07_EPITHELIAL_L2/07G_IFE_ORS_boundary"
TDIR <- file.path(OUT,"tables")

dir.create(TDIR,recursive=TRUE,showWarnings=FALSE)

x <- readRDS(IN)
DefaultAssay(x) <- "RNA"

CC <- "epi_snn_res.0.5"

stopifnot(
  ncol(x)==63893,
  nrow(x)==23412,
  CC %in% colnames(x@meta.data)
)

x$epi_cluster_05 <- as.character(x@meta.data[[CC]])

cat("Cells:",ncol(x),"\n")
cat("Genes:",nrow(x),"\n")
cat("Clusters:",length(unique(x$epi_cluster_05)),"\n")

# ============================================================
# 1. DISCRIMINATING PROGRAMMES
#
# Basal_core is deliberately NOT treated as an anatomical
# discriminator. It identifies basal epithelial phenotype.
#
# Identity calls must integrate:
# - programme evidence
# - de novo cluster markers
# - patient reproducibility
# ============================================================

identity_programmes <- list(

  Basal_core = c(
    "KRT5","KRT14","TP63",
    "ITGA6","ITGB4","COL17A1"
  ),

  IFE_basal_support = c(
    "KRT5","KRT14","TP63",
    "ITGA6","ITGB4"
  ),

  IFE_suprabasal = c(
    "KRT1","KRT10","KRT2",
    "KRTDAP","DSC1","IVL","CALML5"
  ),

  IFE_terminal = c(
    "FLG","FLG2","CDSN",
    "SLURP1","SLURP2","ASPRV1"
  ),

  Infundibular_differentiation = c(
    "KRT1","KRT10","KRT17",
    "KRT79","KRT75","SOX9"
  ),

  ORS_support = c(
    "KRT5","KRT14","KRT15",
    "KRT17","KRT19","SOX9",
    "KRT75","COL17A1"
  ),

  Bulge_Isthmus = c(
    "KRT15","KRT19","SOX9","LHX2",
    "NFATC1","LGR5","LGR6",
    "FGF18","CD34"
  ),

  IRS_companion = c(
    "KRT25","KRT27","KRT28",
    "KRT71","KRT72","KRT73",
    "KRT74","KRT75","TCHH"
  ),

  Hair_shaft = c(
    "KRT31","KRT32","KRT33A","KRT33B",
    "KRT35","KRT36",
    "KRT81","KRT83","KRT85"
  )
)

state_programmes <- list(

  Stress = c(
    "KRT6A","KRT6B","KRT16","KRT17",
    "S100A7","S100A8","S100A9",
    "SERPINB3","SERPINB4"
  ),

  IFN_MHCII = c(
    "CXCL9","CXCL10","CXCL11",
    "ISG15","IFIT2","IFIT3","RSAD2",
    "HLA-DRA","HLA-DRB1",
    "HLA-DQA1","HLA-DPA1","CD74"
  ),

  Cycling = c(
    "MKI67","TOP2A","UBE2C",
    "RRM2","CENPF"
  )
)

programmes <- c(identity_programmes,state_programmes)

# ============================================================
# 2. COVERAGE
# ============================================================

coverage <- bind_rows(
  lapply(names(programmes),function(p){

    requested <- programmes[[p]]
    present <- intersect(requested,rownames(x))

    data.frame(
      programme=p,
      class=ifelse(
        p %in% names(identity_programmes),
        "identity",
        "state"
      ),
      n_present=length(present),
      n_requested=length(requested),
      coverage=length(present)/length(requested),
      missing=paste(
        setdiff(requested,present),
        collapse=";"
      )
    )
  })
)

cat("\n=== MARKER COVERAGE ===\n")
print(as.data.frame(coverage),row.names=FALSE)

write_csv(
  coverage,
  file.path(TDIR,"07G_marker_coverage.csv")
)

# ============================================================
# 3. EXPRESSION MATRIX
# ============================================================

genes <- intersect(
  unique(unlist(programmes)),
  rownames(x)
)

mat <- LayerData(
  x,
  assay="RNA",
  layer="data"
)[genes,,drop=FALSE]

md <- x@meta.data
md$cell <- rownames(md)
md$cluster <- x$epi_cluster_05

# ============================================================
# 4. CELL-LEVEL PROGRAMME SCORES
# ============================================================

score_list <- lapply(
  programmes,
  function(g){

    gg <- intersect(g,rownames(mat))

    Matrix::colMeans(
      mat[gg,,drop=FALSE]
    )
  }
)

score_mat <- do.call(cbind,score_list)
rownames(score_mat) <- colnames(x)

# ============================================================
# 5. CLUSTER-LEVEL EVIDENCE
# ============================================================

cluster_rows <- list()

for(cl in as.character(0:18)){

  cells <- md$cell[md$cluster==cl]

  for(p in names(programmes)){

    v <- score_mat[cells,p]

    cluster_rows[[length(cluster_rows)+1]] <-
      data.frame(
        cluster=cl,
        programme=p,
        class=ifelse(
          p %in% names(identity_programmes),
          "identity",
          "state"
        ),
        n_cells=length(cells),
        mean_programme=mean(v),
        median_programme=median(v),
        pct_positive=100*mean(v>0)
      )
  }
}

cluster_ev <- bind_rows(cluster_rows) %>%
  group_by(programme) %>%
  mutate(
    z={
      s <- sd(mean_programme,na.rm=TRUE)
      if(is.finite(s) && s>0){
        as.numeric(scale(mean_programme))
      } else {
        rep(0,n())
      }
    }
  ) %>%
  ungroup()

write_csv(
  cluster_ev,
  file.path(TDIR,"07G_cluster_programme_evidence.csv")
)

# ============================================================
# 6. GENE-LEVEL EVIDENCE
# ============================================================

gene_rows <- list()

for(cl in as.character(0:18)){

  cells <- md$cell[md$cluster==cl]

  for(g in genes){

    v <- mat[g,cells]

    gene_rows[[length(gene_rows)+1]] <-
      data.frame(
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
  file.path(TDIR,"07G_gene_cluster_evidence.csv")
)

# ============================================================
# 7. PATIENT-LEVEL REPRODUCIBILITY
# ============================================================

patient_rows <- list()

for(cl in as.character(0:18)){

  z <- md[md$cluster==cl,,drop=FALSE]

  for(pt in unique(z$patient_id)){

    cells <- z$cell[z$patient_id==pt]

    for(p in names(programmes)){

      patient_rows[[length(patient_rows)+1]] <-
        data.frame(
          cluster=cl,
          patient_id=pt,
          disease=unique(
            z$disease[z$patient_id==pt]
          )[1],
          batch=unique(
            z$batch[z$patient_id==pt]
          )[1],
          programme=p,
          n_cells=length(cells),
          mean_programme=
            mean(score_mat[cells,p])
        )
    }
  }
}

patient_ev <- bind_rows(patient_rows)

write_csv(
  patient_ev,
  file.path(TDIR,"07G_patient_programme_evidence.csv")
)

# ============================================================
# 8. ANATOMICAL CONTRAST METRICS
#
# These are evidence metrics, NOT automatic labels.
# ============================================================

wide <- cluster_ev %>%
  filter(class=="identity") %>%
  select(
    cluster,
    programme,
    mean_programme,
    z
  ) %>%
  pivot_wider(
    names_from=programme,
    values_from=c(mean_programme,z)
  )

wide <- wide %>%
  mutate(
    ORS_vs_IFE =
      z_ORS_support -
      z_IFE_suprabasal,

    Bulge_vs_ORS =
      z_Bulge_Isthmus -
      z_ORS_support,

    Terminal_vs_ORS =
      z_IFE_terminal -
      z_ORS_support,

    IRS_vs_hair =
      z_IRS_companion -
      z_Hair_shaft
  )

write_csv(
  wide,
  file.path(TDIR,"07G_anatomical_contrasts.csv")
)

# ============================================================
# 9. CLUSTER SUPPORT
# ============================================================

support <- md %>%
  count(
    cluster,
    patient_id,
    disease,
    batch,
    name="n_cells"
  ) %>%
  group_by(cluster) %>%
  summarise(
    total_cells=sum(n_cells),
    n_patients=n_distinct(patient_id),
    n_HC=n_distinct(patient_id[disease=="HC"]),
    n_AA=n_distinct(patient_id[disease=="AA"]),
    n_FFA=n_distinct(patient_id[disease=="FFA"]),
    n_batches=n_distinct(batch),
    dominant_patient=
      patient_id[which.max(n_cells)],
    dominant_patient_fraction=
      max(n_cells)/sum(n_cells),
    .groups="drop"
  )

write_csv(
  support,
  file.path(TDIR,"07G_cluster_support.csv")
)

# ============================================================
# 10. COMPACT REPORT
# ============================================================

cat("\n============================================================\n")
cat("ANATOMICAL BOUNDARY REPORT\n")
cat("============================================================\n")

focus_programmes <- c(
  "Basal_core",
  "IFE_basal_support",
  "IFE_suprabasal",
  "IFE_terminal",
  "Infundibular_differentiation",
  "ORS_support",
  "Bulge_Isthmus",
  "IRS_companion",
  "Hair_shaft"
)

for(cl in as.character(0:18)){

  cat("\n------------------------------------------------------------\n")

  ss <- support %>% filter(cluster==cl)

  cat(
    "CLUSTER",cl,
    "| n =",ss$total_cells,
    "| patients =",ss$n_patients,
    "| max patient =",ss$dominant_patient,
    sprintf("(%.1f%%)",100*ss$dominant_patient_fraction),
    "\n"
  )

  zz <- cluster_ev %>%
    filter(
      cluster==cl,
      programme %in% focus_programmes
    ) %>%
    arrange(desc(z))

  print(
    as.data.frame(
      zz %>%
        select(
          programme,
          z,
          mean_programme,
          pct_positive
        )
    ),
    row.names=FALSE
  )

  cc <- wide %>% filter(cluster==cl)

  cat(
    sprintf(
      "CONTRASTS: ORS-IFE=%.2f | Bulge-ORS=%.2f | Terminal-ORS=%.2f | IRS-Hair=%.2f\n",
      cc$ORS_vs_IFE,
      cc$Bulge_vs_ORS,
      cc$Terminal_vs_ORS,
      cc$IRS_vs_hair
    )
  )
}

# ============================================================
# 11. SPECIFIC CHECKS
# ============================================================

cat("\n============================================================\n")
cat("EXPECTED INTERNAL CONTROLS\n")
cat("============================================================\n")

for(cl in c("10","12","15","16","17","18")){

  cat("\nCLUSTER",cl,"\n")

  print(
    as.data.frame(
      cluster_ev %>%
        filter(
          cluster==cl,
          programme %in% focus_programmes
        ) %>%
        arrange(desc(z)) %>%
        select(
          programme,
          z,
          mean_programme,
          pct_positive
        )
    ),
    row.names=FALSE
  )
}

# ============================================================
# 12. AUDIT
# ============================================================

audit <- c(
  "PIPELINE_V2 — 07G IFE/infundibulum/ORS boundary audit",
  paste0("Generated: ",Sys.time()),
  "Input: epithelial_v2_07B_reclustered.rds",
  "Working resolution: 0.5",
  "Cells: 63893",
  "Purpose: distinguish IFE basal/suprabasal/terminal from infundibular, ORS and bulge/isthmus programmes",
  "Basal_core is treated as shared basal phenotype, not a lineage-specific discriminator",
  "Identity and state programmes remain separate",
  "Disease not used to define anatomical identity",
  "No cells removed",
  "No geometry recalculated",
  "No L2 labels assigned",
  "No L2 freeze performed"
)

writeLines(
  audit,
  file.path(TDIR,"07G_audit.txt")
)

cat("\n============================================================\n")
cat("07G COMPLETE\n")
cat("============================================================\n")
cat("No L2 labels assigned.\n")
cat("No L2 freeze performed.\n")
cat("PASS — anatomical boundary evidence generated.\n")
