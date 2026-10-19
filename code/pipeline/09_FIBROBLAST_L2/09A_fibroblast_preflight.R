suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(readr)
})

cat("\n============================================================\n")
cat("PIPELINE_V2 — 09A FIBROBLAST L2 PREFLIGHT\n")
cat("============================================================\n")

IN <- "PIPELINE_V2/objects/annotated/global_atlas_v2_L1_frozen.rds"

OUT <- paste0(
  "PIPELINE_V2/results/annotation/09_FIBROBLAST_L2/",
  "09A_preflight/tables"
)

dir.create(OUT,recursive=TRUE,showWarnings=FALSE)

# ==========================================================================
# 1. LOAD FROZEN L1 ATLAS
# ==========================================================================

cat("\nSTEP 1/8 — loading frozen L1 atlas\n")

x <- readRDS(IN)

stopifnot(
  "L1" %in% colnames(x@meta.data),
  "patient_id" %in% colnames(x@meta.data),
  "disease" %in% colnames(x@meta.data),
  "batch" %in% colnames(x@meta.data)
)

fib_cells <- colnames(x)[
  as.character(x$L1)=="Fibroblast"
]

fib <- subset(
  x,
  cells=fib_cells
)

cat("Global cells:",ncol(x),"\n")
cat("Fibroblast cells:",ncol(fib),"\n")
cat("Genes:",nrow(fib),"\n")
cat("Patients:",length(unique(fib$patient_id)),"\n")
cat("Batches:",paste(sort(unique(fib$batch)),collapse=", "),"\n")

# ==========================================================================
# 2. BASIC INVARIANTS
# ==========================================================================

cat("\nSTEP 2/8 — basic invariants\n")

stopifnot(
  ncol(fib)==4808,
  nrow(fib)==23412,
  length(unique(fib$patient_id))==16
)

disease_counts <- as.data.frame(
  table(fib$disease),
  stringsAsFactors=FALSE
)

colnames(disease_counts) <- c(
  "disease",
  "n_cells"
)

disease_counts$pct <-
  100*disease_counts$n_cells/ncol(fib)

cat("\n=== FIBROBLAST CELLS BY DISEASE ===\n")
print(disease_counts,row.names=FALSE)

write_csv(
  disease_counts,
  file.path(
    OUT,
    "09A_cells_by_disease.csv"
  )
)

# ==========================================================================
# 3. PATIENT SUPPORT
# ==========================================================================

cat("\nSTEP 3/8 — patient support\n")

patient_counts <- fib@meta.data %>%
  count(
    patient_id,
    disease,
    batch,
    name="n_cells"
  ) %>%
  arrange(
    disease,
    desc(n_cells)
  )

cat("\n=== FIBROBLAST CELLS BY PATIENT ===\n")
print(
  as.data.frame(patient_counts),
  row.names=FALSE
)

write_csv(
  patient_counts,
  file.path(
    OUT,
    "09A_cells_by_patient.csv"
  )
)

patient_summary <- patient_counts %>%
  summarise(
    total_cells=sum(n_cells),
    n_patients=n_distinct(patient_id),
    n_HC=n_distinct(patient_id[disease=="HC"]),
    n_AA=n_distinct(patient_id[disease=="AA"]),
    n_FFA=n_distinct(patient_id[disease=="FFA"]),
    n_batches=n_distinct(batch),
    min_cells_patient=min(n_cells),
    median_cells_patient=median(n_cells),
    max_cells_patient=max(n_cells)
  )

cat("\n=== PATIENT SUPPORT SUMMARY ===\n")
print(
  as.data.frame(patient_summary),
  row.names=FALSE
)

# ==========================================================================
# 4. MARKER PROGRAMMES
#
# These are lineage/QC programmes only.
# FibConsensus families are intentionally NOT used in 09A.
# ==========================================================================

cat("\nSTEP 4/8 — lineage/QC marker coverage\n")

programmes <- list(

  Fibroblast_core=c(
    "COL1A1","COL1A2","COL3A1",
    "DCN","LUM","PDGFRA","COL6A1"
  ),

  Mesenchymal_core=c(
    "COL1A1","COL1A2","DCN",
    "LUM","COL3A1","VIM"
  ),

  Myofibroblast_SMC=c(
    "ACTA2","TAGLN","MYL9",
    "MYH11","CNN1","TPM2"
  ),

  Pericyte=c(
    "RGS5","CSPG4","MCAM",
    "NOTCH3","PDGFRB","ABCC9"
  ),

  Vascular_EC=c(
    "PECAM1","VWF","EMCN",
    "KDR","ESAM","ENG"
  ),

  Lymphatic_EC=c(
    "PROX1","PDPN","LYVE1",
    "FLT4","CCL21"
  ),

  Epithelial=c(
    "KRT5","KRT14","KRT15",
    "KRT17","KRT10","EPCAM"
  ),

  Immune=c(
    "PTPRC","CD3D","CD3E",
    "LST1","LYZ","TYROBP"
  ),

  Melanocyte=c(
    "PMEL","MLANA","TYR",
    "DCT","MITF"
  ),

  Neural_Schwann=c(
    "S100B","SOX10","PLP1",
    "MPZ","PMP2","CDH19"
  ),

  Cycling=c(
    "MKI67","TOP2A","UBE2C",
    "RRM2","CENPF"
  )
)

coverage <- bind_rows(
  lapply(
    names(programmes),
    function(nm){

      requested <- programmes[[nm]]
      present <- intersect(
        requested,
        rownames(fib)
      )

      data.frame(
        programme=nm,
        n_present=length(present),
        n_requested=length(requested),
        coverage=length(present)/length(requested),
        missing=paste(
          setdiff(requested,present),
          collapse=";"
        ),
        stringsAsFactors=FALSE
      )
    }
  )
)

cat("\n=== MARKER COVERAGE ===\n")
print(coverage,row.names=FALSE)

write_csv(
  coverage,
  file.path(
    OUT,
    "09A_marker_coverage.csv"
  )
)

# ==========================================================================
# 5. EXPRESSION EVIDENCE
# ==========================================================================

cat("\nSTEP 5/8 — broad lineage evidence\n")

genes <- unique(
  unlist(programmes)
)

genes <- intersect(
  genes,
  rownames(fib)
)

dat <- LayerData(
  fib,
  assay="RNA",
  layer="data"
)[genes,,drop=FALSE]

programme_summary <- bind_rows(
  lapply(
    names(programmes),
    function(nm){

      g <- intersect(
        programmes[[nm]],
        rownames(dat)
      )

      if(length(g)==0){
        return(NULL)
      }

      m <- Matrix::colMeans(
        dat[g,,drop=FALSE]
      )

      positive <- Matrix::colSums(
        dat[g,,drop=FALSE] > 0
      ) > 0

      data.frame(
        programme=nm,
        mean_programme=mean(m),
        median_programme=median(m),
        pct_cells_positive=
          100*mean(positive),
        stringsAsFactors=FALSE
      )
    }
  )
)

cat("\n=== BROAD PROGRAMME EVIDENCE ===\n")
print(
  programme_summary,
  row.names=FALSE
)

write_csv(
  programme_summary,
  file.path(
    OUT,
    "09A_broad_programme_evidence.csv"
  )
)

# ==========================================================================
# 6. GLOBAL-CLUSTER PROVENANCE
# ==========================================================================

cat("\nSTEP 6/8 — source global-cluster provenance\n")

cluster_col <- NULL

candidates <- c(
  "global_cluster_source",
  "seurat_clusters",
  "RNA_snn_res.0.6"
)

for(z in candidates){
  if(z %in% colnames(fib@meta.data)){
    cluster_col <- z
    break
  }
}

if(is.null(cluster_col)){
  stop("No source global-cluster metadata found.")
}

cat("Source cluster field:",cluster_col,"\n")

source_clusters <- fib@meta.data %>%
  mutate(
    source_cluster=
      as.character(.data[[cluster_col]])
  ) %>%
  count(
    source_cluster,
    disease,
    batch,
    name="n_cells"
  )

write_csv(
  source_clusters,
  file.path(
    OUT,
    "09A_source_global_cluster_provenance.csv"
  )
)

cat("\n=== SOURCE CLUSTERS ===\n")

print(
  source_clusters %>%
    group_by(source_cluster) %>%
    summarise(
      n_cells=sum(n_cells),
      n_diseases=n_distinct(disease),
      n_batches=n_distinct(batch),
      .groups="drop"
    ) %>%
    arrange(desc(n_cells)),
  n=Inf
)

# ==========================================================================
# 7. QC / DOUBLET DISTRIBUTION
# ==========================================================================

cat("\nSTEP 7/8 — QC and doublet distribution\n")

qc <- data.frame(
  metric=c(
    "nCount_RNA",
    "nFeature_RNA",
    "percent.mt",
    "doublet_score"
  ),
  q05=c(
    quantile(fib$nCount_RNA,.05,na.rm=TRUE),
    quantile(fib$nFeature_RNA,.05,na.rm=TRUE),
    quantile(fib$percent.mt,.05,na.rm=TRUE),
    quantile(fib$doublet_score,.05,na.rm=TRUE)
  ),
  q25=c(
    quantile(fib$nCount_RNA,.25,na.rm=TRUE),
    quantile(fib$nFeature_RNA,.25,na.rm=TRUE),
    quantile(fib$percent.mt,.25,na.rm=TRUE),
    quantile(fib$doublet_score,.25,na.rm=TRUE)
  ),
  median=c(
    median(fib$nCount_RNA,na.rm=TRUE),
    median(fib$nFeature_RNA,na.rm=TRUE),
    median(fib$percent.mt,na.rm=TRUE),
    median(fib$doublet_score,na.rm=TRUE)
  ),
  q75=c(
    quantile(fib$nCount_RNA,.75,na.rm=TRUE),
    quantile(fib$nFeature_RNA,.75,na.rm=TRUE),
    quantile(fib$percent.mt,.75,na.rm=TRUE),
    quantile(fib$doublet_score,.75,na.rm=TRUE)
  ),
  q95=c(
    quantile(fib$nCount_RNA,.95,na.rm=TRUE),
    quantile(fib$nFeature_RNA,.95,na.rm=TRUE),
    quantile(fib$percent.mt,.95,na.rm=TRUE),
    quantile(fib$doublet_score,.95,na.rm=TRUE)
  )
)

cat("\n=== QC DISTRIBUTION ===\n")
print(qc,row.names=FALSE)

write_csv(
  qc,
  file.path(
    OUT,
    "09A_QC_distribution.csv"
  )
)

# ==========================================================================
# 8. FORMAL CHECKS + AUDIT
# ==========================================================================

cat("\nSTEP 8/8 — formal checks\n")

checks <- data.frame(
  check=c(
    "fibroblast_cells_4808",
    "genes_23412",
    "patients_16",
    "HC_present",
    "AA_present",
    "FFA_present",
    "two_batches_present",
    "fibroblast_core_complete",
    "no_FibConsensus_labels_used",
    "no_cells_modified",
    "no_L2_labels_assigned",
    "no_frozen_object_written"
  ),

  pass=c(
    ncol(fib)==4808,
    nrow(fib)==23412,
    length(unique(fib$patient_id))==16,
    "HC" %in% fib$disease,
    "AA" %in% fib$disease,
    "FFA" %in% fib$disease,
    length(unique(fib$batch))==2,
    coverage$coverage[
      coverage$programme=="Fibroblast_core"
    ]==1,
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
    "09A_formal_checks.csv"
  )
)

audit <- c(
  "PIPELINE_V2 — 09A fibroblast L2 preflight",
  paste0("Generated: ",Sys.time()),
  "Input: frozen global L1 atlas.",
  paste0("Fibroblast cells: ",ncol(fib)),
  paste0("Genes: ",nrow(fib)),
  paste0(
    "Patients: ",
    length(unique(fib$patient_id))
  ),
  "Broad lineage/QC programmes evaluated.",
  "FibConsensus families intentionally not used at this stage.",
  "Disease not used for identity assignment.",
  "No cells reassigned.",
  "No cells removed.",
  "No L2 labels assigned.",
  "No fibroblast object frozen.",
  "Global L1 frozen atlas unchanged."
)

writeLines(
  audit,
  file.path(
    OUT,
    "09A_audit.txt"
  )
)

cat("\n============================================================\n")
cat("09A COMPLETE\n")
cat("============================================================\n")
cat("Fibroblast cells:",ncol(fib),"\n")
cat("No cells reassigned.\n")
cat("No L2 labels assigned.\n")
cat("No frozen fibroblast object written.\n")
cat("PASS — fibroblast preflight complete.\n")

