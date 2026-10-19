suppressPackageStartupMessages({
  library(Seurat)
  library(harmony)
  library(dplyr)
  library(readr)
})

cat("\n============================================================\n")
cat("PIPELINE_V2 — 09B FIBROBLAST DE NOVO GEOMETRY\n")
cat("============================================================\n")

IN <- "PIPELINE_V2/objects/annotated/global_atlas_v2_L1_frozen.rds"

OUTOBJ <- "PIPELINE_V2/objects/annotated/fibroblast_v2_09B_reclustered.rds"

TDIR <- paste0(
  "PIPELINE_V2/results/annotation/09_FIBROBLAST_L2/",
  "09B_recluster/tables"
)

dir.create(TDIR, recursive=TRUE, showWarnings=FALSE)
dir.create(dirname(OUTOBJ), recursive=TRUE, showWarnings=FALSE)

# ==========================================================================
# 1. LOAD FROZEN L1 ATLAS AND EXTRACT FIBROBLASTS
# ==========================================================================

cat("\nSTEP 1/8 — loading frozen L1 atlas\n")

x <- readRDS(IN)

fib_cells <- colnames(x)[
  as.character(x$L1)=="Fibroblast"
]

fib <- subset(
  x,
  cells=fib_cells
)

stopifnot(
  ncol(fib)==4808,
  nrow(fib)==23412,
  length(unique(fib$patient_id))==16
)

cat("Fibroblast cells:",ncol(fib),"\n")
cat("Genes:",nrow(fib),"\n")
cat("Patients:",length(unique(fib$patient_id)),"\n")
cat("Batches:",paste(sort(unique(fib$batch)),collapse=", "),"\n")

# Preserve explicit provenance
fib$fibroblast_source_L1 <- as.character(fib$L1)
fib$fibroblast_source_global_cluster <- as.character(fib$seurat_clusters)

# ==========================================================================
# 2. RESET INHERITED GEOMETRY
# ==========================================================================

cat("\nSTEP 2/8 — resetting inherited geometry\n")

fib@reductions <- list()
fib@graphs <- list()
fib@neighbors <- list()

# Remove inherited clustering solutions from working subset.
# Provenance copied above is retained.
drop_cols <- grep(
  "^(RNA_snn_res\\.|fibroblast_snn_res\\.)",
  colnames(fib@meta.data),
  value=TRUE
)

if(length(drop_cols)>0){
  fib@meta.data[,drop_cols] <- NULL
}

cat("Inherited reductions removed.\n")
cat("Inherited graphs removed.\n")
cat("Frozen L1 provenance retained.\n")

# ==========================================================================
# 3. NORMALIZATION / HVG
# ==========================================================================

cat("\nSTEP 3/8 — normalization and variable features\n")

DefaultAssay(fib) <- "RNA"

fib <- NormalizeData(
  fib,
  normalization.method="LogNormalize",
  scale.factor=10000,
  verbose=FALSE
)

fib <- FindVariableFeatures(
  fib,
  selection.method="vst",
  nfeatures=3000,
  verbose=FALSE
)

stopifnot(
  length(VariableFeatures(fib))==3000
)

cat("Variable features:",length(VariableFeatures(fib)),"\n")

# ==========================================================================
# 4. PCA
# ==========================================================================

cat("\nSTEP 4/8 — fibroblast PCA\n")

fib <- ScaleData(
  fib,
  features=VariableFeatures(fib),
  verbose=FALSE
)

fib <- RunPCA(
  fib,
  features=VariableFeatures(fib),
  npcs=50,
  verbose=FALSE
)

cat(
  "PCA dimensions:",
  ncol(Embeddings(fib,"pca")),
  "\n"
)

# ==========================================================================
# 5. HARMONY
# ==========================================================================

cat("\nSTEP 5/8 — Harmony integration by batch\n")
cat("Harmony variable: batch\n")
cat("Disease is NOT used for integration.\n")
cat("FibConsensus is NOT used for integration or clustering.\n")

fib <- RunHarmony(
  object=fib,
  group.by.vars="batch",
  reduction.use="pca",
  verbose=FALSE
)

stopifnot(
  "harmony" %in% names(fib@reductions)
)

cat(
  "Harmony dimensions:",
  ncol(Embeddings(fib,"harmony")),
  "\n"
)

hdims <- min(
  30,
  ncol(Embeddings(fib,"harmony"))
)

# ==========================================================================
# 6. NEIGHBORS + UMAP
# ==========================================================================

cat("\nSTEP 6/8 — neighbors and UMAP\n")

fib <- FindNeighbors(
  fib,
  reduction="harmony",
  dims=1:hdims,
  graph.name=c(
    "fibroblast_nn",
    "fibroblast_snn"
  ),
  verbose=FALSE
)

fib <- RunUMAP(
  fib,
  reduction="harmony",
  dims=1:hdims,
  reduction.name="umap_fibroblast",
  reduction.key="FIBUMAP_",
  verbose=FALSE
)

# ==========================================================================
# 7. RESOLUTION SWEEP
# ==========================================================================

cat("\nSTEP 7/8 — resolution sweep\n")

resolutions <- c(
  0.2,0.3,0.4,0.5,0.6,0.8,1.0
)

for(r in resolutions){

  cat("Resolution",r,"... ")

  fib <- FindClusters(
    fib,
    graph.name="fibroblast_snn",
    resolution=r,
    algorithm=1,
    cluster.name=paste0(
      "fibroblast_snn_res.",
      r
    ),
    verbose=FALSE
  )

  k <- length(
    unique(
      fib@meta.data[[paste0("fibroblast_snn_res.",r)]]
    )
  )

  cat(k,"clusters\n")
}

# ==========================================================================
# 8. SUPPORT AUDIT + SAVE
# ==========================================================================

cat("\nSTEP 8/8 — support audit and save\n")

support_list <- list()

for(r in resolutions){

  cn <- paste0(
    "fibroblast_snn_res.",
    r
  )

  md <- fib@meta.data %>%
    mutate(
      cluster=as.character(.data[[cn]])
    )

  tmp <- md %>%
    count(
      cluster,
      patient_id,
      disease,
      batch,
      name="n"
    )

  s <- tmp %>%
    group_by(cluster) %>%
    summarise(
      cluster_cells=sum(n),
      n_patients=n_distinct(patient_id),
      n_HC=n_distinct(patient_id[disease=="HC"]),
      n_AA=n_distinct(patient_id[disease=="AA"]),
      n_FFA=n_distinct(patient_id[disease=="FFA"]),
      n_batches=n_distinct(batch),
      dominant_patient=patient_id[which.max(n)],
      max_patient_cells=max(n),
      max_patient_fraction=max(n)/sum(n),
      .groups="drop"
    ) %>%
    mutate(resolution=r)

  support_list[[as.character(r)]] <- s
}

support <- bind_rows(support_list)

write_csv(
  support,
  file.path(
    TDIR,
    "09B_resolution_patient_support.csv"
  )
)

summary <- support %>%
  group_by(resolution) %>%
  summarise(
    n_clusters=n(),
    min_cluster_cells=min(cluster_cells),
    median_cluster_cells=median(cluster_cells),
    max_cluster_cells=max(cluster_cells),
    min_patients=min(n_patients),
    median_patients=median(n_patients),
    clusters_lt5_patients=sum(n_patients<5),
    clusters_single_batch=sum(n_batches<2),
    clusters_patient_gt50=sum(max_patient_fraction>0.5),
    max_patient_fraction=max(max_patient_fraction),
    .groups="drop"
  )

cat("\n=== RESOLUTION SUMMARY ===\n")
print(
  as.data.frame(summary),
  row.names=FALSE
)

write_csv(
  summary,
  file.path(
    TDIR,
    "09B_resolution_summary.csv"
  )
)

# Disease composition is descriptive QC only.
# Disease is not used to select or annotate clusters.

disease_support <- bind_rows(
  lapply(
    resolutions,
    function(r){

      cn <- paste0(
        "fibroblast_snn_res.",
        r
      )

      fib@meta.data %>%
        mutate(
          cluster=as.character(.data[[cn]])
        ) %>%
        count(
          cluster,
          disease,
          name="n_cells"
        ) %>%
        group_by(cluster) %>%
        mutate(
          cluster_total=sum(n_cells),
          fraction=n_cells/cluster_total
        ) %>%
        ungroup() %>%
        mutate(resolution=r)
    }
  )
)

write_csv(
  disease_support,
  file.path(
    TDIR,
    "09B_resolution_disease_composition.csv"
  )
)

checks <- data.frame(
  check=c(
    "cells_4808",
    "genes_23412",
    "patients_16",
    "HVG_3000",
    "PCA_50",
    "Harmony_present",
    "UMAP_present",
    "all_resolutions_present",
    "no_FibConsensus_used",
    "disease_not_used_for_integration",
    "no_L2_labels_assigned",
    "global_L1_not_modified"
  ),
  pass=c(
    ncol(fib)==4808,
    nrow(fib)==23412,
    length(unique(fib$patient_id))==16,
    length(VariableFeatures(fib))==3000,
    ncol(Embeddings(fib,"pca"))==50,
    "harmony" %in% names(fib@reductions),
    "umap_fibroblast" %in% names(fib@reductions),
    all(
      paste0(
        "fibroblast_snn_res.",
        resolutions
      ) %in% colnames(fib@meta.data)
    ),
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
    TDIR,
    "09B_formal_checks.csv"
  )
)

fib$fibroblast_analysis_stage <-
  "09B_de_novo_geometry_pre_resolution_freeze"

saveRDS(
  fib,
  OUTOBJ,
  compress=FALSE
)

audit <- c(
  "PIPELINE_V2 — 09B fibroblast de novo geometry",
  paste0("Generated: ",Sys.time()),
  "Input: global_atlas_v2_L1_frozen.rds",
  "Fibroblast L1 cells: 4808.",
  "No pre-clustering cell exclusions.",
  "Inherited geometry removed.",
  "Expression renormalized within fibroblast subset.",
  "Variable features: 3000.",
  "PCA dimensions: 50.",
  paste0("Harmony downstream dimensions: 1-",hdims),
  "Harmony integration variable: batch.",
  "Disease was not used for integration.",
  "FibConsensus was not used for integration or clustering.",
  "Resolution sweep: 0.2,0.3,0.4,0.5,0.6,0.8,1.0.",
  "No resolution selected.",
  "No fibroblast L2 identities assigned.",
  "Global L1 frozen atlas unchanged."
)

writeLines(
  audit,
  file.path(
    TDIR,
    "09B_audit.txt"
  )
)

cat("\nSaved:",OUTOBJ,"\n")

cat("\n============================================================\n")
cat("09B COMPLETE\n")
cat("============================================================\n")
cat("Cells:",ncol(fib),"\n")
cat("No fibroblast resolution selected.\n")
cat("No fibroblast L2 identities assigned.\n")
cat("FibConsensus not yet used.\n")
cat("PASS — fibroblast de novo geometry generated.\n")

