suppressPackageStartupMessages({
  library(Seurat)
  library(harmony)
  library(dplyr)
  library(readr)
})

set.seed(1234)

cat("\n============================================================\n")
cat("PIPELINE_V2 — 07B EPITHELIAL DE NOVO GEOMETRY\n")
cat("============================================================\n\n")

IN <- "PIPELINE_V2/objects/annotated/global_atlas_v2_L1_frozen.rds"

OUTROOT <- "PIPELINE_V2/results/annotation/07_EPITHELIAL_L2/07B_recluster"
TDIR <- file.path(OUTROOT,"tables")
OBJDIR <- "PIPELINE_V2/objects/annotated"

dir.create(TDIR,recursive=TRUE,showWarnings=FALSE)
dir.create(OBJDIR,recursive=TRUE,showWarnings=FALSE)

OUTOBJ <- file.path(
  OBJDIR,
  "epithelial_v2_07B_reclustered.rds"
)

# ------------------------------------------------------------
# 1. Load frozen L1 object
# ------------------------------------------------------------

cat("STEP 1/8 — loading L1-frozen atlas\n")

x <- readRDS(IN)
DefaultAssay(x) <- "RNA"

stopifnot(
  ncol(x)==82893,
  nrow(x)==23412,
  all(x$L1_frozen)
)

epi.cells <- colnames(x)[x$L1=="Epithelial"]

stopifnot(length(epi.cells)==63893)

epi <- subset(x,cells=epi.cells)

# Preserve global identity explicitly
epi$global_cluster <- as.character(
  epi$RNA_snn_res.0.6
)

cat("Epithelial cells:",ncol(epi),"\n")
cat("Genes:",nrow(epi),"\n")
cat("Patients:",length(unique(epi$patient_id)),"\n")
cat("Batches:",paste(unique(epi$batch),collapse=", "),"\n\n")

# ------------------------------------------------------------
# 2. Remove inherited reductions/graphs
#
# We retain metadata and expression but construct a genuinely
# epithelial-specific geometry.
# ------------------------------------------------------------

cat("STEP 2/8 — resetting inherited geometry\n")

epi@reductions <- list()
epi@graphs <- list()
epi@neighbors <- list()

# Remove inherited clustering columns from active identity only;
# metadata columns themselves remain preserved.
Idents(epi) <- factor(epi$global_cluster)

cat("Inherited reductions removed from working subset.\n")
cat("Global cluster metadata preserved.\n\n")

# ------------------------------------------------------------
# 3. Normalize + epithelial HVGs
# ------------------------------------------------------------

cat("STEP 3/8 — normalization and variable features\n")

DefaultAssay(epi) <- "RNA"

epi <- NormalizeData(
  epi,
  normalization.method="LogNormalize",
  scale.factor=10000,
  verbose=FALSE
)

epi <- FindVariableFeatures(
  epi,
  selection.method="vst",
  nfeatures=3000,
  verbose=FALSE
)

hvg <- VariableFeatures(epi)

writeLines(
  hvg,
  file.path(TDIR,"07B_variable_features.txt")
)

cat("Variable features:",length(hvg),"\n")

# ------------------------------------------------------------
# 4. Scale + PCA
# ------------------------------------------------------------

cat("\nSTEP 4/8 — epithelial PCA\n")

epi <- ScaleData(
  epi,
  features=hvg,
  verbose=FALSE
)

epi <- RunPCA(
  epi,
  features=hvg,
  npcs=50,
  verbose=FALSE
)

pca.sd <- Stdev(epi,"pca")

pca.tab <- data.frame(
  PC=seq_along(pca.sd),
  stdev=pca.sd,
  variance=pca.sd^2,
  pct_variance=
    100*(pca.sd^2)/sum(pca.sd^2),
  cumulative_pct=
    cumsum(
      100*(pca.sd^2)/sum(pca.sd^2)
    )
)

write_csv(
  pca.tab,
  file.path(TDIR,"07B_PCA_variance.csv")
)

cat("PCA dimensions:",ncol(Embeddings(epi,"pca")),"\n")

# ------------------------------------------------------------
# 5. Harmony by technical batch ONLY
# ------------------------------------------------------------

cat("\nSTEP 5/8 — Harmony integration by batch\n")

stopifnot(
  "batch" %in% colnames(epi@meta.data),
  "disease" %in% colnames(epi@meta.data)
)

cat("Harmony variable: batch\n")
cat("Disease is NOT used for integration.\n")

epi <- RunHarmony(
  epi,
  group.by.vars="batch",
  reduction.use="pca",
  dims.use=1:30,
  assay.use="RNA",
  verbose=FALSE
)

stopifnot(
  "harmony" %in% names(epi@reductions)
)

cat(
  "Harmony dimensions:",
  ncol(Embeddings(epi,"harmony")),
  "\n"
)

# ------------------------------------------------------------
# 6. Neighbors + UMAP
# ------------------------------------------------------------

cat("\nSTEP 6/8 — neighbors and UMAP\n")

epi <- FindNeighbors(
  epi,
  reduction="harmony",
  dims=1:30,
  graph.name=c(
    "epi_nn",
    "epi_snn"
  ),
  verbose=FALSE
)

epi <- RunUMAP(
  epi,
  reduction="harmony",
  dims=1:30,
  reduction.name="epi_umap",
  reduction.key="epiUMAP_",
  seed.use=1234,
  verbose=FALSE
)

# ------------------------------------------------------------
# 7. Resolution sweep
# ------------------------------------------------------------

cat("\nSTEP 7/8 — resolution sweep\n")

resolutions <- c(
  0.2,
  0.3,
  0.4,
  0.5,
  0.6,
  0.8,
  1.0
)

for(res in resolutions){

  cat("Resolution",res,"... ")

  epi <- FindClusters(
    epi,
    graph.name="epi_snn",
    resolution=res,
    algorithm=1,
    random.seed=1234,
    verbose=FALSE
  )

  src <- "seurat_clusters"

  dest <- paste0(
    "epi_snn_res.",
    format(
      res,
      trim=TRUE,
      scientific=FALSE
    )
  )

  epi[[dest]] <- as.character(
    epi@meta.data[[src]]
  )

  cat(
    length(unique(epi@meta.data[[dest]])),
    "clusters\n"
  )
}

# ------------------------------------------------------------
# 8. Audit tables + save
# ------------------------------------------------------------

cat("\nSTEP 8/8 — audit and save\n")

# Cluster number per resolution
resolution_summary <- bind_rows(
  lapply(resolutions,function(res){

    col <- paste0(
      "epi_snn_res.",
      format(res,trim=TRUE,scientific=FALSE)
    )

    tb <- table(epi@meta.data[[col]])

    data.frame(
      resolution=res,
      n_clusters=length(tb),
      min_cluster_cells=min(tb),
      median_cluster_cells=median(tb),
      max_cluster_cells=max(tb)
    )
  })
)

cat("\n=== RESOLUTION SUMMARY ===\n")
print(
  as.data.frame(resolution_summary),
  row.names=FALSE
)

write_csv(
  resolution_summary,
  file.path(
    TDIR,
    "07B_resolution_summary.csv"
  )
)

# Full patient/disease/batch composition for every resolution
composition <- bind_rows(
  lapply(resolutions,function(res){

    col <- paste0(
      "epi_snn_res.",
      format(res,trim=TRUE,scientific=FALSE)
    )

    epi@meta.data %>%
      mutate(
        epithelial_cluster=
          as.character(.data[[col]])
      ) %>%
      count(
        patient_id,
        disease,
        batch,
        epithelial_cluster,
        name="n_cells"
      ) %>%
      mutate(resolution=res)
  })
)

write_csv(
  composition,
  file.path(
    TDIR,
    "07B_resolution_patient_composition.csv"
  )
)

# Cluster support metrics
support <- composition %>%
  group_by(
    resolution,
    epithelial_cluster
  ) %>%
  summarise(
    n_cells=sum(n_cells),
    n_patients=n_distinct(patient_id),
    n_HC=n_distinct(
      patient_id[disease=="HC"]
    ),
    n_AA=n_distinct(
      patient_id[disease=="AA"]
    ),
    n_FFA=n_distinct(
      patient_id[disease=="FFA"]
    ),
    n_batches=n_distinct(batch),
    max_patient_fraction=
      max(n_cells)/sum(n_cells),
    .groups="drop"
  )

write_csv(
  support,
  file.path(
    TDIR,
    "07B_resolution_cluster_support.csv"
  )
)

# Cross-tab against original global clusters
global_mapping <- bind_rows(
  lapply(resolutions,function(res){

    col <- paste0(
      "epi_snn_res.",
      format(res,trim=TRUE,scientific=FALSE)
    )

    epi@meta.data %>%
      mutate(
        epithelial_cluster=
          as.character(.data[[col]])
      ) %>%
      count(
        global_cluster,
        epithelial_cluster,
        name="n_cells"
      ) %>%
      mutate(resolution=res)
  })
)

write_csv(
  global_mapping,
  file.path(
    TDIR,
    "07B_resolution_vs_global_cluster.csv"
  )
)

# Basic QC flags
flags <- support %>%
  mutate(
    flag_low_patient_support =
      n_patients < 3,
    flag_single_batch =
      n_batches < 2,
    flag_patient_dominated =
      max_patient_fraction > 0.50
  )

write_csv(
  flags,
  file.path(
    TDIR,
    "07B_resolution_support_flags.csv"
  )
)

cat("\n=== SUPPORT FLAGS BY RESOLUTION ===\n")

flag_summary <- flags %>%
  group_by(resolution) %>%
  summarise(
    n_clusters=n(),
    low_patient_support=
      sum(flag_low_patient_support),
    single_batch=
      sum(flag_single_batch),
    patient_dominated=
      sum(flag_patient_dominated),
    .groups="drop"
  )

print(
  as.data.frame(flag_summary),
  row.names=FALSE
)

write_csv(
  flag_summary,
  file.path(
    TDIR,
    "07B_resolution_flag_summary.csv"
  )
)

# Final invariants
stopifnot(
  ncol(epi)==63893,
  nrow(epi)==23412,
  length(unique(epi$patient_id))==16,
  identical(
    sort(colnames(epi)),
    sort(epi.cells)
  )
)

epi$L2_analysis_stage <-
  "07B_resolution_sweep_unfrozen"

saveRDS(
  epi,
  OUTOBJ,
  compress=FALSE
)

cat("\nSaved:",OUTOBJ,"\n")

audit <- c(
  "PIPELINE_V2 — 07B epithelial de novo geometry",
  paste0("Generated: ",Sys.time()),
  "Input: global_atlas_v2_L1_frozen.rds",
  "Cells: 63893 epithelial cells",
  "Genes: 23412",
  "Expression: SoupX-corrected RNA",
  "Variable features: 3000",
  "PCA: 50 dimensions",
  "Harmony: batch only",
  "Harmony dimensions used downstream: 1-30",
  "Disease was not used for integration",
  "Neighbors: Harmony 1-30",
  "UMAP: Harmony 1-30",
  "Resolution sweep: 0.2,0.3,0.4,0.5,0.6,0.8,1.0",
  "No L2 resolution selected",
  "No L2 labels assigned",
  "Global cluster metadata preserved",
  "L1 annotation preserved",
  "Output is exploratory L2 object, not frozen"
)

writeLines(
  audit,
  file.path(TDIR,"07B_audit.txt")
)

cat("\n============================================================\n")
cat("07B COMPLETE\n")
cat("============================================================\n")
cat("Cells:",ncol(epi),"\n")
cat("No L2 resolution selected.\n")
cat("No L2 labels assigned.\n")
cat("PASS — epithelial geometry generated; resolution remains unfrozen.\n")
