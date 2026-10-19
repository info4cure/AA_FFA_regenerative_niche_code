suppressPackageStartupMessages({
  library(Seurat)
  library(harmony)
  library(dplyr)
  library(readr)
})

cat("\n============================================================\n")
cat("PIPELINE_V2 — 09F REBUILD REFINED FIBROBLAST UNIVERSE\n")
cat("============================================================\n")

INOBJ <- "PIPELINE_V2/objects/annotated/fibroblast_v2_09B_reclustered.rds"

LEDGER <- paste0(
  "PIPELINE_V2/results/annotation/09_FIBROBLAST_L2/",
  "09E1_final_adjudication/tables/09E1_master_adjudication.csv"
)

OUTOBJ <- "PIPELINE_V2/objects/annotated/fibroblast_v2_09F_refined.rds"

OUT <- paste0(
  "PIPELINE_V2/results/annotation/09_FIBROBLAST_L2/",
  "09F_refined_recluster/tables"
)

dir.create(OUT,recursive=TRUE,showWarnings=FALSE)
dir.create(dirname(OUTOBJ),recursive=TRUE,showWarnings=FALSE)

# ==========================================================================
# 1. LOAD
# ==========================================================================

cat("\nSTEP 1/9 — loading source fibroblast universe and ledger\n")

x <- readRDS(INOBJ)
ledger <- read_csv(LEDGER,show_col_types=FALSE)

stopifnot(
  ncol(x)==4808,
  nrow(x)==23412,
  nrow(ledger)==4808,
  length(unique(ledger$cell))==4808,
  setequal(colnames(x),ledger$cell)
)

cat("Source cells:",ncol(x),"\n")
cat("Genes:",nrow(x),"\n")

# ==========================================================================
# 2. ALIGN IMMUTABLE LEDGER
# ==========================================================================

cat("\nSTEP 2/9 — aligning adjudication ledger\n")

ledger <- ledger[
  match(colnames(x),ledger$cell),
]

stopifnot(
  identical(ledger$cell,colnames(x))
)

x$fibroblast_09E1_final_evidence <-
  ledger$final_evidence

x$fibroblast_09E1_disposition <-
  ledger$proposed_disposition

x$fibroblast_09E1_cluster03 <-
  ledger$cluster03

x$fibroblast_09E1_cluster04 <-
  ledger$cluster04

# ==========================================================================
# 3. FINAL DISPOSITION AUDIT
# ==========================================================================

cat("\nSTEP 3/9 — final disposition audit\n")

disp <- ledger %>%
  count(
    proposed_disposition,
    name="n_cells"
  ) %>%
  mutate(
    pct=100*n_cells/sum(n_cells)
  )

cat("\n=== SOURCE DISPOSITION ===\n")
print(as.data.frame(disp),row.names=FALSE)

stopifnot(
  sum(
    ledger$proposed_disposition ==
      "KEEP_FIBROBLAST_UNIVERSE"
  )==4714,

  sum(
    ledger$proposed_disposition ==
      "EXCLUDE_DOUBLET"
  )==85,

  sum(
    ledger$proposed_disposition ==
      "RETURN_EPITHELIAL"
  )==8,

  sum(
    ledger$proposed_disposition ==
      "RETURN_T"
  )==1
)

# ==========================================================================
# 4. CONSTRUCT REFINED OBJECT
# ==========================================================================

cat("\nSTEP 4/9 — constructing refined fibroblast object\n")

keep_cells <- ledger$cell[
  ledger$proposed_disposition ==
    "KEEP_FIBROBLAST_UNIVERSE"
]

fib <- subset(
  x,
  cells=keep_cells
)

stopifnot(
  ncol(fib)==4714,
  nrow(fib)==23412,
  length(unique(fib$patient_id))==16
)

cat("Refined fibroblast cells:",ncol(fib),"\n")
cat("Patients:",length(unique(fib$patient_id)),"\n")

patient_support <- fib@meta.data %>%
  count(
    patient_id,
    disease,
    batch,
    name="n_cells"
  )

cat("\n=== REFINED CELLS BY PATIENT ===\n")
print(as.data.frame(patient_support),row.names=FALSE)

write_csv(
  patient_support,
  file.path(OUT,"09F_refined_cells_by_patient.csv")
)

# ==========================================================================
# 5. RESET GEOMETRY
# ==========================================================================

cat("\nSTEP 5/9 — resetting inherited geometry\n")

fib@reductions <- list()
fib@graphs <- list()
fib@neighbors <- list()

drop_cols <- grep(
  "^fibroblast_snn_res\\.",
  colnames(fib@meta.data),
  value=TRUE
)

if(length(drop_cols)>0){
  fib@meta.data[,drop_cols] <- NULL
}

cat("Inherited fibroblast reductions removed.\n")
cat("Inherited fibroblast graphs removed.\n")
cat("09E1 provenance retained.\n")

# ==========================================================================
# 6. NORMALIZATION / HVG / PCA
# ==========================================================================

cat("\nSTEP 6/9 — normalization, variable features and PCA\n")

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

cat("Variable features:",length(VariableFeatures(fib)),"\n")

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
# 7. HARMONY + UMAP
# ==========================================================================

cat("\nSTEP 7/9 — Harmony, neighbors and UMAP\n")
cat("Harmony variable: batch\n")
cat("Disease is NOT used for integration.\n")
cat("FibConsensus is NOT used.\n")

fib <- RunHarmony(
  object=fib,
  group.by.vars="batch",
  reduction.use="pca",
  verbose=FALSE
)

hdims <- min(
  30,
  ncol(Embeddings(fib,"harmony"))
)

cat("Harmony dimensions used:",hdims,"\n")

fib <- FindNeighbors(
  fib,
  reduction="harmony",
  dims=1:hdims,
  graph.name=c(
    "fibroblast_refined_nn",
    "fibroblast_refined_snn"
  ),
  verbose=FALSE
)

fib <- RunUMAP(
  fib,
  reduction="harmony",
  dims=1:hdims,
  reduction.name="umap_fibroblast_refined",
  reduction.key="FIBREFUMAP_",
  verbose=FALSE
)

# ==========================================================================
# 8. RESOLUTION SWEEP
# ==========================================================================

cat("\nSTEP 8/9 — refined resolution sweep\n")

resolutions <- c(
  0.2,0.3,0.4,0.5,0.6,0.8,1.0
)

for(r in resolutions){

  cat("Resolution",r,"... ")

  cn <- paste0(
    "fibroblast_refined_snn_res.",
    r
  )

  fib <- FindClusters(
    fib,
    graph.name="fibroblast_refined_snn",
    resolution=r,
    algorithm=1,
    cluster.name=cn,
    verbose=FALSE
  )

  cat(
    length(unique(fib@meta.data[[cn]])),
    "clusters\n"
  )
}

support_list <- list()

for(r in resolutions){

  cn <- paste0(
    "fibroblast_refined_snn_res.",
    r
  )

  md <- fib@meta.data %>%
    mutate(
      cluster=as.character(.data[[cn]])
    )

  bp <- md %>%
    count(
      cluster,
      patient_id,
      disease,
      batch,
      name="n_cells"
    )

  s <- bp %>%
    group_by(cluster) %>%
    summarise(
      cluster_cells=sum(n_cells),
      n_patients=n_distinct(patient_id),
      n_HC=n_distinct(patient_id[disease=="HC"]),
      n_AA=n_distinct(patient_id[disease=="AA"]),
      n_FFA=n_distinct(patient_id[disease=="FFA"]),
      n_batches=n_distinct(batch),
      dominant_patient=
        patient_id[which.max(n_cells)],
      dominant_patient_cells=max(n_cells),
      dominant_patient_fraction=
        max(n_cells)/sum(n_cells),
      .groups="drop"
    ) %>%
    mutate(resolution=r)

  support_list[[as.character(r)]] <- s
}

support <- bind_rows(support_list)

summary_support <- support %>%
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
    clusters_patient_gt50=
      sum(dominant_patient_fraction>0.5),
    max_patient_fraction=
      max(dominant_patient_fraction),
    .groups="drop"
  )

cat("\n=== REFINED RESOLUTION SUMMARY ===\n")
print(as.data.frame(summary_support),row.names=FALSE)

write_csv(
  support,
  file.path(OUT,"09F_resolution_patient_support.csv")
)

write_csv(
  summary_support,
  file.path(OUT,"09F_resolution_summary.csv")
)

# ==========================================================================
# 9. SAVE / FORMAL AUDIT
# ==========================================================================

cat("\nSTEP 9/9 — save and formal audit\n")

fib$fibroblast_analysis_stage <-
  "09F_refined_geometry_pre_identity_projection"

checks <- data.frame(
  check=c(
    "source_cells_4808",
    "refined_cells_4714",
    "genes_23412",
    "patients_16",
    "HVG_3000",
    "PCA_50",
    "Harmony_present",
    "refined_UMAP_present",
    "all_resolutions_present",
    "no_FibConsensus_used",
    "disease_not_used_for_integration",
    "no_L2_labels_assigned"
  ),

  pass=c(
    ncol(x)==4808,
    ncol(fib)==4714,
    nrow(fib)==23412,
    length(unique(fib$patient_id))==16,
    length(VariableFeatures(fib))==3000,
    ncol(Embeddings(fib,"pca"))==50,
    "harmony" %in% names(fib@reductions),
    "umap_fibroblast_refined" %in%
      names(fib@reductions),

    all(
      paste0(
        "fibroblast_refined_snn_res.",
        resolutions
      ) %in% colnames(fib@meta.data)
    ),

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
  file.path(OUT,"09F_formal_checks.csv")
)

saveRDS(
  fib,
  OUTOBJ,
  compress=FALSE
)

audit <- c(
  "PIPELINE_V2 — 09F refined fibroblast universe",
  paste0("Generated: ",Sys.time()),
  "Source fibroblast cells: 4808.",
  "Retained fibroblast-universe cells: 4714.",
  "Excluded mixed high-doublet cells: 85.",
  "Returned epithelial candidates: 8.",
  "Returned T candidate: 1.",
  "Low-information cells retained.",
  "Mixed low-doublet cells retained conservatively.",
  "Geometry rebuilt de novo after adjudication.",
  "Variable features: 3000.",
  "PCA dimensions: 50.",
  paste0("Harmony dimensions used: 1-",hdims),
  "Harmony variable: batch.",
  "Disease not used for integration.",
  "FibConsensus not used.",
  "No fibroblast L2 identity assigned.",
  "No L2 freeze performed."
)

writeLines(
  audit,
  file.path(OUT,"09F_audit.txt")
)

cat("\nSaved:",OUTOBJ,"\n")

cat("\n============================================================\n")
cat("09F COMPLETE\n")
cat("============================================================\n")
cat("Source cells: 4808\n")
cat("Refined cells:",ncol(fib),"\n")
cat("No fibroblast resolution selected.\n")
cat("No fibroblast L2 labels assigned.\n")
cat("FibConsensus not yet used.\n")
cat("PASS — refined fibroblast universe rebuilt.\n")

