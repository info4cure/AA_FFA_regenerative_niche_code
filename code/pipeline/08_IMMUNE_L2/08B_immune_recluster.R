suppressPackageStartupMessages({
  library(Seurat)
  library(harmony)
  library(dplyr)
  library(readr)
})

set.seed(1234)

cat("\n============================================================\n")
cat("PIPELINE_V2 — 08B IMMUNE DE NOVO GEOMETRY\n")
cat("============================================================\n\n")

IN <- "PIPELINE_V2/objects/annotated/global_atlas_v2_L1_frozen.rds"

OUTOBJ <- "PIPELINE_V2/objects/annotated/immune_v2_08B_reclustered.rds"

OUT <- "PIPELINE_V2/results/annotation/08_IMMUNE_L2/08B_recluster"
TDIR <- file.path(OUT,"tables")

dir.create(dirname(OUTOBJ),recursive=TRUE,showWarnings=FALSE)
dir.create(TDIR,recursive=TRUE,showWarnings=FALSE)

# ============================================================
# 1. LOAD + EXTRACT IMMUNE UNIVERSE
# ============================================================

cat("STEP 1/8 — loading L1-frozen atlas\n")

atlas <- readRDS(IN)
DefaultAssay(atlas) <- "RNA"

stopifnot(
  ncol(atlas)==82893,
  nrow(atlas)==23412,
  "L1" %in% colnames(atlas@meta.data)
)

immune_cells <- colnames(atlas)[
  as.character(atlas$L1) %in% c(
    "T_NK",
    "Myeloid_DC"
  )
]

stopifnot(length(immune_cells)==4762)

immune <- subset(
  atlas,
  cells=immune_cells
)

canonical_cells <- colnames(immune)

stopifnot(
  ncol(immune)==4762,
  all(
    as.character(immune$L1) %in%
      c("T_NK","Myeloid_DC")
  )
)

cat("Immune cells:",ncol(immune),"\n")
cat("Genes:",nrow(immune),"\n")
cat("Patients:",length(unique(immune$patient_id)),"\n")
cat("Batches:",
    paste(sort(unique(immune$batch)),collapse=", "),
    "\n")

# Preserve provenance explicitly
immune$immune_source_L1 <- as.character(immune$L1)

if("RNA_snn_res.0.6" %in% colnames(immune@meta.data)){
  immune$global_cluster_source <-
    as.character(immune$RNA_snn_res.0.6)
}

# ============================================================
# 2. RESET INHERITED GEOMETRY
# ============================================================

cat("\nSTEP 2/8 — resetting inherited geometry\n")

for(r in Reductions(immune)){
  immune[[r]] <- NULL
}

# remove inherited graphs
immune@graphs <- list()

cat("Inherited reductions removed.\n")
cat("Inherited graphs removed.\n")
cat("L1/global provenance metadata preserved.\n")

# ============================================================
# 3. NORMALIZATION + VARIABLE FEATURES
# ============================================================

cat("\nSTEP 3/8 — normalization and variable features\n")

DefaultAssay(immune) <- "RNA"

immune <- NormalizeData(
  immune,
  normalization.method="LogNormalize",
  scale.factor=10000,
  verbose=FALSE
)

immune <- FindVariableFeatures(
  immune,
  selection.method="vst",
  nfeatures=3000,
  verbose=FALSE
)

vf <- VariableFeatures(immune)

cat("Variable features:",length(vf),"\n")

writeLines(
  vf,
  file.path(TDIR,"08B_variable_features.txt")
)

# ============================================================
# 4. SCALE + PCA
# ============================================================

cat("\nSTEP 4/8 — immune PCA\n")

immune <- ScaleData(
  immune,
  features=vf,
  verbose=FALSE
)

immune <- RunPCA(
  immune,
  features=vf,
  npcs=50,
  verbose=FALSE
)

cat(
  "PCA dimensions:",
  ncol(Embeddings(immune,"pca")),
  "\n"
)

# ============================================================
# 5. HARMONY — BATCH ONLY
# ============================================================

cat("\nSTEP 5/8 — Harmony integration by batch\n")
cat("Harmony variable: batch\n")
cat("Disease is NOT used for integration.\n")

immune <- RunHarmony(
  object = immune,
  group.by.vars = "batch",
  reduction.use = "pca",
  verbose = FALSE
)

cat(
  "Harmony dimensions:",
  ncol(Embeddings(immune,"harmony")),
  "\n"
)

# ============================================================
# 6. NEIGHBORS + UMAP
# ============================================================

cat("\nSTEP 6/8 — neighbors and UMAP\n")

immune <- FindNeighbors(
  immune,
  reduction="harmony",
  dims=1:30,
  graph.name=c(
    "immune_nn",
    "immune_snn"
  ),
  verbose=FALSE
)

immune <- RunUMAP(
  immune,
  reduction="harmony",
  dims=1:30,
  reduction.name="immune_umap",
  seed.use=1234,
  verbose=FALSE
)

# ============================================================
# 7. RESOLUTION SWEEP
# ============================================================

cat("\nSTEP 7/8 — resolution sweep\n")

resolutions <- c(
  0.2,0.3,0.4,0.5,0.6,0.8,1.0
)

summary_rows <- list()
support_rows <- list()

for(res in resolutions){

  cat("Resolution",res,"... ")

  tmp <- FindClusters(
    immune,
    graph.name="immune_snn",
    resolution=res,
    random.seed=1234,
    verbose=FALSE
  )

  generated_col <- paste0(
    "immune_snn_res.",
    res
  )

  stopifnot(
    generated_col %in%
      colnames(tmp@meta.data)
  )

  cl <- as.character(
    tmp@meta.data[[generated_col]]
  )

  # preserve clustering in main object
  immune@meta.data[[generated_col]] <- cl

  tab <- table(cl)

  summary_rows[[length(summary_rows)+1]] <-
    data.frame(
      resolution=res,
      n_clusters=length(tab),
      min_cluster_cells=min(tab),
      median_cluster_cells=median(tab),
      max_cluster_cells=max(tab)
    )

  md <- tmp@meta.data
  md$cluster_working <- cl

  patient_support <- md %>%
    count(
      cluster_working,
      patient_id,
      disease,
      batch,
      name="n_cells"
    ) %>%
    group_by(cluster_working) %>%
    summarise(
      cluster_cells=sum(n_cells),
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
      max_patient_cells=max(n_cells),
      max_patient_fraction=
        max(n_cells)/sum(n_cells),
      .groups="drop"
    ) %>%
    mutate(
      resolution=res
    )

  support_rows[[length(support_rows)+1]] <-
    patient_support

  cat(length(tab),"clusters\n")
}

resolution_summary <- bind_rows(summary_rows)
support <- bind_rows(support_rows)

cat("\n=== RESOLUTION SUMMARY ===\n")
print(
  resolution_summary,
  row.names=FALSE
)

write_csv(
  resolution_summary,
  file.path(
    TDIR,
    "08B_resolution_summary.csv"
  )
)

write_csv(
  support,
  file.path(
    TDIR,
    "08B_patient_support_all_resolutions.csv"
  )
)

# ============================================================
# 8. COMPACT SUPPORT SUMMARY
# ============================================================

support_summary <- support %>%
  group_by(resolution) %>%
  summarise(
    n_clusters=n(),
    min_patients=min(n_patients),
    median_patients=median(n_patients),
    clusters_lt3_patients=
      sum(n_patients<3),
    clusters_lt5_patients=
      sum(n_patients<5),
    clusters_single_batch=
      sum(n_batches<2),
    clusters_patient_gt50=
      sum(max_patient_fraction>0.50),
    max_patient_fraction=
      max(max_patient_fraction),
    .groups="drop"
  )

cat("\n=== SUPPORT SUMMARY ===\n")
print(
  support_summary,
  row.names=FALSE
)

write_csv(
  support_summary,
  file.path(
    TDIR,
    "08B_support_summary.csv"
  )
)

# ============================================================
# SAVE EXPLORATORY OBJECT
# ============================================================

immune$immune_analysis_stage <-
  "08B_resolution_sweep_unfrozen"

stopifnot(
  ncol(immune)==4762,
  identical(
    colnames(immune),
    canonical_cells
  )
)

saveRDS(
  immune,
  OUTOBJ,
  compress=FALSE
)

cat("\nSaved:",OUTOBJ,"\n")

# ============================================================
# AUDIT
# ============================================================

audit <- c(
  "PIPELINE_V2 — 08B immune de novo geometry",
  paste0("Generated: ",Sys.time()),
  "Input: global_atlas_v2_L1_frozen.rds",
  "Immune universe: frozen T_NK + Myeloid_DC",
  "Cells: 4762",
  "Genes: 23412",
  "Patients: 16",
  "Expression: SoupX-corrected RNA",
  "Variable features: 3000",
  "PCA: 50 dimensions",
  "Harmony: batch only",
  "Harmony dimensions used downstream for neighbors/UMAP: 1-30",
  "Disease was not used for integration",
  "Neighbors: Harmony 1-30",
  "UMAP: Harmony 1-30",
  "Resolution sweep: 0.2,0.3,0.4,0.5,0.6,0.8,1.0",
  "No immune resolution selected",
  "No immune identities assigned",
  "No cells removed",
  "Frozen L1 provenance preserved",
  "Output is exploratory immune object, not frozen"
)

writeLines(
  audit,
  file.path(TDIR,"08B_audit.txt")
)

cat("\n============================================================\n")
cat("08B COMPLETE\n")
cat("============================================================\n")
cat("Cells:",ncol(immune),"\n")
cat("No immune resolution selected.\n")
cat("No immune identities assigned.\n")
cat("PASS — immune geometry generated; resolution remains unfrozen.\n")
