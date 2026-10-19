#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(harmony)
  library(dplyr)
  library(readr)
  library(tibble)
  library(ggplot2)
})

cat("============================================================\n")
cat("PIPELINE_V2 — 09H4B REBUILD FIBROBLAST AFTER C9 RETURN\n")
cat("============================================================\n\n")

# ============================================================
# Paths
# ============================================================

IN <- paste0(
  "PIPELINE_V2/objects/annotated/",
  "fibroblast_v2_09H4A_C9_removed.rds"
)

OUT_OBJECT <- paste0(
  "PIPELINE_V2/objects/annotated/",
  "fibroblast_v2_09H4B_rebuilt_4685.rds"
)

OUTDIR <- paste0(
  "PIPELINE_V2/results/annotation/09_FIBROBLAST_L2/",
  "09H4B_rebuild_after_C9_return"
)

TAB <- file.path(OUTDIR, "tables")
FIG <- file.path(OUTDIR, "figures")

dir.create(TAB, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(IN)) {
  stop("09H4A fibroblast object not found: ", IN)
}

# ============================================================
# STEP 1 — load post-return fibroblast universe
# ============================================================

cat("STEP 1/9 — loading post-return fibroblast universe\n")

fib <- readRDS(IN)

cat("Cells:", ncol(fib), "\n")
cat("Genes:", nrow(fib), "\n")
cat("Patients:", length(unique(fib$patient_id)), "\n\n")

stopifnot(ncol(fib) == 4685)
stopifnot(length(unique(fib$patient_id)) == 16)

source_cells <- colnames(fib)

# ============================================================
# STEP 2 — audit metadata and C9 removal
# ============================================================

cat("STEP 2/9 — auditing metadata and C9 removal\n")

required_meta <- c(
  "patient_id",
  "disease",
  "batch"
)

missing_meta <- setdiff(required_meta, colnames(fib@meta.data))

if (length(missing_meta) > 0) {
  stop(
    "Missing required metadata: ",
    paste(missing_meta, collapse = ", ")
  )
}

old_res03 <- "fibroblast_refined_snn_res.0.3"

if (old_res03 %in% colnames(fib@meta.data)) {
  old_c9_remaining <- sum(
    as.character(fib@meta.data[[old_res03]]) == "9",
    na.rm = TRUE
  )
} else {
  old_c9_remaining <- NA_integer_
}

cat("Old r0.3 C9 cells remaining:", old_c9_remaining, "\n\n")

if (!is.na(old_c9_remaining)) {
  stopifnot(old_c9_remaining == 0)
}

patient_counts <- fib@meta.data %>%
  count(patient_id, disease, batch, name = "n_cells") %>%
  arrange(disease, patient_id)

write_csv(
  patient_counts,
  file.path(TAB, "09H4B_cells_by_patient.csv")
)

print(patient_counts)

# ============================================================
# STEP 3 — reset inherited geometry
# ============================================================

cat("\nSTEP 3/9 — resetting inherited fibroblast geometry\n")

fib@reductions <- list()
fib@graphs <- list()
fib@neighbors <- list()

# Remove inherited clustering identities as active identity.
Idents(fib) <- factor(rep("fibroblast", ncol(fib)))

fib$fibroblast_09H4B_stage <- "POST_C9_RETURN_REBUILD"

cat("Inherited reductions removed.\n")
cat("Inherited graphs removed.\n")
cat("Inherited neighbours removed.\n")
cat("09H4A provenance retained in metadata/misc where available.\n\n")

# ============================================================
# STEP 4 — normalization / HVG / PCA
# ============================================================

cat("STEP 4/9 — normalization, variable features and PCA\n")

DefaultAssay(fib) <- "RNA"

fib <- NormalizeData(
  fib,
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  verbose = FALSE
)

fib <- FindVariableFeatures(
  fib,
  selection.method = "vst",
  nfeatures = 3000,
  verbose = FALSE
)

cat("Variable features:", length(VariableFeatures(fib)), "\n")

fib <- ScaleData(
  fib,
  features = VariableFeatures(fib),
  verbose = FALSE
)

fib <- RunPCA(
  fib,
  features = VariableFeatures(fib),
  npcs = 50,
  verbose = FALSE
)

cat("PCA dimensions:", ncol(Embeddings(fib, "pca")), "\n\n")

# ============================================================
# STEP 5 — Harmony / neighbours / UMAP
# ============================================================

cat("STEP 5/9 — Harmony, neighbours and UMAP\n")

cat("Harmony variable: batch\n")
cat("Disease is NOT used for integration.\n")
cat("FibConsensus is NOT used for integration.\n")

fib <- RunHarmony(
  object = fib,
  group.by.vars = "batch",
  reduction.use = "pca",
  dims.use = 1:30,
  verbose = FALSE
)

fib <- FindNeighbors(
  fib,
  reduction = "harmony",
  dims = 1:30,
  graph.name = c(
    "fibroblast_09H4B_nn",
    "fibroblast_09H4B_snn"
  ),
  verbose = FALSE
)

set.seed(20260929)

fib <- RunUMAP(
  fib,
  reduction = "harmony",
  dims = 1:30,
  reduction.name = "fibroblast_09H4B_umap",
  reduction.key = "FB09H4BUMAP_",
  seed.use = 20260929,
  verbose = FALSE
)

cat("Harmony dimensions used: 30\n")
cat("New fibroblast UMAP generated.\n\n")

# ============================================================
# STEP 6 — resolution sweep
# ============================================================

cat("STEP 6/9 — post-return resolution sweep\n")

resolutions <- c(
  0.2,
  0.3,
  0.4,
  0.5,
  0.6,
  0.8,
  1.0
)

resolution_cols <- character()

for (res in resolutions) {

  cat("Resolution", res, "... ")

  before_cols <- colnames(fib@meta.data)

  fib <- FindClusters(
    fib,
    graph.name = "fibroblast_09H4B_snn",
    resolution = res,
    algorithm = 1,
    random.seed = 20260929,
    verbose = FALSE
  )

  after_cols <- colnames(fib@meta.data)
  new_cols <- setdiff(after_cols, before_cols)

  if (length(new_cols) != 1) {
    stop(
      "Could not uniquely identify clustering column for resolution ",
      res,
      ". New columns: ",
      paste(new_cols, collapse = ", ")
    )
  }

  raw_col <- new_cols

  target_col <- paste0(
    "fibroblast_09H4B_snn_res.",
    res
  )

  fib@meta.data[[target_col]] <- fib@meta.data[[raw_col]]

  if (raw_col != target_col) {
    fib@meta.data[[raw_col]] <- NULL
  }

  resolution_cols <- c(resolution_cols, target_col)

  cat(
    length(unique(fib@meta.data[[target_col]])),
    "clusters\n"
  )
}

# ============================================================
# STEP 7 — structural support audit
# ============================================================

cat("\nSTEP 7/9 — structural support audit\n")

support_all <- list()

for (i in seq_along(resolutions)) {

  res <- resolutions[i]
  col <- resolution_cols[i]

  md <- fib@meta.data %>%
    rownames_to_column("cell") %>%
    mutate(
      cluster = as.character(.data[[col]])
    )

  tmp <- md %>%
    group_by(cluster) %>%
    summarise(
      n_cells = n(),
      n_patients = n_distinct(patient_id),
      n_HC = n_distinct(patient_id[disease == "HC"]),
      n_AA = n_distinct(patient_id[disease == "AA"]),
      n_FFA = n_distinct(patient_id[disease == "FFA"]),
      n_batches = n_distinct(batch),
      .groups = "drop"
    )

  patient_comp <- md %>%
    count(cluster, patient_id, name = "patient_cells") %>%
    group_by(cluster) %>%
    mutate(
      patient_fraction = patient_cells / sum(patient_cells)
    ) %>%
    arrange(cluster, desc(patient_cells)) %>%
    slice_head(n = 1) %>%
    ungroup() %>%
    select(
      cluster,
      dominant_patient = patient_id,
      dominant_patient_cells = patient_cells,
      dominant_patient_fraction = patient_fraction
    )

  tmp <- tmp %>%
    left_join(patient_comp, by = "cluster") %>%
    mutate(resolution = res)

  support_all[[as.character(res)]] <- tmp
}

support <- bind_rows(support_all)

write_csv(
  support,
  file.path(TAB, "09H4B_resolution_patient_support.csv")
)

summary_res <- support %>%
  group_by(resolution) %>%
  summarise(
    n_clusters = n(),
    min_cluster_cells = min(n_cells),
    median_cluster_cells = median(n_cells),
    max_cluster_cells = max(n_cells),
    min_patients = min(n_patients),
    median_patients = median(n_patients),
    clusters_lt5_patients = sum(n_patients < 5),
    clusters_single_batch = sum(n_batches < 2),
    clusters_patient_gt50 =
      sum(dominant_patient_fraction > 0.50),
    max_patient_fraction =
      max(dominant_patient_fraction),
    .groups = "drop"
  )

write_csv(
  summary_res,
  file.path(TAB, "09H4B_resolution_summary.csv")
)

cat("\n=== POST-RETURN RESOLUTION SUMMARY ===\n")
print(summary_res)

# ============================================================
# STEP 8 — UMAP figures
# ============================================================

cat("\nSTEP 8/9 — figures\n")

p03 <- DimPlot(
  fib,
  reduction = "fibroblast_09H4B_umap",
  group.by = "fibroblast_09H4B_snn_res.0.3",
  label = TRUE,
  repel = TRUE,
  pt.size = 0.45
) +
  ggtitle("Fibroblast post-C9-return rebuild: resolution 0.3") +
  theme_classic()

p04 <- DimPlot(
  fib,
  reduction = "fibroblast_09H4B_umap",
  group.by = "fibroblast_09H4B_snn_res.0.4",
  label = TRUE,
  repel = TRUE,
  pt.size = 0.45
) +
  ggtitle("Fibroblast post-C9-return rebuild: resolution 0.4") +
  theme_classic()

ggsave(
  file.path(FIG, "09H4B_UMAP_res03.pdf"),
  p03,
  width = 8,
  height = 7
)

ggsave(
  file.path(FIG, "09H4B_UMAP_res03.png"),
  p03,
  width = 8,
  height = 7,
  dpi = 300
)

ggsave(
  file.path(FIG, "09H4B_UMAP_res04.pdf"),
  p04,
  width = 8,
  height = 7
)

ggsave(
  file.path(FIG, "09H4B_UMAP_res04.png"),
  p04,
  width = 8,
  height = 7,
  dpi = 300
)

# ============================================================
# STEP 9 — save + formal audit
# ============================================================

cat("\nSTEP 9/9 — save and formal audit\n")

fib$fibroblast_analysis_stage <-
  "09H4B_POST_C9_RETURN_REBUILT"

fib@misc[["09H4B"]] <- list(
  source_object = IN,
  source_cells = 4685,
  returned_C9_cells = 29,
  normalization = "LogNormalize",
  variable_features = 3000,
  PCA_dimensions = 50,
  harmony_variable = "batch",
  harmony_dimensions = 30,
  disease_used_for_integration = FALSE,
  FibConsensus_used_for_geometry = FALSE,
  resolution_selected = FALSE,
  L2_labels_assigned = FALSE,
  cells_removed_in_09H4B = 0,
  random_seed = 20260929
)

saveRDS(fib, OUT_OBJECT)

checks <- tibble(
  check = c(
    "source_cells_4685",
    "unique_cells_4685",
    "patients_16",
    "old_C9_absent",
    "HVG_3000",
    "PCA_50",
    "Harmony_present",
    "new_UMAP_present",
    "all_resolutions_present",
    "disease_not_used_for_integration",
    "FibConsensus_not_used_for_geometry",
    "no_resolution_selected",
    "no_L2_labels_assigned",
    "no_cells_removed"
  ),
  pass = c(
    ncol(fib) == 4685,
    length(unique(colnames(fib))) == 4685,
    length(unique(fib$patient_id)) == 16,
    is.na(old_c9_remaining) || old_c9_remaining == 0,
    length(VariableFeatures(fib)) == 3000,
    ncol(Embeddings(fib, "pca")) == 50,
    "harmony" %in% names(fib@reductions),
    "fibroblast_09H4B_umap" %in% names(fib@reductions),
    all(resolution_cols %in% colnames(fib@meta.data)),
    TRUE,
    TRUE,
    TRUE,
    TRUE,
    ncol(fib) == length(source_cells)
  )
)

write_csv(
  checks,
  file.path(TAB, "09H4B_formal_checks.csv")
)

audit <- c(
  "PIPELINE_V2 — 09H4B REBUILD AFTER C9 RETURN",
  "============================================",
  "",
  paste0("Input object: ", IN),
  paste0("Output object: ", OUT_OBJECT),
  "",
  "Source fibroblast cells: 4,685",
  "C9 cells previously returned to epithelial: 29",
  "Patients: 16",
  "",
  "Inherited fibroblast geometry was discarded.",
  "Normalization and variable-feature selection were rerun.",
  "PCA was rerun with 50 components.",
  "Harmony was rerun using batch only.",
  "Disease was not used for integration.",
  "FibConsensus was not used for geometry construction.",
  "Neighbours and UMAP were rebuilt de novo.",
  "Clustering was rerun across resolutions 0.2–1.0.",
  "",
  "No resolution was selected in 09H4B.",
  "No fibroblast L2 identity was assigned.",
  "No additional cells were removed.",
  "",
  "PASS — post-C9-return fibroblast universe rebuilt."
)

writeLines(
  audit,
  file.path(TAB, "09H4B_audit.txt")
)

cat("\n=== FORMAL CHECKS ===\n")
print(checks)

if (!all(checks$pass)) {
  stop("09H4B formal checks failed.")
}

cat("\nSaved:", OUT_OBJECT, "\n")

cat("\n============================================================\n")
cat("09H4B COMPLETE\n")
cat("============================================================\n")
cat("Fibroblast cells: 4685\n")
cat("Patients: 16\n")
cat("C9 remains excluded from fibroblast universe.\n")
cat("Geometry rebuilt de novo.\n")
cat("Disease not used for integration.\n")
cat("FibConsensus not used for geometry.\n")
cat("No resolution selected.\n")
cat("No L2 identities assigned.\n")
cat("No additional cells removed.\n")
cat("PASS — post-C9-return fibroblast universe rebuilt.\n")
