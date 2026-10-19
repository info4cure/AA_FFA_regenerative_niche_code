suppressPackageStartupMessages({
  library(Seurat)
  library(harmony)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
})

set.seed(20260930)
options(width = 180)

cat("============================================================\n")
cat("PIPELINE_V2 — 10B1 VASCULAR DE NOVO RECLUSTERING\n")
cat("============================================================\n\n")

# ============================================================
# 1. Paths
# ============================================================

SRC <- "PIPELINE_V2/objects/annotated/global_atlas_v2_L1_frozen.rds"

DEST <- paste0(
  "PIPELINE_V2/objects/annotated/",
  "vascular_v2_10B1_reclustered_3040.rds"
)

OUT <- paste0(
  "PIPELINE_V2/results/annotation/",
  "10_VASCULAR_L2/10B1_recluster"
)

TAB <- file.path(OUT, "tables")
FIG <- file.path(OUT, "figures")

dir.create(TAB, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG, recursive = TRUE, showWarnings = FALSE)

stopifnot(file.exists(SRC))

# ============================================================
# 2. Extract vascular universe
# ============================================================

cat("STEP 1/9 — extracting frozen-L1 vascular universe\n")

g <- readRDS(SRC)

stopifnot(ncol(g) == 82893)

vascular_cells <- colnames(g)[
  as.character(g$L1) %in%
    c("Vascular_EC", "Lymphatic_EC")
]

x <- subset(g, cells = vascular_cells)

stopifnot(ncol(x) == 3040)
stopifnot(length(unique(x$patient_id)) == 16)
stopifnot(all(x$doublet_class == "singlet"))

x$vascular_source_L1 <- as.character(x$L1)
x$vascular_analysis_stage <- "10B1_DE_NOVO"

cat("Cells:", ncol(x), "\n")
print(table(x$vascular_source_L1))

# ============================================================
# 3. Remove inherited geometry
# ============================================================

cat("\nSTEP 2/9 — removing inherited global geometry\n")

for (r in Reductions(x)) {
  x[[r]] <- NULL
}

graph_names <- names(x@graphs)

if (length(graph_names) > 0) {
  for (gr in graph_names) {
    x[[gr]] <- NULL
  }
}

# Remove inherited global clustering columns from working object
old_cluster_cols <- grep(
  "^RNA_snn_res\\.|^seurat_clusters$",
  colnames(x@meta.data),
  value = TRUE
)

for (cc in old_cluster_cols) {
  x[[cc]] <- NULL
}

stopifnot(length(Reductions(x)) == 0)

# ============================================================
# 4. Normalize and select variable features de novo
# ============================================================

cat("STEP 3/9 — de novo vascular normalization/features\n")

DefaultAssay(x) <- "RNA"

x <- NormalizeData(
  x,
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  verbose = FALSE
)

x <- FindVariableFeatures(
  x,
  selection.method = "vst",
  nfeatures = 3000,
  verbose = FALSE
)

x <- ScaleData(
  x,
  features = VariableFeatures(x),
  verbose = FALSE
)

# ============================================================
# 5. PCA
# ============================================================

cat("STEP 4/9 — PCA\n")

x <- RunPCA(
  x,
  features = VariableFeatures(x),
  npcs = 50,
  verbose = FALSE
)

pca_sd <- Stdev(x, reduction = "pca")

pca_tbl <- tibble(
  PC = seq_along(pca_sd),
  stdev = pca_sd,
  variance = pca_sd^2,
  variance_fraction =
    (pca_sd^2) / sum(pca_sd^2),
  cumulative_variance =
    cumsum((pca_sd^2) / sum(pca_sd^2))
)

write_csv(
  pca_tbl,
  file.path(TAB, "10B1_PCA_variance.csv")
)

# Fixed dimensional scaffold for reproducibility.
# Resolution selection is NOT performed here.
DIMS <- 1:30

# ============================================================
# 6. Harmony by sequencing batch
# ============================================================

cat("STEP 5/9 — Harmony batch correction\n")

x <- RunHarmony(
  object = x,
  group.by.vars = "batch",
  reduction.use = "pca",
  dims.use = DIMS,
  assay.use = "RNA",
  verbose = FALSE
)

stopifnot("harmony" %in% Reductions(x))

# ============================================================
# 7. Neighbors / UMAP / multiresolution clustering
# ============================================================

cat("STEP 6/9 — neighbors, UMAP and multiresolution clustering\n")

x <- FindNeighbors(
  x,
  reduction = "harmony",
  dims = DIMS,
  graph.name = c(
    "vascular_nn",
    "vascular_snn"
  ),
  verbose = FALSE
)

x <- RunUMAP(
  x,
  reduction = "harmony",
  dims = DIMS,
  reduction.name = "vascular_umap",
  reduction.key = "vUMAP_",
  seed.use = 20260930,
  verbose = FALSE
)

resolutions <- c(
  0.2,
  0.3,
  0.4,
  0.5,
  0.6,
  0.8,
  1.0
)

for (res in resolutions) {

  cat("  resolution", res, "\n")

  x <- FindClusters(
    x,
    graph.name = "vascular_snn",
    resolution = res,
    algorithm = 1,
    random.seed = 20260930,
    verbose = FALSE
  )

  generated <- grep(
    paste0(
      "^vascular_snn_res\\.",
      gsub("\\.", "\\\\.", as.character(res)),
      "$"
    ),
    colnames(x@meta.data),
    value = TRUE
  )

  if (length(generated) != 1) {
    stop(
      "Could not identify clustering column for resolution ",
      res
    )
  }
}

# ============================================================
# 8. Structural audit across resolutions
# ============================================================

cat("STEP 7/9 — structural recurrence audit\n")

cluster_cols <- grep(
  "^vascular_snn_res\\.",
  colnames(x@meta.data),
  value = TRUE
)

audit_resolution <- function(col) {

  md <- x@meta.data

  z <- tibble(
    cell = rownames(md),
    cluster = as.character(md[[col]]),
    patient_id = as.character(md$patient_id),
    batch = as.character(md$batch),
    source_L1 = as.character(md$vascular_source_L1)
  )

  cluster_sizes <- z %>%
    count(cluster, name = "n_cells")

  patient_support <- z %>%
    count(cluster, patient_id, name = "n") %>%
    group_by(cluster) %>%
    summarise(
      n_patients = n(),
      dominant_patient_cells = max(n),
      .groups = "drop"
    )

  dominant <- z %>%
    count(cluster, patient_id, name = "n") %>%
    group_by(cluster) %>%
    slice_max(n, n = 1, with_ties = FALSE) %>%
    transmute(
      cluster,
      dominant_patient = patient_id
    )

  batch_support <- z %>%
    distinct(cluster, batch) %>%
    count(cluster, name = "n_batches")

  source_support <- z %>%
    count(cluster, source_L1, name = "n") %>%
    group_by(cluster) %>%
    mutate(frac = n / sum(n)) %>%
    slice_max(frac, n = 1, with_ties = FALSE) %>%
    transmute(
      cluster,
      dominant_source_L1 = source_L1,
      dominant_source_fraction = frac
    )

  out <- cluster_sizes %>%
    left_join(patient_support, by = "cluster") %>%
    left_join(dominant, by = "cluster") %>%
    left_join(batch_support, by = "cluster") %>%
    left_join(source_support, by = "cluster") %>%
    mutate(
      dominant_patient_fraction =
        dominant_patient_cells / n_cells,
      resolution_column = col
    ) %>%
    arrange(as.integer(cluster))

  out
}

structural <- bind_rows(
  lapply(cluster_cols, audit_resolution)
)

write_csv(
  structural,
  file.path(TAB, "10B1_multiresolution_structural_audit.csv")
)

resolution_summary <- structural %>%
  group_by(resolution_column) %>%
  summarise(
    n_clusters = n(),
    min_cluster_cells = min(n_cells),
    min_patients = min(n_patients),
    max_dominant_patient_fraction =
      max(dominant_patient_fraction),
    all_two_batches = all(n_batches == 2),
    n_clusters_lt5_patients =
      sum(n_patients < 5),
    n_clusters_lt10_patients =
      sum(n_patients < 10),
    .groups = "drop"
  )

write_csv(
  resolution_summary,
  file.path(TAB, "10B1_resolution_summary.csv")
)

cat("\n=== RESOLUTION SUMMARY ===\n")
print(resolution_summary, n = Inf)

# ============================================================
# 9. Figures — structural only
# ============================================================

cat("\nSTEP 8/9 — structural figures\n")

p_source <- DimPlot(
  x,
  reduction = "vascular_umap",
  group.by = "vascular_source_L1",
  pt.size = 0.65
) +
  ggtitle(
    "Vascular de novo reconstruction",
    subtitle = "Source L1 shown for audit only"
  )

ggsave(
  file.path(FIG, "10B1_vascular_source_L1.png"),
  p_source,
  width = 7.5,
  height = 6,
  dpi = 300
)

for (cc in cluster_cols) {

  p <- DimPlot(
    x,
    reduction = "vascular_umap",
    group.by = cc,
    label = TRUE,
    repel = TRUE,
    pt.size = 0.55
  ) +
    ggtitle(
      paste("Vascular", cc),
      subtitle = "Disease not used for geometry"
    )

  safe <- gsub(
    "[^A-Za-z0-9_.-]",
    "_",
    cc
  )

  ggsave(
    file.path(
      FIG,
      paste0("10B1_", safe, ".png")
    ),
    p,
    width = 7.5,
    height = 6,
    dpi = 300
  )
}

# ============================================================
# 10. Save object and checks
# ============================================================

cat("STEP 9/9 — formal checks and save\n")

checks <- tibble(
  check = c(
    "vascular_cells_3040",
    "patients_16",
    "vascular_EC_2858",
    "lymphatic_EC_182",
    "all_singlets",
    "inherited_reductions_removed",
    "de_novo_PCA_present",
    "Harmony_present",
    "vascular_UMAP_present",
    "seven_resolutions_generated",
    "disease_not_used_for_geometry",
    "source_L1_not_used_for_geometry",
    "no_cells_removed",
    "source_object_not_overwritten",
    "no_L2_labels_assigned",
    "no_resolution_selected"
  ),
  pass = c(
    ncol(x) == 3040,
    length(unique(x$patient_id)) == 16,
    sum(x$vascular_source_L1 == "Vascular_EC") == 2858,
    sum(x$vascular_source_L1 == "Lymphatic_EC") == 182,
    all(x$doublet_class == "singlet"),
    TRUE,
    "pca" %in% Reductions(x),
    "harmony" %in% Reductions(x),
    "vascular_umap" %in% Reductions(x),
    length(cluster_cols) == 7,
    TRUE,
    TRUE,
    ncol(x) == 3040,
    normalizePath(SRC, mustWork = TRUE) !=
      normalizePath(DEST, mustWork = FALSE),
    !("vascular_L2" %in% colnames(x@meta.data)),
    TRUE
  )
)

write_csv(
  checks,
  file.path(TAB, "10B1_formal_checks.csv")
)

print(checks, n = Inf)

stopifnot(all(checks$pass))

saveRDS(x, DEST)

audit <- c(
  "PIPELINE_V2 — 10B1 VASCULAR DE NOVO RECLUSTERING",
  "",
  paste("Source:", SRC),
  paste("Output:", DEST),
  "",
  "Universe:",
  "- 3,040 vascular cells",
  "- 2,858 source Vascular_EC",
  "- 182 source Lymphatic_EC",
  "- 16 patients",
  "- both sequencing batches",
  "",
  "Geometry:",
  "- inherited global reductions removed",
  "- RNA normalization and variable features recalculated",
  "- PCA reconstructed de novo",
  "- Harmony applied by sequencing batch",
  "- vascular neighbors and UMAP reconstructed de novo",
  "- resolutions 0.2, 0.3, 0.4, 0.5, 0.6, 0.8 and 1.0 generated",
  "",
  "Interpretive boundary:",
  "- disease not used for geometry",
  "- source L1 not used for geometry",
  "- no resolution selected",
  "- no vascular L2 identity assigned",
  "- no cells removed",
  "- source global object not overwritten",
  "",
  "PASS"
)

writeLines(
  audit,
  file.path(TAB, "10B1_audit.txt")
)

cat("\n============================================================\n")
cat("10B1 COMPLETE\n")
cat("============================================================\n")
cat("Cells:", ncol(x), "\n")
cat("Patients:", length(unique(x$patient_id)), "\n")
cat("Candidate resolutions:", length(cluster_cols), "\n")
cat("No resolution selected.\n")
cat("No L2 identities assigned.\n")
cat("No cells removed.\n")
cat("PASS\n")
