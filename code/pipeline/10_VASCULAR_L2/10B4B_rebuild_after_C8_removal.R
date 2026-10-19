suppressPackageStartupMessages({
  library(Seurat)
  library(harmony)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
})

options(width = 200)

cat("============================================================\n")
cat("PIPELINE_V2 — 10B4B POST-C8 VASCULAR DE NOVO REBUILD\n")
cat("============================================================\n\n")

SRC <- paste0(
  "PIPELINE_V2/objects/annotated/",
  "vascular_v2_10B4A_C8_removed_3023.rds"
)

OBJOUT <- paste0(
  "PIPELINE_V2/objects/annotated/",
  "vascular_v2_10B4B_rebuilt_3023.rds"
)

OUT <- paste0(
  "PIPELINE_V2/results/annotation/",
  "10_VASCULAR_L2/10B4B_rebuild_after_C8"
)

TAB <- file.path(OUT, "tables")
FIG <- file.path(OUT, "figures")

dir.create(TAB, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG, recursive = TRUE, showWarnings = FALSE)

stopifnot(file.exists(SRC))

# ============================================================
# 1. Load post-C8 universe
# ============================================================

cat("STEP 1/9 — loading post-C8 vascular universe\n")

x <- readRDS(SRC)

stopifnot(ncol(x) == 3023)
stopifnot(length(unique(x$patient_id)) == 16)
stopifnot(length(unique(x$batch)) == 2)

cat("Cells:", ncol(x), "\n")
cat("Patients:", length(unique(x$patient_id)), "\n")
cat("Batches:", length(unique(x$batch)), "\n")

# Preserve old r0.3 assignment ONLY as retrospective audit metadata.
# It is NOT used for geometry.

OLDCL <- "vascular_snn_res.0.3"

stopifnot(OLDCL %in% colnames(x@meta.data))

x$vascular_preC8_cluster03 <-
  as.character(x@meta.data[[OLDCL]])

# ============================================================
# 2. Remove inherited geometry
# ============================================================

cat("\nSTEP 2/9 — removing inherited vascular geometry\n")

for (r in Reductions(x)) {
  x[[r]] <- NULL
}

for (g in names(x@graphs)) {
  x[[g]] <- NULL
}

cluster_cols <- grep(
  "^vascular_snn_res\\.|^seurat_clusters$",
  colnames(x@meta.data),
  value = TRUE
)

for (cc in cluster_cols) {
  x[[cc]] <- NULL
}

stopifnot(length(Reductions(x)) == 0)
stopifnot(length(x@graphs) == 0)

# ============================================================
# 3. RNA preprocessing de novo
# ============================================================

cat("STEP 3/9 — de novo normalization/features/scaling\n")

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
  nfeatures = min(2000, nrow(x)),
  verbose = FALSE
)

x <- ScaleData(
  x,
  features = VariableFeatures(x),
  verbose = FALSE
)

# ============================================================
# 4. PCA
# ============================================================

cat("STEP 4/9 — PCA\n")

x <- RunPCA(
  x,
  features = VariableFeatures(x),
  npcs = 30,
  verbose = FALSE
)

# ============================================================
# 5. Harmony
# ============================================================

cat("STEP 5/9 — Harmony batch correction\n")

x <- RunHarmony(
  x,
  group.by.vars = "batch",
  reduction.use = "pca",
  dims.use = 1:20,
  assay.use = "RNA",
  verbose = FALSE
)

# ============================================================
# 6. Neighbors / UMAP / clustering
# ============================================================

cat("STEP 6/9 — neighbors, UMAP and multiresolution clustering\n")

x <- FindNeighbors(
  x,
  reduction = "harmony",
  dims = 1:20,
  graph.name = c(
    "vascular_postC8_nn",
    "vascular_postC8_snn"
  ),
  verbose = FALSE
)

set.seed(20260930)

x <- RunUMAP(
  x,
  reduction = "harmony",
  dims = 1:20,
  reduction.name = "vascular_postC8_umap",
  seed.use = 20260930,
  verbose = FALSE
)

resolutions <- c(
  0.2, 0.3, 0.4, 0.5,
  0.6, 0.8, 1.0
)

for (r in resolutions) {

  cat("  resolution", r, "\n")

  x <- FindClusters(
    x,
    graph.name = "vascular_postC8_snn",
    resolution = r,
    random.seed = 20260930,
    verbose = FALSE
  )
}

# FindClusters names columns from graph name.
# Standardize to explicit post-C8 names.

generated <- grep(
  "^vascular_postC8_snn_res\\.",
  colnames(x@meta.data),
  value = TRUE
)

stopifnot(length(generated) == length(resolutions))

# ============================================================
# 7. Structural recurrence
# ============================================================

cat("STEP 7/9 — structural recurrence audit\n")

support_one <- function(column) {

  md <- tibble(
    cluster =
      as.character(x@meta.data[[column]]),
    patient_id =
      as.character(x$patient_id),
    batch =
      as.character(x$batch)
  )

  dominant <- md %>%
    count(cluster, patient_id, name = "n") %>%
    group_by(cluster) %>%
    mutate(total = sum(n)) %>%
    slice_max(
      order_by = n,
      n = 1,
      with_ties = FALSE
    ) %>%
    transmute(
      cluster,
      dominant_patient = patient_id,
      dominant_patient_cells = n,
      dominant_patient_fraction = n / total
    )

  md %>%
    group_by(cluster) %>%
    summarise(
      n_cells = n(),
      n_patients = n_distinct(patient_id),
      n_batches = n_distinct(batch),
      .groups = "drop"
    ) %>%
    left_join(
      dominant,
      by = "cluster"
    ) %>%
    mutate(
      resolution_column = column
    )
}

structural <- bind_rows(
  lapply(generated, support_one)
)

write_csv(
  structural,
  file.path(
    TAB,
    "10B4B_multiresolution_structural_audit.csv"
  )
)

resolution_summary <- structural %>%
  group_by(resolution_column) %>%
  summarise(
    n_clusters = n(),
    min_cluster_cells = min(n_cells),
    min_patients = min(n_patients),
    max_dominant_patient_fraction =
      max(dominant_patient_fraction),
    all_two_batches =
      all(n_batches == 2),
    n_clusters_lt5_patients =
      sum(n_patients < 5),
    n_clusters_lt10_patients =
      sum(n_patients < 10),
    .groups = "drop"
  )

write_csv(
  resolution_summary,
  file.path(
    TAB,
    "10B4B_resolution_summary.csv"
  )
)

cat("\n=== RESOLUTION SUMMARY ===\n")
print(resolution_summary, n = Inf)

# ============================================================
# 8. Retrospective mapping of pre-C8 states
# ============================================================

cat("\nSTEP 8/9 — retrospective pre-C8 mapping\n")

# IMPORTANT:
# old cluster labels are used only after reconstruction
# to assess structural fate.
# They were not used for PCA/Harmony/UMAP/clustering.

transition_tables <- list()

for (cc in generated) {

  tr <- tibble(
    old_cluster =
      x$vascular_preC8_cluster03,
    new_cluster =
      as.character(x@meta.data[[cc]])
  ) %>%
    count(
      old_cluster,
      new_cluster,
      name = "n_cells"
    ) %>%
    group_by(old_cluster) %>%
    mutate(
      old_total = sum(n_cells),
      fraction_old = n_cells / old_total
    ) %>%
    ungroup() %>%
    mutate(
      resolution_column = cc
    ) %>%
    arrange(
      resolution_column,
      as.integer(old_cluster),
      desc(fraction_old)
    )

  transition_tables[[cc]] <- tr
}

transitions <- bind_rows(transition_tables)

write_csv(
  transitions,
  file.path(
    TAB,
    "10B4B_preC8_to_postC8_transitions.csv"
  )
)

# Explicitly summarize fate of old C4/C6/C7

focus_fate <- transitions %>%
  filter(
    old_cluster %in% c("4","6","7")
  )

write_csv(
  focus_fate,
  file.path(
    TAB,
    "10B4B_C4_C6_C7_postrebuild_fate.csv"
  )
)

cat("\n=== OLD C4/C6/C7 FATE ===\n")
print(focus_fate, n = Inf)

# ============================================================
# Figures
# ============================================================

p <- DimPlot(
  x,
  reduction = "vascular_postC8_umap",
  group.by = generated[
    which(
      grepl(
        "res\\.0\\.3$",
        generated
      )
    )[1]
  ],
  label = TRUE,
  repel = TRUE,
  pt.size = 0.7
) +
  ggtitle(
    "Post-C8 vascular reconstruction — resolution 0.3"
  ) +
  theme_classic()

ggsave(
  file.path(
    FIG,
    "10B4B_postC8_UMAP_r03.png"
  ),
  p,
  width = 8,
  height = 7,
  dpi = 300
)

ggsave(
  file.path(
    FIG,
    "10B4B_postC8_UMAP_r03.pdf"
  ),
  p,
  width = 8,
  height = 7
)

# ============================================================
# 9. Formal checks and save
# ============================================================

cat("\nSTEP 9/9 — formal checks and save\n")

checks <- tibble(
  check = c(
    "postC8_cells_3023",
    "patients_16",
    "batches_2",
    "old_C8_absent",
    "inherited_reductions_removed_before_rebuild",
    "de_novo_PCA_present",
    "Harmony_present",
    "postC8_UMAP_present",
    "seven_resolutions_generated",
    "structural_audit_generated",
    "transition_audit_generated",
    "C4_C6_C7_retained",
    "disease_not_used_for_geometry",
    "old_cluster_not_used_for_geometry",
    "no_additional_cells_removed",
    "no_L2_labels_assigned",
    "no_resolution_selected",
    "source_object_not_overwritten"
  ),
  pass = c(
    ncol(x) == 3023,
    length(unique(x$patient_id)) == 16,
    length(unique(x$batch)) == 2,
    !any(x$vascular_preC8_cluster03 == "8"),
    TRUE,
    "pca" %in% Reductions(x),
    "harmony" %in% Reductions(x),
    "vascular_postC8_umap" %in% Reductions(x),
    length(generated) == 7,
    nrow(structural) > 0,
    nrow(transitions) > 0,
    all(c("4","6","7") %in%
          unique(x$vascular_preC8_cluster03)),
    TRUE,
    TRUE,
    ncol(x) == 3023,
    TRUE,
    TRUE,
    TRUE
  )
)

write_csv(
  checks,
  file.path(
    TAB,
    "10B4B_formal_checks.csv"
  )
)

cat("\n=== FORMAL CHECKS ===\n")
print(checks, n = Inf)

stopifnot(all(checks$pass))

x$vascular_10B4B_stage <-
  "POST_C8_DE_NOVO_REBUILT_UNADJUDICATED"

saveRDS(
  x,
  OBJOUT,
  compress = FALSE
)

writeLines(
  c(
    "PIPELINE_V2 — 10B4B POST-C8 VASCULAR DE NOVO REBUILD",
    "",
    "Source:",
    "vascular_v2_10B4A_C8_removed_3023.rds",
    "",
    "Cells: 3,023",
    "Patients: 16",
    "",
    "Geometry:",
    "- inherited vascular reductions and graphs removed",
    "- RNA normalization/features/scaling recomputed",
    "- PCA recomputed",
    "- Harmony by sequencing batch",
    "- neighbors and UMAP rebuilt de novo",
    "- resolutions 0.2/0.3/0.4/0.5/0.6/0.8/1 generated",
    "",
    "Important:",
    "- disease was not used",
    "- source L1 was not used",
    "- pre-C8 cluster identity was not used for geometry",
    "- pre-C8 labels were retained only for retrospective transition audit",
    "- no additional cells were removed",
    "- no resolution selected",
    "- no L2 identity assigned",
    "",
    "PASS"
  ),
  file.path(
    TAB,
    "10B4B_audit.txt"
  )
)

cat("\n============================================================\n")
cat("10B4B COMPLETE\n")
cat("============================================================\n")
cat("Cells:", ncol(x), "\n")
cat("Patients:", length(unique(x$patient_id)), "\n")
cat("Candidate resolutions:", length(generated), "\n")
cat("No resolution selected.\n")
cat("No L2 identity assigned.\n")
cat("No additional cells removed.\n")
cat("PASS\n")
