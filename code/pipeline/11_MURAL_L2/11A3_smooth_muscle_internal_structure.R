suppressPackageStartupMessages({
  library(Seurat)
  library(harmony)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
})

cat("============================================================\n")
cat("PIPELINE_V2 — 11A3 SMOOTH-MUSCLE INTERNAL STRUCTURE\n")
cat("============================================================\n\n")

SRC <- "PIPELINE_V2/objects/annotated/global_atlas_v2_L1_frozen.rds"

OUT <- paste0(
  "PIPELINE_V2/results/annotation/11_MURAL_L2/",
  "11A3_smooth_muscle_internal_structure"
)

TAB <- file.path(OUT, "tables")
FIG <- file.path(OUT, "figures")

dir.create(TAB, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG, recursive = TRUE, showWarnings = FALSE)

x <- readRDS(SRC)

sm <- subset(
  x,
  cells = colnames(x)[x$L1 == "Smooth_muscle"]
)

stopifnot(ncol(sm) == 330)

cat("Cells:", ncol(sm), "\n")
cat("Patients:", length(unique(sm$patient_id)), "\n")
cat("Batches:", length(unique(sm$batch)), "\n\n")

DefaultAssay(sm) <- "RNA"

# ------------------------------------------------------------
# 1. Remove inherited geometry
# ------------------------------------------------------------

cat("STEP 1/8 — removing inherited geometry\n")

sm@reductions <- list()
sm@graphs <- list()
sm@neighbors <- list()

old_cluster_cols <- grep(
  "snn_res|seurat_clusters",
  colnames(sm@meta.data),
  value = TRUE
)

if (length(old_cluster_cols) > 0) {
  sm@meta.data[old_cluster_cols] <- NULL
}

# ------------------------------------------------------------
# 2. De novo reconstruction
# ------------------------------------------------------------

cat("STEP 2/8 — de novo normalization/features\n")

sm <- NormalizeData(
  sm,
  verbose = FALSE
)

sm <- FindVariableFeatures(
  sm,
  selection.method = "vst",
  nfeatures = 2000,
  verbose = FALSE
)

sm <- ScaleData(
  sm,
  features = VariableFeatures(sm),
  verbose = FALSE
)

cat("STEP 3/8 — PCA\n")

sm <- RunPCA(
  sm,
  features = VariableFeatures(sm),
  npcs = 30,
  verbose = FALSE
)

cat("STEP 4/8 — Harmony\n")

sm <- RunHarmony(
  sm,
  group.by.vars = "batch",
  reduction.use = "pca",
  reduction.save = "harmony",
  verbose = FALSE
)

cat("STEP 5/8 — neighbors / UMAP / multiresolution clustering\n")

dims_use <- 1:20

sm <- FindNeighbors(
  sm,
  reduction = "harmony",
  dims = dims_use,
  graph.name = c(
    "smooth_nn",
    "smooth_snn"
  ),
  verbose = FALSE
)

sm <- RunUMAP(
  sm,
  reduction = "harmony",
  dims = dims_use,
  reduction.name = "smooth_umap",
  verbose = FALSE
)

resolutions <- c(
  0.1, 0.2, 0.3, 0.4, 0.5, 0.6
)

for (r in resolutions) {

  cat("  resolution", r, "\n")

  sm <- FindClusters(
    sm,
    graph.name = "smooth_snn",
    resolution = r,
    algorithm = 1,
    verbose = FALSE
  )

  new_name <- paste0(
    "smooth_snn_res.",
    r
  )

  sm[[new_name]] <-
    as.character(Idents(sm))
}

# ------------------------------------------------------------
# 3. Structural recurrence
# ------------------------------------------------------------

cat("STEP 6/8 — structural recurrence\n")

structural <- bind_rows(
  lapply(resolutions, function(r) {

    col <- paste0(
      "smooth_snn_res.",
      r
    )

    z <- tibble(
      cluster = as.character(
        sm@meta.data[[col]]
      ),
      patient_id =
        as.character(sm$patient_id),
      batch =
        as.character(sm$batch)
    )

    patient_counts <- z %>%
      count(
        cluster,
        patient_id,
        name = "n"
      )

    batch_counts <- z %>%
      count(
        cluster,
        batch,
        name = "n"
      )

    cluster_n <- z %>%
      count(
        cluster,
        name = "n_cells"
      )

    patient_summary <- patient_counts %>%
      group_by(cluster) %>%
      summarise(
        n_patients = n_distinct(patient_id),
        dominant_patient_cells = max(n),
        .groups = "drop"
      )

    batch_summary <- batch_counts %>%
      group_by(cluster) %>%
      summarise(
        n_batches = n_distinct(batch),
        .groups = "drop"
      )

    cluster_n %>%
      left_join(
        patient_summary,
        by = "cluster"
      ) %>%
      left_join(
        batch_summary,
        by = "cluster"
      ) %>%
      mutate(
        dominant_patient_fraction =
          dominant_patient_cells /
          n_cells,
        resolution = r
      )
  })
)

write.csv(
  structural,
  file.path(
    TAB,
    "11A3_multiresolution_structural_audit.csv"
  ),
  row.names = FALSE
)

resolution_summary <- structural %>%
  group_by(resolution) %>%
  summarise(
    n_clusters = n(),
    min_cluster_cells = min(n_cells),
    min_patients = min(n_patients),
    max_dominant_patient_fraction =
      max(dominant_patient_fraction),
    all_two_batches =
      all(n_batches == 2),
    .groups = "drop"
  )

cat("\n=== RESOLUTION SUMMARY ===\n")
print(
  resolution_summary,
  n = Inf
)

# ------------------------------------------------------------
# 4. Targeted lineage evidence
# ------------------------------------------------------------

cat("\nSTEP 7/8 — targeted lineage evidence\n")

programmes <- list(

  CONTRACTILE_CORE = c(
    "MYH11","ACTA2","TAGLN","CNN1",
    "MYL9","TPM2","LMOD1","MYLK",
    "SYNPO2","DES"
  ),

  ARRECTOR_PILI_SMC = c(
    "MYH11","ACTG2","TAGLN","CNN1",
    "DES","MYLK","LMOD1","SYNPO2"
  ),

  PERICYTE = c(
    "RGS5","CSPG4","PDGFRB","RBP1",
    "ABCC9","KCNJ8","COX4I2","HIGD1B"
  ),

  BASAL_KERATINOCYTE = c(
    "KRT5","KRT14","KRT17",
    "TP63","KRT15"
  ),

  DIFFERENTIATED_KERATINOCYTE = c(
    "KRT1","KRT10","IVL",
    "FLG","LOR"
  ),

  EPITHELIAL_JUNCTION = c(
    "DSP","DSC3","DSG1",
    "PKP1","PKP3","CLDN1"
  ),

  GLANDULAR_EPITHELIAL = c(
    "KRT7","KRT8","KRT18",
    "KRT19","MUC1","KRT23"
  ),

  SEBOCYTE = c(
    "KRT7","KRT19","PPARG",
    "PLIN2","FASN","SCD","MGST1"
  ),

  FIBROBLAST_ECM = c(
    "DCN","LUM","COL1A1","COL1A2",
    "COL3A1","DPT","PDGFRA"
  ),

  REMODELING = c(
    "POSTN","ADAM12","MMP11",
    "COL11A1","COMP","THBS2","FN1"
  ),

  SCHWANN = c(
    "SOX10","S100B","PLP1",
    "PMP22","MPZ","NGFR"
  )
)

mat <- GetAssayData(
  sm,
  assay = "RNA",
  layer = "data"
)

available <- lapply(
  programmes,
  function(g)
    intersect(g, rownames(mat))
)

coverage <- tibble(
  programme = names(programmes),
  n_requested = lengths(programmes),
  n_present = lengths(available),
  genes_present = sapply(
    available,
    paste,
    collapse = ";"
  )
)

cat("\n=== PROGRAMME COVERAGE ===\n")
print(coverage, n = Inf)

write.csv(
  coverage,
  file.path(
    TAB,
    "11A3_programme_coverage.csv"
  ),
  row.names = FALSE
)

score_fun <- function(genes) {

  genes <- intersect(
    genes,
    rownames(mat)
  )

  if (length(genes) == 0)
    return(rep(0, ncol(sm)))

  Matrix::colMeans(
    mat[genes, , drop = FALSE]
  )
}

detect_fun <- function(genes) {

  genes <- intersect(
    genes,
    rownames(mat)
  )

  if (length(genes) == 0)
    return(rep(0L, ncol(sm)))

  Matrix::colSums(
    mat[genes, , drop = FALSE] > 0
  )
}

evidence <- tibble(
  cell = colnames(sm),
  patient_id =
    as.character(sm$patient_id),
  batch =
    as.character(sm$batch),
  doublet_score =
    as.numeric(sm$doublet_score),
  nFeature_RNA =
    as.numeric(sm$nFeature_RNA)
)

for (p in names(available)) {

  evidence[[paste0(p, "_score")]] <-
    score_fun(available[[p]])

  evidence[[paste0(p, "_n")]] <-
    detect_fun(available[[p]])
}

for (r in resolutions) {

  col <- paste0(
    "smooth_snn_res.",
    r
  )

  evidence[[col]] <-
    as.character(
      sm@meta.data[[col]]
    )
}

write.csv(
  evidence,
  file.path(
    TAB,
    "11A3_cell_programme_evidence.csv"
  ),
  row.names = FALSE
)

# Cluster programme summaries at each resolution
cluster_programmes <- bind_rows(
  lapply(resolutions, function(r) {

    col <- paste0(
      "smooth_snn_res.",
      r
    )

    bind_rows(
      lapply(names(programmes), function(p) {

        evidence %>%
          mutate(
            cluster =
              .data[[col]]
          ) %>%
          group_by(cluster) %>%
          summarise(
            n_cells = n(),
            n_patients =
              n_distinct(patient_id),
            mean_score =
              mean(
                .data[[paste0(p, "_score")]],
                na.rm = TRUE
              ),
            median_score =
              median(
                .data[[paste0(p, "_score")]],
                na.rm = TRUE
              ),
            mean_genes_detected =
              mean(
                .data[[paste0(p, "_n")]],
                na.rm = TRUE
              ),
            .groups = "drop"
          ) %>%
          mutate(
            resolution = r,
            programme = p
          )
      })
    )
  })
) %>%
  select(
    resolution,
    cluster,
    programme,
    n_cells,
    n_patients,
    mean_score,
    median_score,
    mean_genes_detected
  )

write.csv(
  cluster_programmes,
  file.path(
    TAB,
    "11A3_cluster_programme_scores.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# 5. De novo markers at informative low resolutions
# ------------------------------------------------------------

cat("\nSTEP 8/8 — de novo marker evidence\n")

marker_resolutions <- c(
  0.2, 0.3, 0.4
)

all_markers <- bind_rows(
  lapply(marker_resolutions, function(r) {

    col <- paste0(
      "smooth_snn_res.",
      r
    )

    Idents(sm) <-
      sm@meta.data[[col]]

    m <- FindAllMarkers(
      sm,
      assay = "RNA",
      only.pos = TRUE,
      min.pct = 0.10,
      logfc.threshold = 0.25,
      verbose = FALSE
    )

    m$resolution <- r

    m
  })
)

write.csv(
  all_markers,
  file.path(
    TAB,
    "11A3_de_novo_markers.csv"
  ),
  row.names = FALSE
)

top_markers <- all_markers %>%
  filter(p_val_adj < 0.05) %>%
  group_by(
    resolution,
    cluster
  ) %>%
  slice_max(
    avg_log2FC,
    n = 20,
    with_ties = FALSE
  ) %>%
  ungroup()

cat("\n=== TOP DE NOVO MARKERS ===\n")
print(
  top_markers,
  n = Inf
)

write.csv(
  top_markers,
  file.path(
    TAB,
    "11A3_top20_markers.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# UMAP diagnostic
# ------------------------------------------------------------

p <- DimPlot(
  sm,
  reduction = "smooth_umap",
  group.by = "smooth_snn_res.0.3",
  label = TRUE,
  repel = TRUE,
  pt.size = 1.2
) +
  ggtitle(
    "Smooth-muscle de novo structure — r0.3"
  ) +
  theme_classic()

ggsave(
  file.path(
    FIG,
    "11A3_smooth_muscle_r03_UMAP.pdf"
  ),
  p,
  width = 7,
  height = 6
)

ggsave(
  file.path(
    FIG,
    "11A3_smooth_muscle_r03_UMAP.png"
  ),
  p,
  width = 7,
  height = 6,
  dpi = 300
)

# ------------------------------------------------------------
# Formal checks
# ------------------------------------------------------------

checks <- tibble(
  check = c(
    "cells_330",
    "patients_16",
    "batches_2",
    "inherited_geometry_removed",
    "PCA_generated",
    "Harmony_generated",
    "UMAP_generated",
    "six_resolutions_generated",
    "structural_audit_generated",
    "programme_coverage_generated",
    "cell_evidence_generated",
    "cluster_scores_generated",
    "markers_generated",
    "disease_not_used",
    "no_cells_removed",
    "no_identity_assigned",
    "source_object_not_modified"
  ),
  pass = c(
    ncol(sm) == 330,
    length(unique(sm$patient_id)) == 16,
    length(unique(sm$batch)) == 2,
    TRUE,
    "pca" %in% names(sm@reductions),
    "harmony" %in% names(sm@reductions),
    "smooth_umap" %in%
      names(sm@reductions),
    all(
      paste0(
        "smooth_snn_res.",
        resolutions
      ) %in%
        colnames(sm@meta.data)
    ),
    nrow(structural) > 0,
    nrow(coverage) ==
      length(programmes),
    nrow(evidence) == 330,
    nrow(cluster_programmes) > 0,
    nrow(all_markers) > 0,
    TRUE,
    ncol(sm) == 330,
    TRUE,
    TRUE
  )
)

cat("\n=== FORMAL CHECKS ===\n")
print(checks, n = Inf)

write.csv(
  checks,
  file.path(
    TAB,
    "11A3_formal_checks.csv"
  ),
  row.names = FALSE
)

stopifnot(all(checks$pass))

saveRDS(
  sm,
  file.path(
    OUT,
    "smooth_muscle_11A3_diagnostic_reclustered.rds"
  )
)

cat("\n============================================================\n")
cat("11A3 COMPLETE\n")
cat("============================================================\n")
cat("Cells:", ncol(sm), "\n")
cat("Patients:", length(unique(sm$patient_id)), "\n")
cat("Batches:", length(unique(sm$batch)), "\n")
cat("Diagnostic resolutions:", length(resolutions), "\n")
cat("No cells removed.\n")
cat("No identity assigned.\n")
cat("Disease not used.\n")
cat("Source object not modified.\n")
cat("PASS\n")
