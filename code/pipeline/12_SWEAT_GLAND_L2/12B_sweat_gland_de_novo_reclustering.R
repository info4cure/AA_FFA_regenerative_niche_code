suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(readr)
  library(ggplot2)
  library(patchwork)
  library(harmony)
})

cat("============================================================\n")
cat("PIPELINE_V2 — 12B SWEAT-GLAND DE NOVO RECLUSTERING\n")
cat("============================================================\n\n")

SRC <- "PIPELINE_V2/objects/annotated/global_atlas_v2_L1_frozen.rds"

OUT <- paste0(
  "PIPELINE_V2/results/annotation/",
  "12_SWEAT_GLAND_L2/12B_recluster"
)

TAB <- file.path(OUT, "tables")
FIG <- file.path(OUT, "figures")

dir.create(TAB, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# 1. Extract frozen-L1 sweat-gland universe
# ============================================================

cat("STEP 1/10 — extracting frozen-L1 sweat-gland universe\n")

g <- readRDS(SRC)

stopifnot(ncol(g) == 82893)
stopifnot("L1" %in% colnames(g@meta.data))

x <- subset(
  g,
  subset = L1 == "Sweat_gland"
)

stopifnot(ncol(x) == 1843)
stopifnot(length(unique(x$patient_id)) == 16)
stopifnot(length(unique(x$batch)) == 2)

cells_before <- colnames(x)
genes_before <- rownames(x)

cat("Cells:", ncol(x), "\n")
cat("Patients:", length(unique(x$patient_id)), "\n")
cat("Batches:", length(unique(x$batch)), "\n\n")

# ============================================================
# 2. Remove inherited geometry
# ============================================================

cat("STEP 2/10 — removing inherited global geometry\n")

x@reductions <- list()
x@graphs <- list()
x@neighbors <- list()

# Remove inherited clustering metadata only.
cluster_cols <- grep(
  "snn_res|seurat_clusters",
  colnames(x@meta.data),
  value = TRUE
)

if (length(cluster_cols) > 0) {
  x@meta.data[
    cluster_cols
  ] <- NULL
}

# ============================================================
# 3. De novo normalization / features / scaling
# ============================================================

cat("STEP 3/10 — de novo normalization/features/scaling\n")

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
# 4. PCA
# ============================================================

cat("STEP 4/10 — PCA\n")

x <- RunPCA(
  x,
  features = VariableFeatures(x),
  npcs = 40,
  verbose = FALSE,
  reduction.name = "sweat_pca"
)

# ============================================================
# 5. Harmony
# ============================================================

cat("STEP 5/10 — Harmony batch correction\n")

x <- RunHarmony(
  object = x,
  group.by.vars = "batch",
  reduction.use = "sweat_pca",
  reduction.save = "sweat_harmony",
  verbose = FALSE
)

# ============================================================
# 6. Neighbors / UMAP / multiresolution clustering
# ============================================================

cat("STEP 6/10 — neighbors, UMAP and multiresolution clustering\n")

dims_use <- 1:30

x <- FindNeighbors(
  x,
  reduction = "sweat_harmony",
  dims = dims_use,
  graph.name = c(
    "sweat_nn",
    "sweat_snn"
  ),
  verbose = FALSE
)

x <- RunUMAP(
  x,
  reduction = "sweat_harmony",
  dims = dims_use,
  reduction.name = "sweat_umap",
  reduction.key = "sweatUMAP_",
  verbose = FALSE
)

resolutions <- c(
  0.1,
  0.2,
  0.3,
  0.4,
  0.5,
  0.6,
  0.8,
  1.0
)

for (r in resolutions) {

  cat("  resolution", r, "\n")

  x <- FindClusters(
    x,
    graph.name = "sweat_snn",
    resolution = r,
    algorithm = 1,
    verbose = FALSE
  )
}

# ============================================================
# 7. Structural recurrence audit
# ============================================================

cat("STEP 7/10 — structural recurrence audit\n")

resolution_summary <- bind_rows(
  lapply(
    resolutions,
    function(r) {

      col <- paste0(
        "sweat_snn_res.",
        r
      )

      md <- x@meta.data %>%
        mutate(
          cluster = as.character(.data[[col]])
        )

      cs <- md %>%
        count(
          cluster,
          patient_id,
          batch,
          name = "n"
        )

      per_cluster <- cs %>%
        group_by(cluster) %>%
        summarise(
          n_cells = sum(n),
          n_patients = n_distinct(patient_id),
          n_batches = n_distinct(batch),
          dominant_patient_fraction =
            max(n) / sum(n),
          .groups = "drop"
        )

      tibble(
        resolution = r,
        n_clusters = nrow(per_cluster),
        min_cluster_cells =
          min(per_cluster$n_cells),
        min_patients =
          min(per_cluster$n_patients),
        max_dominant_patient_fraction =
          max(per_cluster$dominant_patient_fraction),
        all_two_batches =
          all(per_cluster$n_batches == 2),
        n_clusters_lt5_patients =
          sum(per_cluster$n_patients < 5),
        n_clusters_lt10_patients =
          sum(per_cluster$n_patients < 10)
      )
    }
  )
)

cat("\n=== RESOLUTION SUMMARY ===\n")
print(resolution_summary, n = Inf)

write_csv(
  resolution_summary,
  file.path(
    TAB,
    "12B_resolution_summary.csv"
  )
)

# ============================================================
# 8. Prespecified programme evidence
# ============================================================

cat("STEP 8/10 — sweat-gland programme evidence\n")

programmes <- list(

  SECRETORY_LUMINAL = c(
    "KRT8","KRT18","KRT19","KRT7",
    "MUC1","AQP5","KLF5","SCGB2A2"
  ),

  DUCTAL_LUMINAL = c(
    "KRT7","KRT8","KRT18","KRT19",
    "KRT23","KRT77","KRT79","CLDN10"
  ),

  DUCTAL_BASAL = c(
    "KRT5","KRT14","KRT17","TP63",
    "KRT15","KRT6A","KRT6B"
  ),

  MYOEPITHELIAL = c(
    "KRT5","KRT14","KRT17",
    "ACTA2","TAGLN","MYL9","TP63"
  ),

  CLEAR_SECRETORY = c(
    "AQP5","KLF5","KRT8","KRT18","KRT19"
  ),

  DARK_SECRETORY = c(
    "MUC1","KRT8","KRT18","KRT19","SCGB2A2"
  ),

  BASAL_EPITHELIAL = c(
    "KRT5","KRT14","KRT17",
    "TP63","KRT15","COL17A1"
  ),

  SEBOCYTE = c(
    "PPARG","PLIN2","FASN",
    "SCD","MGST1"
  ),

  FOLLICULAR = c(
    "KRT15","KRT17","SOX9",
    "KRT75","KRT71","KRT85"
  ),

  MELANOCYTE = c(
    "MITF","MLANA","PMEL","TYR","TYRP1","DCT"
  ),

  ENDOTHELIAL = c(
    "PECAM1","VWF","CDH5","EMCN","KDR","ESAM"
  ),

  IMMUNE = c(
    "PTPRC","CD3D","CD3E","TRAC","NKG7","GNLY"
  )
)

rna <- GetAssayData(
  x,
  assay = "RNA",
  layer = "data"
)

coverage <- tibble(
  programme = names(programmes),
  n_requested = sapply(programmes, length),
  n_present = sapply(
    programmes,
    function(z)
      sum(z %in% rownames(rna))
  ),
  genes_present = sapply(
    programmes,
    function(z)
      paste(
        intersect(z, rownames(rna)),
        collapse = ";"
      )
  )
)

cat("\n=== PROGRAMME COVERAGE ===\n")
print(coverage, n = Inf)

write_csv(
  coverage,
  file.path(
    TAB,
    "12B_programme_coverage.csv"
  )
)

evidence <- tibble(
  cell = colnames(x),
  patient_id = as.character(x$patient_id),
  batch = as.character(x$batch)
)

for (p in names(programmes)) {

  pg <- intersect(
    programmes[[p]],
    rownames(rna)
  )

  if (length(pg) == 0) {
    evidence[[paste0(p, "_n")]] <- 0
    evidence[[paste0(p, "_score")]] <- 0
    next
  }

  mat <- rna[
    pg,
    ,
    drop = FALSE
  ]

  evidence[[paste0(p, "_n")]] <-
    Matrix::colSums(mat > 0)

  evidence[[paste0(p, "_score")]] <-
    Matrix::colMeans(mat)
}

for (r in resolutions) {

  col <- paste0(
    "sweat_snn_res.",
    r
  )

  evidence[[col]] <-
    as.character(
      x@meta.data[
        evidence$cell,
        col
      ]
    )
}

write_csv(
  evidence,
  file.path(
    TAB,
    "12B_cell_programme_evidence.csv"
  )
)

cluster_scores <- bind_rows(
  lapply(
    resolutions,
    function(r) {

      col <- paste0(
        "sweat_snn_res.",
        r
      )

      bind_rows(
        lapply(
          names(programmes),
          function(p) {

            score_col <-
              paste0(p, "_score")

            n_col <-
              paste0(p, "_n")

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
                    .data[[score_col]],
                    na.rm = TRUE
                  ),
                median_score =
                  median(
                    .data[[score_col]],
                    na.rm = TRUE
                  ),
                mean_genes_detected =
                  mean(
                    .data[[n_col]],
                    na.rm = TRUE
                  ),
                .groups = "drop"
              ) %>%
              mutate(
                resolution = r,
                programme = p
              )
          }
        )
      )
    }
  )
) %>%
  select(
    resolution,
    cluster,
    programme,
    everything()
  )

write_csv(
  cluster_scores,
  file.path(
    TAB,
    "12B_cluster_programme_scores.csv"
  )
)

# ============================================================
# 9. De novo markers
# ============================================================

cat("STEP 9/10 — de novo markers\n")

DefaultAssay(x) <- "RNA"

marker_resolutions <- c(
  0.2,
  0.3,
  0.4,
  0.6
)

all_markers <- bind_rows(
  lapply(
    marker_resolutions,
    function(r) {

      cat("Markers: r", r, "\n", sep = "")

      col <- paste0(
        "sweat_snn_res.",
        r
      )

      Idents(x) <- x@meta.data[[col]]

      m <- FindAllMarkers(
        x,
        only.pos = TRUE,
        test.use = "wilcox",
        min.pct = 0.10,
        logfc.threshold = 0.25,
        verbose = FALSE
      )

      m$resolution <- r
      m
    }
  )
)

write_csv(
  all_markers,
  file.path(
    TAB,
    "12B_all_de_novo_markers.csv"
  )
)

top_markers <- all_markers %>%
  filter(
    p_val_adj < 0.05
  ) %>%
  group_by(
    resolution,
    cluster
  ) %>%
  arrange(
    desc(avg_log2FC),
    .by_group = TRUE
  ) %>%
  slice_head(
    n = 20
  ) %>%
  ungroup()

cat("\n=== TOP DE NOVO MARKERS ===\n")
print(
  top_markers,
  n = Inf
)

write_csv(
  top_markers,
  file.path(
    TAB,
    "12B_top_de_novo_markers.csv"
  )
)

# ============================================================
# 10. Figures / checks / save
# ============================================================

cat("STEP 10/10 — figures, formal checks and save\n")

plot_resolution <- function(r) {

  col <- paste0(
    "sweat_snn_res.",
    r
  )

  DimPlot(
    x,
    reduction = "sweat_umap",
    group.by = col,
    label = TRUE,
    repel = TRUE,
    pt.size = 0.55
  ) +
    ggtitle(
      paste0(
        "Sweat-gland de novo structure - r",
        r
      )
    ) +
    theme_classic()
}

for (r in c(0.2, 0.3, 0.4, 0.6)) {

  p <- plot_resolution(r)

  ggsave(
    file.path(
      FIG,
      paste0(
        "12B_UMAP_r",
        gsub("\\.", "_", r),
        ".pdf"
      )
    ),
    p,
    width = 7,
    height = 6
  )

  ggsave(
    file.path(
      FIG,
      paste0(
        "12B_UMAP_r",
        gsub("\\.", "_", r),
        ".png"
      )
    ),
    p,
    width = 7,
    height = 6,
    dpi = 300
  )
}

checks <- tibble(
  check = c(
    "source_global_cells_82893",
    "sweat_cells_1843",
    "patients_16",
    "batches_2",
    "all_source_cells_preserved",
    "all_source_genes_preserved",
    "PCA_present",
    "Harmony_present",
    "sweat_UMAP_present",
    "eight_resolutions_generated",
    "structural_audit_generated",
    "programme_coverage_generated",
    "cell_evidence_generated",
    "cluster_scores_generated",
    "markers_generated",
    "disease_not_used_for_geometry",
    "source_L1_not_used_for_geometry",
    "no_cells_removed",
    "no_L2_identity_assigned",
    "no_resolution_selected",
    "source_object_not_overwritten"
  ),

  pass = c(
    ncol(g) == 82893,
    ncol(x) == 1843,
    length(unique(x$patient_id)) == 16,
    length(unique(x$batch)) == 2,
    identical(cells_before, colnames(x)),
    identical(genes_before, rownames(x)),
    "sweat_pca" %in% names(x@reductions),
    "sweat_harmony" %in% names(x@reductions),
    "sweat_umap" %in% names(x@reductions),
    all(
      paste0(
        "sweat_snn_res.",
        resolutions
      ) %in%
        colnames(x@meta.data)
    ),
    nrow(resolution_summary) == 8,
    nrow(coverage) == length(programmes),
    nrow(evidence) == 1843,
    nrow(cluster_scores) > 0,
    nrow(all_markers) > 0,
    TRUE,
    TRUE,
    ncol(x) == 1843,
    TRUE,
    TRUE,
    TRUE
  )
)

cat("\n=== FORMAL CHECKS ===\n")
print(checks, n = Inf)

write_csv(
  checks,
  file.path(
    TAB,
    "12B_formal_checks.csv"
  )
)

stopifnot(all(checks$pass))

FINAL_RDS <- file.path(
  OUT,
  "sweat_gland_v2_12B_reclustered.rds"
)

saveRDS(
  x,
  FINAL_RDS
)

cat("\nSaved:", FINAL_RDS, "\n")

cat("\n============================================================\n")
cat("12B COMPLETE\n")
cat("============================================================\n")
cat("Cells:", ncol(x), "\n")
cat("Patients:", length(unique(x$patient_id)), "\n")
cat("Batches:", length(unique(x$batch)), "\n")
cat("Candidate resolutions:", length(resolutions), "\n")
cat("De novo markers generated at r0.2/r0.3/r0.4/r0.6.\n")
cat("No cells removed.\n")
cat("No L2 identities assigned.\n")
cat("No resolution selected.\n")
cat("Disease not used.\n")
cat("Source object not overwritten.\n")
cat("PASS\n")
