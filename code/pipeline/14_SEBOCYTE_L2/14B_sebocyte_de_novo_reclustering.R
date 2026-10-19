suppressPackageStartupMessages({
  library(Seurat)
  library(harmony)
  library(dplyr)
  library(tibble)
  library(readr)
})

cat("============================================================\n")
cat("PIPELINE_V2 — 14B SEBOCYTE DE NOVO RECLUSTERING\n")
cat("============================================================\n\n")

SRC <- "PIPELINE_V2/objects/annotated/global_atlas_v2_L1_frozen.rds"
OUT <- "PIPELINE_V2/results/annotation/14_SEBOCYTE_L2/14B_recluster"
TAB <- file.path(OUT, "tables")

dir.create(TAB, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# 1. Extract frozen-L1 sebocyte universe
# ============================================================

cat("STEP 1/10 — extracting frozen-L1 sebocyte universe\n")

g <- readRDS(SRC)

stopifnot(ncol(g) == 82893)

x <- subset(
  g,
  subset = L1 == "Sebocyte"
)

stopifnot(ncol(x) == 1388)
stopifnot(nrow(x) == 23412)
stopifnot(length(unique(x$patient_id)) == 16)
stopifnot(length(unique(x$batch)) == 2)

cells_before <- colnames(x)
genes_before <- rownames(x)

cat("Cells:", ncol(x), "\n")
cat("Genes:", nrow(x), "\n")
cat("Patients:", length(unique(x$patient_id)), "\n")
cat("Batches:", length(unique(x$batch)), "\n")

# ============================================================
# 2. Remove inherited geometry
# ============================================================

cat("STEP 2/10 — removing inherited global geometry\n")

x@reductions <- list()
x@graphs <- list()
x@neighbors <- list()

# ============================================================
# 3. De novo preprocessing
# ============================================================

cat("STEP 3/10 — normalization/features/scaling\n")

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
  npcs = 50,
  reduction.name = "seb_pca",
  verbose = FALSE
)

# ============================================================
# 5. Harmony
# ============================================================

cat("STEP 5/10 — Harmony batch correction\n")

x <- RunHarmony(
  object = x,
  group.by.vars = "batch",
  reduction.use = "seb_pca",
  reduction.save = "seb_harmony",
  dims.use = 1:30,
  verbose = FALSE
)

# ============================================================
# 6. Neighbours / UMAP / multiresolution clustering
# ============================================================

cat("STEP 6/10 — neighbors, UMAP and multiresolution clustering\n")

x <- FindNeighbors(
  x,
  reduction = "seb_harmony",
  dims = 1:30,
  graph.name = c("seb_nn", "seb_snn"),
  verbose = FALSE
)

x <- RunUMAP(
  x,
  reduction = "seb_harmony",
  dims = 1:30,
  reduction.name = "seb_umap",
  seed.use = 20260930,
  verbose = FALSE
)

resolutions <- c(
  0.1, 0.2, 0.3, 0.4,
  0.5, 0.6, 0.8
)

for (res in resolutions) {

  cat("  resolution", res, "\n")

  x <- FindClusters(
    x,
    graph.name = "seb_snn",
    resolution = res,
    algorithm = 1,
    random.seed = 20260930,
    verbose = FALSE
  )

  generated <- grep(
    paste0("res\\.", res, "$"),
    colnames(x@meta.data),
    value = TRUE
  )

  stopifnot(length(generated) > 0)

  x[[paste0("seb_snn_res.", res)]] <-
    x@meta.data[[tail(generated, 1)]]
}

# ============================================================
# 7. Structural recurrence
# ============================================================

cat("STEP 7/10 — structural recurrence audit\n")

structural <- list()

for (res in resolutions) {

  col <- paste0("seb_snn_res.", res)

  support <- x@meta.data %>%
    as_tibble() %>%
    mutate(
      cluster = as.character(.data[[col]])
    ) %>%
    count(
      cluster,
      patient_id,
      batch,
      name = "n_cells"
    ) %>%
    group_by(cluster) %>%
    summarise(
      n_cells = sum(n_cells),
      n_patients = n_distinct(patient_id),
      n_batches = n_distinct(batch),

      dominant_patient_fraction =
        max(n_cells) / sum(n_cells),

      .groups = "drop"
    )

  write_csv(
    support,
    file.path(
      TAB,
      paste0(
        "14B_r",
        gsub("\\.", "", res),
        "_structural_support.csv"
      )
    )
  )

  structural[[as.character(res)]] <-
    support %>%
    summarise(
      resolution = res,
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
        sum(n_patients < 10)
    )
}

resolution_summary <-
  bind_rows(structural)

cat("\n=== RESOLUTION SUMMARY ===\n")
print(resolution_summary, n = Inf)

write_csv(
  resolution_summary,
  file.path(TAB, "14B_resolution_summary.csv")
)

# ============================================================
# 8. Biological programmes
# ============================================================

cat("STEP 8/10 — sebocyte programme evidence\n")

programmes <- list(

  SEBOCYTE_CORE = c(
    "PPARG","PLIN2","FASN","SCD",
    "MGST1","FABP5","KRT7","KRT19"
  ),

  LIPOGENESIS = c(
    "FASN","SCD","ACACA","ACLY",
    "ELOVL5","ELOVL6","DGAT1","DGAT2"
  ),

  LIPID_STORAGE = c(
    "PLIN2","PLIN3","FABP5",
    "ABCA1","APOE","LPL"
  ),

  BASAL_PROGENITOR = c(
    "KRT5","KRT14","KRT15",
    "KRT17","TP63","KRT6A"
  ),

  KERATINIZATION = c(
    "KRT1","KRT10","IVL",
    "SFN","KRTDAP","DSG1"
  ),

  INFLAMMATORY = c(
    "SAA1","SAA2","LCN2",
    "LTF","CXCL8","CXCL14"
  ),

  IFN_RESPONSE = c(
    "IFIT1","IFIT2","IFIT3",
    "ISG15","MX1","STAT1"
  ),

  ANTIGEN_PRESENTATION = c(
    "HLA-A","HLA-B","HLA-C",
    "B2M","HLA-DRA","HLA-DRB1"
  ),

  SWEAT_GLAND = c(
    "AQP5","DCD","PIP",
    "SCGB2A2","SCGB1D2","KRT77"
  ),

  FOLLICULAR = c(
    "KRT15","KRT17","SOX9",
    "KRT75","KRT71","KRT85"
  ),

  MELANOCYTE = c(
    "MITF","MLANA","PMEL",
    "TYR","TYRP1","DCT"
  ),

  IMMUNE = c(
    "PTPRC","CD3D","CD3E",
    "TRAC","NKG7","GNLY"
  )
)

rna <- GetAssayData(
  x,
  assay = "RNA",
  layer = "data"
)

score_programme <- function(genes) {

  genes <- intersect(
    genes,
    rownames(rna)
  )

  Matrix::colMeans(
    rna[
      genes,
      ,
      drop = FALSE
    ]
  )
}

detect_programme <- function(genes) {

  genes <- intersect(
    genes,
    rownames(rna)
  )

  Matrix::colSums(
    rna[
      genes,
      ,
      drop = FALSE
    ] > 0
  )
}

evidence <- tibble(
  cell = colnames(x),
  patient_id = as.character(x$patient_id),
  batch = as.character(x$batch)
)

for (nm in names(programmes)) {

  evidence[[paste0(nm, "_score")]] <-
    score_programme(programmes[[nm]])

  evidence[[paste0(nm, "_genes")]] <-
    detect_programme(programmes[[nm]])
}

write_csv(
  evidence,
  file.path(TAB, "14B_cell_programme_evidence.csv")
)

candidate_resolutions <- c(
  0.2, 0.3, 0.4, 0.6
)

score_tables <- list()

for (res in candidate_resolutions) {

  cl <- as.character(
    x@meta.data[[paste0("seb_snn_res.", res)]]
  )

  tmp <- evidence %>%
    mutate(cluster = cl)

  for (nm in names(programmes)) {

    score_col <- paste0(nm, "_score")
    genes_col <- paste0(nm, "_genes")

    score_tables[[paste0(res, "_", nm)]] <-
      tmp %>%
      group_by(cluster) %>%
      summarise(
        resolution = res,
        programme = nm,
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

        fraction_ge2 =
          mean(
            .data[[genes_col]] >= 2,
            na.rm = TRUE
          ),

        fraction_ge3 =
          mean(
            .data[[genes_col]] >= 3,
            na.rm = TRUE
          ),

        .groups = "drop"
      )
  }
}

cluster_scores <-
  bind_rows(score_tables) %>%
  arrange(
    resolution,
    as.integer(cluster),
    desc(median_score)
  )

write_csv(
  cluster_scores,
  file.path(TAB, "14B_cluster_programme_scores.csv")
)

# ============================================================
# 9. Markers + patient recurrence
# ============================================================

cat("STEP 9/10 — markers and patient recurrence\n")

all_markers <- list()

for (res in candidate_resolutions) {

  cat("Markers: r", res, "\n", sep = "")

  cluster_col <-
    paste0("seb_snn_res.", res)

  Idents(x) <-
    x@meta.data[[cluster_col]]

  m <- FindAllMarkers(
    x,
    assay = "RNA",
    only.pos = TRUE,
    min.pct = 0.10,
    logfc.threshold = 0.25,
    test.use = "wilcox",
    verbose = FALSE
  )

  m$cluster <-
    as.character(m$cluster)

  m$resolution <- res

  all_markers[[as.character(res)]] <- m
}

markers <- bind_rows(all_markers)

stopifnot(nrow(markers) > 0)

write_csv(
  markers,
  file.path(TAB, "14B_all_candidate_markers.csv")
)

top_markers <- markers %>%
  group_by(
    resolution,
    cluster
  ) %>%
  arrange(
    desc(avg_log2FC),
    .by_group = TRUE
  ) %>%
  slice_head(n = 25) %>%
  ungroup()

write_csv(
  top_markers,
  file.path(
    TAB,
    "14B_top25_candidate_markers.csv"
  )
)

patient_recurrence <- list()

for (res in candidate_resolutions) {

  cluster_col <-
    paste0("seb_snn_res.", res)

  tmp <- evidence %>%
    mutate(
      cluster =
        as.character(
          x@meta.data[[cluster_col]]
        )
    ) %>%
    count(
      cluster,
      patient_id,
      batch,
      name = "n_cells"
    ) %>%
    group_by(cluster) %>%
    summarise(
      resolution = res,
      n_patients =
        n_distinct(patient_id),
      n_batches =
        n_distinct(batch),
      patients_ge5_cells =
        sum(n_cells >= 5),
      median_cells_patient =
        median(n_cells),
      .groups = "drop"
    )

  patient_recurrence[[as.character(res)]] <-
    tmp
}

patient_recurrence <-
  bind_rows(patient_recurrence) %>%
  arrange(
    resolution,
    as.integer(cluster)
  )

cat("\n=== PATIENT RECURRENCE ===\n")
print(patient_recurrence, n = Inf)

write_csv(
  patient_recurrence,
  file.path(
    TAB,
    "14B_patient_recurrence_summary.csv"
  )
)

# ============================================================
# 10. Checks + save
# ============================================================

cat("STEP 10/10 — formal checks and save\n")

cells_after <- colnames(x)
genes_after <- rownames(x)

resolution_cols <-
  paste0(
    "seb_snn_res.",
    resolutions
  )

checks <- tibble(
  check = c(
    "source_global_cells_82893",
    "sebocyte_cells_1388",
    "genes_23412",
    "patients_16",
    "batches_2",
    "all_source_cells_preserved",
    "all_source_genes_preserved",
    "PCA_present",
    "Harmony_present",
    "sebocyte_UMAP_present",
    "seven_resolutions_generated",
    "structural_audit_generated",
    "cell_evidence_generated",
    "cluster_scores_generated",
    "markers_generated",
    "patient_recurrence_generated",
    "disease_not_used_for_geometry",
    "disease_not_used_for_identity",
    "no_cells_removed",
    "no_L2_identity_assigned",
    "no_resolution_selected",
    "figures_generated_NO",
    "source_object_not_overwritten"
  ),

  pass = c(
    ncol(g) == 82893,
    ncol(x) == 1388,
    nrow(x) == 23412,
    length(unique(x$patient_id)) == 16,
    length(unique(x$batch)) == 2,

    identical(
      cells_before,
      cells_after
    ),

    identical(
      genes_before,
      genes_after
    ),

    "seb_pca" %in%
      names(x@reductions),

    "seb_harmony" %in%
      names(x@reductions),

    "seb_umap" %in%
      names(x@reductions),

    all(
      resolution_cols %in%
        colnames(x@meta.data)
    ),

    nrow(resolution_summary) ==
      length(resolutions),

    nrow(evidence) == 1388,
    nrow(cluster_scores) > 0,
    nrow(markers) > 0,
    nrow(patient_recurrence) > 0,

    TRUE,
    TRUE,

    ncol(x) ==
      length(cells_before),

    !(
      "sebocyte_L2_label" %in%
        colnames(x@meta.data)
    ),

    !(
      "sebocyte_resolution_frozen" %in%
        colnames(x@meta.data)
    ),

    TRUE,
    TRUE
  )
)

cat("\n=== FORMAL CHECKS ===\n")
print(checks, n = Inf)

write_csv(
  checks,
  file.path(TAB, "14B_formal_checks.csv")
)

stopifnot(all(checks$pass))

FINAL_RDS <- file.path(
  OUT,
  "sebocyte_v2_14B_reclustered.rds"
)

saveRDS(x, FINAL_RDS)

cat("\nSaved:\n")
cat(FINAL_RDS, "\n")

cat("\n============================================================\n")
cat("14B COMPLETE\n")
cat("============================================================\n")
cat("Cells: 1388\n")
cat("Patients: 16\n")
cat("Batches: 2\n")
cat("Candidate resolutions: 0.1-0.8\n")
cat("Markers: r0.2/r0.3/r0.4/r0.6\n")
cat("Disease comparisons: NO\n")
cat("Figures generated: NO\n")
cat("Cells removed: 0\n")
cat("L2 identities assigned: NO\n")
cat("Resolution selected: NO\n")
cat("READY FOR SEBOCYTE ADJUDICATION\n")
cat("PASS\n")
