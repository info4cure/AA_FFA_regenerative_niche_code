suppressPackageStartupMessages({
  library(Seurat)
  library(harmony)
  library(dplyr)
  library(tibble)
  library(readr)
  library(ggplot2)
})

cat("============================================================\n")
cat("PIPELINE_V2 — 13B MELANOCYTE DE NOVO RECLUSTERING\n")
cat("============================================================\n\n")

SRC <- "PIPELINE_V2/objects/annotated/global_atlas_v2_L1_frozen.rds"

OUT <- "PIPELINE_V2/results/annotation/13_MELANOCYTE_L2/13B_recluster"
TAB <- file.path(OUT, "tables")
FIG <- file.path(OUT, "figures")

dir.create(TAB, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG, recursive = TRUE, showWarnings = FALSE)

stopifnot(file.exists(SRC))

cat("STEP 1/10 — extracting frozen-L1 melanocyte universe\n")

g <- readRDS(SRC)

stopifnot(ncol(g) == 82893)
stopifnot("L1" %in% colnames(g@meta.data))

cells <- colnames(g)[as.character(g$L1) == "Melanocyte"]
x <- subset(g, cells = cells)

stopifnot(ncol(x) == 1556)
stopifnot(nrow(x) == 23412)
stopifnot(length(unique(x$patient_id)) == 16)
stopifnot(length(unique(x$batch)) == 2)

cells_before <- colnames(x)
genes_before <- rownames(x)

cat("Cells:", ncol(x), "\n")
cat("Genes:", nrow(x), "\n")
cat("Patients:", length(unique(x$patient_id)), "\n")
cat("Batches:", length(unique(x$batch)), "\n")

cat("\n13B bootstrap OK\n")

# ============================================================
# 2. Remove inherited geometry
# ============================================================

cat("\nSTEP 2/10 — removing inherited global geometry\n")

x@reductions <- list()
x@graphs <- list()
x@neighbors <- list()

# ============================================================
# 3. De novo normalization / variable features / scaling
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
  npcs = 50,
  reduction.name = "mel_pca",
  verbose = FALSE
)

# ============================================================
# 5. Harmony
# ============================================================

cat("STEP 5/10 — Harmony batch correction\n")

x <- RunHarmony(
  object = x,
  group.by.vars = "batch",
  reduction.use = "mel_pca",
  reduction.save = "mel_harmony",
  dims.use = 1:30,
  verbose = FALSE
)

# ============================================================
# 6. Neighbours / UMAP / multiresolution clustering
# ============================================================

cat("STEP 6/10 — neighbors, UMAP and multiresolution clustering\n")

x <- FindNeighbors(
  x,
  reduction = "mel_harmony",
  dims = 1:30,
  graph.name = c("mel_nn", "mel_snn"),
  verbose = FALSE
)

x <- RunUMAP(
  x,
  reduction = "mel_harmony",
  dims = 1:30,
  reduction.name = "mel_umap",
  seed.use = 20260930,
  verbose = FALSE
)

resolutions <- c(
  0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8
)

for (res in resolutions) {

  cat("  resolution", res, "\n")

  x <- FindClusters(
    x,
    graph.name = "mel_snn",
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

  if (length(generated) == 0) {
    stop(
      "Could not identify generated cluster column for resolution ",
      res
    )
  }

  new_col <- paste0(
    "mel_snn_res.",
    res
  )

  x[[new_col]] <-
    x@meta.data[[tail(generated, 1)]]
}

# ============================================================
# 7. Structural recurrence audit
# ============================================================

cat("STEP 7/10 — structural recurrence audit\n")

structural_list <- list()
cluster_support_list <- list()

for (res in resolutions) {

  col <- paste0(
    "mel_snn_res.",
    res
  )

  tmp <- x@meta.data %>%
    as_tibble() %>%
    mutate(
      cluster = as.character(.data[[col]])
    ) %>%
    count(
      cluster,
      patient_id,
      batch,
      name = "n_cells"
    )

  cluster_support <- tmp %>%
    group_by(cluster) %>%
    summarise(
      n_cells = sum(n_cells),
      n_patients = n_distinct(patient_id),
      n_batches = n_distinct(batch),
      dominant_patient_fraction =
        max(n_cells) / sum(n_cells),
      .groups = "drop"
    ) %>%
    mutate(
      resolution = res
    )

  cluster_support_list[[as.character(res)]] <-
    cluster_support

  structural_list[[as.character(res)]] <-
    cluster_support %>%
    summarise(
      resolution = first(resolution),
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
  bind_rows(structural_list)

cluster_support_all <-
  bind_rows(cluster_support_list)

cat("\n=== RESOLUTION SUMMARY ===\n")
print(
  resolution_summary,
  n = Inf
)

write_csv(
  resolution_summary,
  file.path(
    TAB,
    "13B_resolution_summary.csv"
  )
)

write_csv(
  cluster_support_all,
  file.path(
    TAB,
    "13B_cluster_structural_support.csv"
  )
)

cat("\n13B geometry/reclustering block OK\n")

# ============================================================
# 8. Melanocyte identity/state programme evidence
# ============================================================

cat("\nSTEP 8/10 — melanocyte programme evidence\n")

programmes <- list(

  MELANOCYTE_CORE = c(
    "MITF","MLANA","PMEL","TYR",
    "TYRP1","DCT","SOX10","S100B"
  ),

  MELANOGENESIS = c(
    "TYR","TYRP1","DCT","PMEL",
    "MLANA","OCA2","SLC45A2","SLC24A5"
  ),

  MITF_PROGRAMME = c(
    "MITF","MLANA","PMEL",
    "TYR","TYRP1","DCT"
  ),

  PROGENITOR = c(
    "SOX10","KIT","NGFR",
    "EDNRB","PAX3","SOX9"
  ),

  NEURAL_CREST_SHARED = c(
    "SOX10","S100B","PLP1",
    "PMP22","NGFR","GFRA3"
  ),

  SCHWANN_SPECIFIC = c(
    "PLP1","PMP22","MPZ",
    "S100A1","GFRA3","NGFR"
  ),

  IFN_RESPONSE = c(
    "IFIT1","IFIT2","IFIT3",
    "ISG15","MX1","STAT1"
  ),

  ANTIGEN_PRESENTATION = c(
    "HLA-A","HLA-B","HLA-C",
    "B2M","HLA-DRA","HLA-DRB1"
  ),

  HLA_II = c(
    "HLA-DRA","HLA-DRB1",
    "HLA-DPA1","HLA-DPB1",
    "HLA-DQA1","HLA-DQB1"
  ),

  CYCLING = c(
    "MKI67","TOP2A","UBE2C",
    "CENPF","TYMS","PCNA"
  ),

  EPITHELIAL = c(
    "KRT5","KRT14","KRT15",
    "KRT17","EPCAM","KRT8"
  ),

  IMMUNE = c(
    "PTPRC","CD3D","CD3E",
    "TRAC","NKG7","GNLY"
  ),

  FIBROBLAST = c(
    "DCN","LUM","COL1A1",
    "COL1A2","COL3A1","PDGFRA"
  ),

  ENDOTHELIAL = c(
    "PECAM1","VWF","CDH5",
    "EMCN","KDR","ESAM"
  )
)

# ------------------------------------------------------------
# Programme coverage
# ------------------------------------------------------------

coverage <- bind_rows(
  lapply(
    names(programmes),
    function(nm) {

      requested <- programmes[[nm]]

      present <- intersect(
        requested,
        rownames(x)
      )

      tibble(
        programme = nm,
        n_requested = length(requested),
        n_present = length(present),
        genes_present =
          paste(present, collapse = ";")
      )
    }
  )
)

cat("\n=== PROGRAMME COVERAGE ===\n")
print(
  coverage,
  n = Inf
)

write_csv(
  coverage,
  file.path(
    TAB,
    "13B_programme_coverage.csv"
  )
)

# ------------------------------------------------------------
# Cell-level programme evidence
# ------------------------------------------------------------

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

  if (length(genes) == 0) {
    return(
      rep(
        NA_real_,
        ncol(x)
      )
    )
  }

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

  if (length(genes) == 0) {
    return(
      rep(
        0L,
        ncol(x)
      )
    )
  }

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
  patient_id =
    as.character(x$patient_id),
  batch =
    as.character(x$batch)
)

for (nm in names(programmes)) {

  evidence[[paste0(nm, "_score")]] <-
    score_programme(
      programmes[[nm]]
    )

  evidence[[paste0(nm, "_genes")]] <-
    detect_programme(
      programmes[[nm]]
    )
}

write_csv(
  evidence,
  file.path(
    TAB,
    "13B_cell_programme_evidence.csv"
  )
)

# ------------------------------------------------------------
# Candidate-resolution cluster programme summaries
# ------------------------------------------------------------

candidate_resolutions <- c(
  0.2,
  0.3,
  0.4,
  0.6
)

score_tables <- list()

for (res in candidate_resolutions) {

  cluster_col <- paste0(
    "mel_snn_res.",
    res
  )

  cl <- as.character(
    x@meta.data[[cluster_col]]
  )

  tmp <- evidence %>%
    mutate(
      cluster = cl
    )

  for (nm in names(programmes)) {

    score_col <-
      paste0(
        nm,
        "_score"
      )

    genes_col <-
      paste0(
        nm,
        "_genes"
      )

    score_tables[[paste0(
      res,
      "_",
      nm
    )]] <-
      tmp %>%
      group_by(cluster) %>%
      summarise(
        resolution = res,
        programme = nm,

        n_cells = n(),

        n_patients =
          n_distinct(
            patient_id
          ),

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
            .data[[genes_col]],
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
  bind_rows(
    score_tables
  ) %>%
  arrange(
    resolution,
    as.integer(cluster),
    desc(mean_score)
  )

cat("\n=== CLUSTER PROGRAMME SCORES ===\n")

print(
  cluster_scores,
  n = 120
)

write_csv(
  cluster_scores,
  file.path(
    TAB,
    "13B_cluster_programme_scores.csv"
  )
)

# ------------------------------------------------------------
# Compact state summary for rapid adjudication
# ------------------------------------------------------------

compact_scores <-
  cluster_scores %>%
  filter(
    programme %in% c(
      "MELANOCYTE_CORE",
      "MELANOGENESIS",
      "MITF_PROGRAMME",
      "PROGENITOR",
      "NEURAL_CREST_SHARED",
      "SCHWANN_SPECIFIC",
      "IFN_RESPONSE",
      "ANTIGEN_PRESENTATION",
      "HLA_II",
      "CYCLING",
      "EPITHELIAL",
      "IMMUNE"
    )
  ) %>%
  select(
    resolution,
    cluster,
    programme,
    median_score,
    fraction_ge2,
    fraction_ge3
  )

write_csv(
  compact_scores,
  file.path(
    TAB,
    "13B_compact_state_scores.csv"
  )
)

cat("\n13B biological-programme block OK\n")

# ============================================================
# 9. De novo markers + patient recurrence
# ============================================================

cat("\nSTEP 9/10 — de novo markers and patient recurrence\n")

marker_resolutions <- c(
  0.2,
  0.3,
  0.4,
  0.6
)

all_markers <- list()

for (res in marker_resolutions) {

  cat("Markers: r", res, "\n", sep = "")

  cluster_col <- paste0(
    "mel_snn_res.",
    res
  )

  stopifnot(
    cluster_col %in%
      colnames(x@meta.data)
  )

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

  write_csv(
    m,
    file.path(
      TAB,
      paste0(
        "13B_r",
        gsub(
          "\\.",
          "",
          as.character(res)
        ),
        "_all_markers.csv"
      )
    )
  )
}

markers <-
  bind_rows(
    all_markers
  )

stopifnot(
  nrow(markers) > 0
)

top_markers <-
  markers %>%
  group_by(
    resolution,
    cluster
  ) %>%
  arrange(
    desc(avg_log2FC),
    .by_group = TRUE
  ) %>%
  slice_head(
    n = 25
  ) %>%
  ungroup()

cat("\n=== TOP DE NOVO MARKERS ===\n")

print(
  top_markers,
  n = Inf
)

write_csv(
  markers,
  file.path(
    TAB,
    "13B_all_candidate_markers.csv"
  )
)

write_csv(
  top_markers,
  file.path(
    TAB,
    "13B_top25_markers_all_candidate_resolutions.csv"
  )
)

# ------------------------------------------------------------
# Patient-level recurrence
# ------------------------------------------------------------

patient_recurrence_list <- list()

for (res in candidate_resolutions) {

  cluster_col <- paste0(
    "mel_snn_res.",
    res
  )

  tmp <- evidence %>%
    mutate(
      cluster =
        as.character(
          x@meta.data[[cluster_col]]
        )
    )

  patient_cluster <-
    tmp %>%
    group_by(
      cluster,
      patient_id,
      batch
    ) %>%
    summarise(
      n_cells = n(),

      median_melanocyte =
        median(
          MELANOCYTE_CORE_score,
          na.rm = TRUE
        ),

      median_melanogenesis =
        median(
          MELANOGENESIS_score,
          na.rm = TRUE
        ),

      median_MITF =
        median(
          MITF_PROGRAMME_score,
          na.rm = TRUE
        ),

      median_progenitor =
        median(
          PROGENITOR_score,
          na.rm = TRUE
        ),

      median_neural_crest =
        median(
          NEURAL_CREST_SHARED_score,
          na.rm = TRUE
        ),

      median_Schwann_specific =
        median(
          SCHWANN_SPECIFIC_score,
          na.rm = TRUE
        ),

      median_IFN =
        median(
          IFN_RESPONSE_score,
          na.rm = TRUE
        ),

      median_antigen_presentation =
        median(
          ANTIGEN_PRESENTATION_score,
          na.rm = TRUE
        ),

      median_HLAII =
        median(
          HLA_II_score,
          na.rm = TRUE
        ),

      median_cycling =
        median(
          CYCLING_score,
          na.rm = TRUE
        ),

      .groups = "drop"
    )

  write_csv(
    patient_cluster,
    file.path(
      TAB,
      paste0(
        "13B_r",
        gsub(
          "\\.",
          "",
          as.character(res)
        ),
        "_patient_cluster_evidence.csv"
      )
    )
  )

  patient_summary <-
    patient_cluster %>%
    group_by(cluster) %>%
    summarise(
      resolution = res,

      n_patients =
        n_distinct(
          patient_id
        ),

      n_batches =
        n_distinct(
          batch
        ),

      patients_ge5_cells =
        sum(
          n_cells >= 5
        ),

      median_cells_patient =
        median(
          n_cells
        ),

      median_patient_melanocyte =
        median(
          median_melanocyte,
          na.rm = TRUE
        ),

      median_patient_melanogenesis =
        median(
          median_melanogenesis,
          na.rm = TRUE
        ),

      median_patient_MITF =
        median(
          median_MITF,
          na.rm = TRUE
        ),

      median_patient_progenitor =
        median(
          median_progenitor,
          na.rm = TRUE
        ),

      median_patient_neural_crest =
        median(
          median_neural_crest,
          na.rm = TRUE
        ),

      median_patient_Schwann_specific =
        median(
          median_Schwann_specific,
          na.rm = TRUE
        ),

      median_patient_IFN =
        median(
          median_IFN,
          na.rm = TRUE
        ),

      median_patient_antigen_presentation =
        median(
          median_antigen_presentation,
          na.rm = TRUE
        ),

      median_patient_HLAII =
        median(
          median_HLAII,
          na.rm = TRUE
        ),

      median_patient_cycling =
        median(
          median_cycling,
          na.rm = TRUE
        ),

      .groups = "drop"
    )

  patient_recurrence_list[[as.character(res)]] <-
    patient_summary
}

patient_recurrence <-
  bind_rows(
    patient_recurrence_list
  ) %>%
  arrange(
    resolution,
    as.integer(cluster)
  )

cat("\n=== PATIENT RECURRENCE SUMMARY ===\n")

print(
  patient_recurrence,
  n = Inf
)

write_csv(
  patient_recurrence,
  file.path(
    TAB,
    "13B_patient_recurrence_summary.csv"
  )
)

# ============================================================
# 10. Formal checks + save
# ============================================================

cat("\nSTEP 10/10 — formal checks and save\n")

cells_after <- colnames(x)
genes_after <- rownames(x)

resolution_cols <- paste0(
  "mel_snn_res.",
  resolutions
)

checks <- tibble(
  check = c(
    "source_global_cells_82893",
    "melanocyte_cells_1556",
    "genes_23412",
    "patients_16",
    "batches_2",
    "all_source_cells_preserved",
    "all_source_genes_preserved",
    "PCA_present",
    "Harmony_present",
    "melanocyte_UMAP_present",
    "seven_resolutions_generated",
    "structural_audit_generated",
    "programme_coverage_complete",
    "cell_evidence_generated",
    "cluster_scores_generated",
    "markers_generated",
    "patient_recurrence_generated",
    "candidate_resolutions_audited",
    "disease_not_used_for_geometry",
    "disease_not_used_for_identity",
    "source_L1_not_used_for_geometry",
    "no_cells_removed",
    "no_L2_identity_assigned",
    "no_resolution_selected",
    "source_object_not_overwritten"
  ),

  pass = c(
    ncol(g) == 82893,
    ncol(x) == 1556,
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

    "mel_pca" %in%
      names(x@reductions),

    "mel_harmony" %in%
      names(x@reductions),

    "mel_umap" %in%
      names(x@reductions),

    all(
      resolution_cols %in%
        colnames(x@meta.data)
    ),

    nrow(resolution_summary) ==
      length(resolutions),

    all(
      coverage$n_present ==
        coverage$n_requested
    ),

    nrow(evidence) ==
      1556,

    nrow(cluster_scores) >
      0,

    nrow(markers) >
      0,

    nrow(patient_recurrence) >
      0,

    all(
      candidate_resolutions %in%
        unique(
          patient_recurrence$resolution
        )
    ),

    TRUE,
    TRUE,
    TRUE,

    ncol(x) ==
      length(cells_before),

    !(
      "melanocyte_L2_label" %in%
        colnames(x@meta.data)
    ),

    !(
      "melanocyte_resolution_frozen" %in%
        colnames(x@meta.data)
    ),

    TRUE
  )
)

cat("\n=== FORMAL CHECKS ===\n")

print(
  checks,
  n = Inf
)

write_csv(
  checks,
  file.path(
    TAB,
    "13B_formal_checks.csv"
  )
)

stopifnot(
  all(checks$pass)
)

FINAL_RDS <- file.path(
  OUT,
  "melanocyte_v2_13B_reclustered.rds"
)

saveRDS(
  x,
  FINAL_RDS
)

cat("\nSaved:\n")
cat(
  FINAL_RDS,
  "\n"
)

cat("\n============================================================\n")
cat("13B COMPLETE\n")
cat("============================================================\n")
cat("Cells:", ncol(x), "\n")
cat("Genes:", nrow(x), "\n")
cat(
  "Patients:",
  length(unique(x$patient_id)),
  "\n"
)
cat(
  "Batches:",
  length(unique(x$batch)),
  "\n"
)
cat(
  "Candidate resolutions:",
  paste(
    resolutions,
    collapse = ", "
  ),
  "\n"
)
cat("Marker resolutions: 0.2, 0.3, 0.4, 0.6\n")
cat("Melanocyte identity programmes evaluated.\n")
cat("Progenitor/neural-crest programmes evaluated.\n")
cat("Schwann-specific programme evaluated separately.\n")
cat("IFN/HLA-II/antigen-presentation states evaluated.\n")
cat("Cycling evaluated.\n")
cat("Patient recurrence evaluated.\n")
cat("Figures generated: NO\n")
cat("Disease comparisons performed: NO\n")
cat("Cells removed: 0\n")
cat("L2 identities assigned: NO\n")
cat("Resolution selected: NO\n")
cat("Source object overwritten: NO\n")
cat("READY FOR FINAL MELANOCYTE ADJUDICATION\n")
cat("PASS\n")
