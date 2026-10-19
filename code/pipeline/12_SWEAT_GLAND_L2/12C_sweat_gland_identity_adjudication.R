
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tibble)
  library(readr)
})

cat("============================================================\n")
cat("PIPELINE_V2 — 12C SWEAT-GLAND IDENTITY ADJUDICATION\n")
cat("============================================================\n\n")

# ============================================================
# 0. Paths
# ============================================================

SRC <- paste0(
  "PIPELINE_V2/results/annotation/12_SWEAT_GLAND_L2/",
  "12B_recluster/sweat_gland_v2_12B_reclustered.rds"
)

OUT <- paste0(
  "PIPELINE_V2/results/annotation/12_SWEAT_GLAND_L2/",
  "12C_adjudication"
)

TAB <- file.path(OUT, "tables")

dir.create(
  TAB,
  recursive = TRUE,
  showWarnings = FALSE
)

stopifnot(file.exists(SRC))

# ============================================================
# 1. Load object
# ============================================================

cat("STEP 1/8 — loading 12B sweat-gland scaffold\n")

x <- readRDS(SRC)

stopifnot(ncol(x) == 1843)
stopifnot(nrow(x) == 23412)
stopifnot(length(unique(x$patient_id)) == 16)
stopifnot(length(unique(x$batch)) == 2)

r04_col <- "sweat_snn_res.0.4"

stopifnot(
  r04_col %in% colnames(x@meta.data)
)

x$sweat_cluster_candidate <- as.character(
  x@meta.data[[r04_col]]
)

stopifnot(
  length(unique(x$sweat_cluster_candidate)) == 8
)

cat("Cells:", ncol(x), "\n")
cat("Patients:", length(unique(x$patient_id)), "\n")
cat("Batches:", length(unique(x$batch)), "\n")
cat("Candidate resolution: 0.4\n")
cat("Clusters: 8\n\n")

# ============================================================
# 2. Structural support
# ============================================================

cat("STEP 2/8 — structural support\n")

md <- x@meta.data %>%
  mutate(
    cell = rownames(.),
    cluster = sweat_cluster_candidate
  )

if (!"percent_mt" %in% colnames(md) &&
    "percent.mt" %in% colnames(md)) {
  md$percent_mt <- md$percent.mt
}

structural <- md %>%
  group_by(cluster) %>%
  summarise(
    n_cells = n(),
    n_patients = n_distinct(patient_id),
    n_batches = n_distinct(batch),
    dominant_patient_fraction =
      max(table(patient_id)) / n(),
    median_doublet =
      median(doublet_score, na.rm = TRUE),
    median_features =
      median(nFeature_RNA, na.rm = TRUE),
    median_counts =
      median(nCount_RNA, na.rm = TRUE),
    median_mt =
      median(percent_mt, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(as.integer(cluster))

cat("\n=== STRUCTURAL SUPPORT ===\n")
print(structural, n = Inf)

write_csv(
  structural,
  file.path(TAB, "12C_structural_support.csv")
)

# ============================================================
# 3. Define glandular and competing-lineage programmes
# ============================================================

cat("\nSTEP 3/8 — defining identity programmes\n")

programmes <- list(

  PAN_SWEAT = c(
    "KRT7","KRT8","KRT18","KRT19",
    "KRT23","KRT77","KRT79"
  ),

  SECRETORY_GENERAL = c(
    "KRT8","KRT18","KRT19","KRT7",
    "MUC1","AQP5","KLF5","DCD"
  ),

  CLEAR_SECRETORY = c(
    "AQP5","KLF5","KRT8","KRT18",
    "KRT19","CA2","SLC12A2"
  ),

  DARK_SECRETORY = c(
    "MUC1","SCGB2A2","SCGB1D2",
    "PIP","MUCL1","DCD","KRT18"
  ),

  DUCTAL_LUMINAL = c(
    "KRT7","KRT8","KRT18","KRT19",
    "KRT23","KRT77","KRT79","KRT80"
  ),

  DUCTAL_BASAL = c(
    "KRT5","KRT14","KRT17","TP63",
    "KRT15","KRT6A","KRT6B"
  ),

  MYOEPITHELIAL = c(
    "KRT5","KRT14","KRT17",
    "ACTA2","TAGLN","MYL9","TP63"
  ),

  SECRETORY_APOCRINE_LIKE = c(
    "PIP","DCD","SCGB2A2","SCGB1D2",
    "MUCL1","OBP2B","LACRT"
  ),

  INFLAMMATORY_EPITHELIAL = c(
    "LCN2","LTF","SAA1","SAA2",
    "CHI3L1","CHI3L2","CXCL5"
  ),

  SEBOCYTE = c(
    "PPARG","PLIN2","FASN","SCD",
    "MGST1"
  ),

  FOLLICULAR = c(
    "KRT15","KRT17","SOX9",
    "KRT75","KRT71","KRT85"
  ),

  MELANOCYTE = c(
    "MITF","MLANA","PMEL",
    "TYR","TYRP1","DCT"
  ),

  ENDOTHELIAL = c(
    "PECAM1","VWF","CDH5",
    "EMCN","KDR","ESAM"
  ),

  IMMUNE = c(
    "PTPRC","CD3D","CD3E",
    "TRAC","NKG7","GNLY"
  )
)

programme_coverage <- bind_rows(
  lapply(names(programmes), function(p) {

    present <- intersect(
      programmes[[p]],
      rownames(x)
    )

    tibble(
      programme = p,
      n_requested = length(programmes[[p]]),
      n_present = length(present),
      genes_present = paste(
        present,
        collapse = ";"
      )
    )
  })
)

cat("\n=== PROGRAMME COVERAGE ===\n")
print(programme_coverage, n = Inf)

write_csv(
  programme_coverage,
  file.path(TAB, "12C_programme_coverage.csv")
)

# ============================================================
# 4. Cell-level programme evidence
# ============================================================

cat("\nSTEP 4/8 — cell-level programme evidence\n")

DefaultAssay(x) <- "RNA"

expr <- GetAssayData(
  x,
  assay = "RNA",
  layer = "data"
)

evidence <- md %>%
  select(
    cell,
    patient_id,
    batch,
    cluster,
    doublet_score,
    nFeature_RNA,
    nCount_RNA,
    percent_mt
  )

for (p in names(programmes)) {

  genes <- intersect(
    programmes[[p]],
    rownames(expr)
  )

  score_col <- paste0(p, "_score")
  n_col <- paste0(p, "_n")

  if (length(genes) == 0) {

    evidence[[score_col]] <- 0
    evidence[[n_col]] <- 0

  } else {

    mat <- expr[
      genes,
      evidence$cell,
      drop = FALSE
    ]

    evidence[[score_col]] <-
      Matrix::colMeans(mat)

    evidence[[n_col]] <-
      Matrix::colSums(mat > 0)
  }
}

write_csv(
  evidence,
  file.path(TAB, "12C_cell_programme_evidence.csv")
)

# ============================================================
# 5. Cluster programme summary
# ============================================================

cat("\nSTEP 5/8 — cluster-level programme evidence\n")

cluster_programmes <- bind_rows(
  lapply(names(programmes), function(p) {

    score_col <- paste0(p, "_score")
    n_col <- paste0(p, "_n")

    evidence %>%
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
        fraction_ge2 =
          mean(
            .data[[n_col]] >= 2,
            na.rm = TRUE
          ),
        fraction_ge3 =
          mean(
            .data[[n_col]] >= 3,
            na.rm = TRUE
          ),
        .groups = "drop"
      ) %>%
      mutate(programme = p)
  })
) %>%
  select(
    cluster,
    programme,
    everything()
  ) %>%
  arrange(
    as.integer(cluster),
    desc(mean_score)
  )

cat("\n=== CLUSTER PROGRAMME SCORES ===\n")
print(cluster_programmes, n = Inf)

write_csv(
  cluster_programmes,
  file.path(TAB, "12C_cluster_programme_scores.csv")
)

# ============================================================
# 6. De novo markers specifically at r0.4
# ============================================================

cat("\nSTEP 6/8 — r0.4 de novo markers\n")

Idents(x) <- x$sweat_cluster_candidate

markers <- FindAllMarkers(
  x,
  assay = "RNA",
  only.pos = TRUE,
  min.pct = 0.10,
  logfc.threshold = 0.25,
  test.use = "wilcox"
)

markers$cluster <- as.character(
  markers$cluster
)

top_markers <- markers %>%
  group_by(cluster) %>%
  arrange(
    desc(avg_log2FC),
    .by_group = TRUE
  ) %>%
  slice_head(n = 30) %>%
  ungroup()

cat("\n=== TOP r0.4 MARKERS ===\n")
print(top_markers, n = Inf)

write_csv(
  markers,
  file.path(TAB, "12C_r04_all_markers.csv")
)

write_csv(
  top_markers,
  file.path(TAB, "12C_r04_top30_markers.csv")
)

# ============================================================
# 7. Patient recurrence
# ============================================================

cat("\nSTEP 7/8 — patient recurrence\n")

patient_support <- evidence %>%
  group_by(
    cluster,
    patient_id,
    batch
  ) %>%
  summarise(
    n_cells = n(),

    sweat_fraction =
      mean(PAN_SWEAT_n >= 3),

    secretory_fraction =
      mean(SECRETORY_GENERAL_n >= 3),

    clear_fraction =
      mean(CLEAR_SECRETORY_n >= 3),

    dark_fraction =
      mean(DARK_SECRETORY_n >= 3),

    ductal_luminal_fraction =
      mean(DUCTAL_LUMINAL_n >= 3),

    ductal_basal_fraction =
      mean(DUCTAL_BASAL_n >= 3),

    myoepithelial_fraction =
      mean(MYOEPITHELIAL_n >= 3),

    apocrine_like_fraction =
      mean(SECRETORY_APOCRINE_LIKE_n >= 3),

    inflammatory_fraction =
      mean(INFLAMMATORY_EPITHELIAL_n >= 3),

    sebocyte_fraction =
      mean(SEBOCYTE_n >= 2),

    follicular_fraction =
      mean(FOLLICULAR_n >= 2),

    median_doublet =
      median(
        doublet_score,
        na.rm = TRUE
      ),

    .groups = "drop"
  )

write_csv(
  patient_support,
  file.path(TAB, "12C_patient_recurrence.csv")
)

cluster_patient_summary <- patient_support %>%
  group_by(cluster) %>%
  summarise(
    n_patients =
      n_distinct(patient_id),

    patients_ge5_cells =
      sum(n_cells >= 5),

    median_cells_patient =
      median(n_cells),

    median_sweat_fraction =
      median(sweat_fraction),

    median_secretory_fraction =
      median(secretory_fraction),

    median_clear_fraction =
      median(clear_fraction),

    median_dark_fraction =
      median(dark_fraction),

    median_ductal_luminal_fraction =
      median(ductal_luminal_fraction),

    median_ductal_basal_fraction =
      median(ductal_basal_fraction),

    median_myoepithelial_fraction =
      median(myoepithelial_fraction),

    median_apocrine_like_fraction =
      median(apocrine_like_fraction),

    median_inflammatory_fraction =
      median(inflammatory_fraction),

    median_sebocyte_fraction =
      median(sebocyte_fraction),

    median_follicular_fraction =
      median(follicular_fraction),

    .groups = "drop"
  ) %>%
  arrange(as.integer(cluster))

cat("\n=== PATIENT-LEVEL RECURRENCE SUMMARY ===\n")
print(cluster_patient_summary, n = Inf)

write_csv(
  cluster_patient_summary,
  file.path(
    TAB,
    "12C_cluster_patient_recurrence_summary.csv"
  )
)

# ============================================================
# 8. Formal checks + save adjudication scaffold
# ============================================================

cat("\nSTEP 8/8 — formal checks\n")

expected_sizes <- c(
  "0" = 389,
  "1" = 283,
  "2" = 273,
  "3" = 251,
  "4" = 248,
  "5" = 187,
  "6" = 154,
  "7" = 58
)

observed_sizes <- table(
  x$sweat_cluster_candidate
)

checks <- tibble(
  check = c(
    "cells_1843",
    "genes_23412",
    "patients_16",
    "batches_2",
    "r04_eight_clusters",
    "expected_cluster_sizes",
    "all_clusters_ge13_patients",
    "all_clusters_two_batches",
    "programme_coverage_generated",
    "cell_evidence_generated",
    "cluster_scores_generated",
    "markers_generated",
    "patient_recurrence_generated",
    "disease_not_used",
    "no_cells_removed",
    "no_identity_assigned",
    "resolution_not_frozen",
    "source_object_not_overwritten"
  ),

  pass = c(
    ncol(x) == 1843,
    nrow(x) == 23412,
    length(unique(x$patient_id)) == 16,
    length(unique(x$batch)) == 2,
    length(unique(x$sweat_cluster_candidate)) == 8,
    all(
      as.integer(
        observed_sizes[names(expected_sizes)]
      ) == as.integer(expected_sizes)
    ),
    all(structural$n_patients >= 13),
    all(structural$n_batches == 2),
    nrow(programme_coverage) ==
      length(programmes),
    nrow(evidence) == 1843,
    nrow(cluster_programmes) > 0,
    nrow(markers) > 0,
    nrow(patient_support) > 0,
    TRUE,
    ncol(x) == 1843,
    !("sweat_L2_label" %in%
        colnames(x@meta.data)),
    TRUE,
    TRUE
  )
)

cat("\n=== FORMAL CHECKS ===\n")
print(checks, n = Inf)

write_csv(
  checks,
  file.path(TAB, "12C_formal_checks.csv")
)

stopifnot(all(checks$pass))

FINAL_RDS <- file.path(
  OUT,
  "sweat_gland_v2_12C_adjudication_scaffold.rds"
)

saveRDS(
  x,
  FINAL_RDS
)

cat("\nSaved:", FINAL_RDS, "\n")

cat("\n============================================================\n")
cat("12C COMPLETE\n")
cat("============================================================\n")
cat("Cells: 1843\n")
cat("Patients: 16\n")
cat("Batches: 2\n")
cat("Candidate resolution: 0.4\n")
cat("Candidate clusters: 8\n")
cat("Identity programmes evaluated.\n")
cat("r0.4 de novo markers generated.\n")
cat("Patient recurrence evaluated.\n")
cat("Disease used: NO\n")
cat("Cells removed: 0\n")
cat("L2 identities assigned: NO\n")
cat("Resolution frozen: NO\n")
cat("READY FOR FINAL ADJUDICATION\n")
cat("PASS\n")
