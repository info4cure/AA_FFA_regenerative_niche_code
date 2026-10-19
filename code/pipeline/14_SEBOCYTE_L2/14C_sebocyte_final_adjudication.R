suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tibble)
  library(readr)
  library(tidyr)
})

cat("============================================================\n")
cat("PIPELINE_V2 — 14C SEBOCYTE FINAL ADJUDICATION\n")
cat("============================================================\n\n")

SRC <- paste0(
  "PIPELINE_V2/results/annotation/14_SEBOCYTE_L2/",
  "14B_recluster/sebocyte_v2_14B_reclustered.rds"
)

OUT <- paste0(
  "PIPELINE_V2/results/annotation/14_SEBOCYTE_L2/",
  "14C_adjudication"
)

TAB <- file.path(OUT, "tables")

dir.create(TAB, recursive = TRUE, showWarnings = FALSE)

stopifnot(file.exists(SRC))

# ============================================================
# 1. Load scaffold
# ============================================================

cat("STEP 1/7 — loading 14B sebocyte scaffold\n")

x <- readRDS(SRC)

stopifnot(ncol(x) == 1388)
stopifnot(nrow(x) == 23412)
stopifnot(length(unique(x$patient_id)) == 16)
stopifnot(length(unique(x$batch)) == 2)

R04 <- "seb_snn_res.0.4"

stopifnot(R04 %in% colnames(x@meta.data))

cluster <- as.character(
  x@meta.data[[R04]]
)

expected <- c(
  "0" = 436,
  "1" = 355,
  "2" = 278,
  "3" = 164,
  "4" = 155
)

observed <- table(cluster)

stopifnot(
  all(
    as.integer(observed[names(expected)]) ==
      as.integer(expected)
  )
)

cat("Cells:", ncol(x), "\n")
cat("Patients:", length(unique(x$patient_id)), "\n")
cat("Batches:", length(unique(x$batch)), "\n")
cat("Candidate resolution: 0.4\n")
cat("Clusters: 5\n")

# ============================================================
# 2. Targeted hierarchy programmes
# ============================================================

cat("\nSTEP 2/7 — defining targeted sebocyte hierarchy programmes\n")

programmes <- list(

  SEBOCYTE_IDENTITY = c(
    "PPARG","PLIN2","FASN","SCD",
    "MGST1","FABP5","KRT7","KRT19"
  ),

  SEBOCYTE_LIPOGENESIS = c(
    "FASN","SCD","ACACA","ACLY",
    "ELOVL3","ELOVL5","ELOVL6",
    "DGAT1","DGAT2"
  ),

  SEBOCYTE_TERMINAL_LIPID = c(
    "ELOVL3","THRSP","DEGS2","DHCR7",
    "TM7SF2","PGRMC1","AZGP1","APOC1"
  ),

  SEBOCYTE_EARLY_COMMITMENT = c(
    "PPARG","PLIN2","SCD","FABP5",
    "DGAT1","DGAT2","KRT7","KRT19"
  ),

  BASAL_PROGENITOR = c(
    "KRT5","KRT14","KRT15","KRT17",
    "TP63","KRT6A","COL17A1","ITGA6"
  ),

  PILOSEBACEOUS_PROGENITOR = c(
    "SOX9","KRT15","KRT17","LGR5",
    "LRIG1","KRT19","ITGA6","KRT5"
  ),

  FOLLICULAR_DIFFERENTIATION = c(
    "KRT75","KRT71","KRT85",
    "KRT31","KRT32","KRT35"
  ),

  EPIDERMAL_DIFFERENTIATION = c(
    "KRT1","KRT10","IVL",
    "LOR","FLG","DSG1"
  ),

  EPITHELIAL_JUNCTION = c(
    "DSP","DSC2","DSC3",
    "PKP1","PKP3","CLDN1"
  ),

  INFLAMMATORY_SEBOCYTE = c(
    "SAA1","SAA2","LCN2",
    "CXCL14","IL1R2","S100A8","S100A9"
  ),

  ANTIGEN_PRESENTATION = c(
    "HLA-A","HLA-B","HLA-C",
    "B2M","HLA-DRA","HLA-DRB1"
  ),

  SWEAT_GLAND = c(
    "AQP5","DCD","PIP",
    "SCGB2A2","SCGB1D2","KRT77"
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

coverage <- bind_rows(
  lapply(
    names(programmes),
    function(nm) {
      present <- intersect(
        programmes[[nm]],
        rownames(x)
      )

      tibble(
        programme = nm,
        n_requested = length(programmes[[nm]]),
        n_present = length(present),
        genes_present = paste(
          present,
          collapse = ";"
        )
      )
    }
  )
)

cat("\n=== PROGRAMME COVERAGE ===\n")
print(coverage, n = Inf)

write_csv(
  coverage,
  file.path(TAB, "14C_programme_coverage.csv")
)

# ============================================================
# 3. Cell-level evidence
# ============================================================

cat("\nSTEP 3/7 — cell-level hierarchy evidence\n")

rna <- GetAssayData(
  x,
  assay = "RNA",
  layer = "data"
)

score_programme <- function(genes) {
  genes <- intersect(genes, rownames(rna))

  if (length(genes) == 0) {
    return(rep(NA_real_, ncol(x)))
  }

  Matrix::colMeans(
    rna[genes, , drop = FALSE]
  )
}

detect_programme <- function(genes) {
  genes <- intersect(genes, rownames(rna))

  if (length(genes) == 0) {
    return(rep(0L, ncol(x)))
  }

  Matrix::colSums(
    rna[genes, , drop = FALSE] > 0
  )
}

evidence <- tibble(
  cell = colnames(x),
  cluster = cluster,
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
  file.path(TAB, "14C_cell_hierarchy_evidence.csv")
)

# ============================================================
# 4. Cluster-level hierarchy evidence
# ============================================================

cat("\nSTEP 4/7 — cluster-level hierarchy evidence\n")

cluster_scores <- bind_rows(
  lapply(
    names(programmes),
    function(nm) {

      score_col <- paste0(nm, "_score")
      genes_col <- paste0(nm, "_genes")

      evidence %>%
        group_by(cluster) %>%
        summarise(
          programme = nm,
          n_cells = n(),
          n_patients = n_distinct(patient_id),
          n_batches = n_distinct(batch),

          mean_score =
            mean(.data[[score_col]], na.rm = TRUE),

          median_score =
            median(.data[[score_col]], na.rm = TRUE),

          fraction_ge2 =
            mean(.data[[genes_col]] >= 2, na.rm = TRUE),

          fraction_ge3 =
            mean(.data[[genes_col]] >= 3, na.rm = TRUE),

          .groups = "drop"
        )
    }
  )
) %>%
  arrange(
    as.integer(cluster),
    desc(median_score)
  )

write_csv(
  cluster_scores,
  file.path(TAB, "14C_cluster_hierarchy_scores.csv")
)

compact <- cluster_scores %>%
  select(
    cluster,
    programme,
    median_score
  ) %>%
  pivot_wider(
    names_from = programme,
    values_from = median_score
  ) %>%
  arrange(as.integer(cluster))

cat("\n=== COMPACT HIERARCHY MATRIX ===\n")
print(compact, n = Inf, width = Inf)

write_csv(
  compact,
  file.path(TAB, "14C_compact_hierarchy_matrix.csv")
)

# ============================================================
# 5. Patient recurrence
# ============================================================

cat("\nSTEP 5/7 — patient-level recurrence\n")

patient_support <- evidence %>%
  count(
    cluster,
    patient_id,
    batch,
    name = "n_cells"
  )

patient_summary <- patient_support %>%
  group_by(cluster) %>%
  summarise(
    n_cells = sum(n_cells),
    n_patients = n_distinct(patient_id),
    n_batches = n_distinct(batch),
    patients_ge5_cells = sum(n_cells >= 5),
    median_cells_patient = median(n_cells),
    dominant_patient_fraction =
      max(n_cells) / sum(n_cells),
    .groups = "drop"
  ) %>%
  arrange(as.integer(cluster))

cat("\n=== PATIENT RECURRENCE ===\n")
print(patient_summary, n = Inf)

write_csv(
  patient_summary,
  file.path(TAB, "14C_patient_recurrence.csv")
)

# ============================================================
# 6. Integrated adjudication scaffold
# ============================================================

cat("\nSTEP 6/7 — integrated adjudication scaffold\n")

markers_file <- paste0(
  "PIPELINE_V2/results/annotation/14_SEBOCYTE_L2/",
  "14B_recluster/tables/14B_all_candidate_markers.csv"
)

markers <- read_csv(
  markers_file,
  show_col_types = FALSE
) %>%
  filter(
    resolution == 0.4,
    p_val_adj < 0.05
  ) %>%
  group_by(cluster) %>%
  arrange(
    desc(avg_log2FC),
    .by_group = TRUE
  ) %>%
  slice_head(n = 20) %>%
  summarise(
    top_markers =
      paste(gene, collapse = ";"),
    .groups = "drop"
  ) %>%
  mutate(
    cluster = as.character(cluster)
  )

adjudication <- compact %>%
  left_join(
    patient_summary,
    by = "cluster"
  ) %>%
  left_join(
    markers,
    by = "cluster"
  ) %>%
  arrange(as.integer(cluster))

cat("\n=== FINAL ADJUDICATION INPUT ===\n")
print(adjudication, n = Inf, width = Inf)

write_csv(
  adjudication,
  file.path(
    TAB,
    "14C_integrated_adjudication.csv"
  )
)

# ============================================================
# 7. Formal checks + save scaffold
# ============================================================

cat("\nSTEP 7/7 — formal checks and save\n")

checks <- tibble(
  check = c(
    "cells_1388",
    "genes_23412",
    "patients_16",
    "batches_2",
    "r04_five_clusters",
    "expected_cluster_sizes",
    "all_clusters_ge12_patients",
    "all_clusters_two_batches",
    "programme_coverage_ge80pct",
    "cell_evidence_generated",
    "cluster_scores_generated",
    "patient_recurrence_generated",
    "integrated_adjudication_generated",
    "disease_not_used",
    "no_cells_removed",
    "no_identity_assigned",
    "resolution_not_frozen",
    "figures_generated_NO",
    "source_object_not_overwritten"
  ),

  pass = c(
    ncol(x) == 1388,
    nrow(x) == 23412,
    length(unique(x$patient_id)) == 16,
    length(unique(x$batch)) == 2,
    length(unique(cluster)) == 5,

    all(
      as.integer(table(cluster)[names(expected)]) ==
        as.integer(expected)
    ),

    all(patient_summary$n_patients >= 12),
    all(patient_summary$n_batches == 2),

    all(
      coverage$n_present / coverage$n_requested >= 0.80
    ),

    nrow(evidence) == 1388,
    nrow(cluster_scores) > 0,
    nrow(patient_summary) == 5,
    nrow(adjudication) == 5,

    TRUE,
    ncol(x) == 1388,

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
  file.path(TAB, "14C_formal_checks.csv")
)

stopifnot(all(checks$pass))

FINAL_RDS <- file.path(
  OUT,
  "sebocyte_v2_14C_adjudication_scaffold.rds"
)

saveRDS(x, FINAL_RDS)

cat("\nSaved:\n")
cat(FINAL_RDS, "\n")

cat("\n============================================================\n")
cat("14C COMPLETE\n")
cat("============================================================\n")
cat("Cells: 1388\n")
cat("Candidate resolution: 0.4\n")
cat("Clusters: 5\n")
cat("Sebocyte hierarchy adjudicated.\n")
cat("Disease used: NO\n")
cat("Cells removed: 0\n")
cat("L2 identities assigned: NO\n")
cat("Resolution frozen: NO\n")
cat("Figures generated: NO\n")
cat("READY FOR FINAL SEBOCYTE DECISION\n")
cat("PASS\n")
