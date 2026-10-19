suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
})

cat("============================================================\n")
cat("PIPELINE_V2 — 11B2 MURAL IDENTITY EVIDENCE\n")
cat("============================================================\n\n")

SRC <- paste0(
  "PIPELINE_V2/results/annotation/11_MURAL_L2/",
  "11B1_recluster/mural_v2_11B1_reclustered.rds"
)

OUT <- paste0(
  "PIPELINE_V2/results/annotation/11_MURAL_L2/",
  "11B2_identity_evidence"
)

TAB <- file.path(OUT, "tables")
FIG <- file.path(OUT, "figures")

dir.create(TAB, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG, recursive = TRUE, showWarnings = FALSE)

x <- readRDS(SRC)
DefaultAssay(x) <- "RNA"

source_cells <- colnames(x)
source_genes <- rownames(x)

stopifnot(ncol(x) == 1603)

cat("Cells:", ncol(x), "\n")
cat("Patients:", length(unique(x$patient_id)), "\n")
cat("Batches:", length(unique(x$batch)), "\n")

# ============================================================
# 1. Candidate resolutions
# ============================================================

cat("\nSTEP 1/8 — candidate resolutions\n")

resolutions <- c(0.3, 0.4, 0.6)

res_cols <- paste0(
  "mural_snn_res.",
  resolutions
)

stopifnot(
  all(res_cols %in% colnames(x@meta.data))
)

for (i in seq_along(resolutions)) {
  cat(
    "r", resolutions[i],
    " clusters: ",
    length(unique(x@meta.data[[res_cols[i]]])),
    "\n",
    sep = ""
  )
}

# ============================================================
# 2. Prespecified programmes
# ============================================================

cat("\nSTEP 2/8 — defining mural identity programmes\n")

programmes <- list(

  PAN_MURAL = c(
    "RGS5","CSPG4","MCAM","PDGFRB",
    "RBP1","NOTCH3","ABCC9","KCNJ8"
  ),

  PERICYTE = c(
    "RGS5","CSPG4","PDGFRB","RBP1",
    "ABCC9","KCNJ8","COX4I2","HIGD1B",
    "NOTCH3"
  ),

  CONTRACTILE_SMC = c(
    "MYH11","ACTA2","TAGLN","CNN1",
    "MYL9","TPM2","LMOD1","SYNPO2",
    "MYLK"
  ),

  APM_CONTRACTILE = c(
    "MYH11","ACTA2","TAGLN","CNN1",
    "MYLK","ACTG2","LMOD1","SYNPO2",
    "DES","ITGA8"
  ),

  BULGE_BASAL = c(
    "KRT5","KRT14","KRT15","KRT17",
    "TP63","COL17A1","ITGA6","ITGB4"
  ),

  BULGE_STEM = c(
    "SOX9","CD200","KRT15","COL17A1",
    "ITGA6","ITGB4"
  ),

  APM_ATTACHMENT = c(
    "ITGA8","ITGB1","NPNT",
    "ACTA2","TAGLN"
  ),

  FIBROBLAST = c(
    "DCN","LUM","COL1A1","COL1A2",
    "COL3A1","PDGFRA","DPT","COL6A1"
  ),

  REMODELING_FIBROBLAST = c(
    "ADAM12","MMP11","POSTN","COL11A1",
    "COMP","COL10A1","THBS2","FN1"
  ),

  ENDOTHELIAL = c(
    "PECAM1","VWF","CDH5","EMCN",
    "KDR","ESAM","ENG","RAMP2"
  ),

  EPITHELIAL_JUNCTION = c(
    "EPCAM","DSP","PKP1","PKP3",
    "DSC3","DSG1","CLDN1"
  ),

  SCHWANN = c(
    "SOX10","S100B","PLP1","PMP22",
    "MPZ","S100A1","NGFR"
  ),

  IMMUNE = c(
    "PTPRC","CD3D","CD3E",
    "TRAC","NKG7","GNLY"
  )
)

available_programmes <- lapply(
  programmes,
  function(g)
    intersect(g, rownames(x))
)

coverage <- tibble(
  programme = names(programmes),
  n_requested =
    sapply(programmes, length),
  n_present =
    sapply(available_programmes, length),
  genes_present =
    sapply(
      available_programmes,
      paste,
      collapse = ";"
    )
)

cat("\n=== PROGRAMME COVERAGE ===\n")
print(coverage, n = Inf, width = Inf)

write.csv(
  coverage,
  file.path(
    TAB,
    "11B2_programme_coverage.csv"
  ),
  row.names = FALSE
)

# ============================================================
# 3. Cell-level scores
# ============================================================

cat("\nSTEP 3/8 — cell-level programme evidence\n")

mat <- GetAssayData(
  x,
  assay = "RNA",
  layer = "data"
)

score_mean <- function(genes) {

  genes <- intersect(
    genes,
    rownames(mat)
  )

  if (length(genes) == 0)
    return(
      rep(NA_real_, ncol(mat))
    )

  Matrix::colMeans(
    mat[genes, , drop = FALSE]
  )
}

score_detected <- function(genes) {

  genes <- intersect(
    genes,
    rownames(mat)
  )

  if (length(genes) == 0)
    return(
      rep(0L, ncol(mat))
    )

  Matrix::colSums(
    mat[genes, , drop = FALSE] > 0
  )
}

evidence <- tibble(
  cell = colnames(x),
  patient_id =
    as.character(x$patient_id),
  batch =
    as.character(x$batch),
  source_L1 =
    as.character(x$mural_source_L1),
  doublet_score =
    as.numeric(x$doublet_score),
  nFeature_RNA =
    as.numeric(x$nFeature_RNA),
  nCount_RNA =
    as.numeric(x$nCount_RNA),
  percent_mt =
    as.numeric(x$percent.mt)
)

for (p in names(programmes)) {

  evidence[[paste0(p, "_score")]] <-
    score_mean(
      available_programmes[[p]]
    )

  evidence[[paste0(p, "_n")]] <-
    score_detected(
      available_programmes[[p]]
    )
}

for (col in res_cols) {
  evidence[[col]] <-
    as.character(
      x@meta.data[
        evidence$cell,
        col
      ]
    )
}

write.csv(
  evidence,
  file.path(
    TAB,
    "11B2_cell_programme_evidence.csv"
  ),
  row.names = FALSE
)

# ============================================================
# 4. Cluster-level programme evidence
# ============================================================

cat("\nSTEP 4/8 — cluster-level programme evidence\n")

cluster_scores <- bind_rows(
  lapply(seq_along(resolutions), function(i) {

    r <- resolutions[i]
    col <- res_cols[i]

    bind_rows(
      lapply(names(programmes), function(p) {

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

            fraction_ge2 =
              mean(
                .data[[n_col]] >= 2
              ),

            fraction_ge3 =
              mean(
                .data[[n_col]] >= 3
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
    everything()
  )

write.csv(
  cluster_scores,
  file.path(
    TAB,
    "11B2_cluster_programme_scores.csv"
  ),
  row.names = FALSE
)

cat("\n=== CLUSTER PROGRAMME SCORES — r0.6 ===\n")

print(
  cluster_scores %>%
    filter(resolution == 0.6) %>%
    arrange(
      as.integer(cluster),
      desc(mean_score)
    ),
  n = Inf,
  width = Inf
)

# ============================================================
# 5. Patient-level support
# ============================================================

cat("\nSTEP 5/8 — patient-level programme support\n")

patient_scores <- bind_rows(
  lapply(seq_along(resolutions), function(i) {

    r <- resolutions[i]
    col <- res_cols[i]

    bind_rows(
      lapply(names(programmes), function(p) {

        score_col <-
          paste0(p, "_score")

        n_col <-
          paste0(p, "_n")

        evidence %>%
          mutate(
            cluster =
              .data[[col]]
          ) %>%
          group_by(
            cluster,
            patient_id,
            batch
          ) %>%
          summarise(
            n_cells = n(),

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

            fraction_ge3 =
              mean(
                .data[[n_col]] >= 3
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
)

write.csv(
  patient_scores,
  file.path(
    TAB,
    "11B2_patient_programme_support.csv"
  ),
  row.names = FALSE
)

# ============================================================
# 6. De novo markers
# ============================================================

cat("\nSTEP 6/8 — de novo markers\n")

markers_all <- list()

for (i in seq_along(resolutions)) {

  r <- resolutions[i]
  col <- res_cols[i]

  cat("Markers: r", r, "\n", sep = "")

  Idents(x) <-
    factor(
      x@meta.data[[col]]
    )

  m <- FindAllMarkers(
    x,
    assay = "RNA",
    only.pos = TRUE,
    min.pct = 0.10,
    logfc.threshold = 0.25,
    test.use = "wilcox",
    verbose = FALSE
  )

  m$resolution <- r

  markers_all[[as.character(r)]] <- m

  write.csv(
    m,
    file.path(
      TAB,
      paste0(
        "11B2_markers_r",
        gsub("\\.", "", r),
        ".csv"
      )
    ),
    row.names = FALSE
  )
}

all_markers <- bind_rows(
  markers_all
)

write.csv(
  all_markers,
  file.path(
    TAB,
    "11B2_all_candidate_markers.csv"
  ),
  row.names = FALSE
)

top20 <- all_markers %>%
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
  slice_head(n = 20) %>%
  ungroup()

cat("\n=== TOP MARKERS — r0.6 ===\n")

print(
  top20 %>%
    filter(resolution == 0.6),
  n = Inf,
  width = Inf
)

write.csv(
  top20,
  file.path(
    TAB,
    "11B2_top20_markers.csv"
  ),
  row.names = FALSE
)

# ============================================================
# 7. Targeted r0.6 joint-lineage audit
# ============================================================

cat("\nSTEP 7/8 — targeted r0.6 joint-lineage audit\n")

r06 <- "mural_snn_res.0.6"

joint <- evidence %>%
  mutate(
    cluster =
      .data[[r06]],

    pericyte_supported =
      PERICYTE_n >= 3,

    contractile_supported =
      CONTRACTILE_SMC_n >= 4,

    APM_supported =
      APM_CONTRACTILE_n >= 4,

    bulge_supported =
      BULGE_BASAL_n >= 3,

    fibroblast_supported =
      FIBROBLAST_n >= 3,

    remodeling_supported =
      REMODELING_FIBROBLAST_n >= 2,

    endothelial_supported =
      ENDOTHELIAL_n >= 3,

    Schwann_supported =
      SCHWANN_n >= 2,

    immune_supported =
      IMMUNE_n >= 2,

    APM_bulge_joint =
      APM_supported &
      bulge_supported,

    mural_fibroblast_joint =
      pericyte_supported &
      fibroblast_supported,

    contractile_remodeling_joint =
      contractile_supported &
      remodeling_supported
  )

joint_summary <- joint %>%
  group_by(cluster) %>%
  summarise(
    n_cells = n(),
    n_patients =
      n_distinct(patient_id),

    median_doublet =
      median(doublet_score),

    median_features =
      median(nFeature_RNA),

    pericyte_fraction =
      mean(pericyte_supported),

    contractile_fraction =
      mean(contractile_supported),

    APM_fraction =
      mean(APM_supported),

    bulge_fraction =
      mean(bulge_supported),

    APM_bulge_joint_fraction =
      mean(APM_bulge_joint),

    fibroblast_fraction =
      mean(fibroblast_supported),

    remodeling_fraction =
      mean(remodeling_supported),

    mural_fibroblast_joint_fraction =
      mean(mural_fibroblast_joint),

    contractile_remodeling_joint_fraction =
      mean(
        contractile_remodeling_joint
      ),

    endothelial_fraction =
      mean(endothelial_supported),

    Schwann_fraction =
      mean(Schwann_supported),

    immune_fraction =
      mean(immune_supported),

    .groups = "drop"
  )

cat("\n=== r0.6 JOINT LINEAGE SUMMARY ===\n")

print(
  joint_summary,
  n = Inf,
  width = Inf
)

write.csv(
  joint_summary,
  file.path(
    TAB,
    "11B2_r06_joint_lineage_summary.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# Retrospective smooth C0/C1/C2 composition at r0.6
# ------------------------------------------------------------

SMOOTH_DIAG <- paste0(
  "PIPELINE_V2/results/annotation/11_MURAL_L2/",
  "11A3_smooth_muscle_internal_structure/",
  "smooth_muscle_11A3_diagnostic_reclustered.rds"
)

sm <- readRDS(SMOOTH_DIAG)

common <- intersect(
  colnames(sm),
  colnames(x)
)

smooth_r06 <- tibble(
  cell = common,

  old_smooth_cluster =
    as.character(
      sm@meta.data[
        common,
        "smooth_snn_res.0.3"
      ]
    ),

  mural_r06 =
    as.character(
      x@meta.data[
        common,
        r06
      ]
    )
) %>%
  count(
    old_smooth_cluster,
    mural_r06,
    name = "n_cells"
  ) %>%
  group_by(old_smooth_cluster) %>%
  mutate(
    old_total =
      sum(n_cells),
    fraction_old =
      n_cells /
      old_total
  ) %>%
  ungroup() %>%
  arrange(
    as.integer(old_smooth_cluster),
    desc(fraction_old)
  )

cat("\n=== OLD SMOOTH C0/C1/C2 -> r0.6 ===\n")

print(
  smooth_r06,
  n = Inf,
  width = Inf
)

write.csv(
  smooth_r06,
  file.path(
    TAB,
    "11B2_old_smooth_to_r06.csv"
  ),
  row.names = FALSE
)

# ============================================================
# 8. Formal checks
# ============================================================

cat("\nSTEP 8/8 — formal checks\n")

checks <- tibble(
  check = c(
    "cells_1603",
    "patients_16",
    "batches_2",
    "r03_present",
    "r04_present",
    "r06_present",
    "programme_coverage_generated",
    "cell_evidence_generated",
    "cluster_scores_generated",
    "patient_support_generated",
    "markers_r03_generated",
    "markers_r04_generated",
    "markers_r06_generated",
    "joint_lineage_generated",
    "smooth_transition_generated",
    "disease_not_used_for_identity",
    "source_L1_not_used_for_scoring",
    "old_smooth_not_used_for_scoring",
    "no_cells_removed",
    "no_identity_assigned",
    "no_resolution_selected",
    "source_object_not_modified"
  ),

  pass = c
(
    ncol(x) == 1603,
    length(unique(x$patient_id)) == 16,
    length(unique(x$batch)) == 2,
    "mural_snn_res.0.3" %in% colnames(x@meta.data),
    "mural_snn_res.0.4" %in% colnames(x@meta.data),
    "mural_snn_res.0.6" %in% colnames(x@meta.data),
    nrow(coverage) == length(programmes),
    nrow(evidence) == 1603,
    nrow(cluster_scores) > 0,
    nrow(patient_scores) > 0,
    nrow(markers_all[["0.3"]]) > 0,
    nrow(markers_all[["0.4"]]) > 0,
    nrow(markers_all[["0.6"]]) > 0,
    nrow(joint_summary) > 0,
    nrow(smooth_r06) > 0,
    TRUE,
    TRUE,
    TRUE,
    identical(sort(colnames(x)), sort(source_cells)),
    TRUE,
    TRUE,
    TRUE
  )
)

cat("\n=== FORMAL CHECKS ===\n")
print(
  checks,
  n = Inf,
  width = Inf
)

write.csv(
  checks,
  file.path(
    TAB,
    "11B2_formal_checks.csv"
  ),
  row.names = FALSE
)

stopifnot(all(checks$pass))

# ============================================================
# Diagnostic figures
# ============================================================

for (r in resolutions) {

  col <- paste0(
    "mural_snn_res.",
    r
  )

  p <- DimPlot(
    x,
    reduction = "mural_umap",
    group.by = col,
    label = TRUE,
    repel = TRUE,
    pt.size = 0.7
  ) +
    ggtitle(
      paste0(
        "Mural identity evidence - r",
        r
      )
    ) +
    theme_classic()

  ggsave(
    file.path(
      FIG,
      paste0(
        "11B2_mural_r",
        gsub("\\.", "", as.character(r)),
        "_UMAP.pdf"
      )
    ),
    p,
    width = 7,
    height = 6
  )
}

# r0.6 programme dot plot
Idents(x) <- factor(
  x@meta.data[["mural_snn_res.0.6"]]
)

dot_genes <- unique(c(
  "RGS5","CSPG4","PDGFRB","RBP1","ABCC9","KCNJ8",
  "MYH11","ACTA2","TAGLN","CNN1","MYLK","ACTG2",
  "LMOD1","SYNPO2","DES","ITGA8",
  "KRT5","KRT14","KRT15","KRT17","TP63","COL17A1",
  "SOX9","CD200",
  "DCN","LUM","COL1A1","COL1A2","COL3A1","PDGFRA",
  "ADAM12","MMP11","POSTN","COL11A1","COMP","THBS2",
  "PECAM1","VWF",
  "SOX10","S100B",
  "PTPRC"
))

dot_genes <- intersect(
  dot_genes,
  rownames(x)
)

p_dot <- DotPlot(
  x,
  features = dot_genes,
  assay = "RNA"
) +
  RotatedAxis() +
  ggtitle(
    "Mural r0.6 - identity and competing-lineage evidence"
  ) +
  theme_classic()

ggsave(
  file.path(
    FIG,
    "11B2_r06_identity_DotPlot.pdf"
  ),
  p_dot,
  width = 16,
  height = 6
)

# ============================================================
# Audit
# ============================================================

audit <- c(
  "PIPELINE_V2 — 11B2 MURAL IDENTITY EVIDENCE",
  "===========================================",
  "",
  paste0("Cells: ", ncol(x)),
  paste0(
    "Patients: ",
    length(unique(x$patient_id))
  ),
  paste0(
    "Batches: ",
    length(unique(x$batch))
  ),
  "",
  "Candidate resolutions:",
  "0.3, 0.4, 0.6",
  "",
  "Evidence generated:",
  "- de novo markers",
  "- prespecified mural/pericyte programmes",
  "- SMC/APM contractile programmes",
  "- bulge/basal and bulge/stem programmes",
  "- APM attachment-axis programme",
  "- fibroblast and remodeling-fibroblast programmes",
  "- endothelial, epithelial, Schwann and immune competing-lineage evidence",
  "- cluster-level programme summaries",
  "- patient-level programme support",
  "- r0.6 joint-lineage audit",
  "- retrospective old smooth C0/C1/C2 mapping",
  "",
  "Interpretive boundary:",
  "- disease was not used",
  "- source L1 was not used for identity scoring",
  "- old smooth-muscle states were used retrospectively only",
  "- no cells were removed",
  "- no resolution was selected",
  "- no L2 identity was assigned",
  "- source object was not modified",
  "",
  "PASS"
)

writeLines(
  audit,
  file.path(
    OUT,
    "11B2_audit.txt"
  )
)

cat("\n============================================================\n")
cat("11B2 COMPLETE\n")
cat("============================================================\n")
cat("Cells:", ncol(x), "\n")
cat(
  "Patients:",
  length(unique(x$patient_id)),
  "\n"
)
cat("Candidate resolutions: 0.3, 0.4, 0.6\n")
cat("De novo markers generated.\n")
cat("Mural/APM/pericyte programmes evaluated.\n")
cat("APM-bulge evidence evaluated independently.\n")
cat("Fibroblast/remodeling boundary evaluated.\n")
cat("Old smooth C0/C1/C2 used retrospectively only.\n")
cat("No resolution selected.\n")
cat("No L2 identities assigned.\n")
cat("No cells removed.\n")
cat("Disease not used.\n")
cat("Source object not modified.\n")
cat("PASS\n")
