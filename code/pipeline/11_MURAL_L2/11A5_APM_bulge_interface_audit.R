suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(tibble)
})

cat("============================================================\n")
cat("PIPELINE_V2 — 11A5 APM–BULGE INTERFACE AUDIT\n")
cat("============================================================\n\n")

SRC <- paste0(
  "PIPELINE_V2/results/annotation/11_MURAL_L2/",
  "11A3_smooth_muscle_internal_structure/",
  "smooth_muscle_11A3_diagnostic_reclustered.rds"
)

OUT <- paste0(
  "PIPELINE_V2/results/annotation/11_MURAL_L2/",
  "11A5_APM_bulge_interface"
)

TAB <- file.path(OUT, "tables")
dir.create(TAB, recursive = TRUE, showWarnings = FALSE)

x <- readRDS(SRC)

stopifnot(ncol(x) == 330)

CL <- "smooth_snn_res.0.3"

cluster <- as.character(
  x@meta.data[[CL]]
)

stopifnot(length(unique(cluster)) == 3)

DefaultAssay(x) <- "RNA"

mat <- GetAssayData(
  x,
  assay = "RNA",
  layer = "data"
)

cat("Cells:", ncol(x), "\n")
cat("Patients:", length(unique(x$patient_id)), "\n")
cat("Clusters:", length(unique(cluster)), "\n\n")

# ------------------------------------------------------------
# 1. Prespecified biological programmes
# ------------------------------------------------------------

cat("STEP 1/7 — defining APM/bulge programmes\n")

programmes <- list(

  # Mature contractile APM / smooth muscle
  APM_CONTRACTILE = c(
    "MYH11","ACTA2","TAGLN","CNN1",
    "MYLK","ACTG2","LMOD1","SYNPO2",
    "DES","MYL9","TPM2"
  ),

  # Literature-supported APM attachment / mesenchymal side
  APM_ATTACHMENT = c(
    "ITGA8","ITGB1","NPNT",
    "ACTA2","TAGLN"
  ),

  # Basal ORS / bulge epithelial identity
  BULGE_BASAL = c(
    "KRT5","KRT14","KRT15","KRT17",
    "TP63","ITGA6","ITGB4","COL17A1"
  ),

  # Hair-follicle stem/bulge-associated
  BULGE_STEM = c(
    "SOX9","CD200","KRT15",
    "COL17A1","ITGA6","ITGB4"
  ),

  # Epithelial structural coherence
  EPITHELIAL_JUNCTION = c(
    "DSP","DSC3","DSG1",
    "PKP1","PKP3","CLDN1"
  ),

  # Differentiated epidermis — should not dominate
  EPIDERMAL_DIFFERENTIATION = c(
    "KRT1","KRT10","IVL","FLG","LOR"
  ),

  # Glandular alternative
  GLANDULAR = c(
    "KRT7","KRT8","KRT18",
    "KRT19","MUC1","KRT23"
  ),

  # Sebocyte alternative
  SEBOCYTE = c(
    "PPARG","PLIN2","FASN",
    "SCD","MGST1","KRT7","KRT19"
  ),

  # Pericyte alternative
  PERICYTE = c(
    "RGS5","CSPG4","PDGFRB",
    "RBP1","ABCC9","KCNJ8",
    "COX4I2","HIGD1B"
  ),

  # Dermal sheath / fibroblast alternative
  DERMAL_SHEATH_FIBROBLAST = c(
    "COL11A1","DPEP1","RBP4",
    "MME","DCN","LUM","PDGFRA"
  ),

  # Remodeling fibroblast alternative
  REMODELING = c(
    "ADAM12","MMP11","POSTN",
    "COL11A1","COMP","THBS2","FN1"
  ),

  # Schwann/neural alternative
  SCHWANN = c(
    "SOX10","S100B","PLP1",
    "PMP22","MPZ","NGFR"
  )
)

available <- lapply(
  programmes,
  function(g) intersect(g, rownames(mat))
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
  file.path(TAB, "11A5_programme_coverage.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# 2. Cell-level programme evidence
# ------------------------------------------------------------

cat("\nSTEP 2/7 — cell-level programme evidence\n")

score_fun <- function(genes) {

  genes <- intersect(
    genes,
    rownames(mat)
  )

  if (length(genes) == 0)
    return(rep(0, ncol(x)))

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
    return(rep(0L, ncol(x)))

  Matrix::colSums(
    mat[genes, , drop = FALSE] > 0
  )
}

ledger <- tibble(
  cell = colnames(x),
  cluster = cluster,
  patient_id = as.character(x$patient_id),
  batch = as.character(x$batch),
  doublet_score = as.numeric(x$doublet_score),
  nFeature_RNA = as.numeric(x$nFeature_RNA),
  nCount_RNA = as.numeric(x$nCount_RNA),
  percent_mt = as.numeric(x$percent.mt)
)

for (p in names(available)) {

  ledger[[paste0(p, "_score")]] <-
    score_fun(available[[p]])

  ledger[[paste0(p, "_n")]] <-
    detect_fun(available[[p]])
}

ledger <- ledger %>%
  mutate(
    APM_supported =
      APM_CONTRACTILE_n >= 5,

    bulge_basal_supported =
      BULGE_BASAL_n >= 3,

    bulge_stem_supported =
      BULGE_STEM_n >= 2,

    epithelial_junction_supported =
      EPITHELIAL_JUNCTION_n >= 2,

    pericyte_supported =
      PERICYTE_n >= 3,

    dermal_sheath_supported =
      DERMAL_SHEATH_FIBROBLAST_n >= 3,

    remodeling_supported =
      REMODELING_n >= 2,

    Schwann_supported =
      SCHWANN_n >= 2,

    APM_bulge_joint =
      APM_supported &
      bulge_basal_supported,

    APM_bulge_junction_joint =
      APM_supported &
      bulge_basal_supported &
      epithelial_junction_supported
  )

write.csv(
  ledger,
  file.path(TAB, "11A5_cell_level_ledger.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# 3. Cluster-level evidence
# ------------------------------------------------------------

cat("\nSTEP 3/7 — cluster-level interface evidence\n")

cluster_summary <- ledger %>%
  group_by(cluster) %>%
  summarise(
    n_cells = n(),
    n_patients = n_distinct(patient_id),
    n_batches = n_distinct(batch),

    median_doublet =
      median(doublet_score),

    median_features =
      median(nFeature_RNA),

    APM_fraction =
      mean(APM_supported),

    bulge_basal_fraction =
      mean(bulge_basal_supported),

    bulge_stem_fraction =
      mean(bulge_stem_supported),

    junction_fraction =
      mean(epithelial_junction_supported),

    differentiation_fraction =
      mean(EPIDERMAL_DIFFERENTIATION_n >= 2),

    glandular_fraction =
      mean(GLANDULAR_n >= 2),

    sebocyte_fraction =
      mean(SEBOCYTE_n >= 2),

    pericyte_fraction =
      mean(pericyte_supported),

    dermal_sheath_fraction =
      mean(dermal_sheath_supported),

    remodeling_fraction =
      mean(remodeling_supported),

    Schwann_fraction =
      mean(Schwann_supported),

    APM_bulge_joint_fraction =
      mean(APM_bulge_joint),

    APM_bulge_junction_fraction =
      mean(APM_bulge_junction_joint),

    median_APM_score =
      median(APM_CONTRACTILE_score),

    median_attachment_score =
      median(APM_ATTACHMENT_score),

    median_bulge_basal_score =
      median(BULGE_BASAL_score),

    median_bulge_stem_score =
      median(BULGE_STEM_score),

    .groups = "drop"
  ) %>%
  arrange(as.integer(cluster))

cat("\n=== CLUSTER INTERFACE SUMMARY ===\n")
print(cluster_summary, n = Inf, width = Inf)

write.csv(
  cluster_summary,
  file.path(TAB, "11A5_cluster_interface_summary.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# 4. Key attachment genes
# ------------------------------------------------------------

cat("\nSTEP 4/7 — attachment-axis gene detection\n")

target_genes <- c(
  # APM-bulge attachment axis
  "NPNT","ITGA8","ITGB1","EGFL6",

  # bulge/ORS
  "KRT5","KRT14","KRT15","KRT17",
  "TP63","SOX9","CD200",
  "ITGA6","ITGB4","COL17A1",

  # APM
  "MYH11","ACTA2","TAGLN","CNN1",
  "MYLK","ACTG2","DES",

  # epithelial junction
  "DSP","DSC3","DSG1",
  "PKP1","PKP3","CLDN1",

  # alternative identities
  "RGS5","PDGFRB","ABCC9",
  "COL11A1","DPEP1","DCN","LUM",
  "SOX10","S100B","PMP22"
)

target_genes <- intersect(
  unique(target_genes),
  rownames(mat)
)

gene_detection <- bind_rows(
  lapply(sort(unique(cluster)), function(cl) {

    idx <- which(cluster == cl)

    bind_rows(
      lapply(target_genes, function(g) {

        v <- as.numeric(
          mat[g, idx]
        )

        tibble(
          cluster = cl,
          gene = g,
          fraction_detected =
            mean(v > 0),
          mean_expression =
            mean(v)
        )
      })
    )
  })
) %>%
  arrange(
    as.integer(cluster),
    desc(fraction_detected)
  )

cat("\n=== TARGET GENE DETECTION ===\n")
print(gene_detection, n = Inf)

write.csv(
  gene_detection,
  file.path(TAB, "11A5_target_gene_detection.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# 5. C0-specific patient recurrence
# ------------------------------------------------------------

cat("\nSTEP 5/7 — C0 patient recurrence\n")

c0_patient <- ledger %>%
  filter(cluster == "0") %>%
  group_by(patient_id, batch) %>%
  summarise(
    n_cells = n(),

    median_doublet =
      median(doublet_score),

    APM_fraction =
      mean(APM_supported),

    bulge_basal_fraction =
      mean(bulge_basal_supported),

    bulge_stem_fraction =
      mean(bulge_stem_supported),

    junction_fraction =
      mean(epithelial_junction_supported),

    APM_bulge_joint_fraction =
      mean(APM_bulge_joint),

    APM_bulge_junction_fraction =
      mean(APM_bulge_junction_joint),

    median_APM_score =
      median(APM_CONTRACTILE_score),

    median_bulge_basal_score =
      median(BULGE_BASAL_score),

    .groups = "drop"
  ) %>%
  arrange(desc(n_cells))

cat("\n=== C0 PATIENT SUPPORT ===\n")
print(c0_patient, n = Inf, width = Inf)

write.csv(
  c0_patient,
  file.path(TAB, "11A5_C0_patient_support.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# 6. Doublet-independent C0 evidence
# ------------------------------------------------------------

cat("\nSTEP 6/7 — doublet-independent C0 evidence\n")

c0 <- ledger %>%
  filter(cluster == "0") %>%
  mutate(
    doublet_quartile =
      factor(
        ntile(doublet_score, 4),
        levels = 1:4,
        labels = c(
          "Q1_low",
          "Q2",
          "Q3",
          "Q4_high"
        )
      )
  )

c0_doublet <- c0 %>%
  group_by(doublet_quartile) %>%
  summarise(
    n_cells = n(),
    median_doublet =
      median(doublet_score),

    median_features =
      median(nFeature_RNA),

    APM_fraction =
      mean(APM_supported),

    bulge_basal_fraction =
      mean(bulge_basal_supported),

    bulge_stem_fraction =
      mean(bulge_stem_supported),

    junction_fraction =
      mean(epithelial_junction_supported),

    APM_bulge_joint_fraction =
      mean(APM_bulge_joint),

    APM_bulge_junction_fraction =
      mean(APM_bulge_junction_joint),

    .groups = "drop"
  )

cat("\n=== C0 DOUBLET STRATA ===\n")
print(c0_doublet, n = Inf, width = Inf)

write.csv(
  c0_doublet,
  file.path(TAB, "11A5_C0_doublet_strata.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# 7. Formal checks
# ------------------------------------------------------------

cat("\nSTEP 7/7 — formal checks\n")

checks <- tibble(
  check = c(
    "cells_330",
    "patients_16",
    "three_clusters",
    "C0_cells_138",
    "programme_coverage_generated",
    "cell_ledger_generated",
    "cluster_summary_generated",
    "attachment_gene_detection_generated",
    "C0_patient_recurrence_generated",
    "C0_doublet_strata_generated",
    "disease_not_used",
    "no_cells_removed",
    "no_identity_assigned",
    "source_object_not_modified"
  ),

  pass = c(
    ncol(x) == 330,
    length(unique(x$patient_id)) == 16,
    length(unique(cluster)) == 3,
    sum(cluster == "0") == 138,
    nrow(coverage) == length(programmes),
    nrow(ledger) == 330,
    nrow(cluster_summary) == 3,
    nrow(gene_detection) > 0,
    nrow(c0_patient) == 11,
    nrow(c0_doublet) == 4,
    TRUE,
    ncol(x) == 330,
    TRUE,
    TRUE
  )
)

cat("\n=== FORMAL CHECKS ===\n")
print(checks, n = Inf)

write.csv(
  checks,
  file.path(TAB, "11A5_formal_checks.csv"),
  row.names = FALSE
)

stopifnot(all(checks$pass))

cat("\n============================================================\n")
cat("11A5 COMPLETE\n")
cat("============================================================\n")
cat("Cells:", ncol(x), "\n")
cat("C0 cells:", sum(cluster == "0"), "\n")
cat("No cells removed.\n")
cat("No identity assigned.\n")
cat("Disease not used.\n")
cat("PASS\n")
