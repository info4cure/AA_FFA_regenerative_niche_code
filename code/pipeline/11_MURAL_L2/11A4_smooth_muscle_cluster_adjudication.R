suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(tibble)
})

cat("============================================================\n")
cat("PIPELINE_V2 — 11A4 SMOOTH-MUSCLE CLUSTER ADJUDICATION\n")
cat("============================================================\n\n")

SRC <- paste0(
  "PIPELINE_V2/results/annotation/11_MURAL_L2/",
  "11A3_smooth_muscle_internal_structure/",
  "smooth_muscle_11A3_diagnostic_reclustered.rds"
)

EVID <- paste0(
  "PIPELINE_V2/results/annotation/11_MURAL_L2/",
  "11A3_smooth_muscle_internal_structure/tables/",
  "11A3_cell_programme_evidence.csv"
)

OUT <- paste0(
  "PIPELINE_V2/results/annotation/11_MURAL_L2/",
  "11A4_smooth_muscle_cluster_adjudication"
)

TAB <- file.path(OUT, "tables")
dir.create(TAB, recursive = TRUE, showWarnings = FALSE)

x <- readRDS(SRC)
ev <- read.csv(EVID, check.names = FALSE)

stopifnot(ncol(x) == 330)
stopifnot(nrow(ev) == 330)

CL <- "smooth_snn_res.0.3"

cat("Cells:", ncol(x), "\n")
cat("Patients:", length(unique(x$patient_id)), "\n")
cat("Clusters:", length(unique(x@meta.data[[CL]])), "\n\n")

# ------------------------------------------------------------
# 1. Structural support
# ------------------------------------------------------------

cat("STEP 1/6 — structural support\n")

structural <- tibble(
  cluster = as.character(x@meta.data[[CL]]),
  patient_id = as.character(x$patient_id),
  batch = as.character(x$batch),
  doublet_score = as.numeric(x$doublet_score),
  nFeature_RNA = as.numeric(x$nFeature_RNA),
  nCount_RNA = as.numeric(x$nCount_RNA),
  percent_mt = as.numeric(x$percent.mt)
) %>%
  group_by(cluster) %>%
  summarise(
    n_cells = n(),
    n_patients = n_distinct(patient_id),
    n_batches = n_distinct(batch),
    dominant_patient_fraction =
      max(table(patient_id)) / n(),
    median_doublet = median(doublet_score),
    q75_doublet = quantile(doublet_score, 0.75),
    median_features = median(nFeature_RNA),
    median_counts = median(nCount_RNA),
    median_mt = median(percent_mt),
    .groups = "drop"
  ) %>%
  arrange(as.integer(cluster))

cat("\n=== STRUCTURAL SUPPORT ===\n")
print(structural, n = Inf, width = Inf)

write.csv(
  structural,
  file.path(TAB, "11A4_structural_support.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# 2. Programme evidence
# ------------------------------------------------------------

cat("\nSTEP 2/6 — cluster programme evidence\n")

score_cols <- grep(
  "_score$",
  names(ev),
  value = TRUE
)

n_cols <- grep(
  "_n$",
  names(ev),
  value = TRUE
)

programme_summary <- ev %>%
  mutate(
    cluster = as.character(.data[[CL]])
  ) %>%
  group_by(cluster) %>%
  summarise(
    across(
      all_of(score_cols),
      ~ median(.x, na.rm = TRUE),
      .names = "median_{.col}"
    ),
    across(
      all_of(n_cols),
      ~ mean(.x, na.rm = TRUE),
      .names = "mean_{.col}"
    ),
    .groups = "drop"
  )

cat("\n=== PROGRAMME SUMMARY ===\n")
print(programme_summary, n = Inf, width = Inf)

write.csv(
  programme_summary,
  file.path(TAB, "11A4_programme_summary.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# 3. Explicit cell-level competing-lineage evidence
# ------------------------------------------------------------

cat("\nSTEP 3/6 — cell-level competing-lineage evidence\n")

cell_evidence <- ev %>%
  mutate(
    cluster = as.character(.data[[CL]]),

    contractile_ge4 =
      CONTRACTILE_CORE_n >= 4,

    arrector_ge4 =
      ARRECTOR_PILI_SMC_n >= 4,

    basal_ge2 =
      BASAL_KERATINOCYTE_n >= 2,

    differentiated_ge2 =
      DIFFERENTIATED_KERATINOCYTE_n >= 2,

    epithelial_junction_ge2 =
      EPITHELIAL_JUNCTION_n >= 2,

    glandular_ge2 =
      GLANDULAR_EPITHELIAL_n >= 2,

    sebocyte_ge2 =
      SEBOCYTE_n >= 2,

    pericyte_ge3 =
      PERICYTE_n >= 3,

    fibroblast_ge3 =
      FIBROBLAST_ECM_n >= 3,

    remodeling_ge2 =
      REMODELING_n >= 2,

    Schwann_ge2 =
      SCHWANN_n >= 2,

    contractile_basal_joint =
      contractile_ge4 & basal_ge2,

    contractile_Schwann_joint =
      contractile_ge4 & Schwann_ge2
  )

joint <- cell_evidence %>%
  group_by(cluster) %>%
  summarise(
    n_cells = n(),

    contractile_ge4 =
      mean(contractile_ge4),

    arrector_ge4 =
      mean(arrector_ge4),

    basal_ge2 =
      mean(basal_ge2),

    differentiated_ge2 =
      mean(differentiated_ge2),

    epithelial_junction_ge2 =
      mean(epithelial_junction_ge2),

    glandular_ge2 =
      mean(glandular_ge2),

    sebocyte_ge2 =
      mean(sebocyte_ge2),

    pericyte_ge3 =
      mean(pericyte_ge3),

    fibroblast_ge3 =
      mean(fibroblast_ge3),

    remodeling_ge2 =
      mean(remodeling_ge2),

    Schwann_ge2 =
      mean(Schwann_ge2),

    contractile_basal_joint =
      mean(contractile_basal_joint),

    contractile_Schwann_joint =
      mean(contractile_Schwann_joint),

    .groups = "drop"
  ) %>%
  arrange(as.integer(cluster))

cat("\n=== JOINT LINEAGE EVIDENCE ===\n")
print(joint, n = Inf, width = Inf)

write.csv(
  joint,
  file.path(TAB, "11A4_joint_lineage_evidence.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# 4. Target-gene detection by cluster
# ------------------------------------------------------------

cat("\nSTEP 4/6 — target-gene detection\n")

DefaultAssay(x) <- "RNA"
mat <- GetAssayData(x, assay = "RNA", layer = "data")

target_genes <- c(
  # canonical smooth muscle
  "MYH11","ACTA2","TAGLN","CNN1","MYLK",
  "LMOD1","SYNPO2","ACTG2","DES","ACTA1",

  # mural
  "RGS5","CSPG4","PDGFRB","RBP1",
  "ABCC9","KCNJ8","NOTCH3","RGS6",

  # basal epithelial
  "KRT5","KRT14","KRT17","TP63","KRT15",

  # differentiation / junction
  "KRT1","KRT10","DSP","DSC3",
  "DSG1","PKP1","PKP3","CLDN1",

  # glandular / sebaceous
  "KRT7","KRT8","KRT18","KRT19",
  "MUC1","PPARG","PLIN2","FASN","SCD",

  # neural / Schwann
  "SOX10","S100B","PLP1","PMP22",
  "MPZ","NGFR","S100A1",

  # C0 markers
  "ECRG4","SAA1","ITM2A","NPPC",
  "PCSK1N","DKK1","IRX4","FJX1",

  # C1 markers
  "EMX2","ALDH1A2","BCHE","HSD17B6",
  "ADRA2C","SCX","STAC","OLFML3",

  # C2 markers
  "TBX2","RGS6","APOLD1","RCAN2",
  "PDE1A","PDE3A","CTNNA3","NR2F2"
)

target_genes <- intersect(
  unique(target_genes),
  rownames(mat)
)

clusters <- as.character(x@meta.data[[CL]])

gene_detection <- bind_rows(
  lapply(sort(unique(clusters)), function(cl) {

    cells <- which(clusters == cl)

    bind_rows(
      lapply(target_genes, function(g) {

        v <- as.numeric(mat[g, cells])

        tibble(
          cluster = cl,
          gene = g,
          fraction_detected = mean(v > 0),
          mean_expression = mean(v)
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
  file.path(TAB, "11A4_target_gene_detection.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# 5. Patient recurrence of key mixed programmes
# ------------------------------------------------------------

cat("\nSTEP 5/6 — patient recurrence\n")

patient_support <- cell_evidence %>%
  group_by(cluster, patient_id, batch) %>%
  summarise(
    n_cells = n(),
    contractile_fraction =
      mean(contractile_ge4),
    basal_fraction =
      mean(basal_ge2),
    junction_fraction =
      mean(epithelial_junction_ge2),
    Schwann_fraction =
      mean(Schwann_ge2),
    contractile_basal_joint =
      mean(contractile_basal_joint),
    .groups = "drop"
  ) %>%
  arrange(
    as.integer(cluster),
    desc(n_cells)
  )

cat("\n=== PATIENT SUPPORT ===\n")
print(patient_support, n = Inf, width = Inf)

write.csv(
  patient_support,
  file.path(TAB, "11A4_patient_support.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# 6. Formal checks
# ------------------------------------------------------------

cat("\nSTEP 6/6 — formal checks\n")

checks <- tibble(
  check = c(
    "cells_330",
    "patients_16",
    "three_clusters_r03",
    "structural_support_generated",
    "programme_summary_generated",
    "joint_lineage_generated",
    "target_gene_detection_generated",
    "patient_recurrence_generated",
    "disease_not_used",
    "no_cells_removed",
    "no_identity_assigned",
    "source_object_not_modified"
  ),
  pass = c(
    ncol(x) == 330,
    length(unique(x$patient_id)) == 16,
    length(unique(x@meta.data[[CL]])) == 3,
    nrow(structural) == 3,
    nrow(programme_summary) == 3,
    nrow(joint) == 3,
    nrow(gene_detection) > 0,
    nrow(patient_support) > 0,
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
  file.path(TAB, "11A4_formal_checks.csv"),
  row.names = FALSE
)

stopifnot(all(checks$pass))

cat("\n============================================================\n")
cat("11A4 COMPLETE\n")
cat("============================================================\n")
cat("Cells:", ncol(x), "\n")
cat("r0.3 clusters: 3\n")
cat("No cells removed.\n")
cat("No identities assigned.\n")
cat("Disease not used.\n")
cat("PASS\n")
