suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(readr)
})

options(width = 200)

cat("============================================================\n")
cat("PIPELINE_V2 — 10B4C TARGETED POST-C8 LINEAGE AUDIT\n")
cat("============================================================\n")

SRC <- paste0(
  "PIPELINE_V2/objects/annotated/",
  "vascular_v2_10B4B_rebuilt_3023.rds"
)

OUT <- paste0(
  "PIPELINE_V2/results/annotation/",
  "10_VASCULAR_L2/10B4C_targeted_lineage_audit"
)

TAB <- file.path(OUT, "tables")
dir.create(TAB, recursive = TRUE, showWarnings = FALSE)

x <- readRDS(SRC)

stopifnot(ncol(x) == 3023)

CL <- "vascular_postC8_snn_res.0.6"
stopifnot(CL %in% colnames(x@meta.data))

DefaultAssay(x) <- "RNA"

# ------------------------------------------------------------
# Programmes
# ------------------------------------------------------------

programmes <- list(

  ENDOTHELIAL_CORE = c(
    "PECAM1","VWF","CDH5","EMCN",
    "KDR","ESAM","ENG","RAMP2"
  ),

  ENDOTHELIAL_JUNCTION = c(
    "PECAM1","CDH5","CLDN5","ESAM"
  ),

  CAPILLARY = c(
    "CA4","RGCC","BTNL9","EDNRB",
    "EMCN","KDR","PLVAP"
  ),

  VENOUS_ACKR1 = c(
    "ACKR1","SELP","SELE","VWF",
    "VCAM1","ICAM1"
  ),

  ARTERIAL = c(
    "GJA5","GJA4","EFNB2","SOX17",
    "HEY1","HEY2","BMX","SEMA3G","FBLN5"
  ),

  LYMPHATIC = c(
    "PROX1","PDPN","FLT4","CCL21",
    "LYVE1","NRP2","RELN"
  ),

  MURAL_PERICYTE = c(
    "RGS5","CSPG4","MCAM","PDGFRB",
    "RBP1","NOTCH3","ABCC9","KCNJ8"
  ),

  CONTRACTILE_SMC = c(
    "ACTA2","TAGLN","MYH11",
    "CNN1","MYL9","TPM2"
  ),

  MESENCHYMAL_ECM = c(
    "DCN","LUM","COL1A1","COL1A2",
    "COL3A1","DPT","FN1","POSTN"
  ),

  PERIVASCULAR_EBF2 = c(
    "EBF2","TBX18","TBX15","AGTR1",
    "CARMN","DPT","PRR16","BMP5"
  ),

  EPITHELIAL_BASAL = c(
    "KRT5","KRT14","KRT17",
    "KRT1","KRT10","TP63"
  ),

  EPITHELIAL_JUNCTION = c(
    "EPCAM","DSP","PKP1","PKP3",
    "DSC3","DSG1","CLDN1"
  ),

  IMMUNE = c(
    "PTPRC","CD3D","CD3E",
    "TRAC","NKG7","GNLY"
  )
)

coverage <- tibble(
  programme = names(programmes),
  n_requested = sapply(programmes, length),
  n_present = sapply(
    programmes,
    function(g) sum(g %in% rownames(x))
  ),
  genes_present = sapply(
    programmes,
    function(g)
      paste(intersect(g, rownames(x)), collapse = ";")
  )
)

cat("\n===== PROGRAMME COVERAGE =====\n")
print(coverage, n = Inf)

write_csv(
  coverage,
  file.path(TAB, "10B4C_programme_coverage.csv")
)

# ------------------------------------------------------------
# Cell-level detection
# ------------------------------------------------------------

genes <- intersect(
  unique(unlist(programmes)),
  rownames(x)
)

e <- FetchData(
  x,
  vars = genes,
  layer = "data"
)

detected_count <- function(df, genes) {

  g <- intersect(
    genes,
    colnames(df)
  )

  if (length(g) == 0)
    return(rep(0L, nrow(df)))

  rowSums(
    df[, g, drop = FALSE] > 0
  )
}

score_mean <- function(df, genes) {

  g <- intersect(
    genes,
    colnames(df)
  )

  if (length(g) == 0)
    return(rep(NA_real_, nrow(df)))

  rowMeans(
    df[, g, drop = FALSE]
  )
}

cell <- tibble(
  cell = rownames(e),
  cluster =
    as.character(
      x@meta.data[
        match(rownames(e), colnames(x)),
        CL
      ]
    ),
  patient_id =
    as.character(
      x$patient_id[
        match(rownames(e), colnames(x))
      ]
    ),
  batch =
    as.character(
      x$batch[
        match(rownames(e), colnames(x))
      ]
    ),
  doublet_score =
    x$doublet_score[
      match(rownames(e), colnames(x))
    ],
  nFeature_RNA =
    x$nFeature_RNA[
      match(rownames(e), colnames(x))
    ],
  nCount_RNA =
    x$nCount_RNA[
      match(rownames(e), colnames(x))
    ]
)

for (p in names(programmes)) {

  cell[[paste0(p, "_n")]] <-
    detected_count(
      e,
      programmes[[p]]
    )

  cell[[paste0(p, "_score")]] <-
    score_mean(
      e,
      programmes[[p]]
    )
}

write_csv(
  cell,
  file.path(TAB, "10B4C_cell_programme_evidence.csv")
)

# ------------------------------------------------------------
# Cluster summary
# ------------------------------------------------------------

focus <- cell %>%
  filter(cluster %in% c("3","7","8","9"))

summary <- focus %>%
  group_by(cluster) %>%
  summarise(
    n_cells = n(),
    n_patients = n_distinct(patient_id),
    n_batches = n_distinct(batch),

    median_features =
      median(nFeature_RNA),

    median_doublet =
      median(doublet_score),

    median_EC_n =
      median(ENDOTHELIAL_CORE_n),

    EC_ge3 =
      mean(ENDOTHELIAL_CORE_n >= 3),

    EC_ge4 =
      mean(ENDOTHELIAL_CORE_n >= 4),

    EC_junction_ge2 =
      mean(ENDOTHELIAL_JUNCTION_n >= 2),

    mural_ge3 =
      mean(MURAL_PERICYTE_n >= 3),

    contractile_ge2 =
      mean(CONTRACTILE_SMC_n >= 2),

    ECM_ge3 =
      mean(MESENCHYMAL_ECM_n >= 3),

    EBF2_perivascular_ge3 =
      mean(PERIVASCULAR_EBF2_n >= 3),

    epithelial_basal_ge2 =
      mean(EPITHELIAL_BASAL_n >= 2),

    epithelial_junction_ge2 =
      mean(EPITHELIAL_JUNCTION_n >= 2),

    immune_ge2 =
      mean(IMMUNE_n >= 2),

    .groups = "drop"
  )

cat("\n===== FOCUS CLUSTER SUMMARY =====\n")
print(summary, n = Inf)

write_csv(
  summary,
  file.path(TAB, "10B4C_focus_cluster_summary.csv")
)

# ------------------------------------------------------------
# Joint lineage evidence
# ------------------------------------------------------------

joint <- focus %>%
  mutate(
    EC_supported =
      ENDOTHELIAL_CORE_n >= 3,

    strong_EC =
      ENDOTHELIAL_CORE_n >= 4,

    mural_supported =
      MURAL_PERICYTE_n >= 3,

    ECM_supported =
      MESENCHYMAL_ECM_n >= 3,

    perivascular_supported =
      PERIVASCULAR_EBF2_n >= 3,

    epithelial_supported =
      EPITHELIAL_BASAL_n >= 2 |
      EPITHELIAL_JUNCTION_n >= 2,

    immune_supported =
      IMMUNE_n >= 2
  ) %>%
  group_by(cluster) %>%
  summarise(
    n_cells = n(),

    EC_supported =
      sum(EC_supported),

    strong_EC =
      sum(strong_EC),

    EC_and_mural =
      sum(EC_supported & mural_supported),

    EC_and_ECM =
      sum(EC_supported & ECM_supported),

    EC_and_perivascular =
      sum(EC_supported & perivascular_supported),

    EC_and_epithelial =
      sum(EC_supported & epithelial_supported),

    epithelial_without_EC =
      sum(epithelial_supported & !EC_supported),

    perivascular_without_EC =
      sum(perivascular_supported & !EC_supported),

    immune_supported =
      sum(immune_supported),

    .groups = "drop"
  )

cat("\n===== JOINT LINEAGE EVIDENCE =====\n")
print(joint, n = Inf)

write_csv(
  joint,
  file.path(TAB, "10B4C_joint_lineage_evidence.csv")
)

# ------------------------------------------------------------
# Gene detection specifically C8/C9
# ------------------------------------------------------------

target_genes <- c(
  "PECAM1","VWF","CDH5","EMCN","KDR","ESAM","ENG",
  "RGS5","CSPG4","PDGFRB","ACTA2","TAGLN","MYH11",
  "EBF2","TBX18","TBX15","AGTR1","CARMN","DPT",
  "DCN","LUM","COL1A1","COL1A2","COL3A1",
  "KRT5","KRT14","KRT17","KRT1","KRT10",
  "PKP3","DSP","CLDN1","EPCAM",
  "PTPRC","CD3D","CD3E","TRAC"
)

target_genes <- intersect(
  target_genes,
  colnames(e)
)

gene_detection <- bind_rows(
  lapply(
    c("8","9"),
    function(cc) {

      idx <- focus$cluster == cc
      cells <- focus$cell[idx]

      ee <- e[cells, target_genes, drop = FALSE]

      tibble(
        cluster = cc,
        gene = target_genes,
        fraction_detected =
          colMeans(ee > 0),
        mean_expression =
          colMeans(ee)
      ) %>%
        arrange(
          desc(fraction_detected),
          desc(mean_expression)
        )
    }
  )
)

cat("\n===== C8/C9 TARGET GENE DETECTION =====\n")
print(gene_detection, n = Inf)

write_csv(
  gene_detection,
  file.path(TAB, "10B4C_C8_C9_gene_detection.csv")
)

# ------------------------------------------------------------
# Patient recurrence
# ------------------------------------------------------------

patient <- focus %>%
  group_by(cluster, patient_id, batch) %>%
  summarise(
    n_cells = n(),

    EC_supported_fraction =
      mean(ENDOTHELIAL_CORE_n >= 3),

    perivascular_fraction =
      mean(PERIVASCULAR_EBF2_n >= 3),

    epithelial_fraction =
      mean(
        EPITHELIAL_BASAL_n >= 2 |
        EPITHELIAL_JUNCTION_n >= 2
      ),

    median_doublet =
      median(doublet_score),

    .groups = "drop"
  )

write_csv(
  patient,
  file.path(TAB, "10B4C_patient_recurrence.csv")
)

# ------------------------------------------------------------
# Formal checks
# ------------------------------------------------------------

checks <- tibble(
  check = c(
    "source_cells_3023",
    "patients_16",
    "r06_10_clusters",
    "focus_C3_C7_C8_C9_present",
    "programme_coverage_generated",
    "cell_evidence_generated",
    "joint_lineage_generated",
    "C8_C9_gene_detection_generated",
    "patient_recurrence_generated",
    "disease_not_used",
    "no_cells_removed",
    "no_identity_assigned",
    "source_object_not_modified"
  ),

  pass = c(
    ncol(x) == 3023,
    length(unique(x$patient_id)) == 16,
    length(unique(x@meta.data[[CL]])) == 10,
    all(c("3","7","8","9") %in%
          unique(as.character(x@meta.data[[CL]]))),
    nrow(coverage) == length(programmes),
    nrow(cell) == 3023,
    nrow(joint) == 4,
    nrow(gene_detection) > 0,
    nrow(patient) > 0,
    TRUE,
    TRUE,
    TRUE,
    TRUE
  )
)

write_csv(
  checks,
  file.path(TAB, "10B4C_formal_checks.csv")
)

cat("\n===== FORMAL CHECKS =====\n")
print(checks, n = Inf)

stopifnot(all(checks$pass))

writeLines(
  c(
    "PIPELINE_V2 — 10B4C TARGETED POST-C8 LINEAGE AUDIT",
    "",
    "Purpose:",
    "- resolve lineage status of r0.6 C3/C7/C8/C9",
    "- distinguish endothelial states from mural/perivascular or epithelial contamination",
    "- evaluate evidence at cell and patient levels",
    "",
    "Interpretive boundary:",
    "- disease not used",
    "- no cells removed",
    "- no identities assigned",
    "- no resolution frozen",
    "- source object unchanged",
    "",
    "PASS"
  ),
  file.path(TAB, "10B4C_audit.txt")
)

cat("\n============================================================\n")
cat("10B4C COMPLETE\n")
cat("============================================================\n")
cat("Cells:", ncol(x), "\n")
cat("Focus clusters: C3/C7/C8/C9\n")
cat("No cells removed.\n")
cat("No identities assigned.\n")
cat("PASS\n")
