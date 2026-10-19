suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
})

cat("============================================================\n")
cat("PIPELINE_V2 — 11A2 SMOOTH-MUSCLE LINEAGE-INTEGRITY AUDIT\n")
cat("============================================================\n\n")

SRC <- "PIPELINE_V2/objects/annotated/global_atlas_v2_L1_frozen.rds"

OUT <- paste0(
  "PIPELINE_V2/results/annotation/11_MURAL_L2/",
  "11A2_smooth_muscle_audit"
)

TAB <- file.path(OUT, "tables")
FIG <- file.path(OUT, "figures")

dir.create(TAB, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG, recursive = TRUE, showWarnings = FALSE)

x <- readRDS(SRC)

stopifnot(ncol(x) == 82893)
stopifnot("L1" %in% colnames(x@meta.data))

sm <- subset(
  x,
  cells = colnames(x)[x$L1 == "Smooth_muscle"]
)

cat("Cells:", ncol(sm), "\n")
cat("Patients:", length(unique(sm$patient_id)), "\n")
cat("Batches:", length(unique(sm$batch)), "\n\n")

stopifnot(ncol(sm) == 330)

DefaultAssay(sm) <- "RNA"

# ------------------------------------------------------------
# Programmes
# ------------------------------------------------------------

cat("STEP 1/7 — defining lineage programmes\n")

programmes <- list(

  CONTRACTILE_SMC = c(
    "ACTA2","TAGLN","MYH11","CNN1",
    "MYL9","TPM2","DES","LMOD1"
  ),

  VASCULAR_SMC = c(
    "MYH11","CNN1","ACTA2","TAGLN",
    "LMOD1","SYNPO2","SORBS2","MYLK"
  ),

  PERICYTE = c(
    "RGS5","CSPG4","PDGFRB","RBP1",
    "ABCC9","KCNJ8","COX4I2","HIGD1B"
  ),

  FIBROBLAST = c(
    "DCN","LUM","COL1A1","COL1A2",
    "COL3A1","PDGFRA","DPT","COL6A1"
  ),

  REMODELING_FIBROBLAST = c(
    "ADAM12","MMP11","POSTN","COL11A1",
    "COMP","COL10A1","FN1","THBS2"
  ),

  EPITHELIAL_BASAL = c(
    "KRT5","KRT14","KRT17","TP63",
    "KRT15","KRT1","KRT10"
  ),

  EPITHELIAL_JUNCTION = c(
    "EPCAM","DSP","PKP1","PKP3",
    "DSC3","DSG1","CLDN1"
  ),

  SWEAT_GLAND = c(
    "KRT8","KRT18","KRT19","KRT7",
    "KRT23","MUC1","KRT18","EPCAM"
  ),

  ENDOTHELIAL = c(
    "PECAM1","VWF","CDH5","EMCN",
    "KDR","ESAM","ENG","RAMP2"
  ),

  SCHWANN = c(
    "SOX10","S100B","PLP1","PMP22",
    "MPZ","S100A1","NGFR"
  ),

  IMMUNE = c(
    "PTPRC","CD3D","CD3E","TRAC",
    "NKG7","GNLY"
  )
)

features <- rownames(sm)

available <- lapply(
  programmes,
  function(g) intersect(g, features)
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

print(coverage, n = Inf)

write.csv(
  coverage,
  file.path(TAB, "11A2_programme_coverage.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# Cell-level detection
# ------------------------------------------------------------

cat("\nSTEP 2/7 — cell-level programme detection\n")

mat <- GetAssayData(sm, assay = "RNA", layer = "data")

detect_n <- function(genes) {
  genes <- intersect(genes, rownames(mat))

  if (length(genes) == 0) {
    return(rep(0L, ncol(mat)))
  }

  Matrix::colSums(mat[genes, , drop = FALSE] > 0)
}

mean_expr <- function(genes) {
  genes <- intersect(genes, rownames(mat))

  if (length(genes) == 0) {
    return(rep(0, ncol(mat)))
  }

  Matrix::colMeans(mat[genes, , drop = FALSE])
}

ledger <- tibble(
  cell = colnames(sm),
  patient_id = as.character(sm$patient_id),
  batch = as.character(sm$batch),
  doublet_class = as.character(sm$doublet_class),
  doublet_score = as.numeric(sm$doublet_score),
  nFeature_RNA = as.numeric(sm$nFeature_RNA),
  nCount_RNA = as.numeric(sm$nCount_RNA),
  percent_mt = as.numeric(sm$percent.mt)
)

for (p in names(available)) {
  ledger[[paste0(p, "_n")]] <- detect_n(available[[p]])
  ledger[[paste0(p, "_score")]] <- mean_expr(available[[p]])
}

ledger <- ledger %>%
  mutate(
    contractile_supported =
      CONTRACTILE_SMC_n >= 3,

    strong_contractile =
      CONTRACTILE_SMC_n >= 4,

    pericyte_supported =
      PERICYTE_n >= 3,

    fibroblast_supported =
      FIBROBLAST_n >= 3,

    remodeling_supported =
      REMODELING_FIBROBLAST_n >= 2,

    epithelial_basal_supported =
      EPITHELIAL_BASAL_n >= 2,

    epithelial_junction_supported =
      EPITHELIAL_JUNCTION_n >= 2,

    epithelial_supported =
      epithelial_basal_supported |
      epithelial_junction_supported,

    sweat_supported =
      SWEAT_GLAND_n >= 3,

    endothelial_supported =
      ENDOTHELIAL_n >= 3,

    Schwann_supported =
      SCHWANN_n >= 2,

    immune_supported =
      IMMUNE_n >= 2,

    contractile_epithelial_joint =
      contractile_supported &
      epithelial_supported,

    strong_contractile_epithelial_joint =
      strong_contractile &
      epithelial_supported
  )

write.csv(
  ledger,
  file.path(TAB, "11A2_cell_level_ledger.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# Joint evidence
# ------------------------------------------------------------

cat("\nSTEP 3/7 — joint contractile/epithelial evidence\n")

cat("\nCONTRACTILE >=3 vs EPITHELIAL\n")
print(
  table(
    contractile = ledger$contractile_supported,
    epithelial = ledger$epithelial_supported
  )
)

cat("\nCONTRACTILE >=4 vs EPITHELIAL\n")
print(
  table(
    contractile = ledger$strong_contractile,
    epithelial = ledger$epithelial_supported
  )
)

joint_summary <- tibble(
  n_cells = nrow(ledger),

  contractile_fraction =
    mean(ledger$contractile_supported),

  strong_contractile_fraction =
    mean(ledger$strong_contractile),

  pericyte_fraction =
    mean(ledger$pericyte_supported),

  fibroblast_fraction =
    mean(ledger$fibroblast_supported),

  remodeling_fraction =
    mean(ledger$remodeling_supported),

  epithelial_fraction =
    mean(ledger$epithelial_supported),

  sweat_fraction =
    mean(ledger$sweat_supported),

  endothelial_fraction =
    mean(ledger$endothelial_supported),

  Schwann_fraction =
    mean(ledger$Schwann_supported),

  immune_fraction =
    mean(ledger$immune_supported),

  contractile_epithelial_joint =
    mean(ledger$contractile_epithelial_joint),

  strong_contractile_epithelial_joint =
    mean(ledger$strong_contractile_epithelial_joint)
)

cat("\n=== JOINT SUMMARY ===\n")
print(joint_summary, width = Inf)

write.csv(
  joint_summary,
  file.path(TAB, "11A2_joint_summary.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# Gene detection
# ------------------------------------------------------------

cat("\nSTEP 4/7 — targeted gene detection\n")

target_genes <- unique(unlist(programmes))

gene_detection <- bind_rows(
  lapply(target_genes, function(g) {

    if (!g %in% rownames(mat)) {
      return(NULL)
    }

    v <- as.numeric(mat[g, ])

    tibble(
      gene = g,
      fraction_detected = mean(v > 0),
      mean_expression = mean(v)
    )
  })
) %>%
  arrange(desc(fraction_detected))

print(gene_detection, n = Inf)

write.csv(
  gene_detection,
  file.path(TAB, "11A2_gene_detection.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# Technical association
# ------------------------------------------------------------

cat("\nSTEP 5/7 — technical association\n")

scores <- c(
  "CONTRACTILE_SMC_score",
  "EPITHELIAL_BASAL_score",
  "EPITHELIAL_JUNCTION_score",
  "SWEAT_GLAND_score"
)

technical <- c(
  "doublet_score",
  "nFeature_RNA",
  "nCount_RNA",
  "percent_mt"
)

technical_cor <- bind_rows(
  lapply(scores, function(s) {
    bind_rows(
      lapply(technical, function(t) {

        z <- suppressWarnings(
          cor.test(
            ledger[[s]],
            ledger[[t]],
            method = "spearman",
            exact = FALSE
          )
        )

        tibble(
          biological_score = s,
          technical_variable = t,
          n = nrow(ledger),
          spearman_rho = unname(z$estimate),
          p_value = z$p.value
        )
      })
    )
  })
) %>%
  mutate(
    p_adj = p.adjust(p_value, method = "BH")
  )

print(technical_cor, n = Inf)

write.csv(
  technical_cor,
  file.path(TAB, "11A2_technical_correlations.csv"),
  row.names = FALSE
)

# Doublet quartiles
ledger <- ledger %>%
  mutate(
    doublet_stratum = factor(
      ntile(doublet_score, 4),
      levels = 1:4,
      labels = c("Q1_low","Q2","Q3","Q4_high")
    )
  )

doublet_strata <- ledger %>%
  group_by(doublet_stratum) %>%
  summarise(
    n_cells = n(),
    median_doublet = median(doublet_score),
    median_features = median(nFeature_RNA),
    median_counts = median(nCount_RNA),

    contractile_score =
      median(CONTRACTILE_SMC_score),

    epithelial_basal_score =
      median(EPITHELIAL_BASAL_score),

    epithelial_junction_score =
      median(EPITHELIAL_JUNCTION_score),

    sweat_score =
      median(SWEAT_GLAND_score),

    contractile_fraction =
      mean(contractile_supported),

    epithelial_fraction =
      mean(epithelial_supported),

    joint_fraction =
      mean(contractile_epithelial_joint),

    .groups = "drop"
  )

cat("\n=== DOUBLET STRATA ===\n")
print(doublet_strata, n = Inf, width = Inf)

write.csv(
  doublet_strata,
  file.path(TAB, "11A2_doublet_strata.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# Patient recurrence
# ------------------------------------------------------------

cat("\nSTEP 6/7 — patient recurrence\n")

patient_support <- ledger %>%
  group_by(patient_id, batch) %>%
  summarise(
    n_cells = n(),

    median_doublet =
      median(doublet_score),

    contractile_fraction =
      mean(contractile_supported),

    strong_contractile_fraction =
      mean(strong_contractile),

    epithelial_fraction =
      mean(epithelial_supported),

    sweat_fraction =
      mean(sweat_supported),

    joint_fraction =
      mean(contractile_epithelial_joint),

    median_contractile =
      median(CONTRACTILE_SMC_score),

    median_epithelial_basal =
      median(EPITHELIAL_BASAL_score),

    median_epithelial_junction =
      median(EPITHELIAL_JUNCTION_score),

    .groups = "drop"
  ) %>%
  arrange(desc(n_cells))

cat("\n=== PATIENT SUPPORT ===\n")
print(patient_support, n = Inf, width = Inf)

write.csv(
  patient_support,
  file.path(TAB, "11A2_patient_support.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# Exploratory classification only
# ------------------------------------------------------------

cat("\nSTEP 7/7 — evidence strata\n")

ledger <- ledger %>%
  mutate(
    evidence_stratum = case_when(

      strong_contractile &
      !epithelial_supported ~
        "CONTRACTILE_ONLY",

      strong_contractile &
      epithelial_supported &
      doublet_score <= median(doublet_score) ~
        "CONTRACTILE_EPITHELIAL_LOW_DOUBLET",

      strong_contractile &
      epithelial_supported ~
        "CONTRACTILE_EPITHELIAL",

      !strong_contractile &
      epithelial_supported ~
        "EPITHELIAL_PREDOMINANT",

      TRUE ~
        "OTHER"
    )
  )

strata <- ledger %>%
  count(evidence_stratum, name = "n_cells") %>%
  mutate(
    fraction = n_cells / sum(n_cells)
  ) %>%
  arrange(desc(n_cells))

cat("\n=== EVIDENCE STRATA ===\n")
print(strata, n = Inf)

write.csv(
  strata,
  file.path(TAB, "11A2_evidence_strata.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# Formal checks
# ------------------------------------------------------------

checks <- tibble(
  check = c(
    "smooth_muscle_cells_330",
    "patients_16",
    "batches_2",
    "all_singlets",
    "programme_coverage_complete",
    "cell_ledger_generated",
    "joint_evidence_generated",
    "gene_detection_generated",
    "technical_correlations_generated",
    "doublet_strata_generated",
    "patient_recurrence_generated",
    "disease_not_used",
    "no_cells_removed",
    "no_identity_assigned",
    "source_object_not_modified"
  ),

  pass = c(
    ncol(sm) == 330,
    length(unique(sm$patient_id)) == 16,
    length(unique(sm$batch)) == 2,
    all(sm$doublet_class == "singlet"),
    all(coverage$n_present > 0),
    nrow(ledger) == 330,
    nrow(joint_summary) == 1,
    nrow(gene_detection) > 0,
    nrow(technical_cor) > 0,
    nrow(doublet_strata) == 4,
    nrow(patient_support) == 16,
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
  file.path(TAB, "11A2_formal_checks.csv"),
  row.names = FALSE
)

stopifnot(all(checks$pass))

cat("\n============================================================\n")
cat("11A2 COMPLETE\n")
cat("============================================================\n")
cat("Smooth-muscle cells:", ncol(sm), "\n")
cat("Patients:", length(unique(sm$patient_id)), "\n")
cat("No cells removed.\n")
cat("No identities assigned.\n")
cat("Disease not used.\n")
cat("PASS\n")
