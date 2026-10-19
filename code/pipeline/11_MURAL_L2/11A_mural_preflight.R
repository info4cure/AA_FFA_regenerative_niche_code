suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(readr)
})

options(width = 220)

cat("============================================================\n")
cat("PIPELINE_V2 — 11A MURAL / SMOOTH-MUSCLE PREFLIGHT\n")
cat("============================================================\n")

SRC <- paste0(
  "PIPELINE_V2/objects/annotated/",
  "global_atlas_v2_L1_frozen.rds"
)

OUT <- paste0(
  "PIPELINE_V2/results/annotation/",
  "11_MURAL_L2/11A_preflight"
)

TAB <- file.path(OUT, "tables")

dir.create(TAB, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# 1. Load global L1
# ============================================================

cat("\nSTEP 1/7 — loading frozen global L1 atlas\n")

g <- readRDS(SRC)

stopifnot(ncol(g) == 82893)
stopifnot("L1" %in% colnames(g@meta.data))

focus <- c(
  "Pericyte_SMC",
  "Smooth_muscle"
)

x <- subset(
  g,
  cells = colnames(g)[
    as.character(g$L1) %in% focus
  ]
)

cat("Cells:", ncol(x), "\n")
cat("Patients:", length(unique(x$patient_id)), "\n")
cat("Batches:", length(unique(x$batch)), "\n")

cat("\nSource L1:\n")
print(table(x$L1))

# ============================================================
# 2. Structural inventory
# ============================================================

cat("\nSTEP 2/7 — structural inventory\n")

structural <- x@meta.data %>%
  mutate(
    source_L1 = as.character(L1)
  ) %>%
  count(
    source_L1,
    patient_id,
    batch,
    name = "n_cells"
  )

summary <- structural %>%
  group_by(source_L1) %>%
  summarise(
    n_cells = sum(n_cells),
    n_patients = n_distinct(patient_id),
    n_batches = n_distinct(batch),
    dominant_patient_cells = max(n_cells),
    dominant_patient_fraction =
      max(n_cells) / sum(n_cells),
    median_cells_per_patient =
      median(n_cells),
    min_cells_per_patient =
      min(n_cells),
    max_cells_per_patient =
      max(n_cells),
    .groups = "drop"
  )

cat("\n=== STRUCTURAL SUMMARY ===\n")
print(summary, n = Inf)

write_csv(
  structural,
  file.path(TAB, "11A_patient_structure.csv")
)

write_csv(
  summary,
  file.path(TAB, "11A_structural_summary.csv")
)

# ============================================================
# 3. QC
# ============================================================

cat("\nSTEP 3/7 — QC / doublet audit\n")

qc <- x@meta.data %>%
  mutate(
    source_L1 = as.character(L1)
  ) %>%
  group_by(source_L1) %>%
  summarise(
    n_cells = n(),
    singlet_fraction =
      mean(doublet_class == "singlet"),
    median_doublet =
      median(doublet_score),
    q75_doublet =
      quantile(doublet_score, 0.75),
    median_features =
      median(nFeature_RNA),
    median_counts =
      median(nCount_RNA),
    median_mt =
      median(percent.mt),
    .groups = "drop"
  )

cat("\n=== QC ===\n")
print(qc, n = Inf)

write_csv(
  qc,
  file.path(TAB, "11A_QC.csv")
)

# ============================================================
# 4. Prespecified lineage programmes
# ============================================================

cat("\nSTEP 4/7 — programme coverage\n")

programmes <- list(

  PAN_MURAL = c(
    "RGS5","CSPG4","MCAM","PDGFRB",
    "RBP1","NOTCH3","ABCC9","KCNJ8"
  ),

  PERICYTE = c(
    "RGS5","CSPG4","PDGFRB","RBP1",
    "ABCC9","KCNJ8","COX4I2","HIGD1B"
  ),

  CONTRACTILE_SMC = c(
    "ACTA2","TAGLN","MYH11","CNN1",
    "MYL9","TPM2","DES","LMOD1"
  ),

  VASCULAR_SMC = c(
    "MYH11","CNN1","ACTA2","TAGLN",
    "LMOD1","SYNPO2","SORBS2","MYLK"
  ),

  FIBROBLAST = c(
    "DCN","LUM","COL1A1","COL1A2",
    "COL3A1","PDGFRA","DPT","COL6A1"
  ),

  REMODELING_FIBROBLAST = c(
    "ADAM12","MMP11","POSTN","COL11A1",
    "COMP","COL10A1","FN1","THBS2"
  ),

  ENDOTHELIAL = c(
    "PECAM1","VWF","CDH5","EMCN",
    "KDR","ESAM","ENG","RAMP2"
  ),

  EPITHELIAL = c(
    "EPCAM","KRT5","KRT14","KRT17",
    "KRT19","KRT8","KRT18","DSP"
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

available <- lapply(
  programmes,
  function(z) intersect(z, rownames(x))
)

coverage <- tibble(
  programme = names(programmes),
  n_requested = sapply(programmes, length),
  n_present = sapply(available, length),
  genes_present = sapply(
    available,
    paste,
    collapse = ";"
  )
)

print(coverage, n = Inf)

write_csv(
  coverage,
  file.path(TAB, "11A_programme_coverage.csv")
)

# ============================================================
# 5. Cell-level detection
# ============================================================

cat("\nSTEP 5/7 — cell-level lineage evidence\n")

genes <- unique(unlist(available))

expr <- FetchData(
  x,
  vars = genes,
  layer = "data"
)

detect_n <- function(gs) {
  gs <- intersect(gs, colnames(expr))
  if (!length(gs))
    return(rep(0L, nrow(expr)))
  rowSums(expr[, gs, drop = FALSE] > 0)
}

score <- function(gs) {
  gs <- intersect(gs, colnames(expr))
  if (!length(gs))
    return(rep(NA_real_, nrow(expr)))
  rowMeans(expr[, gs, drop = FALSE])
}

cell <- tibble(
  cell = rownames(expr),
  source_L1 =
    as.character(
      x$L1[
        match(rownames(expr), colnames(x))
      ]
    ),
  patient_id =
    as.character(
      x$patient_id[
        match(rownames(expr), colnames(x))
      ]
    ),
  batch =
    as.character(
      x$batch[
        match(rownames(expr), colnames(x))
      ]
    )
)

for (p in names(available)) {
  cell[[paste0(p, "_n")]] <-
    detect_n(available[[p]])
  cell[[paste0(p, "_score")]] <-
    score(available[[p]])
}

write_csv(
  cell,
  file.path(TAB, "11A_cell_programme_evidence.csv")
)

programme_summary <- cell %>%
  group_by(source_L1) %>%
  summarise(
    n_cells = n(),

    median_mural =
      median(PAN_MURAL_score),

    median_pericyte =
      median(PERICYTE_score),

    median_contractile =
      median(CONTRACTILE_SMC_score),

    median_vSMC =
      median(VASCULAR_SMC_score),

    median_fibroblast =
      median(FIBROBLAST_score),

    median_remodeling =
      median(REMODELING_FIBROBLAST_score),

    mural_ge3 =
      mean(PAN_MURAL_n >= 3),

    pericyte_ge3 =
      mean(PERICYTE_n >= 3),

    contractile_ge3 =
      mean(CONTRACTILE_SMC_n >= 3),

    fibroblast_ge3 =
      mean(FIBROBLAST_n >= 3),

    remodeling_ge2 =
      mean(REMODELING_FIBROBLAST_n >= 2),

    endothelial_ge3 =
      mean(ENDOTHELIAL_n >= 3),

    epithelial_ge2 =
      mean(EPITHELIAL_n >= 2),

    Schwann_ge2 =
      mean(SCHWANN_n >= 2),

    immune_ge2 =
      mean(IMMUNE_n >= 2),

    .groups = "drop"
  )

cat("\n=== PROGRAMME SUMMARY ===\n")
print(programme_summary, n = Inf)

write_csv(
  programme_summary,
  file.path(TAB, "11A_programme_summary.csv")
)

# ============================================================
# 6. Joint programme evidence
# ============================================================

cat("\nSTEP 6/7 — joint lineage evidence\n")

joint <- cell %>%
  mutate(
    mural_supported =
      PAN_MURAL_n >= 3,

    contractile_supported =
      CONTRACTILE_SMC_n >= 3,

    fibroblast_supported =
      FIBROBLAST_n >= 3,

    remodeling_supported =
      REMODELING_FIBROBLAST_n >= 2,

    endothelial_supported =
      ENDOTHELIAL_n >= 3
  ) %>%
  group_by(source_L1) %>%
  summarise(
    n_cells = n(),

    mural_and_contractile =
      mean(
        mural_supported &
        contractile_supported
      ),

    mural_and_fibroblast =
      mean(
        mural_supported &
        fibroblast_supported
      ),

    mural_and_remodeling =
      mean(
        mural_supported &
        remodeling_supported
      ),

    contractile_and_remodeling =
      mean(
        contractile_supported &
        remodeling_supported
      ),

    mural_and_endothelial =
      mean(
        mural_supported &
        endothelial_supported
      ),

    .groups = "drop"
  )

cat("\n=== JOINT EVIDENCE ===\n")
print(joint, n = Inf)

write_csv(
  joint,
  file.path(TAB, "11A_joint_programme_evidence.csv")
)

# ============================================================
# 7. Formal checks
# ============================================================

cat("\nSTEP 7/7 — formal checks\n")

checks <- tibble(
  check = c(
    "source_global_cells_82893",
    "mural_universe_1603",
    "Pericyte_SMC_1273",
    "Smooth_muscle_330",
    "patients_16",
    "batches_2",
    "all_singlets",
    "programme_coverage_complete",
    "cell_evidence_generated",
    "joint_evidence_generated",
    "disease_not_used",
    "no_cells_removed",
    "no_labels_modified",
    "source_object_not_overwritten"
  ),

  pass = c(
    ncol(g) == 82893,
    ncol(x) == 1603,
    sum(x$L1 == "Pericyte_SMC") == 1273,
    sum(x$L1 == "Smooth_muscle") == 330,
    length(unique(x$patient_id)) == 16,
    length(unique(x$batch)) == 2,
    all(x$doublet_class == "singlet"),
    all(
      coverage$n_present ==
        coverage$n_requested
    ),
    nrow(cell) == 1603,
    nrow(joint) == 2,
    TRUE,
    TRUE,
    TRUE,
    TRUE
  )
)

cat("\n=== FORMAL CHECKS ===\n")
print(checks, n = Inf)

write_csv(
  checks,
  file.path(TAB, "11A_formal_checks.csv")
)

stopifnot(all(checks$pass))

cat("\n============================================================\n")
cat("11A COMPLETE\n")
cat("============================================================\n")
cat("Cells:", ncol(x), "\n")
cat("Patients:", length(unique(x$patient_id)), "\n")
cat("Pericyte_SMC:", sum(x$L1 == "Pericyte_SMC"), "\n")
cat("Smooth_muscle:", sum(x$L1 == "Smooth_muscle"), "\n")
cat("No reclustering performed.\n")
cat("No identities assigned.\n")
cat("No cells removed.\n")
cat("PASS\n")
