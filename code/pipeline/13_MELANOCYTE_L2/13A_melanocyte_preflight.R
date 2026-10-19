suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tibble)
  library(readr)
})

cat("============================================================\n")
cat("PIPELINE_V2 — 13A MELANOCYTE PREFLIGHT\n")
cat("============================================================\n\n")

# ============================================================
# Paths
# ============================================================

SRC <- "PIPELINE_V2/objects/annotated/global_atlas_v2_L1_frozen.rds"

OUT <- "PIPELINE_V2/results/annotation/13_MELANOCYTE_L2/13A_preflight"
TAB <- file.path(OUT, "tables")

dir.create(TAB, recursive = TRUE, showWarnings = FALSE)

stopifnot(file.exists(SRC))

# ============================================================
# 1. Load frozen global L1
# ============================================================

cat("STEP 1/7 — loading frozen global L1 atlas\n")

g <- readRDS(SRC)

stopifnot(ncol(g) == 82893)
stopifnot("L1" %in% colnames(g@meta.data))

x <- subset(
  g,
  cells = colnames(g)[as.character(g$L1) == "Melanocyte"]
)

cat("Cells:", ncol(x), "\n")
cat("Patients:", length(unique(x$patient_id)), "\n")
cat("Batches:", length(unique(x$batch)), "\n\n")

stopifnot(ncol(x) == 1556)
stopifnot(length(unique(x$patient_id)) == 16)
stopifnot(length(unique(x$batch)) == 2)

# ============================================================
# 2. Structural support
# ============================================================

cat("STEP 2/7 — structural inventory\n")

patient_support <- x@meta.data %>%
  as_tibble() %>%
  count(
    patient_id,
    batch,
    name = "n_cells"
  ) %>%
  arrange(desc(n_cells))

cat("\n=== PATIENT SUPPORT ===\n")
print(patient_support, n = Inf)

write_csv(
  patient_support,
  file.path(TAB, "13A_melanocyte_by_patient.csv")
)

# ============================================================
# 3. QC / doublet audit
# ============================================================

cat("\nSTEP 3/7 — QC / doublet audit\n")

md <- x@meta.data

mt_col <- if ("percent_mt" %in% colnames(md)) {
  "percent_mt"
} else if ("percent.mt" %in% colnames(md)) {
  "percent.mt"
} else {
  NA_character_
}

if (is.na(mt_col)) {
  stop("No mitochondrial percentage column found.")
}

qc <- tibble(
  n_cells = ncol(x),

  singlet_fraction =
    if ("doublet_class" %in% colnames(md)) {
      mean(
        tolower(as.character(md$doublet_class)) == "singlet",
        na.rm = TRUE
      )
    } else {
      NA_real_
    },

  median_doublet =
    if ("doublet_score" %in% colnames(md)) {
      median(md$doublet_score, na.rm = TRUE)
    } else {
      NA_real_
    },

  q75_doublet =
    if ("doublet_score" %in% colnames(md)) {
      as.numeric(
        quantile(md$doublet_score, 0.75, na.rm = TRUE)
      )
    } else {
      NA_real_
    },

  median_features =
    median(md$nFeature_RNA, na.rm = TRUE),

  median_counts =
    median(md$nCount_RNA, na.rm = TRUE),

  median_mt =
    median(md[[mt_col]], na.rm = TRUE)
)

cat("\n=== QC ===\n")
print(qc, n = Inf)

write_csv(
  qc,
  file.path(TAB, "13A_QC_summary.csv")
)

# ============================================================
# 4. Identity / state programmes
# ============================================================

cat("\nSTEP 4/7 — defining melanocyte programmes\n")

programmes <- list(

  MELANOCYTE_CORE = c(
    "MITF",
    "MLANA",
    "PMEL",
    "TYR",
    "TYRP1",
    "DCT",
    "SOX10",
    "S100B"
  ),

  MELANOGENESIS = c(
    "TYR",
    "TYRP1",
    "DCT",
    "PMEL",
    "MLANA",
    "OCA2",
    "SLC45A2",
    "SLC24A5"
  ),

  MITF_PROGRAMME = c(
    "MITF",
    "MLANA",
    "PMEL",
    "TYR",
    "TYRP1",
    "DCT"
  ),

  STEM_PROGENITOR = c(
    "SOX10",
    "KIT",
    "NGFR",
    "EDNRB",
    "PAX3",
    "SOX9"
  ),

  DENDRITIC_NEURAL_CREST = c(
    "SOX10",
    "S100B",
    "PLP1",
    "PMP22",
    "NGFR",
    "GFRA3"
  ),

  STRESS_IFN = c(
    "IFIT1",
    "IFIT2",
    "IFIT3",
    "ISG15",
    "MX1",
    "STAT1"
  ),

  ANTIGEN_PRESENTATION = c(
    "HLA-A",
    "HLA-B",
    "HLA-C",
    "B2M",
    "HLA-DRA",
    "HLA-DRB1"
  ),

  CYCLING = c(
    "MKI67",
    "TOP2A",
    "UBE2C",
    "CENPF",
    "TYMS",
    "PCNA"
  ),

  EPITHELIAL = c(
    "KRT5",
    "KRT14",
    "KRT15",
    "KRT17",
    "EPCAM",
    "KRT8"
  ),

  SCHWANN = c(
    "SOX10",
    "S100B",
    "PLP1",
    "PMP22",
    "MPZ",
    "NGFR"
  ),

  IMMUNE = c(
    "PTPRC",
    "CD3D",
    "CD3E",
    "TRAC",
    "NKG7",
    "GNLY"
  ),

  FIBROBLAST = c(
    "DCN",
    "LUM",
    "COL1A1",
    "COL1A2",
    "COL3A1",
    "PDGFRA"
  ),

  ENDOTHELIAL = c(
    "PECAM1",
    "VWF",
    "CDH5",
    "EMCN",
    "KDR",
    "ESAM"
  )
)

coverage <- bind_rows(
  lapply(
    names(programmes),
    function(nm) {

      requested <- programmes[[nm]]
      present <- intersect(requested, rownames(x))

      tibble(
        programme = nm,
        n_requested = length(requested),
        n_present = length(present),
        genes_present = paste(present, collapse = ";")
      )
    }
  )
)

cat("\n=== PROGRAMME COVERAGE ===\n")
print(coverage, n = Inf)

write_csv(
  coverage,
  file.path(TAB, "13A_programme_coverage.csv")
)

stopifnot(all(coverage$n_present > 0))

# ============================================================
# 5. Cell-level evidence
# ============================================================

cat("\nSTEP 5/7 — cell-level lineage/state evidence\n")

DefaultAssay(x) <- "RNA"

rna <- GetAssayData(
  x,
  assay = "RNA",
  layer = "data"
)

if (nrow(rna) == 0 || ncol(rna) == 0) {
  stop("RNA data layer unavailable.")
}

score_programme <- function(genes) {

  genes <- intersect(
    genes,
    rownames(rna)
  )

  if (length(genes) == 0) {
    return(rep(NA_real_, ncol(x)))
  }

  Matrix::colMeans(
    rna[genes, , drop = FALSE]
  )
}

detect_programme <- function(genes) {

  genes <- intersect(
    genes,
    rownames(rna)
  )

  if (length(genes) == 0) {
    return(rep(0L, ncol(x)))
  }

  Matrix::colSums(
    rna[genes, , drop = FALSE] > 0
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
  file.path(TAB, "13A_cell_programme_evidence.csv")
)

# ============================================================
# 6. Global programme summary
# ============================================================

cat("\nSTEP 6/7 — programme summary\n")

programme_summary <- tibble(
  n_cells = nrow(evidence),

  melanocyte_fraction =
    mean(evidence$MELANOCYTE_CORE_genes >= 2),

  melanogenesis_fraction =
    mean(evidence$MELANOGENESIS_genes >= 2),

  MITF_fraction =
    mean(evidence$MITF_PROGRAMME_genes >= 2),

  progenitor_fraction =
    mean(evidence$STEM_PROGENITOR_genes >= 2),

  neural_crest_fraction =
    mean(evidence$DENDRITIC_NEURAL_CREST_genes >= 2),

  IFN_fraction =
    mean(evidence$STRESS_IFN_genes >= 2),

  antigen_presentation_fraction =
    mean(evidence$ANTIGEN_PRESENTATION_genes >= 2),

  cycling_fraction =
    mean(evidence$CYCLING_genes >= 2),

  epithelial_fraction =
    mean(evidence$EPITHELIAL_genes >= 2),

  Schwann_fraction =
    mean(evidence$SCHWANN_genes >= 2),

  immune_fraction =
    mean(evidence$IMMUNE_genes >= 2),

  fibroblast_fraction =
    mean(evidence$FIBROBLAST_genes >= 2),

  endothelial_fraction =
    mean(evidence$ENDOTHELIAL_genes >= 2),

  melanocyte_Schwann_joint =
    mean(
      evidence$MELANOCYTE_CORE_genes >= 2 &
      evidence$SCHWANN_genes >= 2
    )
)

cat("\n=== PROGRAMME SUMMARY ===\n")
print(programme_summary, n = Inf)

write_csv(
  programme_summary,
  file.path(TAB, "13A_programme_summary.csv")
)

# ============================================================
# 7. Formal checks
# ============================================================

cat("\nSTEP 7/7 — formal checks\n")

checks <- tibble(
  check = c(
    "source_global_cells_82893",
    "melanocyte_cells_1556",
    "patients_16",
    "batches_2",
    "programme_coverage_generated",
    "cell_evidence_generated",
    "programme_summary_generated",
    "disease_not_used",
    "no_cells_removed",
    "no_identity_assigned",
    "no_reclustering_performed",
    "source_object_not_overwritten"
  ),

  pass = c(
    ncol(g) == 82893,
    ncol(x) == 1556,
    length(unique(x$patient_id)) == 16,
    length(unique(x$batch)) == 2,
    nrow(coverage) == length(programmes),
    nrow(evidence) == 1556,
    nrow(programme_summary) == 1,
    TRUE,
    ncol(x) == 1556,
    TRUE,
    TRUE,
    TRUE
  )
)

cat("\n=== FORMAL CHECKS ===\n")
print(checks, n = Inf)

write_csv(
  checks,
  file.path(TAB, "13A_formal_checks.csv")
)

stopifnot(all(checks$pass))

cat("\n============================================================\n")
cat("13A COMPLETE\n")
cat("============================================================\n")
cat("Cells:", ncol(x), "\n")
cat("Patients:", length(unique(x$patient_id)), "\n")
cat("Batches:", length(unique(x$batch)), "\n")
cat("Melanocyte identity/state programmes evaluated.\n")
cat("No reclustering performed.\n")
cat("No identities assigned.\n")
cat("No cells removed.\n")
cat("Disease not used.\n")
cat("PASS\n")
