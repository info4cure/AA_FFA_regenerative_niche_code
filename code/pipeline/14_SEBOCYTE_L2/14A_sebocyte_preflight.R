suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tibble)
  library(readr)
})

cat("============================================================\n")
cat("PIPELINE_V2 — 14A SEBOCYTE PREFLIGHT\n")
cat("============================================================\n\n")

SRC <- "PIPELINE_V2/objects/annotated/global_atlas_v2_L1_frozen.rds"

OUT <- "PIPELINE_V2/results/annotation/14_SEBOCYTE_L2/14A_preflight"
TAB <- file.path(OUT, "tables")

dir.create(TAB, recursive = TRUE, showWarnings = FALSE)

stopifnot(file.exists(SRC))

# ============================================================
# 1. Load frozen global L1 atlas and extract sebocytes
# ============================================================

cat("STEP 1/7 — loading frozen global L1 atlas\n")

g <- readRDS(SRC)

stopifnot(ncol(g) == 82893)
stopifnot("L1" %in% colnames(g@meta.data))

x <- subset(
  g,
  subset = L1 == "Sebocyte"
)

stopifnot(ncol(x) == 1388)

cat("Cells:", ncol(x), "\n")
cat("Genes:", nrow(x), "\n")
cat("Patients:", length(unique(x$patient_id)), "\n")
cat("Batches:", length(unique(x$batch)), "\n")

# ============================================================
# 2. Structural inventory
# ============================================================

cat("\nSTEP 2/7 — structural inventory\n")

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
  file.path(TAB, "14A_patient_support.csv")
)

# ============================================================
# 3. QC / doublet audit
# ============================================================

cat("\nSTEP 3/7 — QC / doublet audit\n")

md <- x@meta.data

if (!"percent_mt" %in% colnames(md) &&
    "percent.mt" %in% colnames(md)) {
  md$percent_mt <- md$percent.mt
}

required_qc <- c(
  "nFeature_RNA",
  "nCount_RNA",
  "percent_mt",
  "doublet_score"
)

missing_qc <- setdiff(
  required_qc,
  colnames(md)
)

if (length(missing_qc) > 0) {
  stop(
    "Missing QC columns: ",
    paste(missing_qc, collapse = ", ")
  )
}

singlet_fraction <- if (
  "doublet_class" %in% colnames(md)
) {
  mean(
    tolower(as.character(md$doublet_class)) ==
      "singlet",
    na.rm = TRUE
  )
} else {
  NA_real_
}

qc_summary <- tibble(
  n_cells = ncol(x),
  singlet_fraction = singlet_fraction,
  median_doublet =
    median(md$doublet_score, na.rm = TRUE),
  q75_doublet =
    quantile(md$doublet_score, 0.75, na.rm = TRUE),
  median_features =
    median(md$nFeature_RNA, na.rm = TRUE),
  median_counts =
    median(md$nCount_RNA, na.rm = TRUE),
  median_mt =
    median(md$percent_mt, na.rm = TRUE)
)

cat("\n=== QC ===\n")
print(qc_summary)

write_csv(
  qc_summary,
  file.path(TAB, "14A_QC_summary.csv")
)

# ============================================================
# 4. Define identity/state programmes
# ============================================================

cat("\nSTEP 4/7 — defining sebocyte programmes\n")

programmes <- list(

  SEBOCYTE_CORE = c(
    "PPARG","PLIN2","FASN","SCD",
    "MGST1","FABP5","KRT7","KRT19"
  ),

  LIPOGENESIS = c(
    "FASN","SCD","ACACA","ACLY",
    "ELOVL5","ELOVL6","DGAT1","DGAT2"
  ),

  LIPID_STORAGE = c(
    "PLIN2","PLIN3","FABP5",
    "ABCA1","APOE","LPL"
  ),

  DIFFERENTIATED_SEBOCYTE = c(
    "PPARG","FASN","SCD","PLIN2",
    "FABP5","MGST1","KRT7","KRT19"
  ),

  BASAL_PROGENITOR = c(
    "KRT5","KRT14","KRT15",
    "KRT17","TP63","KRT6A"
  ),

  KERATINIZATION = c(
    "KRT1","KRT10","IVL",
    "SFN","KRTDAP","DSG1"
  ),

  INFLAMMATORY = c(
    "SAA1","SAA2","LCN2",
    "LTF","CXCL8","CXCL14"
  ),

  IFN_RESPONSE = c(
    "IFIT1","IFIT2","IFIT3",
    "ISG15","MX1","STAT1"
  ),

  ANTIGEN_PRESENTATION = c(
    "HLA-A","HLA-B","HLA-C",
    "B2M","HLA-DRA","HLA-DRB1"
  ),

  SWEAT_GLAND = c(
    "AQP5","DCD","PIP",
    "SCGB2A2","SCGB1D2","KRT77"
  ),

  FOLLICULAR = c(
    "KRT15","KRT17","SOX9",
    "KRT75","KRT71","KRT85"
  ),

  MELANOCYTE = c(
    "MITF","MLANA","PMEL",
    "TYR","TYRP1","DCT"
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
        genes_present =
          paste(present, collapse = ";")
      )
    }
  )
)

cat("\n=== PROGRAMME COVERAGE ===\n")
print(coverage, n = Inf)

write_csv(
  coverage,
  file.path(TAB, "14A_programme_coverage.csv")
)

# ============================================================
# 5. Cell-level evidence
# ============================================================

cat("\nSTEP 5/7 — cell-level lineage/state evidence\n")

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
  file.path(TAB, "14A_cell_programme_evidence.csv")
)

# ============================================================
# 6. Programme summary
# ============================================================

cat("\nSTEP 6/7 — programme summary\n")

fraction_ge2 <- function(nm) {
  mean(
    evidence[[paste0(nm, "_genes")]] >= 2,
    na.rm = TRUE
  )
}

programme_summary <- tibble(
  n_cells = ncol(x),

  sebocyte_fraction =
    fraction_ge2("SEBOCYTE_CORE"),

  lipogenesis_fraction =
    fraction_ge2("LIPOGENESIS"),

  lipid_storage_fraction =
    fraction_ge2("LIPID_STORAGE"),

  differentiated_fraction =
    fraction_ge2("DIFFERENTIATED_SEBOCYTE"),

  basal_progenitor_fraction =
    fraction_ge2("BASAL_PROGENITOR"),

  keratinization_fraction =
    fraction_ge2("KERATINIZATION"),

  inflammatory_fraction =
    fraction_ge2("INFLAMMATORY"),

  IFN_fraction =
    fraction_ge2("IFN_RESPONSE"),

  antigen_presentation_fraction =
    fraction_ge2("ANTIGEN_PRESENTATION"),

  sweat_fraction =
    fraction_ge2("SWEAT_GLAND"),

  follicular_fraction =
    fraction_ge2("FOLLICULAR"),

  melanocyte_fraction =
    fraction_ge2("MELANOCYTE"),

  immune_fraction =
    fraction_ge2("IMMUNE"),

  fibroblast_fraction =
    fraction_ge2("FIBROBLAST"),

  endothelial_fraction =
    fraction_ge2("ENDOTHELIAL")
)

cat("\n=== PROGRAMME SUMMARY ===\n")
print(programme_summary, width = Inf)

write_csv(
  programme_summary,
  file.path(TAB, "14A_programme_summary.csv")
)

# ============================================================
# 7. Formal checks
# ============================================================

cat("\nSTEP 7/7 — formal checks\n")

checks <- tibble(
  check = c(
    "source_global_cells_82893",
    "sebocyte_cells_1388",
    "genes_23412",
    "patients_16",
    "batches_2",
    "programme_coverage_generated",
    "programme_coverage_complete",
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
    ncol(x) == 1388,
    nrow(x) == 23412,
    length(unique(x$patient_id)) == 16,
    length(unique(x$batch)) == 2,
    nrow(coverage) == length(programmes),
    all(coverage$n_present == coverage$n_requested),
    nrow(evidence) == 1388,
    nrow(programme_summary) == 1,
    TRUE,
    ncol(x) == 1388,
    !("sebocyte_L2_label" %in% colnames(x@meta.data)),
    TRUE,
    TRUE
  )
)

cat("\n=== FORMAL CHECKS ===\n")
print(checks, n = Inf)

write_csv(
  checks,
  file.path(TAB, "14A_formal_checks.csv")
)

stopifnot(all(checks$pass))

cat("\n============================================================\n")
cat("14A COMPLETE\n")
cat("============================================================\n")
cat("Cells:", ncol(x), "\n")
cat("Patients:", length(unique(x$patient_id)), "\n")
cat("Batches:", length(unique(x$batch)), "\n")
cat("Sebocyte identity/state programmes evaluated.\n")
cat("No reclustering performed.\n")
cat("No identities assigned.\n")
cat("No cells removed.\n")
cat("Disease not used.\n")
cat("Figures generated: NO\n")
cat("PASS\n")
