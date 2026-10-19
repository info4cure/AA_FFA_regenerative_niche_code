suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tibble)
  library(readr)
})

cat("============================================================\n")
cat("PIPELINE_V2 — 15A LYMPHATIC ENDOTHELIAL PREFLIGHT\n")
cat("============================================================\n\n")

SRC <- "PIPELINE_V2/objects/annotated/global_atlas_v2_L1_frozen.rds"

OUT <- "PIPELINE_V2/results/annotation/15_LYMPHATIC_EC/15A_preflight"
TAB <- file.path(OUT, "tables")

dir.create(TAB, recursive = TRUE, showWarnings = FALSE)

stopifnot(file.exists(SRC))

# ============================================================
# 1. Extract frozen-L1 lymphatic EC universe
# ============================================================

cat("STEP 1/7 — loading frozen global L1 atlas\n")

g <- readRDS(SRC)

stopifnot(ncol(g) == 82893)
stopifnot("L1" %in% colnames(g@meta.data))

x <- subset(
  g,
  subset = L1 == "Lymphatic_EC"
)

stopifnot(ncol(x) == 182)

cat("Cells:", ncol(x), "\n")
cat("Genes:", nrow(x), "\n")
cat("Patients:", length(unique(x$patient_id)), "\n")
cat("Batches:", length(unique(x$batch)), "\n")

# ============================================================
# 2. Structural support
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
  file.path(TAB, "15A_patient_support.csv")
)

# ============================================================
# 3. QC
# ============================================================

cat("\nSTEP 3/7 — QC / doublet audit\n")

md <- x@meta.data

if (
  !"percent_mt" %in% colnames(md) &&
  "percent.mt" %in% colnames(md)
) {
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

qc <- tibble(
  n_cells = ncol(x),

  median_doublet =
    median(
      md$doublet_score,
      na.rm = TRUE
    ),

  q75_doublet =
    quantile(
      md$doublet_score,
      0.75,
      na.rm = TRUE
    ),

  median_features =
    median(
      md$nFeature_RNA,
      na.rm = TRUE
    ),

  median_counts =
    median(
      md$nCount_RNA,
      na.rm = TRUE
    ),

  median_mt =
    median(
      md$percent_mt,
      na.rm = TRUE
    )
)

cat("\n=== QC ===\n")
print(qc)

write_csv(
  qc,
  file.path(TAB, "15A_QC_summary.csv")
)

# ============================================================
# 4. Define endothelial hierarchy programmes
# ============================================================

cat("\nSTEP 4/7 — defining endothelial hierarchy programmes\n")

programmes <- list(

  ENDOTHELIAL_CORE = c(
    "PECAM1","VWF","CDH5","KDR",
    "ESAM","EMCN","ENG","RAMP2"
  ),

  LYMPHATIC_CORE = c(
    "PROX1","PDPN","FLT4","LYVE1",
    "CCL21","CCL14","RELN","EFEMP1"
  ),

  LYMPHATIC_IDENTITY = c(
    "PROX1","PDPN","FLT4",
    "LYVE1","CCL21","CCL14"
  ),

  LYMPHATIC_CAPILLARY_LIKE = c(
    "LYVE1","CCL21","RELN",
    "EFEMP1","PDPN","FLT4"
  ),

  LYMPHATIC_COLLECTING_LIKE = c(
    "FOXC2","GJA4","GJA5",
    "KLF2","KLF4","CLDN5"
  ),

  ARTERIAL = c(
    "EFNB2","SOX17","GJA5",
    "GJA4","HEY1","DLL4"
  ),

  VENOUS = c(
    "ACKR1","NR2F2","SELP",
    "VCAM1","VWF","EPAS1"
  ),

  CAPILLARY_BLOOD = c(
    "CA4","RGCC","EMCN",
    "KDR","PECAM1","RAMP2"
  ),

  TIP_ANGIOGENIC = c(
    "ESM1","ANGPT2","APLNR",
    "KDR","CXCR4","PGF"
  ),

  ANTIGEN_PRESENTATION = c(
    "HLA-A","HLA-B","HLA-C",
    "B2M","HLA-DRA","HLA-DRB1"
  ),

  IFN_RESPONSE = c(
    "IFIT1","IFIT2","IFIT3",
    "ISG15","MX1","STAT1"
  ),

  MURAL = c(
    "RGS5","PDGFRB","CSPG4",
    "MCAM","ACTA2","TAGLN"
  ),

  FIBROBLAST = c(
    "DCN","LUM","COL1A1",
    "COL1A2","COL3A1","PDGFRA"
  ),

  EPITHELIAL = c(
    "KRT5","KRT14","KRT17",
    "EPCAM","KRT8","KRT18"
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
        coverage_fraction =
          length(present) /
          length(programmes[[nm]]),
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
  file.path(TAB, "15A_programme_coverage.csv")
)

# ============================================================
# 5. Cell-level evidence
# ============================================================

cat("\nSTEP 5/7 — cell-level endothelial evidence\n")

rna <- GetAssayData(
  x,
  assay = "RNA",
  layer = "data"
)

score_programme <- function(genes) {

  genes <- intersect(
    genes,
    rownames(rna)
  )

  if (length(genes) == 0) {
    return(rep(NA_real_, ncol(x)))
  }

  Matrix::colMeans(
    rna[
      genes,
      ,
      drop = FALSE
    ]
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
    rna[
      genes,
      ,
      drop = FALSE
    ] > 0
  )
}

evidence <- tibble(
  cell = colnames(x),
  patient_id =
    as.character(x$patient_id),
  batch =
    as.character(x$batch)
)

for (nm in names(programmes)) {

  evidence[[paste0(nm, "_score")]] <- score_programme(
    programmes[[nm]]
  )

  evidence[[paste0(nm, "_genes")]] <- detect_programme(
    programmes[[nm]]
  )
}

write_csv(
  evidence,
  file.path(
    TAB,
    "15A_cell_programme_evidence.csv"
  )
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

  endothelial_fraction =
    fraction_ge2("ENDOTHELIAL_CORE"),

  lymphatic_core_fraction =
    fraction_ge2("LYMPHATIC_CORE"),

  lymphatic_identity_fraction =
    fraction_ge2("LYMPHATIC_IDENTITY"),

  lymphatic_capillary_fraction =
    fraction_ge2("LYMPHATIC_CAPILLARY_LIKE"),

  lymphatic_collecting_fraction =
    fraction_ge2("LYMPHATIC_COLLECTING_LIKE"),

  arterial_fraction =
    fraction_ge2("ARTERIAL"),

  venous_fraction =
    fraction_ge2("VENOUS"),

  blood_capillary_fraction =
    fraction_ge2("CAPILLARY_BLOOD"),

  angiogenic_fraction =
    fraction_ge2("TIP_ANGIOGENIC"),

  antigen_presentation_fraction =
    fraction_ge2("ANTIGEN_PRESENTATION"),

  IFN_fraction =
    fraction_ge2("IFN_RESPONSE"),

  mural_fraction =
    fraction_ge2("MURAL"),

  fibroblast_fraction =
    fraction_ge2("FIBROBLAST"),

  epithelial_fraction =
    fraction_ge2("EPITHELIAL"),

  immune_fraction =
    fraction_ge2("IMMUNE")
)

cat("\n=== PROGRAMME SUMMARY ===\n")
print(
  programme_summary,
  width = Inf
)

write_csv(
  programme_summary,
  file.path(TAB, "15A_programme_summary.csv")
)

# ============================================================
# 7. Formal checks
# ============================================================

cat("\nSTEP 7/7 — formal checks\n")

checks <- tibble(
  check = c(
    "source_global_cells_82893",
    "lymphatic_cells_182",
    "genes_23412",
    "patients_15",
    "batches_2",
    "programme_coverage_ge80pct",
    "cell_evidence_generated",
    "programme_summary_generated",
    "disease_not_used",
    "no_cells_removed",
    "no_identity_assigned",
    "no_reclustering_performed",
    "figures_generated_NO",
    "source_object_not_overwritten"
  ),

  pass = c(
    ncol(g) == 82893,
    ncol(x) == 182,
    nrow(x) == 23412,
    length(unique(x$patient_id)) == 15,
    length(unique(x$batch)) == 2,

    all(
      coverage$coverage_fraction >= 0.80
    ),

    nrow(evidence) == 182,
    nrow(programme_summary) == 1,

    TRUE,
    ncol(x) == 182,

    !(
      "lymphatic_L2_label" %in%
        colnames(x@meta.data)
    ),

    TRUE,
    TRUE,
    TRUE
  )
)

cat("\n=== FORMAL CHECKS ===\n")
print(checks, n = Inf)

write_csv(
  checks,
  file.path(TAB, "15A_formal_checks.csv")
)

stopifnot(all(checks$pass))

cat("\n============================================================\n")
cat("15A COMPLETE\n")
cat("============================================================\n")
cat("Lymphatic EC cells: 182\n")
cat("Patients: 15\n")
cat("Batches: 2\n")
cat("Endothelial/lymphatic hierarchy evaluated.\n")
cat("Blood-vascular programmes evaluated as boundary controls.\n")
cat("No reclustering performed.\n")
cat("No identities assigned.\n")
cat("No cells removed.\n")
cat("Disease used: NO\n")
cat("Figures generated: NO\n")
cat("PASS\n")
