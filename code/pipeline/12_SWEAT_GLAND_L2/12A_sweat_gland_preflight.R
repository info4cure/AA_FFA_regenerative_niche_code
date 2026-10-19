suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(readr)
})

cat("============================================================\n")
cat("PIPELINE_V2 — 12A SWEAT-GLAND PREFLIGHT\n")
cat("============================================================\n\n")

SRC <- "PIPELINE_V2/objects/annotated/global_atlas_v2_L1_frozen.rds"

OUT <- paste0(
  "PIPELINE_V2/results/annotation/",
  "12_SWEAT_GLAND_L2/12A_preflight"
)

TAB <- file.path(OUT, "tables")

dir.create(
  TAB,
  recursive = TRUE,
  showWarnings = FALSE
)

# ============================================================
# 1. Load frozen global L1 atlas
# ============================================================

cat("STEP 1/7 — loading frozen global L1 atlas\n")

g <- readRDS(SRC)

stopifnot(ncol(g) == 82893)
stopifnot("L1" %in% colnames(g@meta.data))

x <- subset(
  g,
  subset = L1 == "Sweat_gland"
)

cat("Cells:", ncol(x), "\n")
cat(
  "Patients:",
  length(unique(x$patient_id)),
  "\n"
)
cat(
  "Batches:",
  length(unique(x$batch)),
  "\n\n"
)

stopifnot(ncol(x) == 1843)

# ============================================================
# 2. Structural inventory
# ============================================================

cat("STEP 2/7 — structural inventory\n")

patient_support <- x@meta.data %>%
  as_tibble(
    rownames = "cell"
  ) %>%
  count(
    patient_id,
    batch,
    name = "n_cells"
  ) %>%
  arrange(
    desc(n_cells)
  )

cat("\n=== PATIENT SUPPORT ===\n")
print(
  patient_support,
  n = Inf
)

write_csv(
  patient_support,
  file.path(
    TAB,
    "12A_patient_support.csv"
  )
)

# ============================================================
# 3. QC / doublet audit
# ============================================================

cat("\nSTEP 3/7 — QC / doublet audit\n")

md <- x@meta.data

if (
  !"percent_mt" %in% colnames(md) &&
  "percent.mt" %in% colnames(md)
) {
  md$percent_mt <- md$percent.mt
}

qc <- tibble(
  n_cells = nrow(md),

  singlet_fraction =
    if ("doublet_class" %in% colnames(md))
      mean(
        as.character(md$doublet_class) == "singlet"
      )
    else NA_real_,

  median_doublet =
    if ("doublet_score" %in% colnames(md))
      median(
        md$doublet_score,
        na.rm = TRUE
      )
    else NA_real_,

  q75_doublet =
    if ("doublet_score" %in% colnames(md))
      quantile(
        md$doublet_score,
        0.75,
        na.rm = TRUE,
        names = FALSE
      )
    else NA_real_,

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
    if ("percent_mt" %in% colnames(md))
      median(
        md$percent_mt,
        na.rm = TRUE
      )
    else NA_real_
)

cat("\n=== QC ===\n")
print(qc)

write_csv(
  qc,
  file.path(
    TAB,
    "12A_QC_summary.csv"
  )
)

# ============================================================
# 4. Prespecified lineage programmes
# ============================================================

cat("\nSTEP 4/7 — programme coverage\n")

programmes <- list(

  SWEAT_PAN = c(
    "KRT7",
    "KRT8",
    "KRT18",
    "KRT19",
    "KRT23",
    "KRT77",
    "KRT79"
  ),

  SECRETORY_LUMINAL = c(
    "KRT8",
    "KRT18",
    "KRT19",
    "KRT7",
    "MUC1",
    "KLF5",
    "AQP5",
    "SCGB2A2"
  ),

  DUCTAL_LUMINAL = c(
    "KRT7",
    "KRT8",
    "KRT18",
    "KRT19",
    "KRT23",
    "KRT77",
    "KRT79",
    "CLDN10"
  ),

  DUCTAL_BASAL = c(
    "KRT5",
    "KRT14",
    "KRT17",
    "TP63",
    "KRT15",
    "KRT6A",
    "KRT6B"
  ),

  MYOEPITHELIAL = c(
    "KRT5",
    "KRT14",
    "KRT17",
    "ACTA2",
    "TAGLN",
    "MYL9",
    "TP63"
  ),

  CLEAR_SECRETORY = c(
    "AQP5",
    "KLF5",
    "KRT8",
    "KRT18",
    "KRT19"
  ),

  DARK_SECRETORY = c(
    "MUC1",
    "KRT8",
    "KRT18",
    "KRT19",
    "SCGB2A2"
  ),

  EPITHELIAL_BASAL = c(
    "KRT5",
    "KRT14",
    "KRT17",
    "TP63",
    "KRT15",
    "COL17A1"
  ),

  SEBOCYTE = c(
    "PPARG",
    "PLIN2",
    "FASN",
    "SCD",
    "MGST1",
    "KRT7",
    "KRT19"
  ),

  HAIR_FOLLICLE = c(
    "KRT15",
    "KRT17",
    "SOX9",
    "KRT75",
    "KRT71",
    "KRT85"
  ),

  MELANOCYTE = c(
    "MITF",
    "MLANA",
    "PMEL",
    "TYR",
    "TYRP1",
    "DCT"
  ),

  ENDOTHELIAL = c(
    "PECAM1",
    "VWF",
    "CDH5",
    "EMCN",
    "KDR",
    "ESAM"
  ),

  IMMUNE = c(
    "PTPRC",
    "CD3D",
    "CD3E",
    "TRAC",
    "NKG7",
    "GNLY"
  )
)

genes <- rownames(x)

coverage <- tibble(
  programme = names(programmes),
  n_requested =
    sapply(
      programmes,
      length
    ),
  n_present =
    sapply(
      programmes,
      function(z)
        sum(z %in% genes)
    ),
  genes_present =
    sapply(
      programmes,
      function(z)
        paste(
          z[z %in% genes],
          collapse = ";"
        )
    )
)

print(
  coverage,
  n = Inf
)

write_csv(
  coverage,
  file.path(
    TAB,
    "12A_programme_coverage.csv"
  )
)

# ============================================================
# 5. Cell-level gene detection
# ============================================================

cat("\nSTEP 5/7 — cell-level lineage evidence\n")

DefaultAssay(x) <- "RNA"

rna <- GetAssayData(
  x,
  assay = "RNA",
  layer = "data"
)

evidence <- tibble(
  cell = colnames(x),
  patient_id =
    as.character(x$patient_id),
  batch =
    as.character(x$batch)
)

for (p in names(programmes)) {

  pg <- intersect(
    programmes[[p]],
    rownames(rna)
  )

  if (length(pg) == 0) {
    evidence[[paste0(p, "_n")]] <- 0
    evidence[[paste0(p, "_score")]] <- 0
    next
  }

  mat <- rna[
    pg,
    ,
    drop = FALSE
  ]

  evidence[[paste0(p, "_n")]] <- Matrix::colSums(
    mat > 0
  )

  evidence[[paste0(p, "_score")]] <- Matrix::colMeans(
    mat
  )
}

write_csv(
  evidence,
  file.path(
    TAB,
    "12A_cell_lineage_evidence.csv"
  )
)

# ============================================================
# 6. Programme summary / competing lineage
# ============================================================

cat("\nSTEP 6/7 — programme summary\n")

summary <- tibble(
  n_cells = nrow(evidence),

  sweat_fraction =
    mean(
      evidence$SWEAT_PAN_n >= 3
    ),

  secretory_fraction =
    mean(
      evidence$SECRETORY_LUMINAL_n >= 3
    ),

  ductal_luminal_fraction =
    mean(
      evidence$DUCTAL_LUMINAL_n >= 3
    ),

  ductal_basal_fraction =
    mean(
      evidence$DUCTAL_BASAL_n >= 3
    ),

  myoepithelial_fraction =
    mean(
      evidence$MYOEPITHELIAL_n >= 3
    ),

  sebocyte_fraction =
    mean(
      evidence$SEBOCYTE_n >= 3
    ),

  follicular_fraction =
    mean(
      evidence$HAIR_FOLLICLE_n >= 2
    ),

  melanocyte_fraction =
    mean(
      evidence$MELANOCYTE_n >= 2
    ),

  endothelial_fraction =
    mean(
      evidence$ENDOTHELIAL_n >= 3
    ),

  immune_fraction =
    mean(
      evidence$IMMUNE_n >= 2
    ),

  secretory_ductal_joint =
    mean(
      evidence$SECRETORY_LUMINAL_n >= 3 &
      evidence$DUCTAL_LUMINAL_n >= 3
    ),

  basal_myoepithelial_joint =
    mean(
      evidence$DUCTAL_BASAL_n >= 3 &
      evidence$MYOEPITHELIAL_n >= 3
    )
)

cat("\n=== PROGRAMME SUMMARY ===\n")
print(summary)

write_csv(
  summary,
  file.path(
    TAB,
    "12A_programme_summary.csv"
  )
)

# ============================================================
# 7. Formal checks
# ============================================================

cat("\nSTEP 7/7 — formal checks\n")

checks <- tibble(
  check = c(
    "source_global_cells_82893",
    "sweat_cells_1843",
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
    ncol(x) == 1843,
    length(unique(x$patient_id)) == 16,
    length(unique(x$batch)) == 2,
    nrow(coverage) == length(programmes),
    nrow(evidence) == 1843,
    nrow(summary) == 1,
    TRUE,
    ncol(x) == 1843,
    TRUE,
    TRUE,
    TRUE
  )
)

cat("\n=== FORMAL CHECKS ===\n")
print(
  checks,
  n = Inf
)

write_csv(
  checks,
  file.path(
    TAB,
    "12A_formal_checks.csv"
  )
)

stopifnot(
  all(checks$pass)
)

cat("\n============================================================\n")
cat("12A COMPLETE\n")
cat("============================================================\n")
cat("Cells:", ncol(x), "\n")
cat(
  "Patients:",
  length(unique(x$patient_id)),
  "\n"
)
cat(
  "Batches:",
  length(unique(x$batch)),
  "\n"
)
cat("No reclustering performed.\n")
cat("No identities assigned.\n")
cat("No cells removed.\n")
cat("Disease not used.\n")
cat("PASS\n")
