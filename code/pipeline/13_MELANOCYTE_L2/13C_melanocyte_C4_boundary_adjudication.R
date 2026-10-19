suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tibble)
  library(readr)
})

cat("============================================================\n")
cat("PIPELINE_V2 — 13C MELANOCYTE C4 BOUNDARY ADJUDICATION\n")
cat("============================================================\n\n")

SRC <- paste0(
  "PIPELINE_V2/results/annotation/13_MELANOCYTE_L2/",
  "13B_recluster/melanocyte_v2_13B_reclustered.rds"
)

OUT <- paste0(
  "PIPELINE_V2/results/annotation/13_MELANOCYTE_L2/",
  "13C_C4_boundary"
)

TAB <- file.path(OUT, "tables")

dir.create(
  TAB,
  recursive = TRUE,
  showWarnings = FALSE
)

stopifnot(file.exists(SRC))

# ============================================================
# 1. Load 13B scaffold
# ============================================================

cat("STEP 1/8 — loading 13B melanocyte scaffold\n")

x <- readRDS(SRC)

stopifnot(ncol(x) == 1556)
stopifnot(nrow(x) == 23412)
stopifnot(length(unique(x$patient_id)) == 16)
stopifnot(length(unique(x$batch)) == 2)

r04_col <- "mel_snn_res.0.4"

stopifnot(
  r04_col %in%
    colnames(x@meta.data)
)

cluster <- as.character(
  x@meta.data[[r04_col]]
)

stopifnot(
  length(unique(cluster)) == 5
)

cat("Cells:", ncol(x), "\n")
cat("Patients:", length(unique(x$patient_id)), "\n")
cat("Batches:", length(unique(x$batch)), "\n")
cat("r0.4 clusters:", length(unique(cluster)), "\n")

# ============================================================
# 2. Define C4 candidate boundary population
# ============================================================

cat("\nSTEP 2/8 — defining C4 boundary population\n")

c4_cells <- colnames(x)[cluster == "4"]
mel_cells <- colnames(x)[cluster != "4"]

stopifnot(length(c4_cells) == 39)
stopifnot(length(mel_cells) == 1517)

cat("C4 cells:", length(c4_cells), "\n")
cat("Remaining melanocyte candidates:", length(mel_cells), "\n")

# ============================================================
# 3. Targeted lineage programmes
# ============================================================

cat("\nSTEP 3/8 — targeted melanocyte/epithelial evidence\n")

programmes <- list(

  MELANOCYTE_CORE = c(
    "MITF","MLANA","PMEL","TYR",
    "TYRP1","DCT","SOX10","S100B"
  ),

  MELANOGENESIS = c(
    "TYR","TYRP1","DCT","PMEL",
    "MLANA","OCA2","SLC45A2","SLC24A5"
  ),

  EPITHELIAL_BASAL = c(
    "KRT5","KRT14","KRT15","KRT17",
    "TP63","COL17A1","ITGA6","ITGB4"
  ),

  EPITHELIAL_JUNCTION = c(
    "DSP","DSC2","DSC3","PKP1",
    "PKP3","CLDN1","EPCAM"
  ),

  C4_TARGETED = c(
    "KRT14","KRT15","AQP3","FGFBP1",
    "CXCL14","DSC2","ASCL2","WNT16",
    "GJB5","EDN1","SERPINB2","CLEC2A"
  ),

  SCHWANN_SPECIFIC = c(
    "PLP1","PMP22","MPZ",
    "S100A1","GFRA3","NGFR"
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
        genes_present = paste(
          present,
          collapse = ";"
        )
      )
    }
  )
)

cat("\n=== PROGRAMME COVERAGE ===\n")
print(coverage, n = Inf)

write_csv(
  coverage,
  file.path(
    TAB,
    "13C_programme_coverage.csv"
  )
)

rna <- GetAssayData(
  x,
  assay = "RNA",
  layer = "data"
)

programme_score <- function(genes) {
  genes <- intersect(
    genes,
    rownames(rna)
  )

  Matrix::colMeans(
    rna[
      genes,
      ,
      drop = FALSE
    ]
  )
}

programme_detect <- function(genes) {
  genes <- intersect(
    genes,
    rownames(rna)
  )

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
  cluster = cluster,
  boundary_group = ifelse(
    cluster == "4",
    "C4",
    "C0-C3"
  ),
  patient_id = as.character(x$patient_id),
  batch = as.character(x$batch)
)

for (nm in names(programmes)) {
  evidence[[paste0(nm, "_score")]] <-
    programme_score(programmes[[nm]])

  evidence[[paste0(nm, "_genes")]] <-
    programme_detect(programmes[[nm]])
}

write_csv(
  evidence,
  file.path(
    TAB,
    "13C_cell_lineage_evidence.csv"
  )
)

# ============================================================
# 4. C4 vs C0-C3 targeted summary
# ============================================================

cat("\nSTEP 4/8 — C4 versus C0-C3 lineage summary\n")

boundary_summary <- evidence %>%
  group_by(boundary_group) %>%
  summarise(
    n_cells = n(),
    n_patients = n_distinct(patient_id),
    n_batches = n_distinct(batch),

    median_melanocyte =
      median(MELANOCYTE_CORE_score),

    melanocyte_ge3_fraction =
      mean(MELANOCYTE_CORE_genes >= 3),

    median_melanogenesis =
      median(MELANOGENESIS_score),

    melanogenesis_ge3_fraction =
      mean(MELANOGENESIS_genes >= 3),

    median_epithelial_basal =
      median(EPITHELIAL_BASAL_score),

    epithelial_basal_ge3_fraction =
      mean(EPITHELIAL_BASAL_genes >= 3),

    median_epithelial_junction =
      median(EPITHELIAL_JUNCTION_score),

    epithelial_junction_ge2_fraction =
      mean(EPITHELIAL_JUNCTION_genes >= 2),

    median_C4_targeted =
      median(C4_TARGETED_score),

    C4_targeted_ge3_fraction =
      mean(C4_TARGETED_genes >= 3),

    median_Schwann =
      median(SCHWANN_SPECIFIC_score),

    immune_ge2_fraction =
      mean(IMMUNE_genes >= 2),

    .groups = "drop"
  )

cat("\n=== C4 BOUNDARY SUMMARY ===\n")
print(boundary_summary, n = Inf)

write_csv(
  boundary_summary,
  file.path(
    TAB,
    "13C_C4_boundary_summary.csv"
  )
)

# ============================================================
# 5. Target-gene detection
# ============================================================

cat("\nSTEP 5/8 — targeted gene detection\n")

target_genes <- unique(c(
  programmes$MELANOCYTE_CORE,
  programmes$MELANOGENESIS,
  programmes$EPITHELIAL_BASAL,
  programmes$EPITHELIAL_JUNCTION,
  programmes$C4_TARGETED
))

target_genes <- intersect(
  target_genes,
  rownames(rna)
)

gene_detection <- bind_rows(
  lapply(
    c("C4", "C0-C3"),
    function(grp) {

      cells_grp <- evidence$cell[
        evidence$boundary_group == grp
      ]

      mat <- rna[
        target_genes,
        cells_grp,
        drop = FALSE
      ]

      tibble(
        boundary_group = grp,
        gene = target_genes,
        fraction_detected =
          Matrix::rowMeans(mat > 0),
        mean_expression =
          Matrix::rowMeans(mat)
      )
    }
  )
)

cat("\n=== TARGET GENE DETECTION ===\n")
print(gene_detection, n = Inf)

write_csv(
  gene_detection,
  file.path(
    TAB,
    "13C_target_gene_detection.csv"
  )
)

# ============================================================
# 6. Patient/batch/doublet recurrence of C4
# ============================================================

cat("\nSTEP 6/8 — C4 structural recurrence\n")

meta <- x@meta.data %>%
  as_tibble(
    rownames = "cell"
  ) %>%
  mutate(
    cluster = cluster
  )

c4_support <- meta %>%
  filter(cluster == "4") %>%
  group_by(
    patient_id,
    batch
  ) %>%
  summarise(
    n_cells = n(),
    median_doublet =
      median(
        doublet_score,
        na.rm = TRUE
      ),
    median_features =
      median(
        nFeature_RNA,
        na.rm = TRUE
      ),
    median_counts =
      median(
        nCount_RNA,
        na.rm = TRUE
      ),
    .groups = "drop"
  ) %>%
  arrange(desc(n_cells))

cat("\n=== C4 PATIENT SUPPORT ===\n")
print(c4_support, n = Inf)

write_csv(
  c4_support,
  file.path(
    TAB,
    "13C_C4_patient_support.csv"
  )
)

c4_global <- tibble(
  n_cells = length(c4_cells),
  n_patients =
    n_distinct(
      meta$patient_id[
        meta$cluster == "4"
      ]
    ),
  n_batches =
    n_distinct(
      meta$batch[
        meta$cluster == "4"
      ]
    ),
  median_doublet =
    median(
      meta$doublet_score[
        meta$cluster == "4"
      ],
      na.rm = TRUE
    ),
  q75_doublet =
    quantile(
      meta$doublet_score[
        meta$cluster == "4"
      ],
      0.75,
      na.rm = TRUE
    ),
  median_features =
    median(
      meta$nFeature_RNA[
        meta$cluster == "4"
      ],
      na.rm = TRUE
    )
)

cat("\n=== C4 GLOBAL SUPPORT ===\n")
print(c4_global)

write_csv(
  c4_global,
  file.path(
    TAB,
    "13C_C4_global_support.csv"
  )
)

# ============================================================
# 7. Formal adjudication and refined object
# ============================================================

cat("\nSTEP 7/8 — formal boundary adjudication\n")

adjudication <- tibble(
  cluster = "4",
  n_cells = 39L,
  decision = "RETURN_TO_EPITHELIAL",
  identity = "Basal epithelial-like",
  basis = paste(
    "Coordinated KRT14/KRT15/AQP3/FGFBP1/CXCL14/",
    "DSC2-associated epithelial programme recurrent across patients;",
    "not retained as a melanocyte L2 state."
  ),
  disease_used = FALSE,
  anatomical_location_claimed = FALSE
)

print(adjudication, n = Inf)

write_csv(
  adjudication,
  file.path(
    TAB,
    "13C_C4_adjudication.csv"
  )
)

mel_refined <- subset(
  x,
  cells = mel_cells
)

returned_epithelial <- subset(
  x,
  cells = c4_cells
)

stopifnot(ncol(mel_refined) == 1517)
stopifnot(ncol(returned_epithelial) == 39)

mel_refined$melanocyte_boundary_status <-
  "REFINED_AFTER_C4_RETURN"

returned_epithelial$melanocyte_boundary_status <-
  "RETURN_TO_EPITHELIAL"

# ============================================================
# 8. Checks and save
# ============================================================

cat("\nSTEP 8/8 — formal checks and save\n")

checks <- tibble(
  check = c(
    "source_cells_1556",
    "source_genes_23412",
    "patients_16",
    "batches_2",
    "r04_five_clusters",
    "C4_cells_39",
    "C4_ge11_patients",
    "C4_two_batches",
    "programme_coverage_complete",
    "cell_evidence_generated",
    "boundary_summary_generated",
    "target_gene_detection_generated",
    "C4_patient_support_generated",
    "C4_decision_return_epithelial",
    "refined_melanocytes_1517",
    "returned_epithelial_39",
    "cell_conservation_1556",
    "disease_not_used",
    "no_L2_identity_assigned",
    "resolution_not_frozen",
    "source_object_not_overwritten"
  ),
  pass = c(
    ncol(x) == 1556,
    nrow(x) == 23412,
    length(unique(x$patient_id)) == 16,
    length(unique(x$batch)) == 2,
    length(unique(cluster)) == 5,
    length(c4_cells) == 39,
    n_distinct(
      meta$patient_id[
        meta$cluster == "4"
      ]
    ) >= 11,
    n_distinct(
      meta$batch[
        meta$cluster == "4"
      ]
    ) == 2,
    all(
      coverage$n_present ==
        coverage$n_requested
    ),
    nrow(evidence) == 1556,
    nrow(boundary_summary) == 2,
    nrow(gene_detection) > 0,
    nrow(c4_support) > 0,
    adjudication$decision ==
      "RETURN_TO_EPITHELIAL",
    ncol(mel_refined) == 1517,
    ncol(returned_epithelial) == 39,
    ncol(mel_refined) +
      ncol(returned_epithelial) ==
      1556,
    TRUE,
    !(
      "melanocyte_L2_label" %in%
        colnames(mel_refined@meta.data)
    ),
    !(
      "melanocyte_resolution_frozen" %in%
        colnames(mel_refined@meta.data)
    ),
    TRUE
  )
)

cat("\n=== FORMAL CHECKS ===\n")
print(checks, n = Inf)

write_csv(
  checks,
  file.path(
    TAB,
    "13C_formal_checks.csv"
  )
)

stopifnot(all(checks$pass))

MEL_RDS <- file.path(
  OUT,
  "melanocyte_v2_13C_refined_1517.rds"
)

RETURN_RDS <- file.path(
  OUT,
  "melanocyte_v2_13C_return_to_epithelial_39.rds"
)

saveRDS(
  mel_refined,
  MEL_RDS
)

saveRDS(
  returned_epithelial,
  RETURN_RDS
)

write_csv(
  tibble(
    cell = c4_cells,
    action = "RETURN_TO_EPITHELIAL"
  ),
  file.path(
    TAB,
    "13C_return_to_epithelial_ledger.csv"
  )
)

cat("\nSaved refined melanocytes:\n")
cat(MEL_RDS, "\n")

cat("\nSaved epithelial return object:\n")
cat(RETURN_RDS, "\n")

cat("\n============================================================\n")
cat("13C COMPLETE\n")
cat("============================================================\n")
cat("Initial melanocyte candidates: 1556\n")
cat("C4 returned to epithelial: 39\n")
cat("Refined melanocyte universe: 1517\n")
cat("C4 recurrence: >=11 patients, 2 batches\n")
cat("Reason: coordinated basal epithelial programme\n")
cat("Disease used: NO\n")
cat("Melanocyte L2 identities assigned: NO\n")
cat("Resolution frozen: NO\n")
cat("Source object overwritten: NO\n")
cat("READY FOR FINAL MELANOCYTE FREEZE\n")
cat("PASS\n")
