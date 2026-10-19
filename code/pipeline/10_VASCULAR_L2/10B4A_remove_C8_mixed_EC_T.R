suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(readr)
})

options(width = 200)

cat("============================================================\n")
cat("PIPELINE_V2 — 10B4A C8 MIXED EC–T ADJUDICATION\n")
cat("============================================================\n\n")

SRC <- "PIPELINE_V2/objects/annotated/vascular_v2_10B1_reclustered_3040.rds"

OUT <- paste0(
  "PIPELINE_V2/results/annotation/",
  "10_VASCULAR_L2/10B4A_C8_mixed_EC_T"
)

TAB <- file.path(OUT, "tables")
dir.create(TAB, recursive = TRUE, showWarnings = FALSE)

OBJOUT <- paste0(
  "PIPELINE_V2/objects/annotated/",
  "vascular_v2_10B4A_C8_removed_3023.rds"
)

x <- readRDS(SRC)

CL <- "vascular_snn_res.0.3"
stopifnot(CL %in% colnames(x@meta.data))

cl <- as.character(x@meta.data[[CL]])
c8 <- colnames(x)[cl == "8"]

stopifnot(length(c8) == 17)

cat("Source cells:", ncol(x), "\n")
cat("C8 cells:", length(c8), "\n")

# ------------------------------------------------------------
# Explicit cell-level evidence
# ------------------------------------------------------------

DefaultAssay(x) <- "RNA"

EC <- c(
  "PECAM1","VWF","CDH5","EMCN",
  "KDR","ESAM","ENG"
)

IMM <- c(
  "PTPRC","CD3D","CD3E","TRAC",
  "TRBC1","TRBC2","NKG7","GNLY",
  "PRF1","CD8A","CD8B"
)

genes <- intersect(
  unique(c(EC, IMM)),
  rownames(x)
)

e <- FetchData(
  x,
  vars = genes,
  cells = c8,
  layer = "data"
)

count_detected <- function(df, genes) {
  g <- intersect(genes, colnames(df))
  rowSums(df[, g, drop = FALSE] > 0)
}

ledger <- tibble(
  cell = rownames(e),
  patient_id = as.character(
    x$patient_id[
      match(rownames(e), colnames(x))
    ]
  ),
  batch = as.character(
    x$batch[
      match(rownames(e), colnames(x))
    ]
  ),
  doublet_class = as.character(
    x$doublet_class[
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
    ],
  endothelial_genes_detected =
    count_detected(e, EC),
  immune_TNK_genes_detected =
    count_detected(e, IMM)
) %>%
  mutate(
    endothelial_supported =
      endothelial_genes_detected >= 3,
    immune_supported =
      immune_TNK_genes_detected >= 3,
    joint_EC_immune =
      endothelial_supported & immune_supported,
    adjudication =
      "MIXED_EC_T_PROBABLE_RESIDUAL_DOUBLET",
    disposition =
      "EXCLUDE_FROM_REFINED_VASCULAR_UNIVERSE"
  )

write_csv(
  ledger,
  file.path(TAB, "10B4A_C8_cell_adjudication_ledger.csv")
)

cat("\n=== C8 CELL-LEVEL LEDGER ===\n")
print(ledger, n = Inf)

cat("\nJoint EC + immune:",
    sum(ledger$joint_EC_immune),
    "/", nrow(ledger), "\n")

# ------------------------------------------------------------
# Remove only adjudicated C8
# ------------------------------------------------------------

keep <- setdiff(colnames(x), c8)

y <- subset(x, cells = keep)

stopifnot(ncol(y) == 3023)
stopifnot(length(intersect(colnames(y), c8)) == 0)

y$vascular_10B4A_status <-
  "RETAINED_FOR_POST_C8_REBUILD"

# ------------------------------------------------------------
# Removal ledger
# ------------------------------------------------------------

removal <- tibble(
  source_cells = ncol(x),
  excluded_cells = length(c8),
  retained_cells = ncol(y),
  excluded_cluster_pre_rebuild = "r0.3_C8",
  reason =
    "concurrent strong endothelial and T/NK transcriptional programmes",
  interpretation =
    "mixed EC-T / probable residual doublets",
  returned_to_immune = FALSE,
  disease_used = FALSE
)

write_csv(
  removal,
  file.path(TAB, "10B4A_removal_summary.csv")
)

# ------------------------------------------------------------
# Formal checks
# ------------------------------------------------------------

checks <- tibble(
  check = c(
    "source_cells_3040",
    "C8_cells_17",
    "C8_all_endothelial_supported",
    "C8_joint_EC_immune_16_of_17",
    "retained_cells_3023",
    "all_C8_removed",
    "C8_not_returned_to_immune",
    "disease_not_used",
    "no_other_cells_removed",
    "source_object_not_overwritten"
  ),
  pass = c(
    ncol(x) == 3040,
    length(c8) == 17,
    all(ledger$endothelial_supported),
    sum(ledger$joint_EC_immune) == 16,
    ncol(y) == 3023,
    length(intersect(colnames(y), c8)) == 0,
    TRUE,
    TRUE,
    ncol(x) - ncol(y) == 17,
    TRUE
  )
)

write_csv(
  checks,
  file.path(TAB, "10B4A_formal_checks.csv")
)

cat("\n=== FORMAL CHECKS ===\n")
print(checks, n = Inf)

stopifnot(all(checks$pass))

saveRDS(
  y,
  OBJOUT,
  compress = FALSE
)

writeLines(
  c(
    "PIPELINE_V2 — 10B4A C8 MIXED EC–T ADJUDICATION",
    "",
    "Source vascular universe: 3,040 cells",
    "Excluded: 17 r0.3 C8 cells",
    "Retained: 3,023 cells",
    "",
    "Rationale:",
    "- all 17 cells retained strong endothelial support",
    "- 16/17 simultaneously expressed >=3 immune/T-NK genes",
    "- C8 showed PTPRC/CD3D/CD3E/TRAC/TRBC and cytotoxic signals",
    "- cells were therefore not reassigned as pure immune cells",
    "- pattern adjudicated as mixed EC-T / probable residual doublets",
    "",
    "C4 low-complexity endothelial cells were retained.",
    "C6 endothelial-mesenchymal/contractile cells were retained.",
    "C7 rare bona fide endothelial cells were retained.",
    "",
    "Disease status was not used.",
    "No other cells were removed.",
    "",
    "PASS"
  ),
  file.path(TAB, "10B4A_audit.txt")
)

cat("\n============================================================\n")
cat("10B4A COMPLETE\n")
cat("============================================================\n")
cat("Source:", ncol(x), "\n")
cat("Excluded mixed EC-T:", length(c8), "\n")
cat("Retained:", ncol(y), "\n")
cat("PASS\n")
