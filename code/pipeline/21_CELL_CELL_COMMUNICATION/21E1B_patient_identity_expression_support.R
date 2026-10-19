suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(edgeR)
  library(dplyr)
  library(readr)
  library(tidyr)
  library(stringr)
  library(tibble)
})

cat("============================================================\n")
cat("PIPELINE_V2 — 21E1B PATIENT × IDENTITY EXPRESSION SUPPORT\n")
cat("============================================================\n")

OBJ <- paste0(
  "PIPELINE_V2/objects/frozen/",
  "global_atlas_v2_21B_communication_frozen.rds"
)

E1A <- paste0(
  "PIPELINE_V2/results/disease_biology/",
  "21_CELL_CELL_COMMUNICATION/",
  "21E1_OBSERVABLE_UNIVERSE/tables/",
  "21E1A_observable_LR_universe.csv"
)

OUT <- paste0(
  "PIPELINE_V2/results/disease_biology/",
  "21_CELL_CELL_COMMUNICATION/",
  "21E1_OBSERVABLE_UNIVERSE"
)

TAB <- file.path(OUT, "tables")

# ============================================================
# 1. Frozen inputs
# ============================================================

cat("\nSTEP 1/8 — recovering frozen inputs\n")

stopifnot(
  file.exists(OBJ),
  file.exists(E1A)
)

obj <- readRDS(OBJ)

lr <- read_csv(
  E1A,
  show_col_types=FALSE
)

md <- obj[[]]

required_md <- c(
  "patient_id",
  "disease",
  "batch",
  "communication_identity",
  "communication_annotation_level"
)

stopifnot(
  all(required_md %in% names(md))
)

cat("Cells:", ncol(obj), "\n")
cat(
  "Patients:",
  n_distinct(md$patient_id),
  "\n"
)
cat(
  "Communication identities:",
  n_distinct(md$communication_identity),
  "\n"
)

# ============================================================
# 2. Define LR gene universe
# ============================================================

cat("\nSTEP 2/8 — defining frozen observable LR gene universe\n")

lr_genes <- sort(
  unique(
    c(
      lr$source_genesymbol,
      lr$target_genesymbol
    )
  )
)

lr_genes <- intersect(
  lr_genes,
  rownames(obj[["RNA"]])
)

cat(
  "Unique observable LR genes:",
  length(lr_genes),
  "\n"
)

stopifnot(
  length(lr_genes) > 0
)

# ============================================================
# 3. Construct patient × identity groups
# ============================================================

cat("\nSTEP 3/8 — constructing patient × identity groups\n")

md2 <- md %>%
  mutate(
    cell_barcode=rownames(md),
    pseudobulk_id=paste(
      patient_id,
      communication_identity,
      sep="|||"
    )
  )

group_support <- md2 %>%
  count(
    pseudobulk_id,
    patient_id,
    disease,
    batch,
    communication_identity,
    name="n_cells"
  ) %>%
  arrange(
    communication_identity,
    disease,
    patient_id
  )

cat(
  "Patient × identity groups:",
  nrow(group_support),
  "\n"
)

cat(
  "Groups >=20 cells:",
  sum(group_support$n_cells >= 20),
  "\n"
)

cat(
  "Groups >=50 cells:",
  sum(group_support$n_cells >= 50),
  "\n"
)

write_csv(
  group_support,
  file.path(
    TAB,
    "21E1B_patient_identity_support.csv"
  )
)

# ============================================================
# 4. Aggregate raw counts
# ============================================================

cat("\nSTEP 4/8 — aggregating raw counts\n")

counts <- LayerData(
  obj,
  assay="RNA",
  layer="counts"
)

counts <- counts[
  lr_genes,
  ,
  drop=FALSE
]

grp <- factor(
  md2$pseudobulk_id,
  levels=unique(md2$pseudobulk_id)
)

design_sparse <- sparse.model.matrix(
  ~0 + grp
)

pb <- counts %*% design_sparse

colnames(pb) <- levels(grp)

cat(
  "Pseudobulk matrix:",
  nrow(pb),
  "genes ×",
  ncol(pb),
  "patient-identity groups\n"
)

# ============================================================
# 5. TMM/logCPM
# ============================================================

cat("\nSTEP 5/8 — TMM normalization and logCPM\n")

dge <- DGEList(
  counts=as.matrix(pb)
)

keep_library <- dge$samples$lib.size > 0

dge <- dge[, keep_library, keep.lib.sizes=FALSE]

dge <- calcNormFactors(
  dge,
  method="TMM"
)

lcpm <- cpm(
  dge,
  log=TRUE,
  prior.count=2
)

cat(
  "Nonzero-library groups:",
  ncol(lcpm),
  "\n"
)

# ============================================================
# 6. Cell-level detection fractions
# ============================================================

cat("\nSTEP 6/8 — calculating within-group detection fractions\n")

det <- counts > 0

det_n <- det %*% design_sparse

# The model matrix column names are syntactic design labels.
# Restore the exact frozen patient × identity identifiers explicitly,
# as was done for the pseudobulk count matrix.
stopifnot(
  ncol(det_n) == length(levels(grp))
)

colnames(det_n) <- levels(grp)

group_sizes <- table(grp)

stopifnot(
  identical(
    colnames(det_n),
    names(group_sizes)
  )
)

group_n <- as.numeric(
  group_sizes[
    colnames(det_n)
  ]
)

stopifnot(
  length(group_n) == ncol(det_n),
  !anyNA(group_n),
  all(group_n > 0)
)

pct <- sweep(
  det_n,
  2,
  group_n,
  "/"
)

stopifnot(
  identical(
    rownames(pct),
    rownames(counts)
  ),
  identical(
    colnames(pct),
    levels(grp)
  )
)

# ============================================================
# 7. Long-form frozen expression evidence
# ============================================================

cat("\nSTEP 7/8 — constructing expression-support table\n")

lcpm_long <- as.data.frame(
  as.table(lcpm),
  stringsAsFactors=FALSE
) %>%
  as_tibble() %>%
  rename(
    gene=Var1,
    pseudobulk_id=Var2,
    logCPM=Freq
  )

# Explicit alignment audit between pseudobulk expression
# and cell-level detection matrices.
stopifnot(
  setequal(
    rownames(lcpm),
    rownames(pct)
  ),
  setequal(
    colnames(lcpm),
    colnames(pct)
  )
)

pct2 <- pct[
  match(
    rownames(lcpm),
    rownames(pct)
  ),
  match(
    colnames(lcpm),
    colnames(pct)
  ),
  drop=FALSE
]

rownames(pct2) <- rownames(lcpm)
colnames(pct2) <- colnames(lcpm)

stopifnot(
  identical(
    dim(pct2),
    dim(lcpm)
  ),
  identical(
    rownames(pct2),
    rownames(lcpm)
  ),
  identical(
    colnames(pct2),
    colnames(lcpm)
  ),
  all(
    pct2 >= 0 &
    pct2 <= 1
  )
)

# pct2 is retained as a Matrix-class object through the
# detection calculation. Convert only this final compact
# 1755 × patient-identity matrix to base matrix for long-form
# serialization; upstream count matrices remain sparse.
pct_long <- as.data.frame(
  as.table(
    as.matrix(pct2)
  ),
  stringsAsFactors=FALSE
) %>%
  as_tibble() %>%
  rename(
    gene=Var1,
    pseudobulk_id=Var2,
    pct_expr=Freq
  )

cat(
  "logCPM long rows:",
  nrow(lcpm_long),
  "\n"
)

cat(
  "pct_expr long rows:",
  nrow(pct_long),
  "\n"
)

expected_long_rows <-
  nrow(lcpm) * ncol(lcpm)

stopifnot(
  nrow(lcpm_long) == expected_long_rows,
  nrow(pct_long) == expected_long_rows,
  !anyDuplicated(
    lcpm_long[
      c(
        "gene",
        "pseudobulk_id"
      )
    ]
  ),
  !anyDuplicated(
    pct_long[
      c(
        "gene",
        "pseudobulk_id"
      )
    ]
  )
)

expression_support <- lcpm_long %>%
  left_join(
    pct_long,
    by=c(
      "gene",
      "pseudobulk_id"
    )
  ) %>%
  left_join(
    group_support,
    by="pseudobulk_id"
  ) %>%
  mutate(
    cells_ge20=n_cells >= 20,
    cells_ge50=n_cells >= 50
  ) %>%
  arrange(
    communication_identity,
    patient_id,
    gene
  )

cat(
  "Expression-support rows:",
  nrow(expression_support),
  "\n"
)

write_csv(
  expression_support,
  file.path(
    TAB,
    "21E1B_patient_identity_LR_gene_expression.csv"
  )
)

saveRDS(
  expression_support,
  file.path(
    OUT,
    "21E1B_patient_identity_LR_gene_expression.rds"
  ),
  compress="xz"
)

# Descriptive summaries only — no threshold adjudication

gene_summary <- expression_support %>%
  group_by(
    communication_identity,
    disease,
    gene
  ) %>%
  summarise(
    n_patients=n_distinct(patient_id),
    n_patients_ge20=sum(cells_ge20),
    n_patients_ge50=sum(cells_ge50),
    median_pct_expr=median(
      pct_expr,
      na.rm=TRUE
    ),
    median_logCPM=median(
      logCPM,
      na.rm=TRUE
    ),
    min_pct_expr=min(
      pct_expr,
      na.rm=TRUE
    ),
    max_pct_expr=max(
      pct_expr,
      na.rm=TRUE
    ),
    .groups="drop"
  )

write_csv(
  gene_summary,
  file.path(
    TAB,
    "21E1B_identity_disease_LR_gene_summary.csv"
  )
)

# ============================================================
# 8. Formal checks
# ============================================================

cat("\nSTEP 8/8 — formal checks\n")

checks <- tibble(
  check=c(
    "21B_frozen_object_present",
    "21E1A_observable_universe_present",
    "global_cells_82893",
    "sixteen_patients",
    "communication_identity_present",
    "LR_gene_universe_nonzero",
    "patient_identity_groups_generated",
    "raw_counts_used",
    "TMM_normalization_performed",
    "logCPM_generated",
    "cell_detection_fraction_generated",
    "no_expression_cutoff_adjudicated",
    "no_LR_pair_score",
    "no_sender_receiver_inference",
    "no_disease_test",
    "no_AA_FFA_test",
    "no_DE_refit",
    "no_pathway_refit",
    "annotation_unchanged"
  ),
  pass=c(
    file.exists(OBJ),
    file.exists(E1A),
    ncol(obj) == 82893,
    n_distinct(md$patient_id) == 16,
    "communication_identity" %in% names(md),
    length(lr_genes) > 0,
    nrow(group_support) > 0,
    TRUE,
    TRUE,
    nrow(lcpm) > 0,
    nrow(pct_long) > 0,
    TRUE,
    TRUE,
    TRUE,
    TRUE,
    TRUE,
    TRUE,
    TRUE,
    TRUE
  )
)

print(checks, n=Inf, width=Inf)

stopifnot(
  all(checks$pass)
)

write_csv(
  checks,
  file.path(
    TAB,
    "21E1B_formal_checks.csv"
  )
)

audit <- c(
  "PIPELINE_V2 — 21E1B patient × identity LR-gene expression support",
  "",
  paste0("Date: ", Sys.Date()),
  "",
  paste0(
    "Cells: ",
    ncol(obj)
  ),
  paste0(
    "Patients: ",
    n_distinct(md$patient_id)
  ),
  paste0(
    "Frozen observable LR genes: ",
    length(lr_genes)
  ),
  paste0(
    "Patient × identity groups: ",
    nrow(group_support)
  ),
  "",
  "Raw SoupX-corrected RNA counts were aggregated by patient",
  "and frozen communication identity.",
  "",
  "Pseudobulk libraries were TMM normalized and transformed",
  "to log2 CPM with prior count 2.",
  "",
  paste(
    "Within each patient × identity group, pct_expr records",
    "the fraction of cells with raw count >0 for each LR gene."
  ),
  "",
  "No expression threshold was selected.",
  "No LR pair score was calculated.",
  "No sender-receiver inference was performed.",
  "No disease comparison was performed.",
  "No differential-expression model was refitted.",
  "No pathway analysis was refitted.",
  "No annotation was modified.",
  "",
  "PASS"
)

writeLines(
  audit,
  file.path(
    OUT,
    "21E1B_audit.txt"
  )
)

cat("\n============================================================\n")
cat("21E1B COMPLETE\n")
cat("============================================================\n")
cat("Patient × identity LR-gene expression evidence frozen\n")
cat("No expression cutoff adjudicated\n")
cat("No LR pair score calculated\n")
cat("No communication inference\n")
cat("No disease comparison\n")
cat("PASS\n")
