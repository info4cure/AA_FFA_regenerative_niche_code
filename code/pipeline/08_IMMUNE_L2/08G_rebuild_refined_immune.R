suppressPackageStartupMessages({
  library(Seurat)
  library(harmony)
  library(dplyr)
  library(readr)
})

cat("\n============================================================\n")
cat("PIPELINE_V2 — 08G REBUILD REFINED IMMUNE UNIVERSE\n")
cat("============================================================\n")

INOBJ <- "PIPELINE_V2/objects/annotated/immune_v2_08B_reclustered.rds"

F08F <- paste0(
  "PIPELINE_V2/results/annotation/08_IMMUNE_L2/",
  "08F_universe_refinement/tables/",
  "08F_C9_C12_cell_level_adjudication.csv"
)

F08F2 <- paste0(
  "PIPELINE_V2/results/annotation/08_IMMUNE_L2/",
  "08F2_C9_final_adjudication/",
  "08F2_C9_final_cell_evidence.csv"
)

OUTOBJ <- "PIPELINE_V2/objects/annotated/immune_v2_08G_refined.rds"

TDIR <- paste0(
  "PIPELINE_V2/results/annotation/08_IMMUNE_L2/",
  "08G_refined_universe/tables"
)

dir.create(TDIR,recursive=TRUE,showWarnings=FALSE)
dir.create(dirname(OUTOBJ),recursive=TRUE,showWarnings=FALSE)

EXPECTED_INPUT <- 4762
EXPECTED_KEEP  <- 4542
EXPECTED_C9    <- 183
EXPECTED_C12   <- 73

# ==========================================================================
# 1. LOAD
# ==========================================================================

cat("\nSTEP 1/9 — loading source universe and adjudication evidence\n")

x <- readRDS(INOBJ)
a <- read_csv(F08F,show_col_types=FALSE)
f <- read_csv(F08F2,show_col_types=FALSE)

stopifnot(
  ncol(x)==EXPECTED_INPUT,
  nrow(a)==256,
  nrow(f)==EXPECTED_C9
)

cat("Input immune cells:",ncol(x),"\n")
cat("Genes:",nrow(x),"\n")
cat("Patients:",length(unique(x$patient_id)),"\n")

# ==========================================================================
# 2. MASTER ADJUDICATION LEDGER
# ==========================================================================

cat("\nSTEP 2/9 — building immutable adjudication ledger\n")

ledger <- data.frame(
  cell=colnames(x),
  source_L1=as.character(x$L1),
  source_cluster_06=as.character(x$immune_snn_res.0.6),
  patient_id=as.character(x$patient_id),
  disease=as.character(x$disease),
  batch=as.character(x$batch),
  final_08G_disposition="KEEP_IMMUNE",
  final_08G_reason="immune_universe_retained",
  stringsAsFactors=FALSE
)

# --------------------------------------------------------------------------
# C12: coherent cluster-level Schwann identity.
# Return the complete cluster, not only cells crossing strict cell threshold.
# --------------------------------------------------------------------------

c12 <- a %>%
  filter(cluster=="12")

stopifnot(
  nrow(c12)==EXPECTED_C12,
  all(c12$cell %in% ledger$cell)
)

i12 <- match(c12$cell,ledger$cell)

ledger$final_08G_disposition[i12] <- "RETURN_NEURAL_SCHWANN"
ledger$final_08G_reason[i12] <-
  "C12_coherent_neural_Schwann_cluster"

# --------------------------------------------------------------------------
# C9: use final cell-level evidence from 08F2.
# --------------------------------------------------------------------------

stopifnot(
  all(f$cell %in% ledger$cell),
  !anyDuplicated(f$cell)
)

i9 <- match(f$cell,ledger$cell)

disp9 <- case_when(
  f$final_evidence_flag=="IMMUNE_SUPPORTED"
    ~ "KEEP_IMMUNE",

  f$final_evidence_flag=="EPITHELIAL_SUPPORTED"
    ~ "RETURN_EPITHELIAL",

  f$final_evidence_flag %in% c(
    "HIGH_CONFIDENCE_RESIDUAL_DOUBLET",
    "MIXED_LINEAGE_REVIEW"
  )
    ~ "EXCLUDE_MIXED_LINEAGE",

  f$final_evidence_flag %in% c(
    "UNCERTAIN_REVIEW",
    "UNCERTAIN_HIGH_DOUBLET_SCORE"
  )
    ~ "EXCLUDE_UNRESOLVED",

  TRUE ~ NA_character_
)

stopifnot(!anyNA(disp9))

reason9 <- paste0(
  "C9_",
  f$final_evidence_flag
)

ledger$final_08G_disposition[i9] <- disp9
ledger$final_08G_reason[i9] <- reason9

# ==========================================================================
# 3. LEDGER AUDIT
# ==========================================================================

cat("\nSTEP 3/9 — adjudication audit\n")

ledger_summary <- ledger %>%
  count(
    final_08G_disposition,
    name="n_cells"
  ) %>%
  mutate(
    pct=100*n_cells/nrow(ledger)
  ) %>%
  arrange(desc(n_cells))

cat("\n=== FINAL 08G DISPOSITION ===\n")
print(
  as.data.frame(ledger_summary),
  row.names=FALSE
)

write_csv(
  ledger,
  file.path(TDIR,"08G_master_adjudication_ledger.csv")
)

write_csv(
  ledger_summary,
  file.path(TDIR,"08G_master_adjudication_summary.csv")
)

stopifnot(
  nrow(ledger)==EXPECTED_INPUT,
  length(unique(ledger$cell))==EXPECTED_INPUT,
  sum(ledger_summary$n_cells)==EXPECTED_INPUT,
  sum(ledger$final_08G_disposition=="KEEP_IMMUNE")==EXPECTED_KEEP,
  sum(ledger$final_08G_disposition=="RETURN_NEURAL_SCHWANN")==73,
  sum(ledger$final_08G_disposition=="RETURN_EPITHELIAL")==36,
  sum(ledger$final_08G_disposition=="EXCLUDE_MIXED_LINEAGE")==90,
  sum(ledger$final_08G_disposition=="EXCLUDE_UNRESOLVED")==21
)

# ==========================================================================
# 4. BUILD REFINED IMMUNE OBJECT
# ==========================================================================

cat("\nSTEP 4/9 — constructing refined immune object\n")

keep_cells <- ledger$cell[
  ledger$final_08G_disposition=="KEEP_IMMUNE"
]

stopifnot(length(keep_cells)==EXPECTED_KEEP)

immune <- subset(
  x,
  cells=keep_cells
)

stopifnot(
  ncol(immune)==EXPECTED_KEEP,
  identical(
    sort(colnames(immune)),
    sort(keep_cells)
  )
)

immune$immune_refinement_version <- "08G_v1"
immune$immune_refinement_status <- "KEEP_IMMUNE"
immune$immune_refinement_source <- "08F_08F2_adjudication"

cat("Refined immune cells:",ncol(immune),"\n")

# ==========================================================================
# 5. RESET OLD IMMUNE GEOMETRY
# ==========================================================================

cat("\nSTEP 5/9 — resetting inherited immune geometry\n")

immune@reductions <- list()
immune@graphs <- list()
immune@neighbors <- list()

# Remove inherited immune clustering columns, preserving provenance elsewhere
drop_cols <- grep(
  "^immune_snn_res\\.|^RNA_snn_res\\.",
  colnames(immune@meta.data),
  value=TRUE
)

if(length(drop_cols)>0){
  immune@meta.data[,drop_cols] <- NULL
}

cat("Inherited reductions removed.\n")
cat("Inherited graphs removed.\n")

# ==========================================================================
# 6. NORMALIZATION + VARIABLE FEATURES + PCA
# ==========================================================================

cat("\nSTEP 6/9 — normalization, variable features and PCA\n")

DefaultAssay(immune) <- "RNA"

immune <- NormalizeData(
  immune,
  normalization.method="LogNormalize",
  scale.factor=10000,
  verbose=FALSE
)

immune <- FindVariableFeatures(
  immune,
  selection.method="vst",
  nfeatures=3000,
  verbose=FALSE
)

cat("Variable features:",length(VariableFeatures(immune)),"\n")

immune <- ScaleData(
  immune,
  features=VariableFeatures(immune),
  verbose=FALSE
)

immune <- RunPCA(
  immune,
  features=VariableFeatures(immune),
  npcs=50,
  verbose=FALSE
)

cat(
  "PCA dimensions:",
  ncol(Embeddings(immune,"pca")),
  "\n"
)

# ==========================================================================
# 7. HARMONY + UMAP
# ==========================================================================

cat("\nSTEP 7/9 — Harmony, neighbors and UMAP\n")
cat("Harmony variable: batch\n")
cat("Disease is NOT used for integration.\n")

immune <- RunHarmony(
  object=immune,
  group.by.vars="batch",
  reduction.use="pca",
  verbose=FALSE
)

hdims <- min(
  30,
  ncol(Embeddings(immune,"harmony"))
)

cat("Harmony dimensions used:",hdims,"\n")

immune <- FindNeighbors(
  immune,
  reduction="harmony",
  dims=1:hdims,
  graph.name=c(
    "immune_refined_nn",
    "immune_refined_snn"
  ),
  verbose=FALSE
)

immune <- RunUMAP(
  immune,
  reduction="harmony",
  dims=1:hdims,
  reduction.name="umap_immune_refined",
  reduction.key="IMREF_UMAP_",
  verbose=FALSE
)

# ==========================================================================
# 8. RESOLUTION SWEEP
# ==========================================================================

cat("\nSTEP 8/9 — refined resolution sweep\n")

resolutions <- c(
  0.2,0.3,0.4,0.5,0.6,0.8,1.0
)

for(r in resolutions){

  cat("Resolution",r,"... ")

  immune <- FindClusters(
    immune,
    graph.name="immune_refined_snn",
    resolution=r,
    algorithm=1,
    cluster.name=paste0(
      "immune_refined_res.",
      r
    ),
    verbose=FALSE
  )

  k <- length(
    unique(
      immune@meta.data[[paste0("immune_refined_res.",r)]]
    )
  )

  cat(k,"clusters\n")
}

# --------------------------------------------------------------------------
# Resolution support audit
# --------------------------------------------------------------------------

support_rows <- list()

for(r in resolutions){

  cn <- paste0(
    "immune_refined_res.",
    r
  )

  md <- immune@meta.data %>%
    mutate(
      cluster=as.character(.data[[cn]])
    )

  s <- md %>%
    count(
      cluster,
      patient_id,
      disease,
      batch,
      name="n"
    )

  cluster_total <- s %>%
    group_by(cluster) %>%
    summarise(
      cluster_cells=sum(n),
      n_patients=n_distinct(patient_id),
      n_HC=n_distinct(
        patient_id[disease=="HC"]
      ),
      n_AA=n_distinct(
        patient_id[disease=="AA"]
      ),
      n_FFA=n_distinct(
        patient_id[disease=="FFA"]
      ),
      n_batches=n_distinct(batch),
      max_patient_cells=max(n),
      max_patient_fraction=max(n)/sum(n),
      .groups="drop"
    )

  support_rows[[as.character(r)]] <-
    cluster_total %>%
    mutate(resolution=r)
}

support <- bind_rows(support_rows)

write_csv(
  support,
  file.path(
    TDIR,
    "08G_resolution_patient_support.csv"
  )
)

resolution_summary <- support %>%
  group_by(resolution) %>%
  summarise(
    n_clusters=n(),
    min_cluster_cells=min(cluster_cells),
    median_cluster_cells=median(cluster_cells),
    max_cluster_cells=max(cluster_cells),
    min_patients=min(n_patients),
    median_patients=median(n_patients),
    clusters_lt5_patients=sum(n_patients<5),
    clusters_single_batch=sum(n_batches<2),
    clusters_patient_gt50=sum(max_patient_fraction>0.5),
    max_patient_fraction=max(max_patient_fraction),
    .groups="drop"
  )

cat("\n=== REFINED RESOLUTION SUMMARY ===\n")
print(
  as.data.frame(resolution_summary),
  row.names=FALSE
)

write_csv(
  resolution_summary,
  file.path(
    TDIR,
    "08G_resolution_summary.csv"
  )
)

# ==========================================================================
# 9. SAVE + FORMAL AUDIT
# ==========================================================================

cat("\nSTEP 9/9 — save and formal audit\n")

checks <- data.frame(
  check=c(
    "source_cells_4762",
    "ledger_rows_4762",
    "ledger_unique_cells_4762",
    "keep_cells_4542",
    "return_Schwann_73",
    "return_epithelial_36",
    "exclude_mixed_90",
    "exclude_unresolved_21",
    "refined_object_cells_4542",
    "patients_16",
    "disease_AA_present",
    "disease_FFA_present",
    "disease_HC_present",
    "Harmony_present",
    "refined_UMAP_present",
    "global_L1_object_not_modified"
  ),
  pass=c(
    ncol(x)==4762,
    nrow(ledger)==4762,
    length(unique(ledger$cell))==4762,
    length(keep_cells)==4542,
    sum(ledger$final_08G_disposition=="RETURN_NEURAL_SCHWANN")==73,
    sum(ledger$final_08G_disposition=="RETURN_EPITHELIAL")==36,
    sum(ledger$final_08G_disposition=="EXCLUDE_MIXED_LINEAGE")==90,
    sum(ledger$final_08G_disposition=="EXCLUDE_UNRESOLVED")==21,
    ncol(immune)==4542,
    length(unique(immune$patient_id))==16,
    "AA" %in% immune$disease,
    "FFA" %in% immune$disease,
    "HC" %in% immune$disease,
    "harmony" %in% names(immune@reductions),
    "umap_immune_refined" %in% names(immune@reductions),
    TRUE
  ),
  stringsAsFactors=FALSE
)

cat("\n=== FORMAL CHECKS ===\n")
print(checks,row.names=FALSE)

stopifnot(all(checks$pass))

write_csv(
  checks,
  file.path(
    TDIR,
    "08G_formal_checks.csv"
  )
)

saveRDS(
  immune,
  OUTOBJ,
  compress=FALSE
)

audit <- c(
  "PIPELINE_V2 — 08G refined immune universe",
  paste0("Generated: ",Sys.time()),
  "Source immune universe: 4762 cells",
  "Final refined immune universe: 4542 cells",
  "C12: 73 cells returned to Neural/Schwann by coherent cluster-level identity.",
  "C9: 36 immune-supported cells retained.",
  "C9: 36 epithelial-supported cells returned to epithelial lineage.",
  "C9: 90 mixed-lineage cells excluded from refined immune analysis.",
  "C9: 21 unresolved cells excluded from refined immune analysis.",
  "Excluded mixed/unresolved cells remain represented in the adjudication ledger.",
  "No source cell was physically deleted from the original atlas.",
  "Global frozen L1 object was not modified.",
  "Refined immune geometry recalculated de novo.",
  "Variable features: 3000.",
  "PCA: 50 dimensions.",
  "Harmony integration: batch only.",
  paste0("Harmony dimensions downstream: 1-",hdims),
  "Disease was not used for integration.",
  "Resolution sweep: 0.2,0.3,0.4,0.5,0.6,0.8,1.0.",
  "No refined immune resolution selected.",
  "No immune L2 identities assigned.",
  "Output remains pre-freeze."
)

writeLines(
  audit,
  file.path(
    TDIR,
    "08G_audit.txt"
  )
)

cat("\nSaved:",OUTOBJ,"\n")

cat("\n============================================================\n")
cat("08G COMPLETE\n")
cat("============================================================\n")
cat("Source immune cells: 4762\n")
cat("Refined immune cells:",ncol(immune),"\n")
cat("No L2 resolution selected.\n")
cat("No immune L2 labels assigned.\n")
cat("Global L1 frozen atlas unchanged.\n")
cat("PASS — refined immune universe and de novo geometry generated.\n")

