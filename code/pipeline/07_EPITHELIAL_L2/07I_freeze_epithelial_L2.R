suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(readr)
})

cat("\n============================================================\n")
cat("PIPELINE_V2 — 07I FREEZE EPITHELIAL L2\n")
cat("============================================================\n\n")

IN <- "PIPELINE_V2/objects/annotated/epithelial_v2_07B_reclustered.rds"
MAP <- "PIPELINE_V2/results/annotation/07_EPITHELIAL_L2/07H_prefreeze/tables/07H_provisional_L2_adjudication_map.csv"
CHECKS <- "PIPELINE_V2/results/annotation/07_EPITHELIAL_L2/07H_prefreeze/tables/07H_prefreeze_checks.csv"

OUTOBJ <- "PIPELINE_V2/objects/annotated/epithelial_v2_L2_frozen.rds"
OUT <- "PIPELINE_V2/results/annotation/07_EPITHELIAL_L2/07I_frozen"
TDIR <- file.path(OUT,"tables")

dir.create(dirname(OUTOBJ),recursive=TRUE,showWarnings=FALSE)
dir.create(TDIR,recursive=TRUE,showWarnings=FALSE)

# ============================================================
# 1. LOAD + VERIFY PREFREEZE PASS
# ============================================================

cat("STEP 1/7 — loading canonical epithelial object\n")

x <- readRDS(IN)

stopifnot(
  ncol(x)==63893,
  nrow(x)==23412,
  "epi_snn_res.0.5" %in% colnames(x@meta.data),
  length(unique(x$patient_id))==16
)

canonical_cells <- colnames(x)

checks <- read_csv(
  CHECKS,
  show_col_types=FALSE
)

stopifnot(
  nrow(checks)>0,
  all(checks$pass)
)

cat("07H pre-freeze checks: ALL PASS\n")
cat("Cells:",ncol(x),"\n")
cat("Genes:",nrow(x),"\n")

# ============================================================
# 2. LOAD FROZEN ADJUDICATION MAP
# ============================================================

cat("\nSTEP 2/7 — loading adjudication map\n")

map <- read_csv(
  MAP,
  show_col_types=FALSE,
  col_types=cols(
    cluster=col_character(),
    L2_identity=col_character(),
    L2_state=col_character(),
    L2_confidence=col_character(),
    L2_flag=col_character()
  )
)

stopifnot(
  nrow(map)==19,
  !anyDuplicated(map$cluster),
  setequal(map$cluster,as.character(0:18)),
  !anyNA(map)
)

# Rename provisional terminology in output artifact
frozen_map <- map %>%
  mutate(
    annotation_status="FROZEN",
    annotation_version="V2_EPITHELIAL_L2_20260929",
    clustering_source="epi_snn_res.0.5",
    expression_source="SoupX_corrected"
  )

write_csv(
  frozen_map,
  file.path(TDIR,"07I_epithelial_L2_frozen_map.csv")
)

# ============================================================
# 3. ASSIGN FROZEN METADATA
# ============================================================

cat("STEP 3/7 — assigning frozen L2 metadata\n")

x$EPI_cluster_L2 <- as.character(x$`epi_snn_res.0.5`)

idx <- match(
  x$EPI_cluster_L2,
  map$cluster
)

stopifnot(!anyNA(idx))

x$L2_identity   <- map$L2_identity[idx]
x$L2_state      <- map$L2_state[idx]
x$L2_confidence <- map$L2_confidence[idx]
x$L2_flag       <- map$L2_flag[idx]

x$L2_frozen <- TRUE
x$L2_status <- "FROZEN"
x$L2_annotation_version <- "V2_EPITHELIAL_L2_20260929"
x$L2_clustering_source <- "epi_snn_res.0.5"
x$L2_expression_source <- "SoupX_corrected"

stopifnot(
  !anyNA(x$L2_identity),
  !anyNA(x$L2_state),
  !anyNA(x$L2_confidence),
  !anyNA(x$L2_flag)
)

# ============================================================
# 4. FINAL CELL-COUNT + IDENTITY AUDIT
# ============================================================

cat("STEP 4/7 — final cell-count audit\n")

identity_counts <- x@meta.data %>%
  count(
    L2_identity,
    name="n_cells"
  ) %>%
  mutate(
    pct=100*n_cells/sum(n_cells)
  ) %>%
  arrange(desc(n_cells))

state_counts <- x@meta.data %>%
  count(
    L2_state,
    name="n_cells"
  ) %>%
  mutate(
    pct=100*n_cells/sum(n_cells)
  ) %>%
  arrange(desc(n_cells))

stopifnot(
  sum(identity_counts$n_cells)==63893,
  sum(state_counts$n_cells)==63893,
  length(unique(x$L2_identity))==9,
  length(unique(x$L2_state))==4
)

cat("\n=== FROZEN L2 IDENTITIES ===\n")
print(identity_counts,row.names=FALSE)

cat("\n=== FROZEN L2 STATES ===\n")
print(state_counts,row.names=FALSE)

write_csv(
  identity_counts,
  file.path(TDIR,"07I_L2_identity_counts.csv")
)

write_csv(
  state_counts,
  file.path(TDIR,"07I_L2_state_counts.csv")
)

# ============================================================
# 5. PATIENT SUPPORT
# ============================================================

cat("\nSTEP 5/7 — final patient-support audit\n")

patient_identity <- x@meta.data %>%
  count(
    L2_identity,
    patient_id,
    disease,
    batch,
    name="n_cells"
  )

support <- patient_identity %>%
  group_by(L2_identity) %>%
  summarise(
    total_cells=sum(n_cells),
    n_patients=n_distinct(patient_id),
    n_HC=n_distinct(patient_id[disease=="HC"]),
    n_AA=n_distinct(patient_id[disease=="AA"]),
    n_FFA=n_distinct(patient_id[disease=="FFA"]),
    n_batches=n_distinct(batch),
    dominant_patient=
      patient_id[which.max(n_cells)],
    dominant_patient_fraction=
      max(n_cells)/sum(n_cells),
    .groups="drop"
  ) %>%
  arrange(desc(total_cells))

stopifnot(
  all(support$n_patients>=5),
  all(support$n_batches==2),
  all(support$dominant_patient_fraction<=0.50)
)

cat("\n=== FROZEN IDENTITY SUPPORT ===\n")
print(support,row.names=FALSE)

write_csv(
  patient_identity,
  file.path(TDIR,"07I_L2_identity_by_patient.csv")
)

write_csv(
  support,
  file.path(TDIR,"07I_L2_identity_support.csv")
)

# ============================================================
# 6. VERIFY GEOMETRY/CELLS UNCHANGED + SAVE
# ============================================================

cat("\nSTEP 6/7 — canonical integrity + save\n")

stopifnot(
  ncol(x)==63893,
  nrow(x)==23412,
  identical(
    colnames(x),
    canonical_cells
  )
)

# No reductions, neighbors or clusters recalculated here.
saveRDS(
  x,
  OUTOBJ,
  compress=FALSE
)

cat("Saved:",OUTOBJ,"\n")

# ============================================================
# 7. FINAL AUDIT
# ============================================================

cat("\nSTEP 7/7 — audit trail\n")

audit <- c(
  "PIPELINE_V2 — 07I epithelial L2 freeze",
  paste0("Generated: ",Sys.time()),
  "Input: epithelial_v2_07B_reclustered.rds",
  "Pre-freeze validation: 07H ALL PASS",
  "Cells: 63893",
  "Genes: 23412",
  "Patients: 16",
  "Working/frozen clustering resolution: 0.5",
  "Frozen anatomical identities: 9",
  "Frozen transcriptional states: 4",
  "Identity and state stored separately",
  "Expression source: SoupX-corrected RNA",
  "Disease was not used to define anatomical identity",
  "No cells removed",
  "No PCA recalculated",
  "No Harmony recalculated",
  "No neighbors recalculated",
  "No UMAP recalculated",
  "No clusters recalculated",
  "Cell identity/order preserved exactly",
  "C14 retained within IFE_Infundibulum_suprabasal with patient-dominant QC flag",
  "Ambiguous basal clusters retained as Basal_epithelial rather than overassigned to IFE or ORS",
  "Annotation version: V2_EPITHELIAL_L2_20260929",
  paste0("Frozen object: ",OUTOBJ)
)

writeLines(
  audit,
  file.path(TDIR,"07I_audit.txt")
)

cat("\n============================================================\n")
cat("07I COMPLETE\n")
cat("============================================================\n")
cat("Cells:",ncol(x),"\n")
cat("L2 identities:",length(unique(x$L2_identity)),"\n")
cat("L2 states:",length(unique(x$L2_state)),"\n")
cat("Output:",OUTOBJ,"\n")
cat("PASS — epithelial L2 annotation frozen.\n")
