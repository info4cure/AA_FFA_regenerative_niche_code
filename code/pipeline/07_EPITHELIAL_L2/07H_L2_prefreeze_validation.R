suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(readr)
})

cat("\n============================================================\n")
cat("PIPELINE_V2 — 07H EPITHELIAL L2 PRE-FREEZE VALIDATION\n")
cat("============================================================\n\n")

IN <- "PIPELINE_V2/objects/annotated/epithelial_v2_07B_reclustered.rds"

OUT <- "PIPELINE_V2/results/annotation/07_EPITHELIAL_L2/07H_prefreeze"
TDIR <- file.path(OUT,"tables")
dir.create(TDIR,recursive=TRUE,showWarnings=FALSE)

x <- readRDS(IN)
CC <- "epi_snn_res.0.5"

stopifnot(
  ncol(x) == 63893,
  nrow(x) == 23412,
  CC %in% colnames(x@meta.data),
  length(unique(x$patient_id)) == 16
)

x$epi_cluster_05 <- as.character(x@meta.data[[CC]])

stopifnot(
  setequal(unique(x$epi_cluster_05),as.character(0:18))
)

cat("Cells:",ncol(x),"\n")
cat("Patients:",length(unique(x$patient_id)),"\n")
cat("Working resolution: 0.5\n")
cat("Clusters:",length(unique(x$epi_cluster_05)),"\n")

# ============================================================
# 1. PROVISIONAL ADJUDICATION MAP
#
# Identity = anatomical identity.
# State = orthogonal transcriptional state.
# Flag = QC/adjudication issue, NOT biology.
# ============================================================

map <- data.frame(
  cluster=as.character(0:18),

  L2_identity=c(
    "IFE_Infundibulum_suprabasal",      # 0
    "ORS_basal",                         # 1
    "IFE_Infundibulum_suprabasal",      # 2
    "IFE_Infundibulum_suprabasal",      # 3
    "Basal_epithelial",                  # 4
    "ORS",                               # 5
    "IFE_Infundibulum_differentiating", # 6
    "ORS_basal",                         # 7
    "Basal_epithelial",                  # 8
    "Basal_epithelial",                  # 9
    "Bulge_Isthmus",                     # 10
    "IFE_Infundibulum_suprabasal",      # 11
    "Bulge_Isthmus",                     # 12
    "Basal_epithelial",                  # 13
    "IFE_Infundibulum_suprabasal",      # 14
    "Hair_shaft",                        # 15
    "Bulge_Isthmus",                     # 16
    "IFE_Infundibulum_terminal",         # 17
    "IRS_companion"                      # 18
  ),

  L2_state=c(
    "Baseline",       # 0
    "Baseline",       # 1
    "Baseline",       # 2
    "Stress_high",    # 3
    "Baseline",       # 4
    "Stress_high",    # 5
    "Baseline",       # 6
    "Baseline",       # 7
    "Cycling",        # 8
    "Baseline",       # 9
    "Baseline",       # 10
    "IFN_MHCII_high", # 11
    "Baseline",       # 12
    "Baseline",       # 13
    "Baseline",       # 14
    "Baseline",       # 15
    "Baseline",       # 16
    "Stress_high",    # 17
    "Cycling"         # 18
  ),

  L2_confidence=c(
    "HIGH",     # 0
    "HIGH",     # 1
    "HIGH",     # 2
    "HIGH",     # 3
    "MODERATE", # 4
    "HIGH",     # 5
    "HIGH",     # 6
    "HIGH",     # 7
    "MODERATE", # 8
    "MODERATE", # 9
    "HIGH",     # 10
    "HIGH",     # 11
    "HIGH",     # 12
    "MODERATE", # 13
    "LOW",      # 14
    "HIGH",     # 15
    "HIGH",     # 16
    "HIGH",     # 17
    "HIGH"      # 18
  ),

  L2_flag=c(
    "none", # 0
    "none", # 1
    "none", # 2
    "none", # 3
    "anatomical_boundary_IFE_vs_ORS", # 4
    "none", # 5
    "none", # 6
    "none", # 7
    "anatomical_boundary_IFE_vs_ORS", # 8
    "anatomical_boundary_IFE_vs_ORS", # 9
    "none", # 10
    "patient_enriched_state_cluster", # 11
    "none", # 12
    "anatomical_boundary_IFE_vs_ORS", # 13
    "patient_dominant_cluster_ARE096",# 14
    "none", # 15
    "none", # 16
    "none", # 17
    "none"  # 18
  ),

  stringsAsFactors=FALSE
)

stopifnot(
  nrow(map)==19,
  !anyDuplicated(map$cluster),
  !anyNA(map)
)

write_csv(
  map,
  file.path(TDIR,"07H_provisional_L2_adjudication_map.csv")
)

cat("\n=== PROVISIONAL MAP ===\n")
print(map,row.names=FALSE)

# ============================================================
# 2. APPLY PROVISIONAL MAP IN MEMORY ONLY
# ============================================================

idx <- match(x$epi_cluster_05,map$cluster)

stopifnot(!anyNA(idx))

x$L2_identity_provisional   <- map$L2_identity[idx]
x$L2_state_provisional      <- map$L2_state[idx]
x$L2_confidence_provisional <- map$L2_confidence[idx]
x$L2_flag_provisional       <- map$L2_flag[idx]

stopifnot(
  length(x$L2_identity_provisional)==ncol(x),
  !anyNA(x$L2_identity_provisional),
  !anyNA(x$L2_state_provisional)
)

# ============================================================
# 3. CELL-COUNT INVARIANT
# ============================================================

cat("\nSTEP 1/7 — cell-count invariant\n")

identity_counts <- x@meta.data %>%
  count(
    L2_identity_provisional,
    name="n_cells"
  ) %>%
  mutate(
    pct=100*n_cells/sum(n_cells)
  ) %>%
  arrange(desc(n_cells))

state_counts <- x@meta.data %>%
  count(
    L2_state_provisional,
    name="n_cells"
  ) %>%
  mutate(
    pct=100*n_cells/sum(n_cells)
  ) %>%
  arrange(desc(n_cells))

stopifnot(
  sum(identity_counts$n_cells)==63893,
  sum(state_counts$n_cells)==63893
)

cat("\n=== L2 IDENTITY COUNTS ===\n")
print(identity_counts,row.names=FALSE)

cat("\n=== L2 STATE COUNTS ===\n")
print(state_counts,row.names=FALSE)

write_csv(
  identity_counts,
  file.path(TDIR,"07H_identity_counts.csv")
)

write_csv(
  state_counts,
  file.path(TDIR,"07H_state_counts.csv")
)

# ============================================================
# 4. PATIENT SUPPORT PER FINAL IDENTITY
# ============================================================

cat("\nSTEP 2/7 — patient support by identity\n")

patient_identity <- x@meta.data %>%
  count(
    L2_identity_provisional,
    patient_id,
    disease,
    batch,
    name="n_cells"
  )

identity_support <- patient_identity %>%
  group_by(L2_identity_provisional) %>%
  summarise(
    total_cells=sum(n_cells),
    n_patients=n_distinct(patient_id),
    n_HC=n_distinct(patient_id[disease=="HC"]),
    n_AA=n_distinct(patient_id[disease=="AA"]),
    n_FFA=n_distinct(patient_id[disease=="FFA"]),
    n_batches=n_distinct(batch),

    dominant_patient=
      patient_id[which.max(n_cells)],

    dominant_patient_cells=max(n_cells),

    dominant_patient_fraction=
      max(n_cells)/sum(n_cells),

    .groups="drop"
  ) %>%
  arrange(desc(total_cells))

cat("\n=== FINAL IDENTITY PATIENT SUPPORT ===\n")
print(identity_support,row.names=FALSE)

write_csv(
  patient_identity,
  file.path(TDIR,"07H_identity_by_patient.csv")
)

write_csv(
  identity_support,
  file.path(TDIR,"07H_identity_patient_support.csv")
)

# ============================================================
# 5. DISEASE + BATCH COMPOSITION
#
# Descriptive only. NOT differential abundance.
# ============================================================

cat("\nSTEP 3/7 — disease/batch composition\n")

identity_disease <- x@meta.data %>%
  count(
    L2_identity_provisional,
    disease,
    name="n_cells"
  ) %>%
  group_by(L2_identity_provisional) %>%
  mutate(
    fraction=n_cells/sum(n_cells)
  ) %>%
  ungroup()

identity_batch <- x@meta.data %>%
  count(
    L2_identity_provisional,
    batch,
    name="n_cells"
  ) %>%
  group_by(L2_identity_provisional) %>%
  mutate(
    fraction=n_cells/sum(n_cells)
  ) %>%
  ungroup()

write_csv(
  identity_disease,
  file.path(TDIR,"07H_identity_by_disease.csv")
)

write_csv(
  identity_batch,
  file.path(TDIR,"07H_identity_by_batch.csv")
)

# ============================================================
# 6. STATE × IDENTITY
# ============================================================

cat("\nSTEP 4/7 — state within identity\n")

state_identity <- x@meta.data %>%
  count(
    L2_identity_provisional,
    L2_state_provisional,
    name="n_cells"
  ) %>%
  group_by(L2_identity_provisional) %>%
  mutate(
    fraction_within_identity=
      n_cells/sum(n_cells)
  ) %>%
  ungroup() %>%
  arrange(
    L2_identity_provisional,
    desc(n_cells)
  )

cat("\n=== STATE × IDENTITY ===\n")
print(state_identity,row.names=FALSE)

write_csv(
  state_identity,
  file.path(TDIR,"07H_state_by_identity.csv")
)

# ============================================================
# 7. STATE PATIENT SUPPORT
# ============================================================

cat("\nSTEP 5/7 — patient support by state\n")

patient_state <- x@meta.data %>%
  count(
    L2_state_provisional,
    patient_id,
    disease,
    batch,
    name="n_cells"
  )

state_support <- patient_state %>%
  group_by(L2_state_provisional) %>%
  summarise(
    total_cells=sum(n_cells),
    n_patients=n_distinct(patient_id),
    n_HC=n_distinct(patient_id[disease=="HC"]),
    n_AA=n_distinct(patient_id[disease=="AA"]),
    n_FFA=n_distinct(patient_id[disease=="FFA"]),
    dominant_patient=
      patient_id[which.max(n_cells)],
    dominant_patient_fraction=
      max(n_cells)/sum(n_cells),
    .groups="drop"
  ) %>%
  arrange(desc(total_cells))

cat("\n=== STATE PATIENT SUPPORT ===\n")
print(state_support,row.names=FALSE)

write_csv(
  patient_state,
  file.path(TDIR,"07H_state_by_patient.csv")
)

write_csv(
  state_support,
  file.path(TDIR,"07H_state_patient_support.csv")
)

# ============================================================
# 8. FLAGS / CONFIDENCE
# ============================================================

cat("\nSTEP 6/7 — confidence and flags\n")

flag_counts <- x@meta.data %>%
  count(
    epi_cluster_05,
    L2_identity_provisional,
    L2_state_provisional,
    L2_confidence_provisional,
    L2_flag_provisional,
    name="n_cells"
  ) %>%
  arrange(as.integer(epi_cluster_05))

cat("\n=== CLUSTER ADJUDICATION AUDIT ===\n")
print(flag_counts,row.names=FALSE)

write_csv(
  flag_counts,
  file.path(TDIR,"07H_cluster_adjudication_audit.csv")
)

# ============================================================
# 9. FORMAL PRE-FREEZE CHECKS
# ============================================================

cat("\nSTEP 7/7 — formal pre-freeze checks\n")

checks <- data.frame(
  check=c(
    "all_cells_assigned_identity",
    "all_cells_assigned_state",
    "identity_counts_sum_63893",
    "state_counts_sum_63893",
    "all_identities_multibatch",
    "all_identities_ge5_patients",
    "no_identity_patient_gt50pct",
    "hair_shaft_ge5_patients",
    "IRS_companion_ge5_patients"
  ),

  pass=c(
    !anyNA(x$L2_identity_provisional),

    !anyNA(x$L2_state_provisional),

    sum(identity_counts$n_cells)==63893,

    sum(state_counts$n_cells)==63893,

    all(identity_support$n_batches>=2),

    all(identity_support$n_patients>=5),

    all(
      identity_support$dominant_patient_fraction <= 0.50
    ),

    identity_support$n_patients[
      identity_support$L2_identity_provisional==
        "Hair_shaft"
    ] >= 5,

    identity_support$n_patients[
      identity_support$L2_identity_provisional==
        "IRS_companion"
    ] >= 5
  )
)

cat("\n=== PRE-FREEZE CHECKS ===\n")
print(checks,row.names=FALSE)

write_csv(
  checks,
  file.path(TDIR,"07H_prefreeze_checks.csv")
)

# ============================================================
# 10. SPECIAL C14 ABSORPTION CHECK
# ============================================================

c14_cells <- sum(x$epi_cluster_05=="14")

suprabasal_cells <- sum(
  x$L2_identity_provisional==
    "IFE_Infundibulum_suprabasal"
)

c14_fraction_of_identity <-
  c14_cells/suprabasal_cells

c14_check <- data.frame(
  cluster="14",
  cluster_cells=c14_cells,
  merged_identity=
    "IFE_Infundibulum_suprabasal",
  merged_identity_cells=
    suprabasal_cells,
  cluster_fraction_of_merged_identity=
    c14_fraction_of_identity
)

cat("\n=== CLUSTER 14 ABSORPTION CHECK ===\n")
print(c14_check,row.names=FALSE)

write_csv(
  c14_check,
  file.path(TDIR,"07H_cluster14_absorption_check.csv")
)

# ============================================================
# 11. FINAL STATUS
# ============================================================

all_pass <- all(checks$pass)

audit <- c(
  "PIPELINE_V2 — 07H epithelial L2 pre-freeze validation",
  paste0("Generated: ",Sys.time()),
  "Input: epithelial_v2_07B_reclustered.rds",
  "Working resolution: 0.5",
  "Cells: 63893",
  "Clusters: 19",
  paste0(
    "Provisional anatomical identities: ",
    length(unique(x$L2_identity_provisional))
  ),
  paste0(
    "Provisional states: ",
    length(unique(x$L2_state_provisional))
  ),
  "Identity and state remain separate.",
  "Patient-enrichment retained as QC flag rather than biological state.",
  "Disease/batch composition is descriptive only.",
  "No differential abundance performed.",
  "No geometry recalculated.",
  "No frozen object written.",
  paste0(
    "All formal pre-freeze checks passed: ",
    all_pass
  )
)

writeLines(
  audit,
  file.path(TDIR,"07H_audit.txt")
)

cat("\n============================================================\n")
cat("07H COMPLETE\n")
cat("============================================================\n")

if(all_pass){
  cat("PASS — all formal pre-freeze checks passed.\n")
  cat("No frozen object written.\n")
} else {
  cat("REVIEW REQUIRED — one or more checks failed.\n")
  cat("No frozen object written.\n")
}
