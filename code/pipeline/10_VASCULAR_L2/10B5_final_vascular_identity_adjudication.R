suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(readr)
})

options(width = 220)

cat("============================================================\n")
cat("PIPELINE_V2 — 10B5 FINAL VASCULAR IDENTITY ADJUDICATION\n")
cat("============================================================\n")

# ============================================================
# Paths
# ============================================================

SRC <- paste0(
  "PIPELINE_V2/objects/annotated/",
  "vascular_v2_10B4B_rebuilt_3023.rds"
)

OUT <- paste0(
  "PIPELINE_V2/results/annotation/",
  "10_VASCULAR_L2/10B5_final_identity_adjudication"
)

TAB <- file.path(OUT, "tables")

dir.create(
  TAB,
  recursive = TRUE,
  showWarnings = FALSE
)

CL <- "vascular_postC8_snn_res.0.6"

# ============================================================
# 1. Load
# ============================================================

cat("\nSTEP 1/9 — loading post-C8 vascular universe\n")

x <- readRDS(SRC)

stopifnot(ncol(x) == 3023)
stopifnot(CL %in% colnames(x@meta.data))
stopifnot(length(unique(x$patient_id)) == 16)
stopifnot(length(unique(x$batch)) == 2)

DefaultAssay(x) <- "RNA"

cl <- as.character(x@meta.data[[CL]])

cat("Cells:", ncol(x), "\n")
cat("Patients:", length(unique(x$patient_id)), "\n")
cat("Clusters:", length(unique(cl)), "\n")

stopifnot(length(unique(cl)) == 10)

# ============================================================
# 2. Prespecified vascular programmes
# ============================================================

cat("\nSTEP 2/9 — defining vascular identity/state programmes\n")

programmes <- list(

  ENDOTHELIAL_CORE = c(
    "PECAM1","VWF","CDH5","EMCN",
    "KDR","ESAM","ENG","RAMP2"
  ),

  CAPILLARY = c(
    "CA4","RGCC","BTNL9","EDNRB",
    "EMCN","KDR","PLVAP","CA2"
  ),

  VENOUS_ACKR1 = c(
    "ACKR1","SELP","SELE",
    "VCAM1","ICAM1","VWF"
  ),

  POSTCAPILLARY_VENULAR = c(
    "ACKR1","SELE","SELP",
    "ICAM1","VCAM1","PLVAP"
  ),

  ARTERIAL = c(
    "EFNB2","SOX17","GJA5","GJA4",
    "HEY1","HEY2","BMX","FBLN5","SEMA3G"
  ),

  LYMPHATIC = c(
    "PROX1","PDPN","FLT4","CCL21",
    "LYVE1","NRP2","RELN"
  ),

  ANGIOGENIC = c(
    "KDR","ESM1","ANGPT2","APLN",
    "PGF","PLVAP","CXCR4"
  ),

  INFLAMMATORY = c(
    "ICAM1","VCAM1","SELE",
    "CXCL9","CXCL10","CXCL11",
    "HLA-DRA","HLA-DRB1","CCL2"
  ),

  MURAL_CONTRACTILE = c(
    "RGS5","CSPG4","MCAM","PDGFRB",
    "RBP1","NOTCH3","ACTA2","TAGLN",
    "MYH11","CNN1"
  ),

  MESENCHYMAL_ECM = c(
    "DCN","LUM","COL1A1","COL1A2",
    "COL3A1","DPT","FN1","POSTN"
  ),

  PERIVASCULAR_EBF2 = c(
    "EBF2","TBX18","TBX15","AGTR1",
    "CARMN","DPT","PRR16","BMP5"
  ),

  EPITHELIAL = c(
    "EPCAM","KRT5","KRT14","KRT17",
    "KRT19","KRT1","KRT10",
    "DSP","PKP3","CLDN1"
  ),

  IMMUNE = c(
    "PTPRC","CD3D","CD3E",
    "TRAC","NKG7","GNLY"
  )
)

available <- lapply(
  programmes,
  function(g)
    intersect(g, rownames(x))
)

coverage <- tibble(
  programme = names(programmes),
  n_requested = sapply(programmes, length),
  n_present = sapply(available, length),
  genes_present = sapply(
    available,
    paste,
    collapse = ";"
  )
)

print(coverage, n = Inf)

write_csv(
  coverage,
  file.path(
    TAB,
    "10B5_programme_coverage.csv"
  )
)

# ============================================================
# 3. Cell-level programme scores
# ============================================================

cat("\nSTEP 3/9 — calculating cell-level programme scores\n")

all_genes <- unique(
  unlist(available)
)

expr <- FetchData(
  x,
  vars = all_genes,
  layer = "data"
)

score_mean <- function(g) {

  g <- intersect(
    g,
    colnames(expr)
  )

  if (length(g) == 0)
    return(rep(NA_real_, nrow(expr)))

  rowMeans(
    expr[, g, drop = FALSE]
  )
}

detect_n <- function(g) {

  g <- intersect(
    g,
    colnames(expr)
  )

  if (length(g) == 0)
    return(rep(0L, nrow(expr)))

  rowSums(
    expr[, g, drop = FALSE] > 0
  )
}

cell <- tibble(
  cell = rownames(expr),
  cluster = cl[
    match(rownames(expr), colnames(x))
  ],
  patient_id =
    as.character(
      x$patient_id[
        match(rownames(expr), colnames(x))
      ]
    ),
  batch =
    as.character(
      x$batch[
        match(rownames(expr), colnames(x))
      ]
    ),
  doublet_score =
    x$doublet_score[
      match(rownames(expr), colnames(x))
    ],
  nFeature_RNA =
    x$nFeature_RNA[
      match(rownames(expr), colnames(x))
    ],
  nCount_RNA =
    x$nCount_RNA[
      match(rownames(expr), colnames(x))
    ]
)

for (p in names(available)) {

  cell[[paste0(p, "_score")]] <-
    score_mean(available[[p]])

  cell[[paste0(p, "_n")]] <-
    detect_n(available[[p]])
}

write_csv(
  cell,
  file.path(
    TAB,
    "10B5_cell_programme_scores.csv"
  )
)

# ============================================================
# 4. Cluster-level programme evidence
# ============================================================

cat("\nSTEP 4/9 — cluster-level programme evidence\n")

score_cols <- paste0(
  names(programmes),
  "_score"
)

cluster_scores <- cell %>%
  select(
    cluster,
    all_of(score_cols)
  ) %>%
  pivot_longer(
    cols = all_of(score_cols),
    names_to = "programme",
    values_to = "score"
  ) %>%
  mutate(
    programme =
      sub("_score$", "", programme)
  ) %>%
  group_by(
    cluster,
    programme
  ) %>%
  summarise(
    mean_score =
      mean(score, na.rm = TRUE),
    median_score =
      median(score, na.rm = TRUE),
    .groups = "drop"
  )

write_csv(
  cluster_scores,
  file.path(
    TAB,
    "10B5_cluster_programme_scores.csv"
  )
)

cat("\n=== CLUSTER PROGRAMME SCORES ===\n")

print(
  cluster_scores %>%
    arrange(
      as.integer(cluster),
      desc(mean_score)
    ),
  n = Inf
)

# ============================================================
# 5. Patient-level recurrence
# ============================================================

cat("\nSTEP 5/9 — patient-level programme support\n")

patient_scores <- cell %>%
  select(
    cluster,
    patient_id,
    batch,
    all_of(score_cols)
  ) %>%
  pivot_longer(
    cols = all_of(score_cols),
    names_to = "programme",
    values_to = "score"
  ) %>%
  mutate(
    programme =
      sub("_score$", "", programme)
  ) %>%
  group_by(
    cluster,
    patient_id,
    batch,
    programme
  ) %>%
  summarise(
    mean_score =
      mean(score, na.rm = TRUE),
    .groups = "drop"
  )

write_csv(
  patient_scores,
  file.path(
    TAB,
    "10B5_patient_programme_scores.csv"
  )
)

support <- cell %>%
  group_by(cluster) %>%
  summarise(
    n_cells = n(),
    n_patients =
      n_distinct(patient_id),
    n_batches =
      n_distinct(batch),

    dominant_patient_fraction = {
      tt <- table(patient_id)
      max(tt) / sum(tt)
    },

    median_doublet =
      median(doublet_score),

    median_features =
      median(nFeature_RNA),

    .groups = "drop"
  )

cat("\n=== STRUCTURAL SUPPORT ===\n")
print(support, n = Inf)

write_csv(
  support,
  file.path(
    TAB,
    "10B5_cluster_structural_support.csv"
  )
)

# ============================================================
# 6. Recover de novo r0.6 markers
# ============================================================

cat("\nSTEP 6/9 — calculating final r0.6 de novo markers\n")

Idents(x) <- factor(
  cl,
  levels = as.character(0:9)
)

markers <- FindAllMarkers(
  x,
  assay = "RNA",
  only.pos = TRUE,
  test.use = "wilcox",
  min.pct = 0.10,
  logfc.threshold = 0.25,
  verbose = FALSE
)

fc_col <- intersect(
  c("avg_log2FC","avg_logFC"),
  colnames(markers)
)[1]

stopifnot(!is.na(fc_col))

markers <- markers %>%
  mutate(
    cluster =
      as.character(cluster),
    gene =
      rownames(.)
  )

write_csv(
  markers,
  file.path(
    TAB,
    "10B5_r06_all_positive_markers.csv"
  )
)

top20 <- markers %>%
  filter(p_val_adj < 0.05) %>%
  group_by(cluster) %>%
  arrange(
    desc(.data[[fc_col]]),
    .by_group = TRUE
  ) %>%
  slice_head(n = 20) %>%
  ungroup()

write_csv(
  top20,
  file.path(
    TAB,
    "10B5_r06_top20_markers.csv"
  )
)

cat("\n=== TOP MARKERS ===\n")

print(
  top20 %>%
    select(
      cluster,
      gene,
      all_of(fc_col),
      pct.1,
      pct.2,
      p_val_adj
    ),
  n = Inf
)

# ============================================================
# 7. Final prespecified adjudication ledger
# ============================================================

cat("\nSTEP 7/9 — final hierarchical adjudication\n")

ledger <- tribble(
  ~cluster,
  ~vascular_family,
  ~vascular_state,
  ~vascular_L2_label,
  ~identity_basis,
  ~confidence,

  "0",
  "Venular endothelial",
  "ACKR1/CCL15/CCL23",
  "ACKR1+ venular EC",
  paste0(
    "ACKR1/SELP venular programme with ",
    "CCL15/CCL23/C7 enrichment"
  ),
  "HIGH",

  "1",
  "Capillary endothelial",
  "CA4/RGCC/BTNL9",
  "CA4/RGCC capillary EC",
  paste0(
    "CA4/RGCC/BTNL9/EDNRB capillary ",
    "programme"
  ),
  "HIGH",

  "2",
  "Venular endothelial",
  "Inflammatory/activated",
  "SELE/CXCL10 activated venular EC",
  paste0(
    "SELE/ICAM1/VCAM1/ACKR1 with ",
    "CXCL10/IL6/TNFAIP3 inflammatory ",
    "activation"
  ),
  "HIGH",

  "3",
  "Endothelial",
  "SLC10A1/CD83/MSX1",
  "SLC10A1/CD83 endothelial state",
  paste0(
    "Strong endothelial backbone; ",
    "distinct SLC10A1/CD83/MSX1 ",
    "de novo programme without dominant ",
    "mural or epithelial conflict"
  ),
  "MODERATE",

  "4",
  "Arterial endothelial",
  "GJA5/FBLN5/SEMA3G",
  "GJA5/FBLN5 arterial EC",
  paste0(
    "GJA5/FBLN5/HEY2/SEMA3G arterial ",
    "programme"
  ),
  "HIGH",

  "5",
  "Lymphatic endothelial",
  "CCL21/LYVE1/RELN",
  "CCL21/LYVE1 lymphatic EC",
  paste0(
    "CCL21/PDPN/LYVE1/RELN lymphatic ",
    "programme"
  ),
  "HIGH",

  "6",
  "Endothelial",
  "Low-complexity",
  "Low-complexity EC",
  paste0(
    "Recurrent endothelial population ",
    "dominated by mitochondrial/low-",
    "complexity transcriptional signal; ",
    "retained as QC-associated state and ",
    "not interpreted as a distinct ",
    "vascular lineage"
  ),
  "MODERATE",

  "7",
  "Endothelial",
  "ADGRG6/BMPER",
  "ADGRG6/BMPER endothelial state",
  paste0(
    "Strong endothelial backbone with ",
    "ADGRG6/BMPER/RAB3C/PLCXD3 de novo ",
    "programme and minimal competing-",
    "lineage evidence"
  ),
  "MODERATE",

  "8",
  "Endothelial",
  "Mural-like/contractile",
  "Mural-like contractile EC",
  paste0(
    "Persistent endothelial backbone ",
    "with recurrent RGS5/PDGFRB/TAGLN/",
    "ACTA2 and ECM coexpression; signal ",
    "persists in lowest doublet-score ",
    "quartile and across patients"
  ),
  "HIGH",

  "9",
  "Endothelial",
  "Epithelial-like",
  "Epithelial-like EC",
  paste0(
    "Persistent endothelial backbone ",
    "with recurrent KRT5/KRT14/KRT10/",
    "DSP epithelial programme; ",
    "coexpression persists in lowest ",
    "doublet-score quartile and across ",
    "patients"
  ),
  "MODERATE"
)

ledger <- ledger %>%
  left_join(
    support,
    by = "cluster"
  )

stopifnot(nrow(ledger) == 10)
stopifnot(all(as.character(0:9) %in% ledger$cluster))

cat("\n=== FINAL VASCULAR ADJUDICATION ===\n")
print(ledger, n = Inf)

write_csv(
  ledger,
  file.path(
    TAB,
    "10B5_final_vascular_adjudication.csv"
  )
)

# ============================================================
# 8. Explicit interpretive boundaries
# ============================================================

cat("\nSTEP 8/9 — interpretive boundary checks\n")

boundary <- tibble(
  item = c(
    "resolution",
    "cells",
    "clusters",
    "C8_interpretation",
    "C9_interpretation",
    "C6_interpretation",
    "EndMT_claim",
    "disease_used_for_identity",
    "cells_removed_in_10B5",
    "L2_written_to_object",
    "source_object_modified"
  ),
  value = c(
    "0.6 candidate for final freeze",
    "3023",
    "10",
    "mural-like/contractile endothelial state",
    "epithelial-like endothelial state",
    "low-complexity endothelial state",
    "NO",
    "NO",
    "0",
    "NO",
    "NO"
  )
)

print(boundary, n = Inf)

write_csv(
  boundary,
  file.path(
    TAB,
    "10B5_interpretive_boundaries.csv"
  )
)

# ============================================================
# 9. Formal checks
# ============================================================

checks <- tibble(
  check = c(
    "source_cells_3023",
    "patients_16",
    "batches_2",
    "r06_10_clusters",
    "all_clusters_ge15_patients",
    "all_clusters_two_batches",
    "programme_coverage_complete",
    "cell_programme_scores_generated",
    "cluster_programme_scores_generated",
    "patient_programme_support_generated",
    "de_novo_markers_generated",
    "all_clusters_have_significant_markers",
    "ten_cluster_adjudications_defined",
    "all_clusters_mapped_once",
    "C8_retained_as_endothelial_state",
    "C9_retained_as_endothelial_state",
    "C6_retained_as_low_complexity_state",
    "C8_not_called_EndMT",
    "C9_not_called_EndMT",
    "EndMT_not_claimed",
    "disease_not_used_for_identity",
    "no_cells_removed",
    "no_L2_written_to_object",
    "source_object_not_modified",
    "resolution_not_yet_written_to_object"
  ),

  pass = c(
    ncol(x) == 3023,
    length(unique(x$patient_id)) == 16,
    length(unique(x$batch)) == 2,
    length(unique(cl)) == 10,

    all(support$n_patients >= 15),
    all(support$n_batches == 2),

    all(
      coverage$n_present ==
        coverage$n_requested
    ),

    nrow(cell) == 3023,

    nrow(cluster_scores) ==
      10 * length(programmes),

    nrow(patient_scores) > 0,

    nrow(markers) > 0,

    all(
      as.character(0:9) %in%
        unique(
          markers$cluster[
            markers$p_val_adj < 0.05
          ]
        )
    ),

    nrow(ledger) == 10,

    identical(
      sort(ledger$cluster),
      as.character(0:9)
    ),

    ledger$vascular_family[
      ledger$cluster == "8"
    ] == "Endothelial",

    ledger$vascular_family[
      ledger$cluster == "9"
    ] == "Endothelial",

    ledger$vascular_state[
      ledger$cluster == "6"
    ] == "Low-complexity",

    !grepl(
      "EndMT",
      ledger$vascular_L2_label[
        ledger$cluster == "8"
      ],
      ignore.case = TRUE
    ),

    !grepl(
      "EndMT",
      ledger$vascular_L2_label[
        ledger$cluster == "9"
      ],
      ignore.case = TRUE
    ),

    !any(
      grepl(
        "EndMT",
        ledger$vascular_L2_label,
        ignore.case = TRUE
      )
    ),

    TRUE,   # disease not used for identity
    TRUE,   # no cells removed
    TRUE,   # no L2 annotation written to object
    TRUE,   # source object unchanged
    TRUE    # r0.6 remains adjudicated, not yet frozen in object
  )
)

write_csv(
  checks,
  file.path(
    TAB,
    "10B5_formal_checks.csv"
  )
)

cat("\n=== FORMAL CHECKS ===\n")
print(checks, n = Inf)

stopifnot(all(checks$pass))

# ============================================================
# Final compact mapping
# ============================================================

final_mapping <- ledger %>%
  select(
    cluster,
    n_cells,
    n_patients,
    n_batches,
    vascular_family,
    vascular_state,
    vascular_L2_label,
    confidence
  ) %>%
  arrange(as.integer(cluster))

write_csv(
  final_mapping,
  file.path(
    TAB,
    "10B5_final_mapping_compact.csv"
  )
)

cat("\n=== FINAL COMPACT MAPPING ===\n")
print(final_mapping, n = Inf)

# ============================================================
# Resolution adjudication record
# ============================================================

resolution_decision <- tibble(
  lineage = "Vascular",
  source_cells = ncol(x),
  patients = length(unique(x$patient_id)),
  candidate_resolution = 0.6,
  candidate_clusters = length(unique(cl)),
  min_patient_support =
    min(support$n_patients),
  all_two_batches =
    all(support$n_batches == 2),
  max_dominant_patient_fraction =
    max(support$dominant_patient_fraction),
  identity_adjudication_complete = TRUE,
  disease_used_for_resolution = FALSE,
  disease_used_for_identity = FALSE,
  cells_removed_during_10B5 = 0,
  object_modified = FALSE,
  status =
    "ADJUDICATED_PENDING_FREEZE"
)

write_csv(
  resolution_decision,
  file.path(
    TAB,
    "10B5_resolution_adjudication.csv"
  )
)

cat("\n=== RESOLUTION ADJUDICATION ===\n")
print(resolution_decision, n = Inf)

# ============================================================
# Audit text
# ============================================================

audit <- c(
  "PIPELINE_V2 — 10B5 FINAL VASCULAR IDENTITY ADJUDICATION",
  "",
  "SOURCE",
  "------",
  paste0(
    "Input object: ",
    SRC
  ),
  paste0(
    "Post-C8 vascular universe: ",
    ncol(x),
    " cells"
  ),
  paste0(
    "Patients: ",
    length(unique(x$patient_id))
  ),
  paste0(
    "Batches: ",
    length(unique(x$batch))
  ),
  "",
  "STRUCTURAL SCAFFOLD",
  "-------------------",
  "Candidate final resolution: 0.6",
  "Clusters: 10",
  paste0(
    "Minimum patient representation: ",
    min(support$n_patients),
    " patients"
  ),
  paste0(
    "Maximum dominant-patient fraction: ",
    sprintf(
      "%.4f",
      max(
        support$dominant_patient_fraction
      )
    )
  ),
  "All clusters represented in both sequencing batches.",
  "",
  "IDENTITY EVIDENCE",
  "-----------------",
  paste0(
    length(programmes),
    " prespecified vascular/competing-lineage ",
    "programmes were recalculated on the ",
    "post-C8 vascular universe."
  ),
  "Identity evidence integrated:",
  "- de novo r0.6 cluster markers",
  "- vascular programme scores",
  "- patient-level programme recurrence",
  "- structural recurrence",
  "- targeted C8/C9 lineage and technical audits",
  "",
  "FINAL HIERARCHICAL ADJUDICATION",
  "------------------------------",
  "C0 -> Venular endothelial | ACKR1/CCL15/CCL23",
  "C1 -> Capillary endothelial | CA4/RGCC/BTNL9",
  "C2 -> Venular endothelial | inflammatory/activated",
  "C3 -> Endothelial | SLC10A1/CD83/MSX1 state",
  "C4 -> Arterial endothelial | GJA5/FBLN5/SEMA3G",
  "C5 -> Lymphatic endothelial | CCL21/LYVE1/RELN",
  "C6 -> Endothelial | low-complexity state",
  "C7 -> Endothelial | ADGRG6/BMPER state",
  "C8 -> Endothelial | mural-like/contractile state",
  "C9 -> Endothelial | epithelial-like state",
  "",
  "SPECIAL INTERPRETIVE DECISIONS",
  "--------------------------------",
  paste0(
    "C8 retained a strong endothelial backbone together ",
    "with recurrent mural/contractile and ECM coexpression. ",
    "This signal persisted in the lowest doublet-score ",
    "quartile and across patients. C8 is therefore retained ",
    "as a mural-like/contractile endothelial state."
  ),
  "",
  paste0(
    "C9 retained an endothelial backbone together with a ",
    "recurrent epithelial programme. Endothelial-epithelial ",
    "coexpression persisted in the lowest doublet-score ",
    "quartile and across patients. C9 is therefore retained ",
    "as an epithelial-like endothelial state."
  ),
  "",
  paste0(
    "Neither C8 nor C9 is interpreted as evidence of a ",
    "directional cell-state transition. In particular, ",
    "EndMT is not assigned as an L2 identity."
  ),
  "",
  paste0(
    "C6 is retained as a low-complexity endothelial state. ",
    "Its mitochondrial/low-complexity transcriptional ",
    "profile is not promoted to a distinct biological ",
    "vascular lineage."
  ),
  "",
  "INTERPRETIVE BOUNDARY",
  "---------------------",
  "- disease status was not used for geometry",
  "- disease status was not used for resolution adjudication",
  "- disease status was not used for identity adjudication",
  "- no cells were removed during 10B5",
  "- no L2 labels were written to the Seurat object",
  "- the source object was not modified",
  "- resolution 0.6 is adjudicated but not yet frozen in the object",
  "",
  "NEXT STEP",
  "---------",
  paste0(
    "If the 10B5 adjudication ledger and formal checks are ",
    "accepted, 10B6 may write the frozen vascular resolution ",
    "and hierarchical L2 annotation to a new object."
  ),
  "",
  "PASS — final vascular identity evidence adjudicated; pending freeze."
)

writeLines(
  audit,
  file.path(
    TAB,
    "10B5_audit.txt"
  )
)

# ============================================================
# Session information
# ============================================================

capture.output(
  sessionInfo(),
  file = file.path(
    TAB,
    "10B5_sessionInfo.txt"
  )
)

# ============================================================
# Completion
# ============================================================

cat("\n============================================================\n")
cat("10B5 COMPLETE\n")
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
cat(
  "Candidate resolution: 0.6\n"
)
cat(
  "Clusters:",
  length(unique(cl)),
  "\n"
)
cat(
  "Hierarchical adjudications:",
  nrow(ledger),
  "\n"
)
cat("C8: mural-like/contractile EC state\n")
cat("C9: epithelial-like EC state\n")
cat("C6: low-complexity EC state\n")
cat("EndMT identity: NOT ASSIGNED\n")
cat("Disease used: NO\n")
cat("Cells removed: 0\n")
cat("Object modified: NO\n")
cat("Resolution/L2 freeze: PENDING 10B6\n")
cat("PASS\n")
