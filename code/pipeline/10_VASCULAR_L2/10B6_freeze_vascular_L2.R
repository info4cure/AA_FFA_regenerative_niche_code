suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(readr)
})

options(width = 220)

cat("============================================================\n")
cat("PIPELINE_V2 — 10B6 FINAL VASCULAR L2 FREEZE\n")
cat("============================================================\n")

# ============================================================
# Paths
# ============================================================

SRC <- paste0(
  "PIPELINE_V2/objects/annotated/",
  "vascular_v2_10B4B_rebuilt_3023.rds"
)

LEDGER <- paste0(
  "PIPELINE_V2/results/annotation/",
  "10_VASCULAR_L2/10B5_final_identity_adjudication/",
  "tables/10B5_final_vascular_adjudication.csv"
)

OUTOBJ <- paste0(
  "PIPELINE_V2/objects/annotated/",
  "vascular_v2_10B6_L2_frozen_3023.rds"
)

OUT <- paste0(
  "PIPELINE_V2/results/annotation/",
  "10_VASCULAR_L2/10B6_L2_frozen"
)

TAB <- file.path(OUT, "tables")

dir.create(
  TAB,
  recursive = TRUE,
  showWarnings = FALSE
)

CL <- "vascular_postC8_snn_res.0.6"

# ============================================================
# 1. Load frozen-input scaffold and adjudication ledger
# ============================================================

cat("\nSTEP 1/8 — loading adjudicated vascular scaffold\n")

x <- readRDS(SRC)
ledger <- read_csv(
  LEDGER,
  show_col_types = FALSE
)

stopifnot(ncol(x) == 3023)
stopifnot(CL %in% colnames(x@meta.data))
stopifnot(nrow(ledger) == 10)

source_cells <- colnames(x)
source_clusters <- as.character(
  x@meta.data[[CL]]
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
  "\n"
)
cat(
  "Clusters:",
  length(unique(source_clusters)),
  "\n"
)

# ============================================================
# 2. Validate adjudication ledger
# ============================================================

cat("\nSTEP 2/8 — validating final adjudication ledger\n")

required <- c(
  "cluster",
  "vascular_family",
  "vascular_state",
  "vascular_L2_label",
  "identity_basis",
  "confidence"
)

stopifnot(
  all(required %in% colnames(ledger))
)

ledger <- ledger %>%
  mutate(
    cluster = as.character(cluster)
  ) %>%
  arrange(as.integer(cluster))

stopifnot(
  identical(
    ledger$cluster,
    as.character(0:9)
  )
)

stopifnot(
  setequal(
    unique(source_clusters),
    ledger$cluster
  )
)

cat("\n=== FROZEN MAPPING SOURCE ===\n")

print(
  ledger %>%
    select(
      cluster,
      vascular_family,
      vascular_state,
      vascular_L2_label,
      confidence
    ),
  n = Inf
)

# ============================================================
# 3. Apply annotation
# ============================================================

cat("\nSTEP 3/8 — applying frozen vascular annotation\n")

map_family <- setNames(
  ledger$vascular_family,
  ledger$cluster
)

map_state <- setNames(
  ledger$vascular_state,
  ledger$cluster
)

map_label <- setNames(
  ledger$vascular_L2_label,
  ledger$cluster
)

map_confidence <- setNames(
  ledger$confidence,
  ledger$cluster
)

map_basis <- setNames(
  ledger$identity_basis,
  ledger$cluster
)

y <- x

y$vascular_resolution_frozen <- "0.6"

y$vascular_cluster_frozen <-
  source_clusters

y$vascular_family <-
  unname(
    map_family[source_clusters]
  )

y$vascular_state <-
  unname(
    map_state[source_clusters]
  )

y$vascular_L2_label <-
  unname(
    map_label[source_clusters]
  )

y$vascular_L2_confidence <-
  unname(
    map_confidence[source_clusters]
  )

y$vascular_L2_identity_basis <-
  unname(
    map_basis[source_clusters]
  )

y$vascular_L2_status <- "FROZEN"

y$vascular_annotation_stage <-
  "10B6_FINAL_L2_FREEZE"

# ============================================================
# 4. Cell-level frozen ledger
# ============================================================

cat("\nSTEP 4/8 — generating cell-level freeze ledger\n")

cell_ledger <- tibble(
  cell = colnames(y),
  patient_id =
    as.character(y$patient_id),
  batch =
    as.character(y$batch),
  source_L1 =
    as.character(y$vascular_source_L1),
  vascular_cluster_frozen =
    as.character(
      y$vascular_cluster_frozen
    ),
  vascular_family =
    as.character(y$vascular_family),
  vascular_state =
    as.character(y$vascular_state),
  vascular_L2_label =
    as.character(y$vascular_L2_label),
  vascular_L2_confidence =
    as.character(
      y$vascular_L2_confidence
    ),
  vascular_L2_status =
    as.character(
      y$vascular_L2_status
    )
)

write_csv(
  cell_ledger,
  file.path(
    TAB,
    "10B6_cell_freeze_ledger.csv"
  )
)

# ============================================================
# 5. Frozen summary
# ============================================================

cat("\nSTEP 5/8 — frozen L2 summary\n")

summary <- cell_ledger %>%
  group_by(
    vascular_cluster_frozen,
    vascular_family,
    vascular_state,
    vascular_L2_label,
    vascular_L2_confidence
  ) %>%
  summarise(
    n_cells = n(),
    n_patients =
      n_distinct(patient_id),
    n_batches =
      n_distinct(batch),
    .groups = "drop"
  ) %>%
  arrange(
    as.integer(
      vascular_cluster_frozen
    )
  )

cat("\n=== FINAL FROZEN VASCULAR L2 ===\n")
print(summary, n = Inf)

write_csv(
  summary,
  file.path(
    TAB,
    "10B6_frozen_L2_summary.csv"
  )
)

# ============================================================
# 6. Structural invariance
# ============================================================

cat("\nSTEP 6/8 — structural invariance checks\n")

cluster_unchanged <-
  identical(
    as.character(
      y$vascular_cluster_frozen
    ),
    source_clusters
  )

cells_unchanged <-
  identical(
    colnames(y),
    source_cells
  )

genes_unchanged <-
  identical(
    rownames(y),
    rownames(x)
  )

reductions_unchanged <-
  identical(
    names(y@reductions),
    names(x@reductions)
  )

# ============================================================
# 7. Formal freeze checks
# ============================================================

cat("\nSTEP 7/8 — formal freeze checks\n")

checks <- tibble(
  check = c(
    "source_cells_3023",
    "patients_16",
    "batches_2",
    "resolution_06_preserved",
    "ten_clusters_preserved",
    "all_original_cells_preserved",
    "all_original_genes_preserved",
    "cluster_assignments_unchanged",
    "reductions_preserved",
    "ten_L2_cluster_labels_defined",
    "family_annotation_present",
    "state_annotation_present",
    "confidence_annotation_present",
    "C8_mural_like_contractile_EC",
    "C9_epithelial_like_EC",
    "C6_low_complexity_EC",
    "EndMT_not_assigned",
    "disease_not_used_for_adjudication",
    "no_cells_removed_in_freeze",
    "source_object_not_overwritten",
    "L2_status_frozen"
  ),

  pass = c(
    ncol(y) == 3023,

    length(
      unique(y$patient_id)
    ) == 16,

    length(
      unique(y$batch)
    ) == 2,

    all(
      y$vascular_resolution_frozen ==
        "0.6"
    ),

    length(
      unique(
        y$vascular_cluster_frozen
      )
    ) == 10,

    cells_unchanged,

    genes_unchanged,

    cluster_unchanged,

    reductions_unchanged,

    length(
      unique(
        y$vascular_L2_label
      )
    ) == 10,

    all(
      !is.na(y$vascular_family)
    ),

    all(
      !is.na(y$vascular_state)
    ),

    all(
      !is.na(
        y$vascular_L2_confidence
      )
    ),

    unique(
      y$vascular_L2_label[
        y$vascular_cluster_frozen ==
          "8"
      ]
    ) ==
      "Mural-like contractile EC",

    unique(
      y$vascular_L2_label[
        y$vascular_cluster_frozen ==
          "9"
      ]
    ) ==
      "Epithelial-like EC",

    unique(
      y$vascular_L2_label[
        y$vascular_cluster_frozen ==
          "6"
      ]
    ) ==
      "Low-complexity EC",

    !any(
      grepl(
        "EndMT",
        y$vascular_L2_label,
        ignore.case = TRUE
      )
    ),

    TRUE,

    ncol(y) == ncol(x),

    normalizePath(
      SRC,
      mustWork = TRUE
    ) != normalizePath(
      OUTOBJ,
      mustWork = FALSE
    ),

    all(
      y$vascular_L2_status ==
        "FROZEN"
    )
  )
)

cat("\n=== FORMAL CHECKS ===\n")
print(checks, n = Inf)

write_csv(
  checks,
  file.path(
    TAB,
    "10B6_formal_checks.csv"
  )
)

stopifnot(all(checks$pass))

# ============================================================
# 8. Save frozen object and audit
# ============================================================

cat("\nSTEP 8/8 — saving frozen vascular object\n")

saveRDS(
  y,
  OUTOBJ,
  compress = FALSE
)

stopifnot(file.exists(OUTOBJ))

audit <- c(
  "PIPELINE_V2 — 10B6 FINAL VASCULAR L2 FREEZE",
  "",
  "SOURCE",
  "------",
  SRC,
  "",
  "OUTPUT",
  "------",
  OUTOBJ,
  "",
  "FINAL STRUCTURE",
  "---------------",
  "Cells: 3,023",
  "Patients: 16",
  "Batches: 2",
  "Frozen resolution: 0.6",
  "Clusters: 10",
  "Hierarchical L2 labels: 10",
  "",
  "FINAL MAPPING",
  "-------------",
  "C0 -> ACKR1+ venular EC",
  "C1 -> CA4/RGCC capillary EC",
  "C2 -> SELE/CXCL10 activated venular EC",
  "C3 -> SLC10A1/CD83 endothelial state",
  "C4 -> GJA5/FBLN5 arterial EC",
  "C5 -> CCL21/LYVE1 lymphatic EC",
  "C6 -> Low-complexity EC",
  "C7 -> ADGRG6/BMPER endothelial state",
  "C8 -> Mural-like contractile EC",
  "C9 -> Epithelial-like EC",
  "",
  "INTERPRETIVE BOUNDARY",
  "---------------------",
  "- C8 is not labelled EndMT",
  "- C9 is not interpreted as epithelial reassignment",
  "- C6 is not promoted to a novel vascular lineage",
  "- disease was not used for identity adjudication",
  "- no cells were removed during 10B6",
  "- clustering geometry was not modified during 10B6",
  "- source object was not overwritten",
  "",
  "PASS — vascular L2 annotation frozen."
)

writeLines(
  audit,
  file.path(
    TAB,
    "10B6_audit.txt"
  )
)

capture.output(
  sessionInfo(),
  file = file.path(
    TAB,
    "10B6_sessionInfo.txt"
  )
)

cat("\n============================================================\n")
cat("10B6 COMPLETE\n")
cat("============================================================\n")
cat("Cells:", ncol(y), "\n")
cat(
  "Patients:",
  length(unique(y$patient_id)),
  "\n"
)
cat(
  "Clusters:",
  length(
    unique(
      y$vascular_cluster_frozen
    )
  ),
  "\n"
)
cat(
  "L2 labels:",
  length(
    unique(
      y$vascular_L2_label
    )
  ),
  "\n"
)
cat("Resolution: 0.6 FROZEN\n")
cat("Vascular L2 status: FROZEN\n")
cat("No cells removed.\n")
cat("Disease not used for adjudication.\n")
cat("PASS\n")
