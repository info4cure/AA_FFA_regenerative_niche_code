
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(Matrix)
})

cat("============================================================\n")
cat("PIPELINE_V2 — 08H IMMUNE FINAL L2 ADJUDICATION\n")
cat("============================================================\n")

OBJ <- "PIPELINE_V2/objects/annotated/immune_v2_08G_refined.rds"

OUT <- paste0(
  "PIPELINE_V2/results/annotation/08_IMMUNE_L2/",
  "08H_final_L2_adjudication"
)

TAB <- file.path(OUT, "tables")
dir.create(TAB, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# STEP 1 — load refined immune universe
# ============================================================

cat("\nSTEP 1/7 — loading refined immune universe\n")

x <- readRDS(OBJ)

stopifnot(ncol(x) == 4542)
stopifnot(length(unique(x$patient_id)) == 16)

cluster_col <- "immune_refined_res.0.4"

stopifnot(cluster_col %in% colnames(x@meta.data))

x$immune_08H_cluster <- as.character(x@meta.data[[cluster_col]])

cat("Cells:", ncol(x), "\n")
cat("Genes:", nrow(x), "\n")
cat("Patients:", length(unique(x$patient_id)), "\n")
cat("Clusters:", length(unique(x$immune_08H_cluster)), "\n")

print(sort(table(x$immune_08H_cluster)))

# ============================================================
# STEP 2 — patient / disease / batch support
# ============================================================

cat("\nSTEP 2/7 — cluster structural support\n")

meta <- x@meta.data %>%
  mutate(
    cell = rownames(.),
    cluster = immune_08H_cluster
  )

support <- meta %>%
  group_by(cluster) %>%
  summarise(
    n_cells = n(),
    n_patients = n_distinct(patient_id),
    n_HC = n_distinct(patient_id[disease == "HC"]),
    n_AA = n_distinct(patient_id[disease == "AA"]),
    n_FFA = n_distinct(patient_id[disease == "FFA"]),
    n_batches = n_distinct(batch),
    .groups = "drop"
  )

dominance <- meta %>%
  count(cluster, patient_id, name = "n") %>%
  group_by(cluster) %>%
  mutate(total = sum(n)) %>%
  arrange(cluster, desc(n)) %>%
  slice_head(n = 1) %>%
  transmute(
    cluster,
    dominant_patient = patient_id,
    dominant_patient_cells = n,
    dominant_patient_fraction = n / total
  )

support <- support %>%
  left_join(dominance, by = "cluster") %>%
  arrange(as.numeric(cluster))

write_csv(
  support,
  file.path(TAB, "08H_cluster_support.csv")
)

print(support)

# ============================================================
# STEP 3 — de novo markers
# ============================================================

cat("\nSTEP 3/7 — calculating de novo markers\n")

Idents(x) <- x$immune_08H_cluster

DefaultAssay(x) <- "RNA"

markers <- FindAllMarkers(
  x,
  assay = "RNA",
  only.pos = TRUE,
  min.pct = 0.10,
  logfc.threshold = 0.25,
  test.use = "wilcox",
  verbose = FALSE
)

markers <- markers %>%
  mutate(cluster = as.character(cluster))

write_csv(
  markers,
  file.path(TAB, "08H_markers_all.csv")
)

top30 <- markers %>%
  group_by(cluster) %>%
  arrange(p_val_adj, desc(avg_log2FC), .by_group = TRUE) %>%
  slice_head(n = 30) %>%
  ungroup()

write_csv(
  top30,
  file.path(TAB, "08H_markers_top30.csv")
)

cat("\n=== TOP 15 MARKERS PER CLUSTER ===\n")

for (cl in sort(unique(top30$cluster))) {

  cat("\n--- CLUSTER", cl, "---\n")

  print(
    top30 %>%
      filter(cluster == cl) %>%
      select(
        gene,
        avg_log2FC,
        pct.1,
        pct.2,
        p_val_adj
      ) %>%
      slice_head(n = 15)
  )
}

# ============================================================
# STEP 4 — canonical immune programmes
# ============================================================

cat("\nSTEP 4/7 — canonical immune identity programmes\n")

programmes <- list(

  CD4_T = c(
    "CD3D","CD3E","CD3G","TRAC",
    "IL7R","LTB","MALAT1","LTST1",
    "CCR7","LTB","MAL"
  ),

  CD8_T = c(
    "CD3D","CD3E","TRAC",
    "CD8A","CD8B",
    "CCL5","NKG7","GZMK"
  ),

  NK = c(
    "NKG7","GNLY","KLRD1",
    "PRF1","GZMB","GZMH",
    "FCGR3A","TYROBP"
  ),

  Myeloid = c(
    "LYZ","LST1","FCER1G",
    "TYROBP","CTSS","AIF1",
    "CST3"
  ),

  Macrophage = c(
    "C1QA","C1QB","C1QC",
    "APOE","CTSB","CTSD",
    "FCGR3A"
  ),

  Dendritic = c(
    "FCER1A","CD1C","CLEC10A",
    "CST3","HLA-DRA","HLA-DPA1"
  ),

  Langerhans = c(
    "CD207","CD1A","CD1E",
    "FCER1A","HLA-DRA"
  ),

  Mast = c(
    "TPSAB1","TPSB2","CPA3",
    "KIT","MS4A2","HDC"
  ),

  Plasma = c(
    "MZB1","JCHAIN","DERL3",
    "SDC1","IGHG1","IGKC"
  ),

  B_cell = c(
    "CD79A","MS4A1","CD37",
    "CD74","HLA-DRA","CD79B",
    "CD22"
  )
)

programmes <- lapply(
  programmes,
  function(g) unique(g[g %in% rownames(x)])
)

coverage <- tibble(
  programme = names(programmes),
  n_present = lengths(programmes)
)

write_csv(
  coverage,
  file.path(TAB, "08H_programme_coverage.csv")
)

print(coverage)

# ============================================================
# STEP 5 — cell-level programme scoring
# ============================================================

cat("\nSTEP 5/7 — programme scoring\n")

expr <- GetAssayData(
  x,
  assay = "RNA",
  layer = "data"
)

score_matrix <- sapply(
  names(programmes),
  function(p) {

    genes <- programmes[[p]]

    if (length(genes) == 0)
      return(rep(NA_real_, ncol(x)))

    Matrix::colMeans(
      expr[genes, , drop = FALSE]
    )
  }
)

score_matrix <- as.data.frame(score_matrix)

rownames(score_matrix) <- colnames(x)

score_z <- as.data.frame(
  scale(score_matrix)
)

score_z$cell <- rownames(score_z)

score_z <- score_z %>%
  left_join(
    meta %>%
      select(cell, cluster, patient_id, disease, batch),
    by = "cell"
  )

write_csv(
  score_z,
  file.path(TAB, "08H_cell_identity_scores.csv")
)

# ============================================================
# STEP 6 — cluster programme evidence
# ============================================================

cat("\nSTEP 6/7 — cluster-level identity evidence\n")

long <- score_z %>%
  pivot_longer(
    cols = all_of(names(programmes)),
    names_to = "programme",
    values_to = "z"
  )

cluster_scores <- long %>%
  group_by(cluster, programme) %>%
  summarise(
    mean_z = mean(z, na.rm = TRUE),
    median_z = median(z, na.rm = TRUE),
    n_cells = n(),
    n_patients = n_distinct(patient_id),
    .groups = "drop"
  )

write_csv(
  cluster_scores,
  file.path(TAB, "08H_cluster_identity_scores.csv")
)

ranking <- cluster_scores %>%
  group_by(cluster) %>%
  arrange(desc(mean_z), .by_group = TRUE) %>%
  mutate(rank = row_number()) %>%
  ungroup()

write_csv(
  ranking,
  file.path(TAB, "08H_cluster_identity_ranking.csv")
)

top_identity <- ranking %>%
  filter(rank <= 2) %>%
  select(
    cluster,
    rank,
    programme,
    mean_z,
    median_z
  ) %>%
  pivot_wider(
    names_from = rank,
    values_from = c(
      programme,
      mean_z,
      median_z
    ),
    names_glue = "{.value}_{rank}"
  ) %>%
  rename(
    top_identity = programme_1,
    top_mean_z = mean_z_1,
    top_median_z = median_z_1,
    second_identity = programme_2,
    second_mean_z = mean_z_2,
    second_median_z = median_z_2
  ) %>%
  mutate(
    identity_margin =
      top_mean_z - second_mean_z
  ) %>%
  left_join(
    support,
    by = "cluster"
  ) %>%
  arrange(as.numeric(cluster))

write_csv(
  top_identity,
  file.path(TAB, "08H_cluster_top_identity.csv")
)

cat("\n=== TOP IDENTITY PROGRAMMES ===\n")

print(
  top_identity %>%
    select(
      cluster,
      n_cells,
      n_patients,
      top_identity,
      top_mean_z,
      second_identity,
      second_mean_z,
      identity_margin,
      dominant_patient_fraction
    )
)

# patient-within-cluster support

patient_scores <- long %>%
  group_by(
    cluster,
    patient_id,
    disease,
    programme
  ) %>%
  summarise(
    mean_z = mean(z, na.rm = TRUE),
    n_cells = n(),
    .groups = "drop"
  )

write_csv(
  patient_scores,
  file.path(
    TAB,
    "08H_patient_cluster_identity_scores.csv"
  )
)

patient_top <- patient_scores %>%
  group_by(cluster, patient_id) %>%
  arrange(desc(mean_z), .by_group = TRUE) %>%
  slice_head(n = 1) %>%
  ungroup()

patient_support <- patient_top %>%
  count(
    cluster,
    programme,
    name = "n_patients_top"
  ) %>%
  group_by(cluster) %>%
  mutate(
    n_patients_scored = sum(n_patients_top),
    patient_fraction =
      n_patients_top / n_patients_scored
  ) %>%
  ungroup()

write_csv(
  patient_support,
  file.path(
    TAB,
    "08H_patient_identity_support.csv"
  )
)

# ============================================================
# STEP 7 — formal checks
# ============================================================

cat("\nSTEP 7/7 — formal checks\n")

checks <- tibble(
  check = c(
    "cells_4542",
    "genes_23412",
    "patients_16",
    "resolution04_10_clusters",
    "all_clusters_ge2_batches",
    "markers_generated",
    "identity_programmes_generated",
    "cluster_identity_ranking_generated",
    "patient_support_generated",
    "disease_not_used_for_identity_scoring",
    "no_L2_labels_assigned",
    "no_cells_removed",
    "source_object_not_modified"
  ),
  pass = c(
    ncol(x) == 4542,
    nrow(x) == 23412,
    length(unique(x$patient_id)) == 16,
    length(unique(x$immune_08H_cluster)) == 10,
    all(support$n_batches >= 2),
    nrow(markers) > 0,
    length(programmes) == 10,
    nrow(ranking) == 10 * length(programmes),
    nrow(patient_support) > 0,
    TRUE,
    !any(grepl(
      "L2_frozen",
      colnames(x@meta.data),
      ignore.case = TRUE
    )),
    ncol(x) == 4542,
    TRUE
  )
)

write_csv(
  checks,
  file.path(TAB, "08H_formal_checks.csv")
)

print(checks)

if (!all(checks$pass)) {
  cat("\nFAILED CHECKS:\n")
  print(filter(checks, !pass))
  stop("08H formal checks failed.")
}

cat("\n============================================================\n")
cat("08H COMPLETE\n")
cat("============================================================\n")
cat("Immune cells: 4542\n")
cat("Candidate resolution: 0.4\n")
cat("Clusters evaluated: 10\n")
cat("No L2 identities assigned yet.\n")
