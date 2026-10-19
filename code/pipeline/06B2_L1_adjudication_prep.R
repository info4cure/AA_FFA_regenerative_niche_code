suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
})

cat("\n============================================================\n")
cat("PIPELINE_V2 — 06B2 L1 ADJUDICATION PREP\n")
cat("============================================================\n\n")

indir <- "PIPELINE_V2/results/annotation/06_GLOBAL_L1/06B1_corrected/tables"

markers <- read_csv(
  file.path(indir, "06B1_all_positive_markers.csv"),
  show_col_types = FALSE
)

ranking <- read_csv(
  file.path(indir, "06B1_cluster_programme_ranking.csv"),
  show_col_types = FALSE
)

evidence <- read_csv(
  file.path(indir, "06B1_cluster_programme_evidence.csv"),
  show_col_types = FALSE
)

patient <- read_csv(
  file.path(indir, "06B1_patient_cluster_programmes.csv"),
  show_col_types = FALSE
)

# ----------------------------------------------------------
# 1. Top discriminating markers per cluster
# ----------------------------------------------------------

top_markers <- markers %>%
  mutate(cluster = as.character(cluster)) %>%
  filter(
    avg_log2FC > 0,
    pct.1 >= 0.10
  ) %>%
  arrange(cluster, desc(avg_log2FC)) %>%
  group_by(cluster) %>%
  slice_head(n = 15) %>%
  summarise(
    top_markers = paste(gene, collapse = ", "),
    .groups = "drop"
  )

# ----------------------------------------------------------
# 2. Programme evidence: top 3 raw programmes
# ----------------------------------------------------------

top_programmes_raw <- evidence %>%
  mutate(cluster = as.character(cluster)) %>%
  group_by(cluster) %>%
  arrange(desc(mean_programme), .by_group = TRUE) %>%
  slice_head(n = 3) %>%
  mutate(
    txt = paste0(
      programme,
      "=",
      sprintf("%.3f", mean_programme),
      " [",
      sprintf("%.1f", mean_gene_detection_pct),
      "%]"
    )
  ) %>%
  summarise(
    raw_top3 = paste(txt, collapse = " | "),
    .groups = "drop"
  )

# ----------------------------------------------------------
# 3. Patient-level reproducibility of ranked top programme
# ----------------------------------------------------------

ranking2 <- ranking %>%
  mutate(cluster = as.character(cluster))

patient2 <- patient %>%
  mutate(cluster = as.character(cluster)) %>%
  left_join(
    ranking2 %>% select(cluster, top_programme),
    by = "cluster"
  )

patient_support <- patient2 %>%
  group_by(cluster, patient_id) %>%
  mutate(
    patient_rank = rank(
      -mean_programme,
      ties.method = "min"
    )
  ) %>%
  ungroup() %>%
  filter(programme == top_programme) %>%
  group_by(cluster) %>%
  summarise(
    n_patients = n_distinct(patient_id),
    n_patient_top1 = sum(patient_rank == 1),
    n_patient_top2 = sum(patient_rank <= 2),
    patient_top1_frac = n_patient_top1 / n_patients,
    patient_top2_frac = n_patient_top2 / n_patients,
    .groups = "drop"
  )

# ----------------------------------------------------------
# 4. Cluster cell counts
# ----------------------------------------------------------

ncells <- evidence %>%
  mutate(cluster = as.character(cluster)) %>%
  distinct(cluster, n_cells)

# ----------------------------------------------------------
# 5. Combined adjudication evidence table
# ----------------------------------------------------------

summary_tbl <- ranking2 %>%
  left_join(ncells, by = "cluster") %>%
  left_join(patient_support, by = "cluster") %>%
  left_join(top_programmes_raw, by = "cluster") %>%
  left_join(top_markers, by = "cluster") %>%
  mutate(cluster_num = suppressWarnings(as.integer(cluster))) %>%
  arrange(cluster_num) %>%
  select(
    cluster,
    n_cells,
    top_programme,
    second_programme,
    top_z,
    second_z,
    margin_z,
    n_patients,
    n_patient_top1,
    patient_top1_frac,
    n_patient_top2,
    patient_top2_frac,
    raw_top3,
    top_markers
  )

# ----------------------------------------------------------
# 6. Print compact cluster-by-cluster report
# ----------------------------------------------------------

for (i in seq_len(nrow(summary_tbl))) {

  x <- summary_tbl[i, ]

  cat("------------------------------------------------------------\n")
  cat("CLUSTER", x$cluster, " | n =", x$n_cells, "\n")

  cat(
    "RANK:",
    x$top_programme,
    sprintf("(z=%.2f)", x$top_z),
    "vs",
    x$second_programme,
    sprintf("(z=%.2f)", x$second_z),
    "| margin =",
    sprintf("%.2f", x$margin_z),
    "\n"
  )

  cat(
    "PATIENT SUPPORT:",
    x$n_patient_top1, "/", x$n_patients,
    "top1;",
    x$n_patient_top2, "/", x$n_patients,
    "top2\n"
  )

  cat("RAW PROGRAMMES:", x$raw_top3, "\n")
  cat("MARKERS:", x$top_markers, "\n\n")
}

cat("============================================================\n")
cat("06B2 PREP COMPLETE\n")
cat("============================================================\n")
cat("No labels assigned.\n")
cat("No object modified.\n")
cat("No cells removed.\n")

