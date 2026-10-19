suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
})

cat("============================================================\n")
cat("PIPELINE_V2 — 18D2 AA BATCH-DIRECTION AUDIT\n")
cat("============================================================\n")

BASE18C <- paste0(
  "PIPELINE_V2/results/disease_biology/",
  "18_PSEUDOBULK_DE/18C_DE"
)

OUT <- paste0(
  "PIPELINE_V2/results/disease_biology/",
  "18_PSEUDOBULK_DE/18D_batch_audit/tables"
)

sig <- read_csv(
  file.path(BASE18C, "tables", "18C_all_FDR005_genes.csv"),
  show_col_types = FALSE
)

aa <- sig %>%
  filter(
    contrast == "AA_vs_HC",
    robust_ge50
  )

cat("\nAA-HC robust significant calls:", nrow(aa), "\n")
cat("L2:", n_distinct(aa$final_L2_label), "\n")

# ------------------------------------------------------------
# Recurrent genes across AA L2
# ------------------------------------------------------------

rec <- aa %>%
  mutate(
    direction = if_else(logFC > 0, "UP", "DOWN")
  ) %>%
  group_by(gene) %>%
  summarise(
    n_L2 = n_distinct(final_L2_label),
    n_up = sum(direction == "UP"),
    n_down = sum(direction == "DOWN"),
    median_logFC = median(logFC),
    max_abs_logFC = max(abs(logFC)),
    L2 = paste(
      sort(unique(final_L2_label)),
      collapse = ";"
    ),
    .groups = "drop"
  ) %>%
  arrange(
    desc(n_L2),
    desc(max_abs_logFC)
  )

write_csv(
  rec,
  file.path(
    OUT,
    "18D2_AA_HC_recurrent_DE_genes_across_L2.csv"
  )
)

cat("\n===== RECURRENT AA-HC GENES >=3 L2 =====\n")

print(
  as.data.frame(
    rec %>%
      filter(n_L2 >= 3)
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# Directional consistency
# ------------------------------------------------------------

direction_summary <- rec %>%
  mutate(
    direction_consistent =
      n_up == n_L2 |
      n_down == n_L2
  ) %>%
  summarise(
    genes_total = n(),
    genes_ge2_L2 = sum(n_L2 >= 2),
    genes_ge3_L2 = sum(n_L2 >= 3),
    genes_ge4_L2 = sum(n_L2 >= 4),
    genes_ge3_direction_consistent =
      sum(
        n_L2 >= 3 &
        direction_consistent
      )
  )

cat("\n===== RECURRENCE SUMMARY =====\n")
print(direction_summary)

write_csv(
  direction_summary,
  file.path(
    OUT,
    "18D2_AA_HC_recurrence_summary.csv"
  )
)

# ------------------------------------------------------------
# Candidate biological signals of special interest
# ------------------------------------------------------------

genes_interest <- c(
  "CXCL9", "CXCL10", "CXCL11",
  "GBP1", "GBP2", "STAT1",
  "IRF1", "IFIT1", "IFIT2",
  "IFIT3", "ISG15",
  "HLA-DRA", "HLA-DRB1",
  "HLA-DPA1", "HLA-DPB1",
  "B2M", "TAP1",
  "CD200",
  "KRT75", "KRT14",
  "KRT17", "KRT6A",
  "KRT6B", "KRT6C"
)

interest <- aa %>%
  filter(gene %in% genes_interest) %>%
  arrange(
    gene,
    final_L2_label
  ) %>%
  select(
    final_L2_label,
    gene,
    logFC,
    logCPM,
    FDR_BH
  )

cat("\n===== AA BIOLOGICAL ANCHOR GENES =====\n")
print(as.data.frame(interest), row.names = FALSE)

write_csv(
  interest,
  file.path(
    OUT,
    "18D2_AA_HC_biological_anchor_genes.csv"
  )
)

cat("\n============================================================\n")
cat("18D2 PHASE A COMPLETE\n")
cat("No new DE performed\n")
cat("No annotation modified\n")
cat("PASS\n")
cat("============================================================\n")
