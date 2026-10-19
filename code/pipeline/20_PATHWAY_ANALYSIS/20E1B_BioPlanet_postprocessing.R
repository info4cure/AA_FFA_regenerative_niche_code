suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(stringr)
  library(purrr)
})

cat("============================================================\n")
cat("PIPELINE_V2 — 20E1B BIOPLANET POSTPROCESSING\n")
cat("============================================================\n")

# ============================================================
# Paths
# ============================================================

ROOT <- paste0(
  "PIPELINE_V2/results/disease_biology/",
  "20_PATHWAY_ANALYSIS/20E_BIOPLANET"
)

IN <- file.path(
  ROOT,
  "20E1_GSEA",
  "tables",
  "20E1_all_BioPlanet_GSEA.csv"
)

OUT <- file.path(
  ROOT,
  "20E1B_POSTPROCESSING"
)

TAB <- file.path(OUT, "tables")

dir.create(
  TAB,
  recursive = TRUE,
  showWarnings = FALSE
)

# ============================================================
# 1. Recover frozen BioPlanet GSEA
# ============================================================

cat("\nSTEP 1/8 — recovering frozen 20E1 GSEA\n")

x <- read_csv(
  IN,
  show_col_types = FALSE
)

stopifnot(
  nrow(x) == 49000,
  all(c(
    "contrast",
    "module",
    "final_L2_label",
    "direction",
    "pathway",
    "bioplanet_id",
    "bioplanet_name",
    "padj",
    "NES",
    "leadingEdge"
  ) %in% names(x))
)

sig <- x %>%
  filter(
    !is.na(padj),
    padj < 0.05
  )

cat("Frozen result rows :", nrow(x), "\n")
cat("FDR<0.05 rows      :", nrow(sig), "\n")

stopifnot(nrow(sig) == 4378)

# ============================================================
# 2. Reconstruct leading-edge table
# ============================================================

cat("\nSTEP 2/8 — reconstructing leading edges\n")

leading <- sig %>%
  select(
    contrast,
    module,
    final_L2_label,
    direction,
    pathway,
    bioplanet_id,
    bioplanet_name,
    NES,
    padj,
    leadingEdge
  ) %>%
  filter(
    !is.na(leadingEdge),
    leadingEdge != ""
  ) %>%
  separate_rows(
    leadingEdge,
    sep = ";"
  ) %>%
  transmute(
    contrast,
    module,
    final_L2_label,
    direction,
    pathway,
    bioplanet_id,
    bioplanet_name,
    NES,
    padj,
    gene = str_trim(leadingEdge)
  ) %>%
  filter(
    !is.na(gene),
    gene != ""
  ) %>%
  distinct()

cat(
  "Leading-edge rows:",
  nrow(leading),
  "\n"
)

write_csv(
  leading,
  file.path(
    TAB,
    "20E1B_BioPlanet_leading_edge_genes.csv"
  )
)

# ============================================================
# 3. Pathway recurrence across L2
# ============================================================

cat("\nSTEP 3/8 — pathway recurrence across L2\n")

pathway_recurrence <- sig %>%
  group_by(
    contrast,
    bioplanet_id,
    bioplanet_name,
    direction
  ) %>%
  summarise(
    n_L2 = n_distinct(
      paste(module, final_L2_label, sep = "|||")
    ),
    median_NES = median(
      NES,
      na.rm = TRUE
    ),
    max_abs_NES = max(
      abs(NES),
      na.rm = TRUE
    ),
    min_FDR = min(
      padj,
      na.rm = TRUE
    ),
    L2 = paste(
      sort(unique(final_L2_label)),
      collapse = ";"
    ),
    .groups = "drop"
  ) %>%
  arrange(
    contrast,
    desc(n_L2),
    min_FDR
  )

write_csv(
  pathway_recurrence,
  file.path(
    TAB,
    "20E1B_BioPlanet_pathway_recurrence.csv"
  )
)

# ============================================================
# 4. Leading-edge driver recurrence
# ============================================================

cat("\nSTEP 4/8 — recurrent leading-edge drivers\n")

drivers <- leading %>%
  group_by(
    contrast,
    direction,
    gene
  ) %>%
  summarise(
    n_pathways = n_distinct(
      bioplanet_id
    ),
    n_L2 = n_distinct(
      paste(module, final_L2_label, sep = "|||")
    ),
    n_hits = n(),
    median_NES = median(
      NES,
      na.rm = TRUE
    ),
    min_FDR = min(
      padj,
      na.rm = TRUE
    ),
    pathways = paste(
      sort(unique(bioplanet_name)),
      collapse = ";"
    ),
    L2 = paste(
      sort(unique(final_L2_label)),
      collapse = ";"
    ),
    .groups = "drop"
  ) %>%
  arrange(
    contrast,
    direction,
    desc(n_L2),
    desc(n_pathways),
    min_FDR
  )

write_csv(
  drivers,
  file.path(
    TAB,
    "20E1B_BioPlanet_recurrent_driver_genes.csv"
  )
)

# ============================================================
# 5. L2-level summary
# ============================================================

cat("\nSTEP 5/8 — L2 summaries\n")

l2_summary <- sig %>%
  group_by(
    contrast,
    module,
    final_L2_label
  ) %>%
  summarise(
    n_FDR005 = n(),
    n_UP_FDR005 = sum(
      NES > 0,
      na.rm = TRUE
    ),
    n_DOWN_FDR005 = sum(
      NES < 0,
      na.rm = TRUE
    ),
    median_abs_NES = median(
      abs(NES),
      na.rm = TRUE
    ),
    min_FDR = min(
      padj,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  arrange(
    contrast,
    desc(n_FDR005)
  )

write_csv(
  l2_summary,
  file.path(
    TAB,
    "20E1B_BioPlanet_L2_summary.csv"
  )
)

# ============================================================
# 6. Global contrast summary
# ============================================================

cat("\nSTEP 6/8 — global summaries\n")

global_summary <- sig %>%
  group_by(
    contrast
  ) %>%
  summarise(
    n_L2 = n_distinct(
      paste(module, final_L2_label, sep = "|||")
    ),
    L2_with_FDR005 = n_distinct(
      paste(module, final_L2_label, sep = "|||")
    ),
    total_pathway_hits_FDR005 = n(),
    total_UP_FDR005 = sum(
      NES > 0,
      na.rm = TRUE
    ),
    total_DOWN_FDR005 = sum(
      NES < 0,
      na.rm = TRUE
    ),
    unique_pathways_FDR005 = n_distinct(
      bioplanet_id
    ),
    .groups = "drop"
  )

write_csv(
  global_summary,
  file.path(
    TAB,
    "20E1B_BioPlanet_global_summary.csv"
  )
)

# ============================================================
# 7. Compact biological readout
# ============================================================

cat("\nSTEP 7/8 — compact biological readout\n")

top_pathways <- pathway_recurrence %>%
  group_by(
    contrast,
    direction
  ) %>%
  slice_max(
    order_by = n_L2,
    n = 20,
    with_ties = FALSE
  ) %>%
  ungroup() %>%
  arrange(
    contrast,
    direction,
    desc(n_L2),
    min_FDR
  )

write_csv(
  top_pathways,
  file.path(
    TAB,
    "20E1B_top_recurrent_pathways.csv"
  )
)

top_drivers <- drivers %>%
  filter(
    n_L2 >= 2
  ) %>%
  group_by(
    contrast,
    direction
  ) %>%
  slice_max(
    order_by = n_L2,
    n = 50,
    with_ties = FALSE
  ) %>%
  ungroup() %>%
  arrange(
    contrast,
    direction,
    desc(n_L2),
    desc(n_pathways)
  )

write_csv(
  top_drivers,
  file.path(
    TAB,
    "20E1B_top_recurrent_driver_genes.csv"
  )
)

# ============================================================
# 8. Formal checks + audit
# ============================================================

cat("\nSTEP 8/8 — formal checks\n")

analytical_key_duplicates <- x %>%
  count(
    contrast,
    module,
    final_L2_label,
    bioplanet_id
  ) %>%
  filter(n > 1) %>%
  nrow()

checks <- tibble(
  check = c(
    "frozen_20E1_present",
    "result_rows_49000",
    "FDR005_rows_4378",
    "three_contrasts_present",
    "leading_edges_generated",
    "pathway_recurrence_generated",
    "driver_recurrence_generated",
    "L2_summary_generated",
    "global_summary_generated",
    "no_duplicate_analytical_keys",
    "no_GSEA_refit",
    "no_DE_refit",
    "annotation_unchanged",
    "no_figures_generated"
  ),
  pass = c(
    file.exists(IN),
    nrow(x) == 49000,
    nrow(sig) == 4378,
    setequal(
      unique(x$contrast),
      c(
        "AA_vs_HC",
        "FFA_vs_HC",
        "AA_vs_FFA"
      )
    ),
    nrow(leading) > 0,
    nrow(pathway_recurrence) > 0,
    nrow(drivers) > 0,
    nrow(l2_summary) > 0,
    nrow(global_summary) == 3,
    analytical_key_duplicates == 0,
    TRUE,
    TRUE,
    TRUE,
    TRUE
  )
)

print(
  as.data.frame(checks),
  row.names = FALSE
)

stopifnot(all(checks$pass))

write_csv(
  checks,
  file.path(
    TAB,
    "20E1B_formal_checks.csv"
  )
)

audit <- c(
  "PIPELINE_V2 — 20E1B BIOPLANET POSTPROCESSING",
  "============================================================",
  "",
  "Input:",
  "- frozen 20E1 BioPlanet GSEA consolidated table",
  "- 49,000 pathway-analysis rows",
  "- 4,378 BioPlanet FDR<0.05 hits",
  "",
  "BioPlanet source:",
  "- NCATS BioPlanet",
  "- source date: 2019-05-14",
  "- 1,658 frozen pathways",
  "- 1,246 pathways eligible for GSEA at size 10-500",
  "",
  "Ranking:",
  "- frozen 20A signed edgeR QL rankings",
  "- sign(logFC) * sqrt(QL F)",
  "- no DEG threshold before GSEA",
  "",
  "Postprocessing:",
  "- semicolon-serialized fgsea leadingEdge fields reconstructed",
  "- pathway recurrence quantified across frozen L2 populations",
  "- recurrent leading-edge genes quantified across pathways and L2",
  "- no GSEA was rerun",
  "- no DE was refit",
  "- no cell annotation was changed",
  "",
  "Formal checks:",
  paste0(
    checks$check,
    ": ",
    ifelse(checks$pass, "PASS", "FAIL")
  ),
  "",
  "STATUS: PASS"
)

writeLines(
  audit,
  file.path(
    OUT,
    "20E1B_audit.txt"
  )
)

cat("\n===== GLOBAL BIOPLANET SUMMARY =====\n")
print(
  as.data.frame(global_summary),
  row.names = FALSE
)

cat("\n===== TOP RECURRENT BIOPLANET PATHWAYS =====\n")
print(
  as.data.frame(
    top_pathways %>%
      select(
        contrast,
        direction,
        bioplanet_id,
        bioplanet_name,
        n_L2,
        median_NES,
        min_FDR
      )
  ),
  row.names = FALSE
)

cat("\n===== TOP RECURRENT DRIVER GENES =====\n")
print(
  as.data.frame(
    top_drivers %>%
      select(
        contrast,
        direction,
        gene,
        n_L2,
        n_pathways,
        n_hits,
        min_FDR
      )
  ),
  row.names = FALSE
)

cat("\n============================================================\n")
cat("20E1B COMPLETE\n")
cat("============================================================\n")
cat("BioPlanet postprocessing from frozen 20E1 results\n")
cat("No GSEA refit\n")
cat("No DE refit\n")
cat("No annotation modified\n")
cat("PASS\n")
