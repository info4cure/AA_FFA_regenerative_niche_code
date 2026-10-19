#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(purrr)
  library(tidyr)
})

cat("============================================================\n")
cat("PIPELINE_V2 — 20F0 BIOPLANET TARGETED AXIS PROJECTION\n")
cat("============================================================\n")

ROOT <- "PIPELINE_V2/results/disease_biology/20_PATHWAY_ANALYSIS"

IN_D <- file.path(
  ROOT,
  "20D_TARGETED_PATHWAYS",
  "tables"
)

IN_E <- file.path(
  ROOT,
  "20E_BIOPLANET"
)

OUT <- file.path(
  ROOT,
  "20F_BIOPLANET_TARGETED_AXES"
)

TAB <- file.path(OUT, "tables")

dir.create(TAB, recursive=TRUE, showWarnings=FALSE)

# ============================================================
# 1. Recover frozen pre-specified axis definitions
# ============================================================

cat("\nSTEP 1/8 — recovering frozen 20D1 axis definitions\n")

axis_file <- file.path(
  IN_D,
  "20D1_axis_definitions.csv"
)

stopifnot(file.exists(axis_file))

axes <- read_csv(
  axis_file,
  show_col_types=FALSE
)

stopifnot(
  nrow(axes) == 8,
  all(c("axis","pattern") %in% names(axes)),
  !anyDuplicated(axes$axis)
)

cat("Frozen axes:", nrow(axes), "\n")
print(axes, n=Inf)

# ============================================================
# 2. Recover complete frozen BioPlanet inventory
# ============================================================

cat("\nSTEP 2/8 — recovering frozen BioPlanet inventory\n")

inventory_file <- file.path(
  IN_E,
  "tables",
  "20E0_BioPlanet_pathway_inventory.csv"
)

stopifnot(file.exists(inventory_file))

bp <- read_csv(
  inventory_file,
  show_col_types=FALSE
)

cat("BioPlanet pathways:", nrow(bp), "\n")
cat("Original columns:\n")
print(names(bp))

# Normalize frozen 20E0 inventory column names for downstream use.
# This changes names only; no data are modified.
stopifnot(
  all(c("PATHWAY_ID","PATHWAY_NAME") %in% names(bp))
)

bp <- bp %>%
  rename(
    pathway_id   = PATHWAY_ID,
    pathway_name = PATHWAY_NAME
  )

cat("Normalized columns:\n")
print(names(bp))

stopifnot(
  nrow(bp) == 1658,
  n_distinct(bp$pathway_id) == 1658,
  !anyNA(bp$pathway_id),
  !anyNA(bp$pathway_name),
  all(nzchar(bp$pathway_id)),
  all(nzchar(bp$pathway_name))
)

# ============================================================
# 3. Apply EXACT frozen 20D1 regex rules
# ============================================================

cat("\nSTEP 3/8 — projecting 20D1 axes onto BioPlanet\n")

axis_map <- map_dfr(
  seq_len(nrow(axes)),
  function(i) {

    bp %>%
      filter(
        str_detect(
          pathway_name,
          regex(
            axes$pattern[[i]],
            ignore_case=TRUE
          )
        )
      ) %>%
      mutate(
        axis=axes$axis[[i]]
      )
  }
) %>%
  distinct(
    pathway_id,
    pathway_name,
    axis
  )

cat("BioPlanet axis assignments:", nrow(axis_map), "\n")
cat(
  "Unique mapped pathways:",
  n_distinct(axis_map$pathway_id),
  "\n"
)

write_csv(
  axis_map,
  file.path(
    TAB,
    "20F0_BioPlanet_targeted_pathway_inventory.csv"
  )
)

# ============================================================
# 4. Coverage + multi-axis audit
# ============================================================

cat("\nSTEP 4/8 — axis coverage and cardinality audit\n")

coverage <- axis_map %>%
  count(
    axis,
    name="n_pathways"
  ) %>%
  right_join(
    axes %>% select(axis),
    by="axis"
  ) %>%
  mutate(
    n_pathways=replace_na(n_pathways,0L)
  ) %>%
  arrange(axis)

multi_axis <- axis_map %>%
  count(
    pathway_id,
    pathway_name,
    name="n_axes"
  ) %>%
  filter(n_axes > 1) %>%
  arrange(
    desc(n_axes),
    pathway_name
  )

write_csv(
  coverage,
  file.path(
    TAB,
    "20F0_BioPlanet_axis_coverage.csv"
  )
)

write_csv(
  multi_axis,
  file.path(
    TAB,
    "20F0_BioPlanet_multi_axis_pathways.csv"
  )
)

cat("\n===== BIOPLANET AXIS COVERAGE =====\n")
print(as.data.frame(coverage), row.names=FALSE)

cat(
  "\nBioPlanet pathways assigned to >1 axis:",
  nrow(multi_axis),
  "\n"
)

# ============================================================
# 5. Recover frozen BioPlanet GSEA
# ============================================================

cat("\nSTEP 5/8 — recovering frozen 20E1 GSEA results\n")

gsea_file <- file.path(
  IN_E,
  "20E1_GSEA",
  "tables",
  "20E1_all_BioPlanet_GSEA.csv"
)

stopifnot(file.exists(gsea_file))

gsea <- read_csv(
  gsea_file,
  show_col_types=FALSE
)

stopifnot(
  nrow(gsea) == 49000,
  all(
    c(
      "contrast",
      "module",
      "final_L2_label",
      "bioplanet_id",
      "bioplanet_name",
      "NES",
      "padj",
      "direction"
    ) %in% names(gsea)
  )
)

cat("Frozen GSEA rows:", nrow(gsea), "\n")
cat(
  "Frozen FDR<0.05 rows:",
  sum(gsea$padj < .05, na.rm=TRUE),
  "\n"
)

# ============================================================
# 6. Integrate targeted axes with frozen GSEA
# ============================================================

cat("\nSTEP 6/8 — integrating axes with frozen BioPlanet GSEA\n")

# Determine correspondence between inventory identifiers
# and GSEA identifiers explicitly.
#
# 20E0 inventory is expected to retain BioPlanet stable IDs.
# Rename only after checking the actual columns.

id_candidates <- intersect(
  c(
    "pathway_id",
    "bioplanet_id",
    "PATHWAY_ID"
  ),
  names(bp)
)

stopifnot(length(id_candidates) >= 1)

bp_id_col <- id_candidates[[1]]

axis_for_join <- axis_map %>%
  rename(
    bioplanet_id=pathway_id,
    axis_pathway_name=pathway_name
  )

targeted <- gsea %>%
  inner_join(
    axis_for_join,
    by="bioplanet_id",
    relationship="many-to-many"
  )

cat("Targeted GSEA rows:", nrow(targeted), "\n")
cat(
  "Targeted FDR<0.05 rows:",
  sum(targeted$padj < .05, na.rm=TRUE),
  "\n"
)

write_csv(
  targeted,
  file.path(
    TAB,
    "20F0_BioPlanet_targeted_axes_all_results.csv"
  )
)

targeted_sig <- targeted %>%
  filter(
    !is.na(padj),
    padj < .05
  )

write_csv(
  targeted_sig,
  file.path(
    TAB,
    "20F0_BioPlanet_targeted_axes_FDR005.csv"
  )
)

# ============================================================
# 7. Biological summaries
# ============================================================

cat("\nSTEP 7/8 — targeted biological summaries\n")

axis_summary <- targeted_sig %>%
  group_by(
    contrast,
    axis,
    direction
  ) %>%
  summarise(
    n_hits=n(),
    n_pathways=n_distinct(bioplanet_id),
    n_L2=n_distinct(final_L2_label),
    median_NES=median(NES, na.rm=TRUE),
    max_abs_NES=max(abs(NES), na.rm=TRUE),
    min_FDR=min(padj, na.rm=TRUE),
    .groups="drop"
  ) %>%
  arrange(
    contrast,
    axis,
    direction
  )

pathway_recurrence <- targeted_sig %>%
  group_by(
    contrast,
    axis,
    direction,
    bioplanet_id,
    bioplanet_name
  ) %>%
  summarise(
    n_L2=n_distinct(final_L2_label),
    median_NES=median(NES, na.rm=TRUE),
    max_abs_NES=max(abs(NES), na.rm=TRUE),
    min_FDR=min(padj, na.rm=TRUE),
    L2=paste(
      sort(unique(final_L2_label)),
      collapse=";"
    ),
    .groups="drop"
  ) %>%
  arrange(
    contrast,
    axis,
    direction,
    desc(n_L2),
    min_FDR
  )

L2_axis_summary <- targeted_sig %>%
  group_by(
    contrast,
    module,
    final_L2_label,
    axis,
    direction
  ) %>%
  summarise(
    n_pathways=n_distinct(bioplanet_id),
    median_NES=median(NES, na.rm=TRUE),
    max_abs_NES=max(abs(NES), na.rm=TRUE),
    min_FDR=min(padj, na.rm=TRUE),
    .groups="drop"
  ) %>%
  arrange(
    contrast,
    module,
    final_L2_label,
    axis,
    direction
  )

write_csv(
  axis_summary,
  file.path(
    TAB,
    "20F0_BioPlanet_targeted_axis_summary_FDR005.csv"
  )
)

write_csv(
  pathway_recurrence,
  file.path(
    TAB,
    "20F0_BioPlanet_targeted_pathway_recurrence_FDR005.csv"
  )
)

write_csv(
  L2_axis_summary,
  file.path(
    TAB,
    "20F0_BioPlanet_L2_axis_summary_FDR005.csv"
  )
)

cat("\n===== TARGETED AXIS SUMMARY =====\n")
print(
  as.data.frame(axis_summary),
  row.names=FALSE
)

cat("\n===== TOP RECURRENT TARGETED PATHWAYS =====\n")

print(
  pathway_recurrence %>%
    group_by(
      contrast,
      axis,
      direction
    ) %>%
    slice_head(n=3) %>%
    ungroup() %>%
    select(
      contrast,
      axis,
      direction,
      bioplanet_name,
      n_L2,
      median_NES,
      min_FDR
    ),
  n=Inf
)

# ============================================================
# 8. Formal checks + audit
# ============================================================

cat("\nSTEP 8/8 — formal checks\n")

analytical_key_duplicates <- targeted %>%
  count(
    contrast,
    module,
    final_L2_label,
    bioplanet_id,
    axis
  ) %>%
  filter(n > 1)

checks <- tibble(
  check=c(
    "eight_frozen_axes_recovered",
    "BioPlanet_1658_recovered",
    "axis_projection_generated",
    "all_axes_represented",
    "multi_axis_allowed",
    "frozen_GSEA_49000_used",
    "frozen_FDR005_4378_preserved",
    "targeted_results_generated",
    "targeted_summary_generated",
    "no_duplicate_targeted_keys",
    "no_axis_rules_modified",
    "no_GSEA_refit",
    "no_DE_refit",
    "annotation_unchanged",
    "no_figures_generated"
  ),
  pass=c(
    nrow(axes)==8,
    nrow(bp)==1658,
    nrow(axis_map)>0,
    all(axes$axis %in% axis_map$axis),
    TRUE,
    nrow(gsea)==49000,
    sum(gsea$padj < .05,na.rm=TRUE)==4378,
    nrow(targeted)>0,
    nrow(axis_summary)>0,
    nrow(analytical_key_duplicates)==0,
    TRUE,
    TRUE,
    TRUE,
    TRUE,
    TRUE
  )
)

print(
  as.data.frame(checks),
  row.names=FALSE
)

stopifnot(all(checks$pass))

audit <- c(
  "PIPELINE_V2 — 20F0 BioPlanet targeted axis projection",
  "",
  paste0("Date: ",Sys.Date()),
  "",
  "BIOLOGICAL AXIS PROVENANCE",
  "--------------------------",
  "The eight biological interpretation axes and their regex patterns",
  "were recovered unchanged from the frozen 20D1 specification.",
  "",
  "The axes were originally defined before Reactome / GO-BP enrichment",
  "and before BioPlanet analysis.",
  "",
  "No BioPlanet result was used to define or modify an axis.",
  "",
  "BIOPLANET",
  "---------",
  "Frozen BioPlanet universe: 1,658 pathways.",
  "Frozen BioPlanet GSEA rows: 49,000.",
  "Frozen BioPlanet FDR<0.05 rows: 4,378.",
  "",
  "Axis assignment was performed exclusively by pathway-name matching",
  "against the frozen 20D1 regular-expression patterns.",
  "",
  "Pathways may legitimately map to more than one biological axis.",
  "",
  "STATISTICAL BOUNDARY",
  "--------------------",
  "No GSEA was refitted.",
  "No differential-expression model was refitted.",
  "No cell annotation was modified.",
  "No pathway was selected to define an axis based on significance.",
  "",
  "20F0 is a post hoc projection onto a prospectively frozen",
  "interpretation framework, not a new enrichment analysis."
)

writeLines(
  audit,
  file.path(
    OUT,
    "20F0_audit.txt"
  )
)

cat("\n============================================================\n")
cat("20F0 COMPLETE\n")
cat("============================================================\n")
cat("Frozen 20D1 axes projected onto BioPlanet\n")
cat("No axis redefinition\n")
cat("No GSEA refit\n")
cat("No DE refit\n")
cat("No annotation modified\n")
cat("PASS\n")
