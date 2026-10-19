suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(readr)
  library(tibble)
})

cat("============================================================\n")
cat("PIPELINE_V2 — 18B BUILD FROZEN-L2 PSEUDOBULK COUNTS\n")
cat("============================================================\n")

BASE <- paste0(
  "PIPELINE_V2/results/disease_biology/",
  "18_PSEUDOBULK_DE/18B_pseudobulk"
)

TAB <- file.path(BASE, "tables")
OBJ <- file.path(BASE, "objects")

dir.create(TAB, recursive = TRUE, showWarnings = FALSE)
dir.create(OBJ, recursive = TRUE, showWarnings = FALSE)

MASTER <- readRDS(
  paste0(
    "PIPELINE_V2/results/disease_biology/",
    "17A_frozen_analysis_atlas/tables/",
    "17A_master_cell_ledger_FINAL.rds"
  )
)

ELIG <- read_csv(
  paste0(
    "PIPELINE_V2/results/disease_biology/",
    "18_PSEUDOBULK_DE/18A_preflight/tables/",
    "18A_pseudobulk_contrast_eligibility.csv"
  ),
  show_col_types = FALSE
)

objects <- tribble(
  ~module, ~path,

  "07_EPITHELIAL_L2",
  "PIPELINE_V2/objects/annotated/epithelial_v2_L2_frozen.rds",

  "08_IMMUNE_L2",
  "PIPELINE_V2/objects/annotated/immune_v2_08I_L2_frozen.rds",

  "09_FIBROBLAST_L2",
  paste0(
    "PIPELINE_V2/objects/annotated/",
    "fibroblast_v2_09H5B_L2_frozen_4685.rds"
  ),

  "10_VASCULAR_L2",
  paste0(
    "PIPELINE_V2/objects/annotated/",
    "vascular_v2_10B6_L2_frozen_3023.rds"
  ),

  "11_MURAL_L2",
  paste0(
    "PIPELINE_V2/results/annotation/11_MURAL_L2/",
    "11B5_freeze/mural_v2_11B5_L2_frozen.rds"
  ),

  "12_SWEAT_GLAND_L2",
  paste0(
    "PIPELINE_V2/results/annotation/12_SWEAT_GLAND_L2/",
    "12D_freeze/sweat_gland_v2_12D_L2_frozen.rds"
  ),

  "13_MELANOCYTE_L2",
  paste0(
    "PIPELINE_V2/results/annotation/13_MELANOCYTE_L2/",
    "13D_freeze/melanocyte_v2_13D_L2_frozen.rds"
  ),

  "14_SEBOCYTE_L2",
  paste0(
    "PIPELINE_V2/results/annotation/14_SEBOCYTE_L2/",
    "14D_freeze/sebocyte_v2_14D_L2_frozen.rds"
  )
)

# ============================================================
# 1. Eligible frozen-L2 cells
# ============================================================

cat("\nSTEP 1/8 — frozen-L2 master universe\n")

ledger <- MASTER %>%
  filter(
    analysis_eligible_global,
    analysis_eligible_L2,
    final_L2_status == "FROZEN"
  ) %>%
  select(
    cell,
    patient_id,
    disease,
    batch,
    final_compartment,
    final_L2_label,
    annotation_module
  )

stopifnot(
  nrow(ledger) == 82494,
  n_distinct(ledger$cell) == 82494
)

cat("Cells:", nrow(ledger), "\n")
cat("L2:", n_distinct(ledger$final_L2_label), "\n")


# ============================================================
# 2. Recover RNA count matrices from frozen objects
# ============================================================

cat("\nSTEP 2/8 — recovering RNA counts\n")

count_parts <- list()
cell_seen <- character()

for (i in seq_len(nrow(objects))) {

  module <- objects$module[i]
  f <- objects$path[i]

  cat(module, "\n")

  o <- readRDS(f)

  DefaultAssay(o) <- "RNA"

  cnt <- LayerData(
    o,
    assay = "RNA",
    layer = "counts"
  )

  cells <- colnames(cnt)

  stopifnot(
    !anyDuplicated(cells),
    all(cells %in% ledger$cell)
  )

  if (length(intersect(cell_seen, cells)) > 0) {
    stop("Cells duplicated across frozen L2 objects")
  }

  cell_seen <- c(cell_seen, cells)

  count_parts[[module]] <- cnt

  rm(o, cnt)
  gc()
}

stopifnot(
  length(cell_seen) == 82494,
  setequal(cell_seen, ledger$cell)
)

cat("All 82,494 frozen-L2 cells recovered exactly once.\n")


# ============================================================
# 3. Build pseudobulk independently within each module
# ============================================================

cat("\nSTEP 3/8 — aggregating patient × L2 counts\n")

pb_list <- list()
meta_list <- list()

for (module in names(count_parts)) {

  cnt <- count_parts[[module]]

  md <- ledger %>%
    filter(annotation_module == module) %>%
    filter(cell %in% colnames(cnt))

  stopifnot(
    nrow(md) == ncol(cnt)
  )

  md <- md[
    match(colnames(cnt), md$cell),
    ,
    drop = FALSE
  ]

  stopifnot(
    identical(md$cell, colnames(cnt))
  )

  groups <- paste(
    md$patient_id,
    md$final_L2_label,
    sep = "|||"
  )

  ug <- unique(groups)

  # Sparse incidence matrix:
  # cells × pseudobulk samples
  G <- sparseMatrix(
    i = seq_along(groups),
    j = match(groups, ug),
    x = 1,
    dims = c(
      length(groups),
      length(ug)
    )
  )

  colnames(G) <- ug

  pb <- cnt %*% G

  colnames(pb) <- paste0(
    module,
    "|||",
    ug
  )

  parts <- strsplit(
    ug,
    "\\|\\|\\|"
  )

  pm <- tibble(
    module = module,
    pseudobulk_id = colnames(pb),
    patient_id = vapply(parts, `[`, "", 1),
    final_L2_label = vapply(parts, `[`, "", 2)
  ) %>%
    left_join(
      ledger %>%
        distinct(
          patient_id,
          disease,
          batch
        ),
      by = "patient_id"
    )

  ncell <- table(groups)

  pm$n_cells <- as.integer(
    ncell[
      paste(
        pm$patient_id,
        pm$final_L2_label,
        sep = "|||"
      )
    ]
  )

  pm$library_size <- Matrix::colSums(pb)
  pm$genes_detected <- Matrix::colSums(pb > 0)

  pb_list[[module]] <- pb
  meta_list[[module]] <- pm
}

pb_meta <- bind_rows(meta_list)

stopifnot(
  !anyDuplicated(pb_meta$pseudobulk_id)
)

write_csv(
  pb_meta,
  file.path(
    TAB,
    "18B_pseudobulk_metadata_ALL.csv"
  )
)


# ============================================================
# 4. Conservation audit
# ============================================================

cat("\nSTEP 4/8 — pseudobulk conservation audit\n")

cell_conservation <- pb_meta %>%
  summarise(
    pseudobulk_samples = n(),
    summed_cells = sum(n_cells),
    min_cells = min(n_cells),
    median_cells = median(n_cells),
    max_cells = max(n_cells)
  )

print(
  as.data.frame(cell_conservation),
  row.names = FALSE
)

stopifnot(
  cell_conservation$summed_cells == 82494
)

module_conservation <- pb_meta %>%
  group_by(module) %>%
  summarise(
    n_pseudobulks = n(),
    summed_cells = sum(n_cells),
    median_cells = median(n_cells),
    median_library_size = median(library_size),
    median_genes_detected = median(genes_detected),
    .groups = "drop"
  )

print(
  as.data.frame(module_conservation),
  row.names = FALSE
)

write_csv(
  module_conservation,
  file.path(
    TAB,
    "18B_module_conservation.csv"
  )
)


# ============================================================
# 5. Attach contrast eligibility
# ============================================================

cat("\nSTEP 5/8 — contrast-specific pseudobulk manifests\n")

manifest_list <- list()

for (cc in unique(ELIG$contrast)) {

  e <- ELIG %>%
    filter(
      contrast == cc,
      eligible_ge20
    )

  groups <- unique(
    c(
      e$numerator,
      e$denominator
    )
  )

  m <- pb_meta %>%
    filter(
      disease %in% groups
    ) %>%
    inner_join(
      e %>%
        select(
          contrast,
          primary,
          numerator,
          denominator,
          final_compartment,
          final_L2_label,
          eligible_ge20,
          robust_ge50
        ),
      by = "final_L2_label",
      relationship = "many-to-many"
    ) %>%
    filter(
      n_cells >= 20
    )

  manifest_list[[cc]] <- m

  write_csv(
    m,
    file.path(
      TAB,
      paste0(
        "18B_manifest_",
        cc,
        ".csv"
      )
    )
  )

  cat(
    cc,
    ":",
    n_distinct(m$final_L2_label),
    "L2;",
    nrow(m),
    "pseudobulk samples\n"
  )
}

manifest <- bind_rows(manifest_list)

write_csv(
  manifest,
  file.path(
    TAB,
    "18B_contrast_manifest_ALL.csv"
  )
)


# ============================================================
# 6. QC summaries
# ============================================================

cat("\nSTEP 6/8 — pseudobulk QC\n")

qc <- pb_meta %>%
  group_by(
    module,
    final_L2_label,
    disease
  ) %>%
  summarise(
    n_pseudobulks = n(),
    median_cells = median(n_cells),
    min_cells = min(n_cells),
    max_cells = max(n_cells),
    median_library_size = median(library_size),
    min_library_size = min(library_size),
    median_genes_detected = median(genes_detected),
    .groups = "drop"
  )

write_csv(
  qc,
  file.path(
    TAB,
    "18B_pseudobulk_QC_summary.csv"
  )
)

cat("\nLowest library-size pseudobulks:\n")

print(
  as.data.frame(
    pb_meta %>%
      arrange(library_size) %>%
      select(
        module,
        patient_id,
        disease,
        final_L2_label,
        n_cells,
        library_size,
        genes_detected
      ) %>%
      slice_head(n = 20)
  ),
  row.names = FALSE
)


# ============================================================
# 7. Save sparse count matrices
# ============================================================

cat("\nSTEP 7/8 — saving pseudobulk matrices\n")

for (module in names(pb_list)) {

  saveRDS(
    pb_list[[module]],
    file.path(
      OBJ,
      paste0(
        "18B_",
        module,
        "_pseudobulk_counts.rds"
      )
    )
  )
}

saveRDS(
  pb_meta,
  file.path(
    OBJ,
    "18B_pseudobulk_metadata_ALL.rds"
  )
)


# ============================================================
# 8. Formal checks
# ============================================================

cat("\nSTEP 8/8 — formal checks\n")

manifest_counts <- manifest %>%
  distinct(
    contrast,
    final_L2_label
  ) %>%
  count(
    contrast,
    name = "n_L2"
  )

expected <- tibble(
  contrast = c(
    "AA_vs_HC",
    "FFA_vs_HC",
    "AA_vs_FFA"
  ),
  expected_n = c(
    24,
    18,
    18
  )
)

mc <- expected %>%
  left_join(
    manifest_counts,
    by = "contrast"
  )

checks <- tibble(
  check = c(
    "frozen_L2_cells_82494",
    "all_frozen_cells_recovered_once",
    "eight_modules",
    "counts_are_nonnegative",
    "counts_are_integer_like",
    "pseudobulk_cell_conservation_82494",
    "AA_HC_24_L2",
    "FFA_HC_18_L2",
    "AA_FFA_18_L2",
    "all_manifest_samples_ge20_cells",
    "metadata_complete",
    "no_DE_performed",
    "no_gene_filtering",
    "no_figures_generated",
    "annotation_unchanged"
  ),

  pass = c(
    nrow(ledger) == 82494,

    length(cell_seen) == 82494 &&
      !anyDuplicated(cell_seen),

    length(pb_list) == 8,

    all(
      vapply(
        pb_list,
        function(m)
          all(m@x >= 0),
        logical(1)
      )
    ),

    all(
      vapply(
        pb_list,
        function(m)
          all(
            abs(m@x - round(m@x)) <
              1e-8
          ),
        logical(1)
      )
    ),

    sum(pb_meta$n_cells) == 82494,

    mc$n_L2[
      mc$contrast == "AA_vs_HC"
    ] == 24,

    mc$n_L2[
      mc$contrast == "FFA_vs_HC"
    ] == 18,

    mc$n_L2[
      mc$contrast == "AA_vs_FFA"
    ] == 18,

    all(manifest$n_cells >= 20),

    all(
      !is.na(pb_meta$patient_id) &
      !is.na(pb_meta$disease) &
      !is.na(pb_meta$batch) &
      !is.na(pb_meta$final_L2_label)
    ),

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

write_csv(
  checks,
  file.path(
    TAB,
    "18B_formal_checks.csv"
  )
)

stopifnot(all(checks$pass))

cat("\n============================================================\n")
cat("18B COMPLETE\n")
cat("============================================================\n")
cat("Frozen-L2 cells aggregated: 82,494\n")
cat("Pseudobulk unit: patient × frozen L2\n")
cat("Expression source: RNA counts\n")
cat("Eligibility carried forward from 18A\n")
cat("Minimum cells per DE pseudobulk: 20\n")
cat("Sparse count matrices saved by module\n")
cat("DE performed: NO\n")
cat("Gene filtering: NO\n")
cat("Figures generated: NO\n")
cat("Annotation modified: NO\n")
cat("READY FOR 18C PSEUDOBULK DE\n")
cat("PASS\n")
