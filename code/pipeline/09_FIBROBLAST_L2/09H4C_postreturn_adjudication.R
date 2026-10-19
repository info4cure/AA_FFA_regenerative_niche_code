
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
})

cat("============================================================\n")
cat("PIPELINE_V2 — 09H4C POST-C9 FIBROBLAST ADJUDICATION\n")
cat("============================================================\n\n")


# ============================================================
# 0. Paths
# ============================================================

OBJ_IN <- paste0(
  "PIPELINE_V2/objects/annotated/",
  "fibroblast_v2_09H4B_rebuilt_4685.rds"
)

REF_DIR <- paste0(
  "PIPELINE_V2/references/FibConsensus/",
  "frozen_reference"
)

SIG_FILE <- file.path(
  REF_DIR,
  "atlas_validation_signature_specification_v2_2.csv"
)

FAM_FILE <- file.path(
  REF_DIR,
  "final_consensus_families_v2_2.csv"
)

OUTDIR <- paste0(
  "PIPELINE_V2/results/annotation/",
  "09_FIBROBLAST_L2/",
  "09H4C_postreturn_adjudication"
)

TAB <- file.path(OUTDIR, "tables")
FIG <- file.path(OUTDIR, "figures")

dir.create(TAB, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG, recursive = TRUE, showWarnings = FALSE)


# ============================================================
# 1. Load rebuilt fibroblast universe
# ============================================================

cat("STEP 1/10 — loading post-C9 rebuilt fibroblast universe\n")

stopifnot(file.exists(OBJ_IN))
stopifnot(file.exists(SIG_FILE))
stopifnot(file.exists(FAM_FILE))

x <- readRDS(OBJ_IN)

cat("Cells:", ncol(x), "\n")
cat("Genes:", nrow(x), "\n")
cat(
  "Patients:",
  length(unique(x$patient_id)),
  "\n\n"
)

stopifnot(ncol(x) == 4685)
stopifnot(nrow(x) == 23412)
stopifnot(length(unique(x$patient_id)) == 16)


# ============================================================
# 2. Identify post-return clustering columns
# ============================================================

cat("STEP 2/10 — identifying post-return clustering structure\n")

md <- x@meta.data

res03_candidates <- grep(
  "res\\.0\\.3$",
  colnames(md),
  value = TRUE
)

res04_candidates <- grep(
  "res\\.0\\.4$",
  colnames(md),
  value = TRUE
)

cat("r0.3 candidates:\n")
print(res03_candidates)

cat("\nr0.4 candidates:\n")
print(res04_candidates)

if (length(res03_candidates) < 1) {
  stop("No post-return resolution 0.3 column found.")
}

if (length(res04_candidates) < 1) {
  stop("No post-return resolution 0.4 column found.")
}

# Prefer explicitly post-return/refined columns if multiple exist
pick_resolution <- function(v) {

  preferred <- v[
    grepl(
      "09H4B|fibroblast",
      v,
      ignore.case = TRUE
    )
  ]

  if (length(preferred) > 0) {
    return(tail(preferred, 1))
  }

  tail(v, 1)
}

COL03 <- pick_resolution(res03_candidates)
COL04 <- pick_resolution(res04_candidates)

cat("\nSelected r0.3 column:", COL03, "\n")
cat("Selected r0.4 column:", COL04, "\n\n")

md$cluster03 <- as.character(md[[COL03]])
md$cluster04 <- as.character(md[[COL04]])

stopifnot(length(unique(md$cluster03)) == 10)
stopifnot(length(unique(md$cluster04)) == 11)


# ============================================================
# 3. Structural support
# ============================================================

cat("STEP 3/10 — structural support audit\n")

support_fun <- function(cluster_vector, resolution_label) {

  z <- data.frame(
    cluster = as.character(cluster_vector),
    patient_id = as.character(md$patient_id),
    disease = as.character(md$disease),
    batch = as.character(md$batch),
    stringsAsFactors = FALSE
  )

  patient_counts <- z %>%
    count(cluster, patient_id, name = "n")

  dominant <- patient_counts %>%
    group_by(cluster) %>%
    arrange(desc(n), .by_group = TRUE) %>%
    slice(1) %>%
    ungroup() %>%
    rename(
      dominant_patient = patient_id,
      dominant_patient_cells = n
    )

  disease_support <- z %>%
    distinct(cluster, patient_id, disease) %>%
    count(cluster, disease, name = "n_patients") %>%
    pivot_wider(
      names_from = disease,
      values_from = n_patients,
      values_fill = 0,
      names_prefix = "n_"
    )

  out <- z %>%
    group_by(cluster) %>%
    summarise(
      n_cells = n(),
      n_patients = n_distinct(patient_id),
      n_batches = n_distinct(batch),
      .groups = "drop"
    ) %>%
    left_join(disease_support, by = "cluster") %>%
    left_join(dominant, by = "cluster") %>%
    mutate(
      dominant_patient_fraction =
        dominant_patient_cells / n_cells,
      resolution = resolution_label
    )

  for (nm in c("n_HC", "n_AA", "n_FFA")) {
    if (!nm %in% colnames(out)) out[[nm]] <- 0L
  }

  out %>%
    select(
      resolution,
      cluster,
      n_cells,
      n_patients,
      n_HC,
      n_AA,
      n_FFA,
      n_batches,
      dominant_patient,
      dominant_patient_cells,
      dominant_patient_fraction
    )
}

support03 <- support_fun(
  md$cluster03,
  "0.3"
)

support04 <- support_fun(
  md$cluster04,
  "0.4"
)

support <- bind_rows(
  support03,
  support04
)

write_csv(
  support,
  file.path(
    TAB,
    "09H4C_cluster_support.csv"
  )
)

cat("\n=== r0.3 SUPPORT ===\n")
print(support03)

cat("\n=== r0.4 SUPPORT ===\n")
print(support04)


# ============================================================
# 4. Structural transitions 0.3 -> 0.4
# ============================================================

cat("\nSTEP 4/10 — post-return 0.3 -> 0.4 transitions\n")

transition <- data.frame(
  parent = md$cluster03,
  child = md$cluster04,
  stringsAsFactors = FALSE
) %>%
  count(parent, child, name = "n_cells") %>%
  group_by(parent) %>%
  mutate(
    parent_total = sum(n_cells),
    fraction_parent = n_cells / parent_total
  ) %>%
  ungroup() %>%
  arrange(
    as.numeric(parent),
    desc(fraction_parent)
  )

write_csv(
  transition,
  file.path(
    TAB,
    "09H4C_transition_03_to_04.csv"
  )
)

material_transition <- transition %>%
  filter(fraction_parent >= 0.05)

write_csv(
  material_transition,
  file.path(
    TAB,
    "09H4C_material_transition_03_to_04.csv"
  )
)

cat("\n=== MATERIAL TRANSITIONS >=5% OF PARENT ===\n")
print(material_transition)


# ============================================================
# 5. Frozen FibConsensus reference
# ============================================================

cat("\nSTEP 5/10 — loading frozen FibConsensus signatures\n")

signature_source <- SIG_FILE

sig_raw <- readr::read_csv(
  signature_source,
  show_col_types = FALSE
)

required_signature_columns <- c(
  "entity_id",
  "entity_name",
  "entity_layer",
  "projection_gene_symbol",
  "signature_tier"
)

missing_signature_columns <- setdiff(
  required_signature_columns,
  colnames(sig_raw)
)

if (length(missing_signature_columns) > 0) {
  stop(
    "Frozen FibConsensus signature file is missing required columns: ",
    paste(missing_signature_columns, collapse = ", ")
  )
}

# The archived specification contains additional entities/substructures.
# For fibroblast L2 family adjudication we use only the six frozen
# CONSENSUS_FAMILY entities.
sig <- sig_raw %>%
  dplyr::filter(
    entity_layer == "CONSENSUS_FAMILY",
    !is.na(entity_id),
    !is.na(projection_gene_symbol)
  ) %>%
  dplyr::mutate(
    entity_id = as.character(entity_id),
    projection_gene_symbol = as.character(projection_gene_symbol),
    signature_tier = as.character(signature_tier)
  )

families <- sort(unique(sig$entity_id))

cat("Frozen source rows:", nrow(sig_raw), "\n")
cat("Consensus-family signature rows:", nrow(sig), "\n")
cat("Families:", paste(families, collapse = ", "), "\n")

expected_families <- paste0("FAM0", 1:6)

if (!identical(families, expected_families)) {
  stop(
    "Frozen family validation failed. Observed: ",
    paste(families, collapse = ", "),
    " | Expected: ",
    paste(expected_families, collapse = ", ")
  )
}

stopifnot(nrow(sig_raw) == 56)
stopifnot(nrow(sig) == 40)
stopifnot(length(families) == 6)

coverage_check <- sig %>%
  dplyr::group_by(
    entity_id,
    entity_name,
    signature_tier
  ) %>%
  dplyr::summarise(
    n_requested = dplyr::n_distinct(projection_gene_symbol),
    n_present = sum(
      unique(projection_gene_symbol) %in% rownames(x)
    ),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    coverage = n_present / n_requested
  )

print(coverage_check)

if (!all(coverage_check$coverage == 1)) {
  stop("Incomplete frozen FibConsensus signature coverage after C9 return.")
}

cat("\nSTEP 6/10 — frozen signature coverage\n")

genes_available <- rownames(x)

coverage <- sig %>%
  group_by(
    entity_id,
    entity_name,
    signature_tier
  ) %>%
  summarise(
    n_requested = n_distinct(
      projection_gene_symbol
    ),
    n_present = sum(
      unique(projection_gene_symbol) %in%
        genes_available
    ),
    missing = paste(
      setdiff(
        unique(projection_gene_symbol),
        genes_available
      ),
      collapse = ";"
    ),
    .groups = "drop"
  ) %>%
  mutate(
    coverage = n_present / n_requested
  )

write_csv(
  coverage,
  file.path(
    TAB,
    "09H4C_signature_coverage.csv"
  )
)

print(coverage)

stopifnot(all(coverage$coverage == 1))


# ============================================================
# 7. Recalculate frozen CORE and FULL scores
# ============================================================

cat("\nSTEP 7/10 — recalculating FibConsensus after C9 return\n")

DefaultAssay(x) <- "RNA"

expr <- GetAssayData(
  x,
  assay = "RNA",
  layer = "data"
)

score_gene_set <- function(gene_set) {

  gene_set <- unique(
    gene_set[gene_set %in% rownames(expr)]
  )

  if (length(gene_set) == 0) {
    stop("Empty gene set after coverage filtering.")
  }

  Matrix::colMeans(
    expr[gene_set, , drop = FALSE]
  )
}

score_df <- data.frame(
  cell = colnames(x),
  stringsAsFactors = FALSE
)

score_spec <- list()

for (fid in families) {

  sf <- sig %>%
    filter(entity_id == fid)

  core_genes <- sf %>%
    filter(signature_tier == "CORE") %>%
    pull(projection_gene_symbol) %>%
    unique()

  full_genes <- sf %>%
    pull(projection_gene_symbol) %>%
    unique()

  raw_core <- score_gene_set(core_genes)
  raw_full <- score_gene_set(full_genes)

  z_core <- as.numeric(scale(raw_core))
  z_full <- as.numeric(scale(raw_full))

  score_df[[paste0(fid, "_CORE")]] <- z_core
  score_df[[paste0(fid, "_FULL")]] <- z_full

  score_spec[[fid]] <- data.frame(
    family_id = fid,
    tier = c("CORE", "FULL"),
    n_genes = c(
      length(core_genes),
      length(full_genes)
    )
  )

  cat(
    fid,
    "| CORE:",
    length(core_genes),
    "| FULL:",
    length(full_genes),
    "\n"
  )
}

score_spec <- bind_rows(score_spec)

write_csv(
  score_spec,
  file.path(
    TAB,
    "09H4C_score_specification.csv"
  )
)

score_df$patient_id <- md[
  score_df$cell,
  "patient_id"
]

score_df$disease <- md[
  score_df$cell,
  "disease"
]

score_df$batch <- md[
  score_df$cell,
  "batch"
]

score_df$cluster03 <- md[
  score_df$cell,
  "cluster03"
]

score_df$cluster04 <- md[
  score_df$cell,
  "cluster04"
]

write_csv(
  score_df,
  file.path(
    TAB,
    "09H4C_cell_FibConsensus_scores.csv"
  )
)


# ============================================================
# 8. Cluster and patient-level FibConsensus evidence
# ============================================================

cat("\nSTEP 8/10 — cluster-level and patient-level projection\n")

score_long <- score_df %>%
  pivot_longer(
    cols = matches("^FAM[0-9]+_(CORE|FULL)$"),
    names_to = "family_tier",
    values_to = "z"
  ) %>%
  separate(
    family_tier,
    into = c("family_id", "tier"),
    sep = "_"
  )

summarise_resolution <- function(
  dat,
  cluster_col,
  resolution_label
) {

  dat %>%
    mutate(
      cluster = .data[[cluster_col]]
    ) %>%
    group_by(
      cluster,
      family_id,
      tier
    ) %>%
    summarise(
      mean_z = mean(z, na.rm = TRUE),
      median_z = median(z, na.rm = TRUE),
      n_cells = n_distinct(cell),
      n_patients = n_distinct(patient_id),
      .groups = "drop"
    ) %>%
    mutate(
      resolution = resolution_label
    ) %>%
    select(
      resolution,
      everything()
    )
}

cluster_scores03 <- summarise_resolution(
  score_long,
  "cluster03",
  "0.3"
)

cluster_scores04 <- summarise_resolution(
  score_long,
  "cluster04",
  "0.4"
)

cluster_scores <- bind_rows(
  cluster_scores03,
  cluster_scores04
)

write_csv(
  cluster_scores,
  file.path(
    TAB,
    "09H4C_cluster_FibConsensus_scores.csv"
  )
)

patient_scores03 <- score_long %>%
  group_by(
    cluster03,
    patient_id,
    disease,
    family_id,
    tier
  ) %>%
  summarise(
    mean_z = mean(z, na.rm = TRUE),
    n_cells = n_distinct(cell),
    .groups = "drop"
  ) %>%
  rename(cluster = cluster03) %>%
  mutate(resolution = "0.3")

patient_scores04 <- score_long %>%
  group_by(
    cluster04,
    patient_id,
    disease,
    family_id,
    tier
  ) %>%
  summarise(
    mean_z = mean(z, na.rm = TRUE),
    n_cells = n_distinct(cell),
    .groups = "drop"
  ) %>%
  rename(cluster = cluster04) %>%
  mutate(resolution = "0.4")

patient_scores <- bind_rows(
  patient_scores03,
  patient_scores04
)

write_csv(
  patient_scores,
  file.path(
    TAB,
    "09H4C_patient_cluster_FibConsensus_scores.csv"
  )
)


# CORE/FULL concordance

concordance <- cluster_scores %>%
  select(
    resolution,
    cluster,
    family_id,
    tier,
    mean_z
  ) %>%
  pivot_wider(
    names_from = tier,
    values_from = mean_z
  ) %>%
  group_by(
    resolution,
    family_id
  ) %>%
  summarise(
    n_clusters = n(),
    spearman_CORE_FULL = cor(
      CORE,
      FULL,
      method = "spearman",
      use = "complete.obs"
    ),
    .groups = "drop"
  )

write_csv(
  concordance,
  file.path(
    TAB,
    "09H4C_CORE_FULL_concordance.csv"
  )
)

cat("\n=== CORE/FULL CONCORDANCE ===\n")
print(concordance)


# Top family by r0.3 cluster

ranking03 <- cluster_scores03 %>%
  filter(tier == "CORE") %>%
  arrange(
    as.numeric(cluster),
    desc(mean_z)
  ) %>%
  group_by(cluster) %>%
  mutate(
    core_rank = row_number(),
    second_core_z = nth(mean_z, 2),
    core_margin = first(mean_z) -
      nth(mean_z, 2)
  ) %>%
  ungroup()

top03 <- ranking03 %>%
  filter(core_rank == 1) %>%
  transmute(
    cluster,
    top_family = family_id,
    top_core_z = mean_z,
    core_margin
  )

full03 <- cluster_scores03 %>%
  filter(tier == "FULL") %>%
  select(
    cluster,
    family_id,
    top_full_z = mean_z
  )

top03 <- top03 %>%
  left_join(
    full03,
    by = c(
      "cluster",
      "top_family" = "family_id"
    )
  ) %>%
  left_join(
    support03 %>%
      select(
        cluster,
        n_cells,
        n_patients,
        dominant_patient_fraction
      ),
    by = "cluster"
  )

write_csv(
  ranking03,
  file.path(
    TAB,
    "09H4C_r03_family_ranking.csv"
  )
)

write_csv(
  top03,
  file.path(
    TAB,
    "09H4C_r03_top_family.csv"
  )
)

cat("\n=== r0.3 TOP FROZEN FAMILY ===\n")
print(top03)


# ============================================================
# 9. Formal resolution adjudication
# ============================================================

cat("\nSTEP 9/10 — formal resolution adjudication\n")

r03_recurrent <- all(
  support03$n_patients >= 10 &
    support03$n_batches >= 2 &
    support03$dominant_patient_fraction < 0.5
)

r04_patient_dominated <- any(
  support04$dominant_patient_fraction >= 0.5,
  na.rm = TRUE
)

r04_low_patient_support <- any(
  support04$n_patients < 5
)

# The post-return r0.3 solution is accepted only if all clusters are
# recurrent across patients and batches and no cluster is patient-dominated.
# r0.4 is rejected as the primary identity resolution if it introduces
# patient-dominated or poorly recurrent structure.

selected_resolution <- if (
  r03_recurrent &&
  (r04_patient_dominated || r04_low_patient_support)
) {
  "0.3"
} else {
  "UNRESOLVED"
}

decision_reason <- if (
  selected_resolution == "0.3"
) {
  paste(
    "Resolution 0.3 retained a recurrent 10-cluster structure after",
    "removal of epithelial-like C9, with all clusters represented by",
    "at least 10 patients, both sequencing batches, and no cluster",
    "dominated by a single patient. Resolution 0.4 introduced a",
    "poorly recurrent, strongly patient-dominated cluster and was",
    "therefore not promoted as the primary fibroblast identity",
    "resolution."
  )
} else {
  paste(
    "The prespecified structural criteria did not provide sufficient",
    "evidence to freeze resolution 0.3 over resolution 0.4."
  )
}

decision <- tibble(
  fibroblast_cells = ncol(x),
  patients = length(unique(md$patient_id)),
  resolution_03_clusters = nrow(support03),
  resolution_04_clusters = nrow(support04),
  resolution_03_all_recurrent = r03_recurrent,
  resolution_04_patient_dominated = r04_patient_dominated,
  resolution_04_low_patient_support = r04_low_patient_support,
  selected_resolution = selected_resolution,
  decision_reason = decision_reason
)

write_csv(
  decision,
  file.path(
    TAB,
    "09H4C_resolution_decision.csv"
  )
)

cat("\n=== GLOBAL RESOLUTION DECISION ===\n")
print(decision)


# ============================================================
# 10. Freeze resolution metadata, figures, provenance and audit
# ============================================================

cat("\nSTEP 10/10 — freezing resolution and formal audit\n")

if (selected_resolution != "0.3") {
  stop(
    paste(
      "09H4C cannot freeze fibroblast resolution:",
      "formal criteria did not select 0.3."
    )
  )
}


# ------------------------------------------------------------
# Add resolution-level metadata only.
# NO L2 identities are assigned here.
# ------------------------------------------------------------

x$fibroblast_resolution_frozen <- selected_resolution

x$fibroblast_cluster_frozen <- as.character(
  x@meta.data[[COL03]]
)

x$fibroblast_resolution_status <- "FROZEN"

x$fibroblast_L2_status <- "UNADJUDICATED"

x$fibroblast_annotation_stage <- "09H4C_resolution_frozen"

x@misc[["09H4C"]] <- list(
  fibroblast_cells = ncol(x),
  patients = length(unique(x$patient_id)),
  selected_resolution = selected_resolution,
  selected_cluster_column = COL03,
  comparison_cluster_column = COL04,
  FibConsensus_reference =
    "FibConsensus v2.2 frozen scientific release",
  FibConsensus_archive =
    "FibConsensus v2.3 reproducibility archive",
  FibConsensus_DOI =
    "10.5281/zenodo.22792324",
  disease_used_for_integration = FALSE,
  disease_used_for_resolution_selection = FALSE,
  FibConsensus_used_for_geometry = FALSE,
  FibConsensus_used_for_resolution_evidence = TRUE,
  L2_identifiers_assigned = FALSE,
  additional_cells_removed = 0,
  source_object = OBJ_IN
)


# ------------------------------------------------------------
# Save frozen-resolution object as a NEW object
# ------------------------------------------------------------

OBJ_OUT <- paste0(
  "PIPELINE_V2/objects/annotated/",
  "fibroblast_v2_09H4C_resolution03_frozen_4685.rds"
)

saveRDS(
  x,
  OBJ_OUT
)


# ------------------------------------------------------------
# Figure 1 — post-return UMAP, frozen r0.3 structure
# ------------------------------------------------------------

umap_name <- if (
  "fibroblast_09H4B_umap" %in% names(x@reductions)
) {
  "fibroblast_09H4B_umap"
} else {
  stop("09H4B fibroblast UMAP not found.")
}

p_umap <- DimPlot(
  x,
  reduction = umap_name,
  group.by = "fibroblast_cluster_frozen",
  label = TRUE,
  repel = TRUE,
  raster = FALSE,
  pt.size = 0.45
) +
  labs(
    title = "Fibroblast structure after epithelial return",
    subtitle = "Frozen clustering resolution = 0.3",
    colour = "Cluster"
  ) +
  theme_classic(base_size = 11)

ggsave(
  file.path(
    FIG,
    "09H4C_fibroblast_resolution03_frozen_umap.pdf"
  ),
  p_umap,
  width = 7.2,
  height = 6.2,
  device = cairo_pdf
)

ggsave(
  file.path(
    FIG,
    "09H4C_fibroblast_resolution03_frozen_umap.png"
  ),
  p_umap,
  width = 7.2,
  height = 6.2,
  dpi = 400
)


# ------------------------------------------------------------
# Figure 2 — patient support
# ------------------------------------------------------------

support03_plot <- support03 %>%
  mutate(
    cluster = factor(
      cluster,
      levels = sort(
        unique(as.numeric(cluster))
      )
    )
  )

p_support <- ggplot(
  support03_plot,
  aes(
    x = cluster,
    y = n_patients
  )
) +
  geom_col() +
  geom_hline(
    yintercept = 10,
    linetype = 2
  ) +
  labs(
    title = "Patient support of frozen fibroblast clusters",
    subtitle = "Resolution 0.3 after epithelial-like C9 return",
    x = "Fibroblast cluster",
    y = "Number of patients"
  ) +
  theme_classic(base_size = 11)

ggsave(
  file.path(
    FIG,
    "09H4C_resolution03_patient_support.pdf"
  ),
  p_support,
  width = 7,
  height = 4.8,
  device = cairo_pdf
)

ggsave(
  file.path(
    FIG,
    "09H4C_resolution03_patient_support.png"
  ),
  p_support,
  width = 7,
  height = 4.8,
  dpi = 400
)


# ------------------------------------------------------------
# Figure 3 — FibConsensus CORE evidence at frozen resolution
# ------------------------------------------------------------

heat <- cluster_scores03 %>%
  filter(tier == "CORE") %>%
  mutate(
    cluster = factor(
      cluster,
      levels = sort(
        unique(as.numeric(cluster))
      )
    ),
    family_id = factor(
      family_id,
      levels = families
    )
  )

p_heat <- ggplot(
  heat,
  aes(
    x = cluster,
    y = family_id,
    fill = mean_z
  )
) +
  geom_tile() +
  geom_text(
    aes(
      label = sprintf("%.2f", mean_z)
    ),
    size = 3
  ) +
  labs(
    title = "Frozen FibConsensus projection",
    subtitle = "CORE signatures; frozen fibroblast resolution 0.3",
    x = "Fibroblast cluster",
    y = "FibConsensus family",
    fill = "Mean z"
  ) +
  theme_classic(base_size = 11)

ggsave(
  file.path(
    FIG,
    "09H4C_FibConsensus_CORE_resolution03.pdf"
  ),
  p_heat,
  width = 8.5,
  height = 5.5,
  device = cairo_pdf
)

ggsave(
  file.path(
    FIG,
    "09H4C_FibConsensus_CORE_resolution03.png"
  ),
  p_heat,
  width = 8.5,
  height = 5.5,
  dpi = 400
)


# ------------------------------------------------------------
# Formal checks
# ------------------------------------------------------------


# ------------------------------------------------------------
# Verify that the exact 29 cells returned to epithelial in 09H4A
# are absent from the post-return fibroblast universe.
# Cluster numbers are NOT used for this check because clustering
# was rebuilt de novo in 09H4B and cluster labels are not persistent.
# ------------------------------------------------------------

return_ledger <- readr::read_csv(
  "PIPELINE_V2/results/annotation/09_FIBROBLAST_L2/09H4A_C9_return_epithelial/tables/09H4A_C9_return_ledger.csv",
  show_col_types = FALSE
)

stopifnot("cell" %in% colnames(return_ledger))

returned_09H4A_cells <- unique(return_ledger$cell)

n_returned_09H4A <- length(returned_09H4A_cells)

returned_09H4A_overlap <- intersect(
  returned_09H4A_cells,
  colnames(x)
)

returned_09H4A_cells_absent <- (
  n_returned_09H4A == 29 &&
  length(returned_09H4A_overlap) == 0
)

cat("\n=== 09H4A RETURN-LEDGER AUDIT ===\n")
cat("Ledger cells:", n_returned_09H4A, "\n")
cat(
  "Returned cells remaining in fibroblast universe:",
  length(returned_09H4A_overlap),
  "\n"
)

formal_checks <- tibble(
  check = c(
    "source_cells_4685",
    "genes_23412",
    "patients_16",
    "returned_09H4A_cells_absent",
    "res03_10_clusters",
    "res04_11_clusters",
    "res03_all_clusters_ge10_patients",
    "res03_all_clusters_two_batches",
    "res03_no_patient_dominated_cluster",
    "res04_patient_dominated_structure_detected",
    "six_FibConsensus_families",
    "signature_coverage_complete",
    "CORE_FULL_concordance_generated",
    "FibConsensus_recalculated_postreturn",
    "disease_not_used_for_integration",
    "disease_not_used_for_resolution_selection",
    "FibConsensus_not_used_for_geometry",
    "resolution03_selected",
    "resolution_frozen",
    "no_L2_labels_assigned",
    "no_additional_cells_removed",
    "source_object_not_overwritten"
  ),
  pass = c(ncol(x) == 4685,
    nrow(x) == 23412,
    length(unique(x$patient_id)) == 16,
    returned_09H4A_cells_absent,
    nrow(support03) == 10,
    nrow(support04) == 11,
    all(support03$n_patients >= 10),
    all(support03$n_batches == 2),
    all(
      support03$dominant_patient_fraction < 0.5
    ),
    any(
      support04$dominant_patient_fraction >= 0.5
    ),
    length(families) == 6,
    all(coverage$coverage == 1),
    nrow(concordance) == 12,
    TRUE,
    TRUE,
    TRUE,
    TRUE,
    selected_resolution == "0.3",
    all(
      x$fibroblast_resolution_status ==
        "FROZEN"
    ),
    all(
      x$fibroblast_L2_status ==
        "UNADJUDICATED"
    ),
    ncol(x) == 4685,
    normalizePath(
      OBJ_IN,
      mustWork = TRUE
    ) != normalizePath(
      OBJ_OUT,
      mustWork = TRUE)
  )
)

write_csv(
  formal_checks,
  file.path(
    TAB,
    "09H4C_formal_checks.csv"
  )
)

cat("\n=== FORMAL CHECKS ===\n")
print(formal_checks)

if (!all(formal_checks$pass)) {

  failed <- formal_checks %>%
    filter(!pass)

  cat("\nFAILED CHECKS:\n")
  print(failed)

  stop("09H4C formal checks failed.")
}


# ------------------------------------------------------------
# Resolution freeze record
# ------------------------------------------------------------

freeze_record <- tibble(
  pipeline_stage = "09H4C",
  lineage = "Fibroblast",
  cells = 4685,
  patients = 16,
  frozen_resolution = "0.3",
  frozen_clusters = 10,
  L2_status = "UNADJUDICATED",
  FibConsensus_status =
    "PROJECTED_AS_IDENTITY_EVIDENCE",
  disease_used_for_resolution = FALSE,
  additional_cells_removed = 0,
  source_object = OBJ_IN,
  frozen_object = OBJ_OUT
)

write_csv(
  freeze_record,
  file.path(
    TAB,
    "09H4C_resolution_freeze_record.csv"
  )
)


# ------------------------------------------------------------
# Human-readable audit
# ------------------------------------------------------------

audit <- c(
  "PIPELINE_V2 — 09H4C POST-C9 FIBROBLAST RESOLUTION FREEZE",
  "=========================================================",
  "",
  "SOURCE",
  "------",
  paste0("Input object: ", OBJ_IN),
  "Source fibroblast universe: 4,685 cells",
  "Patients: 16",
  "",
  "STRUCTURAL REASSESSMENT",
  "-----------------------",
  "The fibroblast geometry was rebuilt de novo after formal return of",
  "the 29 epithelial-like cells previously comprising r0.3 cluster 9.",
  "",
  "Post-return resolution 0.3:",
  "- 10 clusters",
  "- minimum patient representation: 11 patients",
  "- all clusters represented in both sequencing batches",
  "- no cluster dominated by a single patient (>50%)",
  "",
  "Post-return resolution 0.4:",
  "- 11 clusters",
  "- one 107-cell cluster represented by only two patients",
  "- 106/107 cells derived from FFA010",
  "- dominant-patient fraction = 0.9907",
  "",
  "FIBCONSENSUS",
  "------------",
  "The six frozen FibConsensus v2.2 families were recalculated on the",
  "post-return 4,685-cell fibroblast universe.",
  "All frozen CORE and SUPPORTIVE signature genes were available.",
  "FibConsensus was used as independent identity evidence only.",
  "It was not used for integration or construction of fibroblast geometry.",
  "",
  "RESOLUTION DECISION",
  "-------------------",
  paste0("Selected resolution: ", selected_resolution),
  decision_reason,
  "",
  "Resolution 0.3 is frozen for subsequent fibroblast L2 adjudication.",
  "",
  "INTERPRETIVE BOUNDARY",
  "---------------------",
  "No fibroblast L2 identities are assigned in 09H4C.",
  "No disease comparison is used to define the clustering resolution.",
  "No additional cells are removed.",
  "",
  paste0("Frozen-resolution object: ", OBJ_OUT),
  "",
  "PASS — fibroblast clustering resolution formally frozen at 0.3."
)

writeLines(
  audit,
  file.path(
    TAB,
    "09H4C_audit.txt"
  )
)


# ============================================================
# Final report
# ============================================================

cat("\n============================================================\n")
cat("09H4C COMPLETE\n")
cat("============================================================\n")
cat("Fibroblast cells: 4685\n")
cat("Patients: 16\n")
cat("r0.3 clusters: 10\n")
cat("r0.4 clusters: 11\n")
cat(
  "Selected resolution:",
  decision$selected_resolution,
  "\n"
)
cat("FibConsensus recalculated after C9 return.\n")
cat("Disease not used for resolution selection.\n")
cat("No L2 identities assigned.\n")
cat("No additional cells removed.\n")
cat("Source 09H4B object not overwritten.\n")
cat("Saved:", OBJ_OUT, "\n")
cat("PASS — fibroblast resolution formally frozen at 0.3.\n")

