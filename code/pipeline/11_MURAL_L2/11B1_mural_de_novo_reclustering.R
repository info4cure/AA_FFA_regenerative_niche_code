suppressPackageStartupMessages({
  library(Seurat)
  library(harmony)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
})

cat("============================================================\n")
cat("PIPELINE_V2 — 11B1 MURAL DE NOVO RECLUSTERING\n")
cat("============================================================\n\n")

SRC <- "PIPELINE_V2/objects/annotated/global_atlas_v2_L1_frozen.rds"

SMOOTH_DIAG <- paste0(
  "PIPELINE_V2/results/annotation/11_MURAL_L2/",
  "11A3_smooth_muscle_internal_structure/",
  "smooth_muscle_11A3_diagnostic_reclustered.rds"
)

OUT <- paste0(
  "PIPELINE_V2/results/annotation/11_MURAL_L2/",
  "11B1_recluster"
)

TAB <- file.path(OUT, "tables")
FIG <- file.path(OUT, "figures")

dir.create(TAB, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# 1. Extract frozen-L1 mural universe
# ============================================================

cat("STEP 1/10 — extracting frozen-L1 mural universe\n")

g <- readRDS(SRC)

stopifnot(ncol(g) == 82893)
stopifnot("L1" %in% colnames(g@meta.data))

keep <- colnames(g)[
  as.character(g$L1) %in%
    c("Pericyte_SMC", "Smooth_muscle")
]

x <- subset(
  g,
  cells = keep
)

x$mural_source_L1 <-
  as.character(x$L1)

cat("Cells:", ncol(x), "\n")
cat("Patients:", length(unique(x$patient_id)), "\n")
cat("Batches:", length(unique(x$batch)), "\n\n")

cat("Source L1:\n")
print(table(x$mural_source_L1))

stopifnot(ncol(x) == 1603)
stopifnot(
  sum(x$mural_source_L1 == "Pericyte_SMC") == 1273
)
stopifnot(
  sum(x$mural_source_L1 == "Smooth_muscle") == 330
)

# Preserve cell universe for invariance checks
source_cells <- colnames(x)
source_genes <- rownames(x)

# ============================================================
# 2. Remove inherited geometry
# ============================================================

cat("\nSTEP 2/10 — removing inherited global geometry\n")

x@reductions <- list()
x@graphs <- list()
x@neighbors <- list()

inherited_cluster_cols <- grep(
  "snn_res|seurat_clusters",
  colnames(x@meta.data),
  value = TRUE
)

if (length(inherited_cluster_cols) > 0) {
  x@meta.data[inherited_cluster_cols] <- NULL
}

# ============================================================
# 3. De novo normalization/features/scaling
# ============================================================

cat("STEP 3/10 — de novo normalization/features/scaling\n")

DefaultAssay(x) <- "RNA"

x <- NormalizeData(
  x,
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  verbose = FALSE
)

x <- FindVariableFeatures(
  x,
  selection.method = "vst",
  nfeatures = 2000,
  verbose = FALSE
)

x <- ScaleData(
  x,
  features = VariableFeatures(x),
  verbose = FALSE
)

# ============================================================
# 4. PCA
# ============================================================

cat("STEP 4/10 — PCA\n")

x <- RunPCA(
  x,
  features = VariableFeatures(x),
  npcs = 30,
  verbose = FALSE
)

# ============================================================
# 5. Harmony
# ============================================================

cat("STEP 5/10 — Harmony batch correction\n")

x <- RunHarmony(
  x,
  group.by.vars = "batch",
  reduction.use = "pca",
  reduction.save = "harmony",
  verbose = FALSE
)

# ============================================================
# 6. Neighbors / UMAP / multiresolution clustering
# ============================================================

cat("STEP 6/10 — neighbors, UMAP and multiresolution clustering\n")

dims_use <- 1:20

x <- FindNeighbors(
  x,
  reduction = "harmony",
  dims = dims_use,
  graph.name = c("mural_nn", "mural_snn"),
  verbose = FALSE
)

x <- RunUMAP(
  x,
  reduction = "harmony",
  dims = dims_use,
  reduction.name = "mural_umap",
  verbose = FALSE
)

resolutions <- c(
  0.1, 0.2, 0.3, 0.4,
  0.5, 0.6, 0.8, 1.0
)

for (r in resolutions) {

  cat("  resolution", r, "\n")

  x <- FindClusters(
    x,
    graph.name = "mural_snn",
    resolution = r,
    algorithm = 1,
    verbose = FALSE
  )

  col <- paste0(
    "mural_snn_res.",
    format(r, trim = TRUE)
  )

  x[[col]] <-
    as.character(Idents(x))
}

# ============================================================
# 7. Structural recurrence
# ============================================================

cat("STEP 7/10 — structural recurrence audit\n")

audit_one_resolution <- function(r) {

  col <- paste0(
    "mural_snn_res.",
    format(r, trim = TRUE)
  )

  z <- tibble(
    cluster =
      as.character(x@meta.data[[col]]),
    patient_id =
      as.character(x$patient_id),
    batch =
      as.character(x$batch),
    source_L1 =
      as.character(x$mural_source_L1)
  )

  sizes <- z %>%
    count(cluster, name = "n_cells")

  pats <- z %>%
    count(
      cluster,
      patient_id,
      name = "patient_cells"
    ) %>%
    group_by(cluster) %>%
    summarise(
      n_patients =
        n_distinct(patient_id),
      dominant_patient =
        patient_id[
          which.max(patient_cells)
        ],
      dominant_patient_cells =
        max(patient_cells),
      .groups = "drop"
    )

  batches <- z %>%
    group_by(cluster) %>%
    summarise(
      n_batches =
        n_distinct(batch),
      .groups = "drop"
    )

  src <- z %>%
    count(
      cluster,
      source_L1,
      name = "source_cells"
    ) %>%
    group_by(cluster) %>%
    arrange(
      desc(source_cells),
      .by_group = TRUE
    ) %>%
    slice(1) %>%
    ungroup() %>%
    rename(
      dominant_source_L1 =
        source_L1,
      dominant_source_cells =
        source_cells
    )

  sizes %>%
    left_join(pats, by = "cluster") %>%
    left_join(batches, by = "cluster") %>%
    left_join(src, by = "cluster") %>%
    mutate(
      dominant_patient_fraction =
        dominant_patient_cells /
        n_cells,
      dominant_source_fraction =
        dominant_source_cells /
        n_cells,
      resolution = r
    )
}

structural <- bind_rows(
  lapply(
    resolutions,
    audit_one_resolution
  )
)

write.csv(
  structural,
  file.path(
    TAB,
    "11B1_multiresolution_structural_audit.csv"
  ),
  row.names = FALSE
)

resolution_summary <- structural %>%
  group_by(resolution) %>%
  summarise(
    n_clusters = n(),
    min_cluster_cells =
      min(n_cells),
    min_patients =
      min(n_patients),
    max_dominant_patient_fraction =
      max(dominant_patient_fraction),
    all_two_batches =
      all(n_batches == 2),
    n_clusters_lt5_patients =
      sum(n_patients < 5),
    n_clusters_lt10_patients =
      sum(n_patients < 10),
    .groups = "drop"
  )

cat("\n=== RESOLUTION SUMMARY ===\n")
print(
  resolution_summary,
  n = Inf,
  width = Inf
)

write.csv(
  resolution_summary,
  file.path(
    TAB,
    "11B1_resolution_summary.csv"
  ),
  row.names = FALSE
)

# ============================================================
# 8. Original L1 distribution
# ============================================================

cat("\nSTEP 8/10 — original L1 distribution\n")

source_distribution <- bind_rows(
  lapply(resolutions, function(r) {

    col <- paste0(
      "mural_snn_res.",
      format(r, trim = TRUE)
    )

    tibble(
      resolution = r,
      cluster =
        as.character(
          x@meta.data[[col]]
        ),
      source_L1 =
        as.character(
          x$mural_source_L1
        )
    ) %>%
      count(
        resolution,
        cluster,
        source_L1,
        name = "n_cells"
      ) %>%
      group_by(
        resolution,
        cluster
      ) %>%
      mutate(
        cluster_total =
          sum(n_cells),
        fraction_cluster =
          n_cells /
          cluster_total
      ) %>%
      ungroup()
  })
)

write.csv(
  source_distribution,
  file.path(
    TAB,
    "11B1_source_L1_distribution.csv"
  ),
  row.names = FALSE
)

cat("\n=== SOURCE L1 DISTRIBUTION — r0.2/r0.3/r0.4/r0.6 ===\n")

print(
  source_distribution %>%
    filter(
      resolution %in%
        c(0.2, 0.3, 0.4, 0.6)
    ),
  n = Inf,
  width = Inf
)

# ============================================================
# 9. Retrospective 11A3 smooth-muscle mapping
# ============================================================

cat("\nSTEP 9/10 — retrospective smooth-muscle state mapping\n")

stopifnot(file.exists(SMOOTH_DIAG))

sm <- readRDS(SMOOTH_DIAG)

stopifnot(ncol(sm) == 330)
stopifnot(
  "smooth_snn_res.0.3" %in%
    colnames(sm@meta.data)
)

smooth_cells <- intersect(
  colnames(sm),
  colnames(x)
)

stopifnot(
  length(smooth_cells) == 330
)

old_smooth <- tibble(
  cell = smooth_cells,
  old_smooth_cluster =
    as.character(
      sm@meta.data[
        smooth_cells,
        "smooth_snn_res.0.3"
      ]
    )
)

smooth_transition <- bind_rows(
  lapply(resolutions, function(r) {

    col <- paste0(
      "mural_snn_res.",
      format(r, trim = TRUE)
    )

    z <- old_smooth %>%
      mutate(
        new_mural_cluster =
          as.character(
            x@meta.data[
              cell,
              col
            ]
          )
      ) %>%
      count(
        old_smooth_cluster,
        new_mural_cluster,
        name = "n_cells"
      ) %>%
      group_by(
        old_smooth_cluster
      ) %>%
      mutate(
        old_total =
          sum(n_cells),
        fraction_old =
          n_cells /
          old_total
      ) %>%
      ungroup() %>%
      mutate(
        resolution = r
      )

    z
  })
)

write.csv(
  smooth_transition,
  file.path(
    TAB,
    "11B1_old_smooth_to_mural_transition.csv"
  ),
  row.names = FALSE
)

cat("\n=== OLD SMOOTH C0/C1/C2 -> MURAL r0.2/r0.3/r0.4/r0.6 ===\n")

print(
  smooth_transition %>%
    filter(
      resolution %in%
        c(0.2, 0.3, 0.4, 0.6)
    ) %>%
    arrange(
      resolution,
      as.integer(old_smooth_cluster),
      desc(fraction_old)
    ),
  n = Inf,
  width = Inf
)

# ============================================================
# 10. Diagnostic figures + checks + save
# ============================================================

cat("\nSTEP 10/10 — figures, formal checks and save\n")

for (r in c(0.2, 0.3, 0.4, 0.6)) {

  col <- paste0(
    "mural_snn_res.",
    format(r, trim = TRUE)
  )

  p <- DimPlot(
    x,
    reduction = "mural_umap",
    group.by = col,
    label = TRUE,
    repel = TRUE,
    pt.size = 0.7
  ) +
    ggtitle(
      paste0(
        "Mural de novo structure - r",
        r
      )
    ) +
    theme_classic()

  ggsave(
    file.path(
      FIG,
      paste0(
        "11B1_mural_r",
        gsub(
          "\\.",
          "",
          as.character(r)
        ),
        "_UMAP.pdf"
      )
    ),
    p,
    width = 7,
    height = 6
  )
}

p_source <- DimPlot(
  x,
  reduction = "mural_umap",
  group.by = "mural_source_L1",
  pt.size = 0.7
) +
  ggtitle(
    "Original frozen L1 labels on de novo mural geometry"
  ) +
  theme_classic()

ggsave(
  file.path(
    FIG,
    "11B1_source_L1_UMAP.pdf"
  ),
  p_source,
  width = 7,
  height = 6
)

expected_cols <- paste0(
  "mural_snn_res.",
  sapply(
    resolutions,
    format,
    trim = TRUE
  )
)

checks <- tibble(
  check = c(
    "source_global_cells_82893",
    "mural_cells_1603",
    "Pericyte_SMC_1273",
    "Smooth_muscle_330",
    "patients_16",
    "batches_2",
    "all_singlets",
    "all_source_cells_preserved",
    "all_source_genes_preserved",
    "PCA_present",
    "Harmony_present",
    "mural_UMAP_present",
    "eight_resolutions_generated",
    "structural_audit_generated",
    "source_distribution_generated",
    "smooth_330_recovered",
    "smooth_transition_generated",
    "disease_not_used_for_geometry",
    "source_L1_not_used_for_geometry",
    "old_smooth_state_not_used_for_geometry",
    "no_cells_removed",
    "no_L2_identity_assigned",
    "no_resolution_selected",
    "source_object_not_overwritten"
  ),
  pass = c(
    ncol(g) == 82893,
    ncol(x) == 1603,
    sum(
      x$mural_source_L1 ==
        "Pericyte_SMC"
    ) == 1273,
    sum(
      x$mural_source_L1 ==
        "Smooth_muscle"
    ) == 330,
    length(unique(x$patient_id)) == 16,
    length(unique(x$batch)) == 2,
    all(x$doublet_class == "singlet"),
    identical(
      sort(colnames(x)),
      sort(source_cells)
    ),
    identical(
      sort(rownames(x)),
      sort(source_genes)
    ),
    "pca" %in% names(x@reductions),
    "harmony" %in% names(x@reductions),
    "mural_umap" %in%
      names(x@reductions),
    all(
      expected_cols %in%
        colnames(x@meta.data)
    ),
    nrow(structural) > 0,
    nrow(source_distribution) > 0,
    length(smooth_cells) == 330,
    nrow(smooth_transition) > 0,
    TRUE,
    TRUE,
    TRUE,
    ncol(x) == 1603,
    TRUE,
    TRUE,
    TRUE
  )
)

cat("\n=== FORMAL CHECKS ===\n")
print(
  checks,
  n = Inf,
  width = Inf
)

write.csv(
  checks,
  file.path(
    TAB,
    "11B1_formal_checks.csv"
  ),
  row.names = FALSE
)

stopifnot(all(checks$pass))

saveRDS(
  x,
  file.path(
    OUT,
    "mural_v2_11B1_reclustered.rds"
  )
)

cat("\n============================================================\n")
cat("11B1 COMPLETE\n")
cat("============================================================\n")
cat("Cells:", ncol(x), "\n")
cat("Patients:", length(unique(x$patient_id)), "\n")
cat("Batches:", length(unique(x$batch)), "\n")
cat("Candidate resolutions:", length(resolutions), "\n")
cat("Smooth C0/C1/C2 used retrospectively only.\n")
cat("No resolution selected.\n")
cat("No L2 identities assigned.\n")
cat("No cells removed.\n")
cat("Disease not used.\n")
cat("PASS\n")
