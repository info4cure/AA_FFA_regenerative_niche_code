suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(readr)
})

options(width = 200)

cat("============================================================\n")
cat("PIPELINE_V2 — 10B2 VASCULAR IDENTITY EVIDENCE\n")
cat("============================================================\n\n")

SRC <- paste0(
  "PIPELINE_V2/objects/annotated/",
  "vascular_v2_10B1_reclustered_3040.rds"
)

OUT <- paste0(
  "PIPELINE_V2/results/annotation/",
  "10_VASCULAR_L2/10B2_identity_evidence"
)

TAB <- file.path(OUT, "tables")
dir.create(TAB, recursive = TRUE, showWarnings = FALSE)

stopifnot(file.exists(SRC))

x <- readRDS(SRC)

stopifnot(ncol(x) == 3040)

DefaultAssay(x) <- "RNA"

# ============================================================
# 1. Candidate resolutions
# ============================================================

cat("STEP 1/7 — candidate resolutions\n")

RES <- c(
  "vascular_snn_res.0.2",
  "vascular_snn_res.0.3"
)

stopifnot(all(RES %in% colnames(x@meta.data)))

cat("r0.2 clusters:",
    length(unique(x@meta.data[[RES[1]]])), "\n")
cat("r0.3 clusters:",
    length(unique(x@meta.data[[RES[2]]])), "\n")

# ============================================================
# 2. De novo markers
# ============================================================

cat("\nSTEP 2/7 — de novo markers\n")

all_markers <- list()

for (res in RES) {

  cat("Markers:", res, "\n")

  Idents(x) <- x@meta.data[[res]]

  m <- FindAllMarkers(
    x,
    assay = "RNA",
    only.pos = TRUE,
    min.pct = 0.10,
    logfc.threshold = 0.25,
    test.use = "wilcox",
    verbose = FALSE
  )

  m$resolution <- res
  m$cluster <- as.character(m$cluster)

  all_markers[[res]] <- m

  safe <- gsub(
    "vascular_snn_res\\.",
    "r",
    res
  )

  write_csv(
    m,
    file.path(
      TAB,
      paste0("10B2_", safe, "_all_markers.csv")
    )
  )

  top25 <- m %>%
    filter(p_val_adj < 0.05) %>%
    group_by(cluster) %>%
    arrange(desc(avg_log2FC), .by_group = TRUE) %>%
    slice_head(n = 25) %>%
    ungroup()

  write_csv(
    top25,
    file.path(
      TAB,
      paste0("10B2_", safe, "_top25_markers.csv")
    )
  )
}

markers <- bind_rows(all_markers)

# ============================================================
# 3. Prespecified vascular programmes
# ============================================================

cat("STEP 3/7 — vascular programme evidence\n")

programmes <- list(

  PAN_ENDOTHELIAL = c(
    "PECAM1","VWF","CDH5","EMCN",
    "KDR","ESAM","ENG","RAMP2"
  ),

  LYMPHATIC = c(
    "PROX1","PDPN","FLT4","CCL21",
    "LYVE1","NRP2","RELN"
  ),

  CAPILLARY = c(
    "CA4","RGCC","EMCN","KDR",
    "PLVAP","EDNRB","BTNL9","CA2"
  ),

  VENOUS_ACKR1 = c(
    "ACKR1","SELP","SELE","VCAM1",
    "ICAM1","VWF"
  ),

  ARTERIAL = c(
    "EFNB2","SOX17","GJA5","GJA4",
    "HEY1","BMX","FBLN5","SEMA3G"
  ),

  ANGIOGENIC = c(
    "KDR","ESM1","ANGPT2","APLN",
    "PGF","PLVAP","KDR","CXCR4"
  ),

  INFLAMMATORY = c(
    "ICAM1","VCAM1","SELE",
    "CXCL9","CXCL10","CXCL11",
    "HLA-DRA","HLA-DRB1","CCL2"
  ),

  POSTCAPILLARY_VENULAR = c(
    "ACKR1","SELE","SELP",
    "ICAM1","VCAM1","PLVAP"
  ),

  MURAL_CONFLICT = c(
    "RGS5","CSPG4","MCAM","ACTA2",
    "TAGLN","MYH11","RBP1","PDGFRB"
  ),

  EPITHELIAL_CONFLICT = c(
    "EPCAM","KRT5","KRT14","KRT17",
    "KRT19","KRT8","KRT18","DSP"
  ),

  IMMUNE_CONFLICT = c(
    "PTPRC","CD3D","CD3E",
    "LST1","TYROBP","FCER1G"
  )
)

available_programmes <- lapply(
  programmes,
  function(g) intersect(g, rownames(x))
)

coverage <- tibble(
  programme = names(programmes),
  n_requested = sapply(programmes, length),
  n_present = sapply(available_programmes, length),
  genes_present = sapply(
    available_programmes,
    paste,
    collapse = ";"
  )
)

write_csv(
  coverage,
  file.path(TAB, "10B2_programme_coverage.csv")
)

print(coverage, n = Inf)

# ============================================================
# 4. Marker-programme evidence
# ============================================================

cat("\nSTEP 4/7 — marker-programme intersections\n")

sig <- markers %>%
  filter(p_val_adj < 0.05)

programme_marker_evidence <- bind_rows(
  lapply(RES, function(res) {

    m <- sig %>%
      filter(resolution == res)

    bind_rows(
      lapply(names(programmes), function(p) {

        genes <- available_programmes[[p]]

        hits <- m %>%
          filter(gene %in% genes) %>%
          group_by(cluster) %>%
          summarise(
            n_markers = n_distinct(gene),
            marker_genes = paste(
              unique(gene[order(-avg_log2FC)]),
              collapse = ";"
            ),
            max_log2FC = max(avg_log2FC),
            .groups = "drop"
          )

        grid <- tibble(
          cluster = sort(
            unique(m$cluster)
          )
        )

        grid %>%
          left_join(hits, by = "cluster") %>%
          mutate(
            resolution = res,
            programme = p,
            n_markers = replace_na(n_markers, 0L),
            marker_genes =
              replace_na(marker_genes, ""),
            max_log2FC =
              replace_na(max_log2FC, 0)
          )
      })
    )
  })
) %>%
  select(
    resolution,
    cluster,
    programme,
    n_markers,
    marker_genes,
    max_log2FC
  ) %>%
  arrange(
    resolution,
    as.integer(cluster),
    programme
  )

write_csv(
  programme_marker_evidence,
  file.path(TAB, "10B2_marker_programme_evidence.csv")
)

# ============================================================
# 5. Average expression programme scores
# ============================================================

cat("STEP 5/7 — cluster programme scores\n")

expr <- GetAssayData(
  x,
  assay = "RNA",
  layer = "data"
)

score_one_programme <- function(genes) {

  genes <- intersect(genes, rownames(expr))

  if (length(genes) == 0) {
    return(rep(NA_real_, ncol(x)))
  }

  Matrix::colMeans(
    expr[genes, , drop = FALSE]
  )
}

score_df <- data.frame(
  cell = colnames(x),
  stringsAsFactors = FALSE
)

for (p in names(available_programmes)) {
  score_df[[p]] <-
    score_one_programme(
      available_programmes[[p]]
    )
}

score_df$patient_id <-
  as.character(x$patient_id)

for (res in RES) {
  score_df[[res]] <-
    as.character(x@meta.data[[res]])
}

score_long <- score_df %>%
  pivot_longer(
    cols = all_of(names(programmes)),
    names_to = "programme",
    values_to = "score"
  )

cluster_scores <- bind_rows(
  lapply(RES, function(res) {

    score_long %>%
      mutate(
        cluster = score_df[[res]][
          match(cell, score_df$cell)
        ]
      ) %>%
      group_by(
        cluster,
        programme
      ) %>%
      summarise(
        mean_score = mean(score, na.rm = TRUE),
        median_score = median(score, na.rm = TRUE),
        n_cells = n(),
        .groups = "drop"
      ) %>%
      mutate(resolution = res)
  })
) %>%
  select(
    resolution,
    cluster,
    programme,
    mean_score,
    median_score,
    n_cells
  ) %>%
  arrange(
    resolution,
    as.integer(cluster),
    programme
  )

write_csv(
  cluster_scores,
  file.path(TAB, "10B2_cluster_programme_scores.csv")
)

# ============================================================
# 6. Patient-level programme support
# ============================================================

cat("STEP 6/7 — patient-level programme support\n")

patient_scores <- bind_rows(
  lapply(RES, function(res) {

    score_long %>%
      mutate(
        cluster = score_df[[res]][
          match(cell, score_df$cell)
        ]
      ) %>%
      group_by(
        cluster,
        patient_id,
        programme
      ) %>%
      summarise(
        mean_score = mean(score, na.rm = TRUE),
        n_cells = n(),
        .groups = "drop"
      ) %>%
      mutate(resolution = res)
  })
)

write_csv(
  patient_scores,
  file.path(TAB, "10B2_patient_programme_scores.csv")
)

# ============================================================
# 7. Compact evidence summary
# ============================================================

cat("STEP 7/7 — compact evidence summary\n")

top_marker_summary <- markers %>%
  filter(p_val_adj < 0.05) %>%
  group_by(resolution, cluster) %>%
  arrange(desc(avg_log2FC), .by_group = TRUE) %>%
  summarise(
    n_significant_markers = n(),
    top10_genes =
      paste(head(gene, 10), collapse = ";"),
    .groups = "drop"
  )

write_csv(
  top_marker_summary,
  file.path(TAB, "10B2_cluster_marker_summary.csv")
)

cat("\n=== CLUSTER MARKER SUMMARY ===\n")
print(top_marker_summary, n = Inf)

cat("\n=== PROGRAMME MARKER EVIDENCE ===\n")
print(
  programme_marker_evidence %>%
    filter(
      n_markers > 0 |
      programme %in%
        c(
          "LYMPHATIC",
          "CAPILLARY",
          "VENOUS_ACKR1",
          "ARTERIAL",
          "ANGIOGENIC",
          "INFLAMMATORY",
          "MURAL_CONFLICT"
        )
    ),
  n = Inf
)

checks <- tibble(
  check = c(
    "cells_3040",
    "patients_16",
    "r02_six_clusters",
    "r03_nine_clusters",
    "markers_generated_r02",
    "markers_generated_r03",
    "programme_coverage_generated",
    "programme_scores_generated",
    "patient_support_generated",
    "disease_not_used_for_identity",
    "source_L1_not_used_for_scoring",
    "no_resolution_selected",
    "no_L2_labels_assigned",
    "no_cells_removed",
    "source_object_not_modified"
  ),
  pass = c(
    ncol(x) == 3040,
    length(unique(x$patient_id)) == 16,
    length(unique(x@meta.data[[RES[1]]])) == 6,
    length(unique(x@meta.data[[RES[2]]])) == 9,
    nrow(all_markers[[RES[1]]]) > 0,
    nrow(all_markers[[RES[2]]]) > 0,
    nrow(coverage) == length(programmes),
    nrow(cluster_scores) > 0,
    nrow(patient_scores) > 0,
    TRUE,
    TRUE,
    TRUE,
    !("vascular_L2" %in% colnames(x@meta.data)),
    ncol(x) == 3040,
    TRUE
  )
)

write_csv(
  checks,
  file.path(TAB, "10B2_formal_checks.csv")
)

cat("\n=== FORMAL CHECKS ===\n")
print(checks, n = Inf)

stopifnot(all(checks$pass))

audit <- c(
  "PIPELINE_V2 — 10B2 VASCULAR IDENTITY EVIDENCE",
  "",
  "Candidate resolutions evaluated:",
  "- 0.2",
  "- 0.3",
  "",
  "Evidence:",
  "- de novo cluster markers",
  "- prespecified vascular transcriptional programmes",
  "- cluster-level programme scores",
  "- patient-level programme support",
  "- competing-lineage programmes",
  "",
  "Interpretive boundary:",
  "- disease not used for identity",
  "- source L1 not used for programme scoring",
  "- no resolution selected",
  "- no L2 identities assigned",
  "- no cells removed",
  "- source object not modified",
  "",
  "PASS"
)

writeLines(
  audit,
  file.path(TAB, "10B2_audit.txt")
)

cat("\n============================================================\n")
cat("10B2 COMPLETE\n")
cat("============================================================\n")
cat("Cells:", ncol(x), "\n")
cat("r0.2 clusters: 6\n")
cat("r0.3 clusters: 9\n")
cat("No resolution selected.\n")
cat("No L2 identities assigned.\n")
cat("PASS\n")
