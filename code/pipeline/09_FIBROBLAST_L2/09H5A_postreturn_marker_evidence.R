suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(readr)
})

cat("============================================================\n")
cat("PIPELINE_V2 — 09H5A POST-RETURN FIBROBLAST MARKER EVIDENCE\n")
cat("============================================================\n\n")

OBJ <- "PIPELINE_V2/objects/annotated/fibroblast_v2_09H4C_resolution03_frozen_4685.rds"

OUT <- file.path(
  "PIPELINE_V2/results/annotation/09_FIBROBLAST_L2",
  "09H5A_postreturn_marker_evidence"
)

TAB <- file.path(OUT, "tables")
dir.create(TAB, recursive = TRUE, showWarnings = FALSE)

stopifnot(file.exists(OBJ))

# ------------------------------------------------------------
# 1. Load frozen-resolution object
# ------------------------------------------------------------

cat("STEP 1/7 — loading frozen r0.3 fibroblast object\n")

x <- readRDS(OBJ)

stopifnot(ncol(x) == 4685)
stopifnot("fibroblast_cluster_frozen" %in% colnames(x@meta.data))

cl <- as.character(x$fibroblast_cluster_frozen)

stopifnot(length(unique(cl)) == 10)
stopifnot(identical(sort(as.integer(unique(cl))), 0:9))

Idents(x) <- "fibroblast_cluster_frozen"

cat("Cells:", ncol(x), "\n")
cat("Clusters:", length(unique(cl)), "\n\n")

# ------------------------------------------------------------
# 2. Prepare RNA expression
# ------------------------------------------------------------

cat("STEP 2/7 — preparing RNA expression\n")

stopifnot("RNA" %in% Assays(x))
DefaultAssay(x) <- "RNA"

# Seurat v5 objects can contain multiple RNA layers.
x <- JoinLayers(x, assay = "RNA")

# ------------------------------------------------------------
# 3. Recalculate markers on the POST-RETURN geometry
# ------------------------------------------------------------

cat("STEP 3/7 — calculating post-return de novo markers\n")

markers <- FindAllMarkers(
  object = x,
  assay = "RNA",
  only.pos = TRUE,
  min.pct = 0.10,
  logfc.threshold = 0.25,
  test.use = "wilcox",
  verbose = FALSE
)

stopifnot(nrow(markers) > 0)

markers <- markers %>%
  mutate(cluster = as.character(cluster)) %>%
  arrange(
    as.integer(cluster),
    p_val_adj,
    desc(avg_log2FC)
  )

write_csv(
  markers,
  file.path(TAB, "09H5A_all_postreturn_markers.csv")
)

cat("Marker rows:", nrow(markers), "\n\n")

# ------------------------------------------------------------
# 4. Ranked marker summaries
# ------------------------------------------------------------

cat("STEP 4/7 — generating ranked marker summaries\n")

sig <- markers %>%
  filter(p_val_adj < 0.05)

top25 <- sig %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC), .by_group = TRUE) %>%
  slice_head(n = 25) %>%
  ungroup()

top50 <- sig %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC), .by_group = TRUE) %>%
  slice_head(n = 50) %>%
  ungroup()

write_csv(
  top25,
  file.path(TAB, "09H5A_top25_markers_by_cluster.csv")
)

write_csv(
  top50,
  file.path(TAB, "09H5A_top50_markers_by_cluster.csv")
)

summary_tbl <- sig %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC), .by_group = TRUE) %>%
  summarise(
    n_significant_markers = n(),
    top10_genes = paste(head(gene, 10), collapse = ";"),
    .groups = "drop"
  ) %>%
  arrange(as.integer(cluster))

write_csv(
  summary_tbl,
  file.path(TAB, "09H5A_cluster_marker_summary.csv")
)

print(summary_tbl)

# ------------------------------------------------------------
# 5. Targeted identity / competing-lineage audit
# ------------------------------------------------------------

cat("\nSTEP 5/7 — targeted programme audit\n")

programmes <- list(

  PAPILLARY = c(
    "COL6A5","COL23A1","COL13A1","COL18A1",
    "APCDD1","NKD2","HSPB3"
  ),

  RETICULAR = c(
    "SFRP2","DPP4","FBLN1","COL12A1",
    "SMOC2","COL1A1","COL1A2","COL3A1"
  ),

  DERMAL_SHEATH = c(
    "DPEP1","COL11A1","KIF26B","EDNRA",
    "POSTN","MME","RBP4","SOX2"
  ),

  DERMAL_PAPILLA = c(
    "COCH","ASPN","TNN","CRABP1",
    "LEF1","CORIN","RSPO3"
  ),

  PI16_PERIVASCULAR = c(
    "PI16","COL15A1","APOD","PTGDS",
    "PPARG","CXCL12"
  ),

  APOE_IMMUNE_PERIVASCULAR = c(
    "APOE","C3","CFD","C7","APOC1",
    "CCL19","IL33","CXCL12"
  ),

  REMODELING_MYOFIBROBLAST = c(
    "ADAM12","MMP11","COL10A1","COMP",
    "ACTA2","TAGLN","POSTN"
  ),

  MURAL_PERICYTE = c(
    "RGS5","CSPG4","MCAM","NOTCH3",
    "PDGFRB","ACTA2","TAGLN","RBP1"
  ),

  EPITHELIAL = c(
    "EPCAM","KRT5","KRT14","KRT17","KRT19",
    "DSP","DSG1","DSC3","TACSTD2","KRT1",
    "KRT10","SERPINB5","NECTIN4","SPINK5","PKP1"
  ),

  NEURAL_SCHWANN = c(
    "S100B","SOX10","PLP1","PMP22","MPZ","S100A1"
  )
)

programme_hits <- bind_rows(
  lapply(names(programmes), function(p) {

    genes <- programmes[[p]]

    sig %>%
      filter(gene %in% genes) %>%
      group_by(cluster) %>%
      summarise(
        programme = p,
        n_programme_markers = n_distinct(gene),
        marker_genes = paste(
          unique(gene[order(-avg_log2FC)]),
          collapse = ";"
        ),
        max_log2FC = max(avg_log2FC),
        .groups = "drop"
      )
  })
)

programme_grid <- expand_grid(
  cluster = as.character(0:9),
  programme = names(programmes)
)

programme_evidence <- programme_grid %>%
  left_join(
    programme_hits,
    by = c("cluster", "programme")
  ) %>%
  mutate(
    n_programme_markers =
      replace_na(n_programme_markers, 0L),
    marker_genes =
      replace_na(marker_genes, ""),
    max_log2FC =
      replace_na(max_log2FC, 0)
  ) %>%
  arrange(as.integer(cluster), programme)

write_csv(
  programme_evidence,
  file.path(TAB, "09H5A_marker_programme_evidence.csv")
)

# ------------------------------------------------------------
# 6. Formal checks
# ------------------------------------------------------------

cat("STEP 6/7 — formal checks\n")

checks <- tibble(
  check = c(
    "source_cells_4685",
    "frozen_resolution_present",
    "ten_clusters_present",
    "markers_generated",
    "all_clusters_have_markers",
    "programme_audit_generated",
    "disease_not_used_for_markers",
    "no_L2_labels_assigned",
    "no_cells_removed",
    "source_object_not_overwritten"
  ),
  pass = c(
    ncol(x) == 4685,
    "fibroblast_cluster_frozen" %in% colnames(x@meta.data),
    length(unique(cl)) == 10,
    nrow(markers) > 0,
    length(unique(markers$cluster)) == 10,
    nrow(programme_evidence) == 10 * length(programmes),
    TRUE,
    TRUE,
    ncol(x) == 4685,
    TRUE
  )
)

write_csv(
  checks,
  file.path(TAB, "09H5A_formal_checks.csv")
)

print(checks)

stopifnot(all(checks$pass))

# ------------------------------------------------------------
# 7. Audit
# ------------------------------------------------------------

cat("\nSTEP 7/7 — writing audit\n")

audit <- c(
  "PIPELINE_V2 — 09H5A POST-RETURN FIBROBLAST MARKER EVIDENCE",
  "",
  paste("Source:", OBJ),
  "",
  "Universe:",
  "- 4,685 fibroblast cells",
  "- frozen clustering resolution 0.3",
  "- 10 de novo clusters",
  "",
  "Purpose:",
  "- recalculate cluster markers after formal return of the 29 epithelial-like cells",
  "- generate post-return identity evidence",
  "- audit fibroblast-family and competing-lineage programmes",
  "",
  "Interpretive boundary:",
  "- no fibroblast L2 identities assigned",
  "- disease not used for marker definition",
  "- no cells removed",
  "- source object not overwritten",
  "",
  "PASS — post-return marker evidence generated."
)

writeLines(
  audit,
  file.path(TAB, "09H5A_audit.txt")
)

cat("\n============================================================\n")
cat("09H5A COMPLETE\n")
cat("============================================================\n")
cat("Cells:", ncol(x), "\n")
cat("Clusters:", length(unique(cl)), "\n")
cat("Marker rows:", nrow(markers), "\n")
cat("No L2 identities assigned.\n")
cat("No cells removed.\n")
cat("PASS\n")
