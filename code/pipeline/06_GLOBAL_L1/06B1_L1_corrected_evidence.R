suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(ggplot2)
  library(patchwork)
})

options(stringsAsFactors = FALSE)
set.seed(20260924)

source("PIPELINE_V2/config/paths.R")

cat("\n============================================================\n")
cat("PIPELINE_V2 — 06B1 L1 CORRECTED EVIDENCE\n")
cat("============================================================\n\n")

OBJ_FILE <- file.path(
  V2_ROOT,
  "objects",
  "ambient_corrected",
  "global_atlas_v2_SoupX_expression.rds"
)

OUT <- file.path(
  V2_ROOT,
  "results",
  "annotation",
  "06_GLOBAL_L1",
  "06B1_corrected"
)

TD <- file.path(OUT, "tables")
FD <- file.path(OUT, "figures")

dir.create(TD, recursive = TRUE, showWarnings = FALSE)
dir.create(FD, recursive = TRUE, showWarnings = FALSE)

# -------------------------------------------------------------------------
# LOAD CORRECTED EXPRESSION OBJECT
# -------------------------------------------------------------------------

x <- readRDS(OBJ_FILE)

DefaultAssay(x) <- "RNA"

Idents(x) <- factor(
  as.character(x$RNA_snn_res.0.6),
  levels = as.character(0:26)
)

stopifnot(
  ncol(x) == 82893,
  nrow(x) == 23412,
  length(levels(Idents(x))) == 27
)

cat("Cells:", ncol(x), "\n")
cat("Genes:", nrow(x), "\n")
cat("Clusters:", length(levels(Idents(x))), "\n\n")

# -------------------------------------------------------------------------
# BROAD L1 EVIDENCE PROGRAMMES
#
# These are evidence programmes.
# They are NOT final/frozen annotations.
# -------------------------------------------------------------------------

programs <- list(

  Keratinocyte = c(
    "KRT5",
    "KRT14",
    "KRT15",
    "KRT17",
    "KRT1",
    "KRT10",
    "KRTDAP",
    "DSG1"
  ),

  Fibroblast = c(
    "COL1A1",
    "COL1A2",
    "DCN",
    "LUM",
    "COL3A1",
    "PDGFRA"
  ),

  T_NK = c(
    "CD3D",
    "CD3E",
    "TRAC",
    "NKG7",
    "GNLY",
    "GZMK",
    "GZMH"
  ),

  Myeloid_DC = c(
    "LYZ",
    "FCER1G",
    "TYROBP",
    "CD1C",
    "CLEC10A",
    "C1QC",
    "LILRB2"
  ),

  Melanocyte = c(
    "MLANA",
    "PMEL",
    "TYR",
    "TYRP1",
    "DCT"
  ),

  Vascular_EC = c(
    "PECAM1",
    "VWF",
    "EMCN",
    "KDR",
    "PLVAP",
    "ACKR1"
  ),

  Lymphatic_EC = c(
    "PROX1",
    "PDPN",
    "LYVE1",
    "CCL21",
    "FLT4"
  ),

  Pericyte = c(
    "RGS5",
    "CSPG4",
    "PDGFRB",
    "NOTCH3",
    "ABCC9",
    "KCNJ8"
  ),

  Smooth_muscle = c(
    "ACTG2",
    "TAGLN",
    "CNN1",
    "MYH11",
    "ACTA2",
    "MYL9"
  ),

  Sebaceous = c(
    "KRT7",
    "PPARG",
    "PLIN2",
    "MGST1",
    "AWAT2",
    "MOGAT1"
  ),

  Sweat_gland = c(
    "DCD",
    "AQP5",
    "KRT19",
    "KRT8",
    "KRT18",
    "FOXA1",
    "MUC7"
  ),

  Neural_Schwann = c(
    "S100B",
    "SOX10",
    "S100A1",
    "PLP1",
    "SLC1A3"
  ),

  Cycling = c(
    "MKI67",
    "TOP2A",
    "UBE2C",
    "AURKB",
    "RRM2",
    "ASPM"
  )
)

requested_n <- vapply(
  programs,
  length,
  integer(1)
)

programs <- lapply(
  programs,
  function(g) intersect(g, rownames(x))
)

coverage <- data.frame(
  programme = names(programs),
  n_present = vapply(programs, length, integer(1)),
  n_requested = requested_n,
  stringsAsFactors = FALSE
)

coverage$coverage <- (
  coverage$n_present /
  coverage$n_requested
)

write.csv(
  coverage,
  file.path(
    TD,
    "06B1_marker_programme_coverage.csv"
  ),
  row.names = FALSE
)

cat("Marker coverage:\n")
print(
  coverage,
  row.names = FALSE
)

if (any(coverage$n_present == 0)) {
  stop(
    "At least one L1 programme has zero genes present."
  )
}

# -------------------------------------------------------------------------
# EXPRESSION MATRICES
# -------------------------------------------------------------------------

counts <- LayerData(
  x,
  assay = "RNA",
  layer = "counts"
)

dat <- LayerData(
  x,
  assay = "RNA",
  layer = "data"
)

stopifnot(
  nrow(counts) == 23412,
  nrow(dat) == 23412,
  identical(
    rownames(counts),
    rownames(dat)
  )
)

clusters <- as.character(
  Idents(x)
)

# -------------------------------------------------------------------------
# CLUSTER-LEVEL PROGRAMME EVIDENCE
# -------------------------------------------------------------------------

cat(
  "\nCalculating cluster-level programme evidence...\n"
)

ev <- list()

for (k in 0:26) {

  ii <- which(
    clusters == as.character(k)
  )

  for (nm in names(programs)) {

    g <- programs[[nm]]

    cell_programme <- Matrix::colMeans(
      dat[
        g,
        ii,
        drop = FALSE
      ]
    )

    gene_detection <- 100 *
      Matrix::rowMeans(
        counts[
          g,
          ii,
          drop = FALSE
        ] > 0
      )

    ev[[length(ev) + 1]] <- data.frame(

      cluster = k,

      programme = nm,

      n_cells = length(ii),

      mean_programme = mean(
        cell_programme
      ),

      median_programme = median(
        cell_programme
      ),

      mean_gene_detection_pct = mean(
        gene_detection
      ),

      stringsAsFactors = FALSE
    )
  }
}

ev <- do.call(
  rbind,
  ev
)

write.csv(
  ev,
  file.path(
    TD,
    "06B1_cluster_programme_evidence.csv"
  ),
  row.names = FALSE
)

# -------------------------------------------------------------------------
# PATIENT × CLUSTER SUPPORT
# -------------------------------------------------------------------------

cat(
  "Calculating patient-level support...\n"
)

pc <- list()

for (k in 0:26) {

  cc <- colnames(x)[
    clusters == as.character(k)
  ]

  pts <- unique(
    as.character(
      x$patient_id[cc]
    )
  )

  for (pt in pts) {

    jj <- cc[
      x$patient_id[cc] == pt
    ]

    for (nm in names(programs)) {

      g <- programs[[nm]]

      pc[[length(pc) + 1]] <- data.frame(

        cluster = k,

        patient_id = pt,

        disease = unique(
          as.character(
            x$disease[jj]
          )
        ),

        batch = unique(
          as.character(
            x$batch[jj]
          )
        ),

        programme = nm,

        n_cells = length(jj),

        mean_programme = mean(
          Matrix::colMeans(
            dat[
              g,
              jj,
              drop = FALSE
            ]
          )
        ),

        stringsAsFactors = FALSE
      )
    }
  }
}

pc <- do.call(
  rbind,
  pc
)

write.csv(
  pc,
  file.path(
    TD,
    "06B1_patient_cluster_programmes.csv"
  ),
  row.names = FALSE
)

# -------------------------------------------------------------------------
# PROGRAMME RANKING
#
# Standardize each programme across clusters.
# This is descriptive evidence only.
# -------------------------------------------------------------------------

wide <- reshape(
  ev[
    ,
    c(
      "cluster",
      "programme",
      "mean_programme"
    )
  ],
  idvar = "cluster",
  timevar = "programme",
  direction = "wide"
)

rownames(wide) <- wide$cluster

M <- as.matrix(
  wide[
    ,
    -1,
    drop = FALSE
  ]
)

colnames(M) <- sub(
  "^mean_programme\\.",
  "",
  colnames(M)
)

Z <- scale(M)

ranktab <- do.call(
  rbind,
  lapply(
    seq_len(nrow(Z)),
    function(i) {

      o <- order(
        Z[i, ],
        decreasing = TRUE
      )

      data.frame(

        cluster = as.integer(
          rownames(Z)[i]
        ),

        top_programme =
          colnames(Z)[o[1]],

        second_programme =
          colnames(Z)[o[2]],

        top_z =
          Z[i, o[1]],

        second_z =
          Z[i, o[2]],

        margin_z =
          Z[i, o[1]] -
          Z[i, o[2]],

        stringsAsFactors = FALSE
      )
    }
  )
)

ranktab <- ranktab[
  order(ranktab$cluster),
]

rownames(ranktab) <- NULL

write.csv(
  ranktab,
  file.path(
    TD,
    "06B1_cluster_programme_ranking.csv"
  ),
  row.names = FALSE
)

# -------------------------------------------------------------------------
# FULL POSITIVE MARKERS ON CORRECTED EXPRESSION
# -------------------------------------------------------------------------

cat(
  "\nFinding corrected-expression positive markers across 27 clusters...\n"
)

markers <- FindAllMarkers(
  x,
  assay = "RNA",
  slot = "data",
  only.pos = TRUE,
  min.pct = 0.10,
  logfc.threshold = 0.25,
  test.use = "wilcox",
  verbose = TRUE
)

write.csv(
  markers,
  file.path(
    TD,
    "06B1_all_positive_markers.csv"
  ),
  row.names = FALSE
)

top20 <- do.call(
  rbind,
  lapply(
    split(
      markers,
      markers$cluster
    ),
    function(z) {

      z <- z[
        order(
          z$p_val_adj,
          -z$avg_log2FC
        ),
      ]

      head(
        z,
        20
      )
    }
  )
)

rownames(top20) <- NULL

write.csv(
  top20,
  file.path(
    TD,
    "06B1_top20_markers.csv"
  ),
  row.names = FALSE
)

# -------------------------------------------------------------------------
# HEATMAP OF BROAD PROGRAMMES
# -------------------------------------------------------------------------

hm <- as.data.frame(
  as.table(Z)
)

colnames(hm) <- c(
  "cluster",
  "programme",
  "z"
)

hm$cluster <- factor(
  hm$cluster,
  levels = rev(
    as.character(
      0:26
    )
  )
)

p1 <- ggplot(
  hm,
  aes(
    programme,
    cluster,
    fill = z
  )
) +
  geom_tile() +
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0
  ) +
  theme_classic(
    base_size = 9
  ) +
  theme(
    axis.text.x =
      element_text(
        angle = 45,
        hjust = 1
      )
  ) +
  labs(
    x = NULL,
    y = "Global cluster",
    fill = "z-score",
    title =
      "L1 programme evidence - SoupX-corrected expression"
  )

ggsave(
  file.path(
    FD,
    "06B1_L1_programme_heatmap.pdf"
  ),
  p1,
  width = 9,
  height = 8,
  bg = "white"
)

ggsave(
  file.path(
    FD,
    "06B1_L1_programme_heatmap.png"
  ),
  p1,
  width = 9,
  height = 8,
  dpi = 600,
  bg = "white"
)

# -------------------------------------------------------------------------
# STACKED VIOLIN
# -------------------------------------------------------------------------

plot_genes <- unique(
  unlist(
    programs
  )
)

plot_genes <- plot_genes[
  plot_genes %in% rownames(x)
]

cat(
  "\nGenerating corrected stacked violin...\n"
)

p2 <- VlnPlot(
  x,
  features = plot_genes,
  group.by = "RNA_snn_res.0.6",
  stack = TRUE,
  flip = TRUE,
  pt.size = 0,
  combine = TRUE
)

ggsave(
  file.path(
    FD,
    "06B1_L1_stacked_violin_corrected.pdf"
  ),
  p2,
  width = 13,
  height = max(
    12,
    length(plot_genes) * 0.22
  ),
  bg = "white",
  limitsize = FALSE
)

# -------------------------------------------------------------------------
# DOTPLOT
#
# Quantitative companion only.
# It is not used as the sole annotation criterion.
# -------------------------------------------------------------------------

cat(
  "Generating corrected DotPlot...\n"
)

p3 <- DotPlot(
  x,
  features = programs,
  group.by = "RNA_snn_res.0.6",
  assay = "RNA"
) +
  RotatedAxis() +
  theme_classic(
    base_size = 9
  ) +
  labs(
    title =
      "L1 marker evidence - SoupX-corrected"
  )

ggsave(
  file.path(
    FD,
    "06B1_L1_dotplot_corrected.pdf"
  ),
  p3,
  width = 14,
  height = 8,
  bg = "white"
)

ggsave(
  file.path(
    FD,
    "06B1_L1_dotplot_corrected.png"
  ),
  p3,
  width = 14,
  height = 8,
  dpi = 600,
  bg = "white"
)

# -------------------------------------------------------------------------
# KNOWN FLAG FROM 05B3-05B6
# -------------------------------------------------------------------------

flags <- data.frame(

  cluster = 25,

  flag =
    "structured_KRT_high_contractile_subpopulation",

  evidence =
    paste0(
      "SoupX-resistant KRT-high subset; ",
      "enriched scDblFinder scores but ",
      "insufficient evidence for exclusion; ",
      "retain all cells pending biological ",
      "L1/L2 adjudication"
    ),

  stringsAsFactors = FALSE
)

write.csv(
  flags,
  file.path(
    TD,
    "06B1_annotation_flags.csv"
  ),
  row.names = FALSE
)

# -------------------------------------------------------------------------
# AUDIT
# -------------------------------------------------------------------------

writeLines(
  c(
    "PIPELINE_V2 — 06B1 corrected L1 evidence",
    paste0(
      "Generated: ",
      format(
        Sys.time(),
        "%Y-%m-%d %H:%M:%S"
      )
    ),
    paste0(
      "Cells: ",
      ncol(x)
    ),
    paste0(
      "Genes: ",
      nrow(x)
    ),
    "Expression source: SoupX-corrected counts/data.",
    "Clusters/geometry: canonical frozen V2 atlas.",
    "Broad programmes are annotation evidence, not frozen labels.",
    "Cell-level marker testing is exploratory annotation evidence.",
    "Patient-cluster programme table is provided to assess cross-individual support.",
    "Cluster 25 is retained and explicitly flagged.",
    "No cells removed; no identities frozen; no geometry recalculated."
  ),
  file.path(
    TD,
    "06B1_audit.txt"
  )
)

# -------------------------------------------------------------------------
# FINAL OUTPUT
# -------------------------------------------------------------------------

cat(
  "\n=== PROGRAMME RANKING ===\n"
)

print(
  ranktab,
  row.names = FALSE
)

cat(
  "\n============================================================\n"
)

cat(
  "06B1 COMPLETE\n"
)

cat(
  "============================================================\n"
)

cat(
  "Markers:",
  nrow(markers),
  "\n"
)

cat(
  "Tables:",
  TD,
  "\n"
)

cat(
  "Figures:",
  FD,
  "\n"
)

cat(
  "\nPASS — evidence only; no L1 labels frozen.\n"
)
