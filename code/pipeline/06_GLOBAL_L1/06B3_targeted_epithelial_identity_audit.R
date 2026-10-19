suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
})

cat("\n============================================================\n")
cat("PIPELINE_V2 — 06B3 TARGETED EPITHELIAL / FOLLICULAR AUDIT\n")
cat("============================================================\n\n")

OBJ <- "PIPELINE_V2/objects/ambient_corrected/global_atlas_v2_SoupX_expression.rds"
OUT <- "PIPELINE_V2/results/annotation/06_GLOBAL_L1/06B3_epithelial_audit"
TDIR <- file.path(OUT, "tables")
FDIR <- file.path(OUT, "figures")
CLUSTER_COL <- "RNA_snn_res.0.6"

dir.create(TDIR, recursive=TRUE, showWarnings=FALSE)
dir.create(FDIR, recursive=TRUE, showWarnings=FALSE)

# ============================================================
# 1. LOAD CANONICAL SOUPX-CORRECTED OBJECT
# ============================================================

x <- readRDS(OBJ)
DefaultAssay(x) <- "RNA"

stopifnot(
  ncol(x) == 82893,
  nrow(x) == 23412,
  CLUSTER_COL %in% colnames(x@meta.data)
)

x$global_cluster <- as.character(x@meta.data[[CLUSTER_COL]])

observed_clusters <- sort(unique(as.integer(x$global_cluster)))

stopifnot(
  identical(observed_clusters, 0:26)
)

cat("Cells:", ncol(x), "\n")
cat("Genes:", nrow(x), "\n")
cat("Clusters:", length(observed_clusters), "\n")
cat("Expression source:", paste(unique(x$expression_source), collapse=", "), "\n\n")

# ============================================================
# 2. TARGETED PROGRAMMES
#
# PURPOSE:
#   - establish epithelial vs non-epithelial L1 evidence
#   - retain anatomical follicular evidence for later L2
#
# IMPORTANT:
#   programme scores are evidence, NOT automatic labels.
#   Cycling is a state, NOT an L1 lineage.
# ============================================================

programmes <- list(

  # ---- Broad epithelial identity ----

  Epithelial_core = c(
    "EPCAM","KRT8","KRT18","KRT19",
    "KRT5","KRT14","KRT17","DSP"
  ),

  Keratinocyte_core = c(
    "KRT5","KRT14","KRT1","KRT10",
    "KRT2","KRT15","KRT17","KRTDAP"
  ),

  # ---- Epidermal differentiation ----

  Epidermal_basal = c(
    "KRT5","KRT14","KRT15","KRT17",
    "TP63","COL17A1","ITGA6","KRT19"
  ),

  Epidermal_suprabasal = c(
    "KRT1","KRT10","KRT2","IVL",
    "KRTDAP","DSC1","DSP","SFN"
  ),

  Epidermal_terminal = c(
    "FLG","FLG2","LOR","CDSN",
    "TGM1","CASP14","KLK7","KPRP"
  ),

  Epithelial_stress = c(
    "KRT6A","KRT6B","KRT6C","KRT16",
    "KRT17","S100A7","S100A8","S100A9"
  ),

  # ---- Hair follicle: ORS / stem-associated ----

  HF_ORS = c(
    "KRT5","KRT14","KRT15","KRT17",
    "KRT19","SOX9","KRT75","KRT79"
  ),

  HF_stem_bulge = c(
    "SOX9","LHX2","KRT15","KRT19",
    "CD200","NFATC1","KRT5","KRT14"
  ),

  # ---- Hair follicle: differentiated compartments ----

  HF_IRS_companion = c(
    "KRT25","KRT27","KRT28","KRT71",
    "KRT72","KRT73","KRT74","TCHH",
    "KRT75"
  ),

  HF_hair_shaft = c(
    "KRT31","KRT32","KRT33A","KRT33B",
    "KRT34","KRT35","KRT36",
    "KRT81","KRT83","KRT85"
  ),

  HF_matrix_differentiation = c(
    "LEF1","SHH","HOXC13","FOXN1",
    "MSX2","DLX3","WNT3","WNT10B"
  ),

  # ---- Glandular / appendage ----

  Sebaceous = c(
    "PPARG","FASN","MOGAT1","MOGAT2",
    "AWAT1","AWAT2","PLIN2","THRSP"
  ),

  Sweat_secretory = c(
    "AQP5","DCD","PIP","MUC7",
    "KRT7","KRT8","KRT18","KRT19"
  ),

  Sweat_ductal = c(
    "KRT5","KRT14","KRT19","KRT77",
    "CFTR","SLC12A2","KRT8","KRT18"
  ),

  Myoepithelial = c(
    "KRT5","KRT14","TP63","KRT17",
    "ACTA2","TAGLN","CNN1","MYL9"
  ),

  # ---- State ----

  Cycling = c(
    "MKI67","TOP2A","UBE2C","BIRC5",
    "CENPF","TYMS","RRM2","PCNA"
  ),

  # ---- Non-epithelial controls ----

  Melanocyte = c(
    "MITF","MLANA","PMEL",
    "TYR","TYRP1","DCT"
  ),

  Fibroblast = c(
    "COL1A1","COL1A2","COL3A1",
    "DCN","LUM","PDGFRA"
  ),

  Pericyte = c(
    "RGS5","CSPG4","MCAM",
    "PDGFRB","NOTCH3","KCNJ8"
  ),

  Smooth_muscle = c(
    "ACTA2","TAGLN","MYH11",
    "CNN1","MYOCD","DES"
  ),

  Vascular_EC = c(
    "PECAM1","VWF","EMCN",
    "KDR","CDH5","PLVAP"
  ),

  Lymphatic_EC = c(
    "PROX1","PDPN","LYVE1",
    "FLT4","CCL21"
  ),

  T_NK = c(
    "CD3D","CD3E","TRBC1",
    "NKG7","GNLY","KLRD1"
  ),

  Myeloid = c(
    "LYZ","FCER1G","TYROBP",
    "CTSS","LST1","C1QC"
  )
)

# ============================================================
# 3. MARKER COVERAGE
# ============================================================

coverage <- bind_rows(
  lapply(names(programmes), function(p) {

    requested <- programmes[[p]]
    present <- intersect(requested, rownames(x))
    missing <- setdiff(requested, rownames(x))

    data.frame(
      programme = p,
      n_present = length(present),
      n_requested = length(requested),
      coverage = length(present) / length(requested),
      present_genes = paste(present, collapse=";"),
      missing_genes = paste(missing, collapse=";")
    )
  })
)

write_csv(
  coverage,
  file.path(TDIR, "06B3_marker_programme_coverage.csv")
)

cat("=== MARKER COVERAGE ===\n")

print(
  as.data.frame(
    coverage %>%
      select(programme, n_present, n_requested, coverage)
  ),
  row.names=FALSE
)

if(any(coverage$coverage < 0.50)){
  warning(
    "At least one programme has <50% marker coverage; interpret cautiously."
  )
}

# ============================================================
# 4. SOUPX-CORRECTED EXPRESSION MATRIX
# ============================================================

genes_use <- unique(unlist(programmes))
genes_use <- intersect(genes_use, rownames(x))

mat <- LayerData(
  x,
  assay="RNA",
  layer="data"
)[genes_use, , drop=FALSE]

cat("\nGenes used:", length(genes_use), "\n")
cat("Calculating programme scores...\n")

# ============================================================
# 5. CELL-LEVEL PROGRAMME SCORES
# ============================================================

score_list <- list()

for(p in names(programmes)){

  gg <- intersect(
    programmes[[p]],
    rownames(mat)
  )

  if(length(gg) == 0){
    stop("No markers available for programme: ", p)
  }

  score_list[[p]] <- Matrix::colMeans(
    mat[gg, , drop=FALSE]
  )
}

score_mat <- do.call(cbind, score_list)
rownames(score_mat) <- colnames(x)

# ============================================================
# 6. CLUSTER-LEVEL EVIDENCE
# ============================================================

meta <- x@meta.data %>%
  mutate(
    cell = rownames(.),
    global_cluster = as.character(.data[[CLUSTER_COL]])
  )

cluster_rows <- list()

for(cl in as.character(0:26)){

  cells <- meta$cell[
    meta$global_cluster == cl
  ]

  for(p in colnames(score_mat)){

    vals <- score_mat[cells, p]

    cluster_rows[[length(cluster_rows)+1]] <- data.frame(
      cluster = cl,
      programme = p,
      n_cells = length(cells),
      mean_programme = mean(vals),
      median_programme = median(vals),
      pct_positive = mean(vals > 0) * 100
    )
  }
}

cluster_ev <- bind_rows(cluster_rows)

write_csv(
  cluster_ev,
  file.path(TDIR, "06B3_cluster_programme_evidence.csv")
)

# ============================================================
# 7. RELATIVE ENRICHMENT ACROSS CLUSTERS
#
# z-score is calculated WITHIN programme across clusters.
# It measures relative enrichment, not cell-type truth.
# ============================================================

cluster_z <- cluster_ev %>%
  group_by(programme) %>%
  mutate(
    programme_z = {
      s <- sd(mean_programme, na.rm=TRUE)
      if(is.finite(s) && s > 0) {
        as.numeric(scale(mean_programme))
      } else {
        rep(0, n())
      }
    }
  ) %>%
  ungroup()

write_csv(
  cluster_z,
  file.path(TDIR, "06B3_cluster_programme_zscores.csv")
)

# ============================================================
# 8. TOP-5 RELATIVE PROGRAMMES PER CLUSTER
# ============================================================

rank_tbl <- cluster_z %>%
  group_by(cluster) %>%
  arrange(desc(programme_z), .by_group=TRUE) %>%
  mutate(rank = row_number()) %>%
  filter(rank <= 5) %>%
  ungroup()

write_csv(
  rank_tbl,
  file.path(TDIR, "06B3_top5_programmes_by_cluster.csv")
)

# ============================================================
# 9. PATIENT-LEVEL EVIDENCE
# ============================================================

needed_meta <- c(
  "patient_id",
  "disease",
  "batch"
)

missing_meta <- setdiff(
  needed_meta,
  colnames(meta)
)

if(length(missing_meta) > 0){
  stop(
    "Missing metadata: ",
    paste(missing_meta, collapse=", ")
  )
}

patient_rows <- list()

for(cl in as.character(0:26)){

  mm <- meta %>%
    filter(global_cluster == cl)

  for(pt in unique(mm$patient_id)){

    cells <- mm$cell[
      mm$patient_id == pt
    ]

    for(p in colnames(score_mat)){

      patient_rows[[length(patient_rows)+1]] <- data.frame(
        cluster = cl,
        patient_id = pt,
        disease = unique(
          mm$disease[
            mm$patient_id == pt
          ]
        )[1],
        batch = unique(
          mm$batch[
            mm$patient_id == pt
          ]
        )[1],
        programme = p,
        n_cells = length(cells),
        mean_programme = mean(
          score_mat[cells, p]
        )
      )
    }
  }
}

patient_ev <- bind_rows(patient_rows)

write_csv(
  patient_ev,
  file.path(TDIR, "06B3_patient_cluster_programmes.csv")
)

# ============================================================
# 10. FOCUSED EPITHELIAL / APPENDAGE CLUSTERS
# ============================================================

focus_clusters <- c(
  0,1,2,3,5,6,7,8,
  10,12,13,14,15,16,17,
  20,22,23,24,25
)

focus <- rank_tbl %>%
  filter(
    as.integer(cluster) %in% focus_clusters
  ) %>%
  mutate(
    cluster_num = as.integer(cluster)
  ) %>%
  arrange(cluster_num, rank)

write_csv(
  focus,
  file.path(TDIR, "06B3_epithelial_focus_top5.csv")
)

# ============================================================
# 11. PRINT CLUSTER-BY-CLUSTER REPORT
# ============================================================

cat("\n============================================================\n")
cat("TARGETED PROGRAMME TOP 5\n")
cat("============================================================\n")

for(cl in focus_clusters){

  zz <- focus %>%
    filter(
      as.integer(cluster) == cl
    ) %>%
    arrange(rank)

  cat("\n------------------------------------------------------------\n")
  cat("CLUSTER", cl, "\n")

  for(i in seq_len(nrow(zz))){

    cat(
      sprintf(
        "%d. %-25s z=%6.2f  mean=%7.3f  positive=%5.1f%%\n",
        zz$rank[i],
        zz$programme[i],
        zz$programme_z[i],
        zz$mean_programme[i],
        zz$pct_positive[i]
      )
    )
  }
}

# ============================================================
# 12. INTERNAL POSITIVE CONTROLS
#
# C16: differentiated follicular epithelium expected.
# C23/C24: sweat gland programmes expected.
# ============================================================

cat("\n============================================================\n")
cat("INTERNAL POSITIVE CONTROLS\n")
cat("============================================================\n")

for(cl in c(16,23,24)){

  zz <- cluster_z %>%
    filter(cluster == as.character(cl)) %>%
    arrange(desc(programme_z)) %>%
    slice_head(n=10)

  cat("\nCLUSTER", cl, "\n")

  print(
    as.data.frame(
      zz %>%
        select(
          programme,
          programme_z,
          mean_programme,
          pct_positive
        )
    ),
    row.names=FALSE
  )
}

# ============================================================
# 13. SPECIAL AMBIGUOUS-CLUSTER TABLE
# ============================================================

ambiguous_clusters <- c(
  2,6,8,13,14,15,17,22
)

ambiguous <- cluster_z %>%
  filter(
    as.integer(cluster) %in% ambiguous_clusters
  ) %>%
  mutate(
    cluster_num = as.integer(cluster)
  ) %>%
  arrange(
    cluster_num,
    desc(programme_z)
  )

write_csv(
  ambiguous,
  file.path(
    TDIR,
    "06B3_ambiguous_clusters_all_programmes.csv"
  )
)

cat("\n============================================================\n")
cat("AMBIGUOUS CLUSTERS — TOP 8\n")
cat("============================================================\n")

for(cl in ambiguous_clusters){

  zz <- ambiguous %>%
    filter(cluster_num == cl) %>%
    slice_head(n=8)

  cat("\nCLUSTER", cl, "\n")

  print(
    as.data.frame(
      zz %>%
        select(
          programme,
          programme_z,
          mean_programme,
          pct_positive
        )
    ),
    row.names=FALSE
  )
}

# ============================================================
# 14. HEATMAP
# ============================================================

plot_df <- cluster_z %>%
  filter(
    as.integer(cluster) %in% focus_clusters
  ) %>%
  mutate(
    cluster = factor(
      cluster,
      levels=as.character(focus_clusters)
    ),
    programme = factor(
      programme,
      levels=rev(names(programmes))
    )
  )

p <- ggplot(
  plot_df,
  aes(
    x=cluster,
    y=programme,
    fill=programme_z
  )
) +
  geom_tile() +
  scale_fill_gradient2(
    midpoint=0,
    low="#2166AC",
    mid="white",
    high="#B2182B"
  ) +
  labs(
    x="Global cluster",
    y=NULL,
    fill="Programme\nz-score",
    title="Targeted epithelial and follicular identity audit",
    subtitle="SoupX-corrected expression; canonical frozen global clusters"
  ) +
  theme_classic(base_size=10) +
  theme(
    axis.text.x=element_text(size=9),
    axis.text.y=element_text(size=8)
  )

ggsave(
  file.path(
    FDIR,
    "06B3_targeted_programme_heatmap.pdf"
  ),
  p,
  width=10,
  height=8.5
)

ggsave(
  file.path(
    FDIR,
    "06B3_targeted_programme_heatmap.png"
  ),
  p,
  width=10,
  height=8.5,
  dpi=300
)

# ============================================================
# 15. AUDIT TRAIL
# ============================================================

audit <- c(
  "PIPELINE_V2 — 06B3 targeted epithelial/follicular identity audit",
  paste0("Generated: ", Sys.time()),
  paste0("Cells: ", ncol(x)),
  paste0("Genes: ", nrow(x)),
  "Expression source: SoupX-corrected RNA/data.",
  "Clusters: canonical frozen RNA_snn_res.0 .6 (0-26).",
  "Programmes are annotation evidence, not automatic labels.",
  "Follicular programmes retained for subsequent L2 adjudication.",
  "Cycling treated as state rather than L1 lineage.",
  "Cluster 25 retained in full and remains flagged from prior audit.",
  "No labels assigned.",
  "No cells removed.",
  "No object modified.",
  "No PCA/Harmony/UMAP/neighbors/clusters recalculated."
)

writeLines(
  audit,
  file.path(TDIR, "06B3_audit.txt")
)

cat("\n============================================================\n")
cat("06B3 COMPLETE\n")
cat("============================================================\n")
cat("Tables:", TDIR, "\n")
cat("Figures:", FDIR, "\n")
cat("PASS — targeted evidence only; no labels frozen.\n")

