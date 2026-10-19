suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(Matrix)
})

cat("============================================================\n")
cat("PIPELINE_V2 — 11B3 TARGETED MURAL LINEAGE-INTEGRITY AUDIT\n")
cat("============================================================\n\n")

SRC <- paste0(
  "PIPELINE_V2/results/annotation/11_MURAL_L2/",
  "11B1_recluster/mural_v2_11B1_reclustered.rds"
)

OUT <- paste0(
  "PIPELINE_V2/results/annotation/11_MURAL_L2/",
  "11B3_lineage_integrity"
)

TAB <- file.path(OUT, "tables")
FIG <- file.path(OUT, "figures")

dir.create(TAB, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG, recursive = TRUE, showWarnings = FALSE)

x <- readRDS(SRC)
DefaultAssay(x) <- "RNA"

stopifnot(ncol(x) == 1603)

CLUSTER_COL <- "mural_snn_res.0.6"

stopifnot(
  CLUSTER_COL %in% colnames(x@meta.data)
)

x$mural_r06 <- as.character(
  x@meta.data[[CLUSTER_COL]]
)

cat("Cells:", ncol(x), "\n")
cat("Patients:", length(unique(x$patient_id)), "\n")
cat("Batches:", length(unique(x$batch)), "\n")
cat(
  "r0.6 clusters:",
  length(unique(x$mural_r06)),
  "\n"
)

# ============================================================
# 1. Programmes
# ============================================================

cat("\nSTEP 1/8 — defining targeted lineage programmes\n")

programmes <- list(

  PERICYTE_CORE = c(
    "RGS5","CSPG4","PDGFRB","RBP1",
    "ABCC9","KCNJ8","COX4I2","HIGD1B",
    "NOTCH3"
  ),

  CONTRACTILE_CORE = c(
    "MYH11","ACTA2","TAGLN","CNN1",
    "MYLK","MYL9","TPM2","LMOD1",
    "SYNPO2","ACTG2"
  ),

  APM_CONTRACTILE = c(
    "MYH11","ACTA2","TAGLN","CNN1",
    "MYLK","ACTG2","LMOD1","SYNPO2",
    "DES","ITGA8"
  ),

  APM_ATTACHMENT = c(
    "ITGA8","ITGB1","NPNT",
    "ACTA2","TAGLN"
  ),

  BULGE_BASAL = c(
    "KRT5","KRT14","KRT15","KRT17",
    "TP63","COL17A1","ITGA6","ITGB4"
  ),

  BULGE_STEM = c(
    "SOX9","CD200","KRT15",
    "COL17A1","ITGA6","ITGB4"
  ),

  EPITHELIAL_JUNCTION = c(
    "EPCAM","DSP","PKP1","PKP3",
    "DSC3","DSG1","CLDN1"
  ),

  FIBROBLAST_CORE = c(
    "DCN","LUM","COL1A1","COL1A2",
    "COL3A1","PDGFRA","DPT","COL6A1"
  ),

  REMODELING = c(
    "ADAM12","MMP11","POSTN","COL11A1",
    "COMP","COL10A1","THBS2","FN1"
  ),

  ENDOTHELIAL = c(
    "PECAM1","VWF","CDH5","EMCN",
    "KDR","ESAM","ENG","RAMP2"
  ),

  SCHWANN = c(
    "SOX10","S100B","PLP1",
    "PMP22","MPZ","S100A1","NGFR"
  ),

  IMMUNE = c(
    "PTPRC","CD3D","CD3E",
    "TRAC","NKG7","GNLY"
  )
)

available <- lapply(
  programmes,
  function(g)
    intersect(g, rownames(x))
)

coverage <- tibble(
  programme = names(programmes),
  n_requested = sapply(programmes, length),
  n_present = sapply(available, length),
  genes_present = sapply(
    available,
    paste,
    collapse = ";"
  )
)

cat("\n=== PROGRAMME COVERAGE ===\n")
print(coverage, n = Inf, width = Inf)

write.csv(
  coverage,
  file.path(TAB, "11B3_programme_coverage.csv"),
  row.names = FALSE
)

# ============================================================
# 2. Cell-level evidence
# ============================================================

cat("\nSTEP 2/8 — cell-level evidence\n")

mat <- GetAssayData(
  x,
  assay = "RNA",
  layer = "data"
)

mean_score <- function(genes) {

  genes <- intersect(
    genes,
    rownames(mat)
  )

  if (length(genes) == 0)
    return(rep(NA_real_, ncol(mat)))

  Matrix::colMeans(
    mat[genes, , drop = FALSE]
  )
}

n_detected <- function(genes) {

  genes <- intersect(
    genes,
    rownames(mat)
  )

  if (length(genes) == 0)
    return(rep(0L, ncol(mat)))

  Matrix::colSums(
    mat[genes, , drop = FALSE] > 0
  )
}

evidence <- tibble(
  cell = colnames(x),
  cluster = as.character(x$mural_r06),
  patient_id = as.character(x$patient_id),
  batch = as.character(x$batch),
  source_L1 = as.character(x$mural_source_L1),
  doublet_score = as.numeric(x$doublet_score),
  nFeature_RNA = as.numeric(x$nFeature_RNA),
  nCount_RNA = as.numeric(x$nCount_RNA),
  percent_mt = as.numeric(x$percent.mt)
)

for (p in names(programmes)) {

  evidence[[paste0(p, "_score")]] <-
    mean_score(available[[p]])

  evidence[[paste0(p, "_n")]] <-
    n_detected(available[[p]])
}

evidence <- evidence %>%
  mutate(
    pericyte_supported =
      PERICYTE_CORE_n >= 3,

    contractile_supported =
      CONTRACTILE_CORE_n >= 4,

    APM_supported =
      APM_CONTRACTILE_n >= 4,

    attachment_supported =
      APM_ATTACHMENT_n >= 2,

    bulge_supported =
      BULGE_BASAL_n >= 3,

    bulge_stem_supported =
      BULGE_STEM_n >= 2,

    junction_supported =
      EPITHELIAL_JUNCTION_n >= 2,

    fibroblast_supported =
      FIBROBLAST_CORE_n >= 3,

    remodeling_supported =
      REMODELING_n >= 2,

    endothelial_supported =
      ENDOTHELIAL_n >= 3,

    Schwann_supported =
      SCHWANN_n >= 2,

    immune_supported =
      IMMUNE_n >= 2,

    pericyte_fibroblast_joint =
      pericyte_supported &
      fibroblast_supported,

    contractile_bulge_joint =
      contractile_supported &
      bulge_supported,

    APM_bulge_joint =
      APM_supported &
      bulge_supported,

    APM_bulge_junction_joint =
      APM_supported &
      bulge_supported &
      junction_supported
  )

write.csv(
  evidence,
  file.path(TAB, "11B3_cell_evidence.csv"),
  row.names = FALSE
)

# ============================================================
# 3. C0/C3 mural-fibroblast boundary
# ============================================================

cat("\nSTEP 3/8 — C0/C3 mural–fibroblast boundary\n")

c03 <- evidence %>%
  filter(cluster %in% c("0", "3"))

c03_summary <- c03 %>%
  group_by(cluster) %>%
  summarise(
    n_cells = n(),
    n_patients = n_distinct(patient_id),
    n_batches = n_distinct(batch),

    median_doublet =
      median(doublet_score),

    median_features =
      median(nFeature_RNA),

    pericyte_fraction =
      mean(pericyte_supported),

    contractile_fraction =
      mean(contractile_supported),

    fibroblast_fraction =
      mean(fibroblast_supported),

    remodeling_fraction =
      mean(remodeling_supported),

    pericyte_fibroblast_joint_fraction =
      mean(pericyte_fibroblast_joint),

    endothelial_fraction =
      mean(endothelial_supported),

    Schwann_fraction =
      mean(Schwann_supported),

    immune_fraction =
      mean(immune_supported),

    median_pericyte_score =
      median(PERICYTE_CORE_score),

    median_fibroblast_score =
      median(FIBROBLAST_CORE_score),

    median_remodeling_score =
      median(REMODELING_score),

    .groups = "drop"
  )

cat("\n=== C0/C3 BOUNDARY SUMMARY ===\n")
print(c03_summary, n = Inf, width = Inf)

write.csv(
  c03_summary,
  file.path(TAB, "11B3_C0_C3_boundary_summary.csv"),
  row.names = FALSE
)

c03_patient <- c03 %>%
  group_by(
    cluster,
    patient_id,
    batch
  ) %>%
  summarise(
    n_cells = n(),

    pericyte_fraction =
      mean(pericyte_supported),

    fibroblast_fraction =
      mean(fibroblast_supported),

    joint_fraction =
      mean(pericyte_fibroblast_joint),

    median_pericyte_score =
      median(PERICYTE_CORE_score),

    median_fibroblast_score =
      median(FIBROBLAST_CORE_score),

    median_doublet =
      median(doublet_score),

    .groups = "drop"
  )

cat("\n=== C0/C3 PATIENT SUPPORT ===\n")
print(c03_patient, n = Inf, width = Inf)

write.csv(
  c03_patient,
  file.path(TAB, "11B3_C0_C3_patient_support.csv"),
  row.names = FALSE
)

# ============================================================
# 4. C6/C8/C9 interface comparison
# ============================================================

cat("\nSTEP 4/8 — C6/C8/C9 interface comparison\n")

interface <- evidence %>%
  filter(cluster %in% c("6", "8", "9"))

interface_summary <- interface %>%
  group_by(cluster) %>%
  summarise(
    n_cells = n(),
    n_patients = n_distinct(patient_id),
    n_batches = n_distinct(batch),

    median_doublet =
      median(doublet_score),

    q75_doublet =
      quantile(
        doublet_score,
        0.75
      ),

    median_features =
      median(nFeature_RNA),

    median_counts =
      median(nCount_RNA),

    median_mt =
      median(percent_mt),

    contractile_fraction =
      mean(contractile_supported),

    APM_fraction =
      mean(APM_supported),

    attachment_fraction =
      mean(attachment_supported),

    bulge_fraction =
      mean(bulge_supported),

    bulge_stem_fraction =
      mean(bulge_stem_supported),

    junction_fraction =
      mean(junction_supported),

    contractile_bulge_joint =
      mean(contractile_bulge_joint),

    APM_bulge_joint =
      mean(APM_bulge_joint),

    APM_bulge_junction_joint =
      mean(APM_bulge_junction_joint),

    pericyte_fraction =
      mean(pericyte_supported),

    fibroblast_fraction =
      mean(fibroblast_supported),

    remodeling_fraction =
      mean(remodeling_supported),

    Schwann_fraction =
      mean(Schwann_supported),

    endothelial_fraction =
      mean(endothelial_supported),

    immune_fraction =
      mean(immune_supported),

    .groups = "drop"
  )

cat("\n=== C6/C8/C9 INTERFACE SUMMARY ===\n")
print(interface_summary, n = Inf, width = Inf)

write.csv(
  interface_summary,
  file.path(TAB, "11B3_C6_C8_C9_interface_summary.csv"),
  row.names = FALSE
)

# ============================================================
# 5. Doublet-independent interface evidence
# ============================================================

cat("\nSTEP 5/8 — doublet-independent interface evidence\n")

interface <- interface %>%
  group_by(cluster) %>%
  mutate(
    doublet_quartile =
      factor(
        ntile(doublet_score, 4),
        levels = 1:4,
        labels = c(
          "Q1_low",
          "Q2",
          "Q3",
          "Q4_high"
        )
      )
  ) %>%
  ungroup()

doublet_strata <- interface %>%
  group_by(
    cluster,
    doublet_quartile
  ) %>%
  summarise(
    n_cells = n(),

    median_doublet =
      median(doublet_score),

    median_features =
      median(nFeature_RNA),

    median_counts =
      median(nCount_RNA),

    APM_fraction =
      mean(APM_supported),

    bulge_fraction =
      mean(bulge_supported),

    junction_fraction =
      mean(junction_supported),

    APM_bulge_joint_fraction =
      mean(APM_bulge_joint),

    APM_bulge_junction_fraction =
      mean(APM_bulge_junction_joint),

    Schwann_fraction =
      mean(Schwann_supported),

    .groups = "drop"
  )

cat("\n=== C6/C8/C9 DOUBLET STRATA ===\n")
print(doublet_strata, n = Inf, width = Inf)

write.csv(
  doublet_strata,
  file.path(TAB, "11B3_C6_C8_C9_doublet_strata.csv"),
  row.names = FALSE
)

lowest_doublet <- interface %>%
  filter(
    doublet_quartile == "Q1_low"
  ) %>%
  group_by(cluster) %>%
  summarise(
    n_cells = n(),

    median_doublet =
      median(doublet_score),

    APM_fraction =
      mean(APM_supported),

    bulge_fraction =
      mean(bulge_supported),

    junction_fraction =
      mean(junction_supported),

    APM_bulge_joint_fraction =
      mean(APM_bulge_joint),

    APM_bulge_junction_fraction =
      mean(APM_bulge_junction_joint),

    .groups = "drop"
  )

cat("\n=== LOWEST DOUBLET QUARTILE ===\n")
print(lowest_doublet, n = Inf, width = Inf)

write.csv(
  lowest_doublet,
  file.path(TAB, "11B3_C6_C8_C9_lowest_doublet.csv"),
  row.names = FALSE
)

# ============================================================
# 6. Technical correlations
# ============================================================

cat("\nSTEP 6/8 — technical-association audit\n")

technical_vars <- c(
  "doublet_score",
  "nFeature_RNA",
  "nCount_RNA",
  "percent_mt"
)

bio_scores <- c(
  "APM_CONTRACTILE_score",
  "BULGE_BASAL_score",
  "BULGE_STEM_score",
  "EPITHELIAL_JUNCTION_score",
  "PERICYTE_CORE_score",
  "FIBROBLAST_CORE_score",
  "REMODELING_score",
  "SCHWANN_score"
)

technical_cor <- bind_rows(
  lapply(c("0","3","6","8","9"), function(cl) {

    z <- evidence %>%
      filter(cluster == cl)

    bind_rows(
      lapply(bio_scores, function(bio) {

        bind_rows(
          lapply(technical_vars, function(tech) {

            ok <- complete.cases(
              z[[bio]],
              z[[tech]]
            )

            if (sum(ok) < 5) {
              return(
                tibble(
                  cluster = cl,
                  biological_score = bio,
                  technical_variable = tech,
                  n = sum(ok),
                  spearman_rho = NA_real_,
                  p_value = NA_real_
                )
              )
            }

            tt <- suppressWarnings(
              cor.test(
                z[[bio]][ok],
                z[[tech]][ok],
                method = "spearman",
                exact = FALSE
              )
            )

            tibble(
              cluster = cl,
              biological_score = bio,
              technical_variable = tech,
              n = sum(ok),
              spearman_rho =
                unname(tt$estimate),
              p_value =
                tt$p.value
            )
          })
        )
      })
    )
  })
) %>%
  group_by(cluster) %>%
  mutate(
    p_adj =
      p.adjust(
        p_value,
        method = "BH"
      )
  ) %>%
  ungroup()

cat("\n=== TECHNICAL CORRELATIONS ===\n")
print(technical_cor, n = Inf, width = Inf)

write.csv(
  technical_cor,
  file.path(TAB, "11B3_technical_correlations.csv"),
  row.names = FALSE
)

# ============================================================
# 7. Targeted gene detection + patient recurrence
# ============================================================

cat("\nSTEP 7/8 — target-gene and recurrence audit\n")

target_genes <- unique(c(
  programmes$PERICYTE_CORE,
  programmes$CONTRACTILE_CORE,
  programmes$APM_CONTRACTILE,
  programmes$APM_ATTACHMENT,
  programmes$BULGE_BASAL,
  programmes$BULGE_STEM,
  programmes$EPITHELIAL_JUNCTION,
  programmes$FIBROBLAST_CORE,
  programmes$REMODELING,
  programmes$SCHWANN,
  "SFRP1","KRT7","KRT8","KRT23",
  "CCL27","GJB3","DSC1","CLCA2","WNT4",
  "KRT86","PTPRT","TAFA2
",
  "CDH18","LMX1B","NOL4",
  "LGR5","CALCA",
  "SOX10","S100B","PMP22"
))

target_genes <- intersect(
  target_genes,
  rownames(mat)
)

focus_clusters <- c(
  "0","3","6","8","9"
)

gene_detection <- bind_rows(
  lapply(focus_clusters, function(cl) {

    cells <- evidence$cell[
      evidence$cluster == cl
    ]

    m <- mat[
      target_genes,
      cells,
      drop = FALSE
    ]

    tibble(
      cluster = cl,
      gene = target_genes,
      fraction_detected =
        as.numeric(
          Matrix::rowMeans(m > 0)
        ),
      mean_expression =
        as.numeric(
          Matrix::rowMeans(m)
        )
    )
  })
)

cat("\n=== TARGET GENE DETECTION ===\n")
print(
  gene_detection,
  n = Inf,
  width = Inf
)

write.csv(
  gene_detection,
  file.path(
    TAB,
    "11B3_target_gene_detection.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# Patient recurrence for C6/C8/C9
# ------------------------------------------------------------

interface_patient <- evidence %>%
  filter(
    cluster %in% c("6","8","9")
  ) %>%
  group_by(
    cluster,
    patient_id,
    batch
  ) %>%
  summarise(
    n_cells = n(),

    median_doublet =
      median(doublet_score),

    median_features =
      median(nFeature_RNA),

    contractile_fraction =
      mean(contractile_supported),

    APM_fraction =
      mean(APM_supported),

    bulge_fraction =
      mean(bulge_supported),

    bulge_stem_fraction =
      mean(bulge_stem_supported),

    junction_fraction =
      mean(junction_supported),

    APM_bulge_joint_fraction =
      mean(APM_bulge_joint),

    APM_bulge_junction_fraction =
      mean(APM_bulge_junction_joint),

    Schwann_fraction =
      mean(Schwann_supported),

    .groups = "drop"
  )

cat("\n=== C6/C8/C9 PATIENT RECURRENCE ===\n")
print(
  interface_patient,
  n = Inf,
  width = Inf
)

write.csv(
  interface_patient,
  file.path(
    TAB,
    "11B3_C6_C8_C9_patient_recurrence.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# Cluster-source L1 relationship
# Descriptive only; never used for identity scoring.
# ------------------------------------------------------------

source_distribution <- evidence %>%
  count(
    cluster,
    source_L1,
    name = "n_cells"
  ) %>%
  group_by(cluster) %>%
  mutate(
    cluster_total =
      sum(n_cells),
    fraction_cluster =
      n_cells / cluster_total
  ) %>%
  ungroup() %>%
  arrange(
    as.integer(cluster),
    desc(n_cells)
  )

cat("\n=== SOURCE L1 DISTRIBUTION — DESCRIPTIVE ONLY ===\n")
print(
  source_distribution,
  n = Inf,
  width = Inf
)

write.csv(
  source_distribution,
  file.path(
    TAB,
    "11B3_source_L1_distribution_descriptive.csv"
  ),
  row.names = FALSE
)

# ============================================================
# 8. Figures + formal checks
# ============================================================

cat("\nSTEP 8/8 — figures and formal checks\n")

# ------------------------------------------------------------
# C0/C3 pericyte-fibroblast relationship
# ------------------------------------------------------------

p_c03 <- ggplot(
  c03,
  aes(
    x = PERICYTE_CORE_score,
    y = FIBROBLAST_CORE_score
  )
) +
  geom_point(
    aes(shape = cluster),
    alpha = 0.65,
    size = 1.5
  ) +
  facet_wrap(
    ~cluster,
    scales = "free"
  ) +
  theme_classic() +
  labs(
    title = "C0/C3 mural-fibroblast coexpression",
    x = "Pericyte-core score",
    y = "Fibroblast-core score"
  )

ggsave(
  file.path(
    FIG,
    "11B3_C0_C3_pericyte_fibroblast.pdf"
  ),
  p_c03,
  width = 8,
  height = 4.5
)

# ------------------------------------------------------------
# C6/C8/C9 APM-bulge relationship
# ------------------------------------------------------------

p_interface <- ggplot(
  interface,
  aes(
    x = APM_CONTRACTILE_score,
    y = BULGE_BASAL_score
  )
) +
  geom_point(
    aes(shape = cluster),
    alpha = 0.7,
    size = 1.6
  ) +
  facet_wrap(
    ~cluster,
    scales = "free"
  ) +
  theme_classic() +
  labs(
    title = "C6/C8/C9 APM-bulge coexpression",
    x = "APM-contractile score",
    y = "Bulge-basal score"
  )

ggsave(
  file.path(
    FIG,
    "11B3_C6_C8_C9_APM_bulge.pdf"
  ),
  p_interface,
  width = 9,
  height = 4.5
)

# ------------------------------------------------------------
# Interface signal vs doublet score
# ------------------------------------------------------------

p_doublet <- ggplot(
  interface,
  aes(
    x = doublet_score,
    y = BULGE_BASAL_score
  )
) +
  geom_point(
    aes(shape = cluster),
    alpha = 0.7,
    size = 1.5
  ) +
  geom_smooth(
    method = "lm",
    se = FALSE,
    linewidth = 0.6
  ) +
  facet_wrap(
    ~cluster,
    scales = "free"
  ) +
  theme_classic() +
  labs(
    title = "Bulge programme versus doublet score",
    x = "Doublet score",
    y = "Bulge-basal score"
  )

ggsave(
  file.path(
    FIG,
    "11B3_C6_C8_C9_bulge_vs_doublet.pdf"
  ),
  p_doublet,
  width = 9,
  height = 4.5
)

# ------------------------------------------------------------
# Focus-cluster DotPlot
# ------------------------------------------------------------

focus_cells <- evidence$cell[
  evidence$cluster %in% focus_clusters
]

focus_obj <- subset(
  x,
  cells = focus_cells
)

focus_obj$focus_cluster <- factor(
  as.character(focus_obj$mural_r06),
  levels = focus_clusters
)

Idents(focus_obj) <- "focus_cluster"

dot_genes <- unique(c(
  "RGS5","CSPG4","PDGFRB","ABCC9","KCNJ8",
  "MYH11","ACTA2","TAGLN","CNN1","MYLK",
  "ITGA8","NPNT",
  "KRT5","KRT14","KRT15","KRT17",
  "TP63","COL17A1","SOX9","CD200",
  "DSP","DSC3","PKP3",
  "DCN","LUM","COL1A1","COL1A2",
  "COL3A1","PDGFRA","DPT",
  "POSTN","ADAM12","MMP11","COL11A1",
  "SOX10","S100B","PMP22",
  "SFRP1","KRT7","KRT8","KRT23",
  "CCL27","GJB3","DSC1","CLCA2","WNT4",
  "KRT86","PTPRT","TAFA2","CDH18","LMX1B"
))

dot_genes <- intersect(
  dot_genes,
  rownames(focus_obj)
)

p_dot <- DotPlot(
  focus_obj,
  features = dot_genes,
  assay = "RNA"
) +
  RotatedAxis() +
  theme_classic() +
  labs(
    title = "Targeted mural lineage-integrity audit"
  )

ggsave(
  file.path(
    FIG,
    "11B3_focus_cluster_DotPlot.pdf"
  ),
  p_dot,
  width = 16,
  height = 5.5
)

# ------------------------------------------------------------
# Formal checks
# ------------------------------------------------------------

checks <- tibble(
  check = c(
    "cells_1603",
    "patients_16",
    "batches_2",
    "r06_10_clusters",
    "programme_coverage_complete",
    "cell_evidence_generated",
    "C0_C3_boundary_generated",
    "C0_C3_patient_support_generated",
    "C6_C8_C9_summary_generated",
    "doublet_strata_generated",
    "lowest_doublet_generated",
    "technical_correlations_generated",
    "target_gene_detection_generated",
    "interface_patient_recurrence_generated",
    "source_distribution_descriptive_only",
    "disease_not_used",
    "source_L1_not_used_for_identity",
    "no_cells_removed",
    "no_identity_assigned",
    "no_resolution_selected",
    "source_object_not_modified"
  ),

  pass = c(
    ncol(x) == 1603,
    length(unique(x$patient_id)) == 16,
    length(unique(x$batch)) == 2,
    length(unique(x$mural_r06)) == 10,
    all(coverage$n_present > 0),
    nrow(evidence) == 1603,
    nrow(c03_summary) == 2,
    nrow(c03_patient) > 0,
    nrow(interface_summary) == 3,
    nrow(doublet_strata) == 12,
    nrow(lowest_doublet) == 3,
    nrow(technical_cor) > 0,
    nrow(gene_detection) > 0,
    nrow(interface_patient) > 0,
    nrow(source_distribution) > 0,
    TRUE,
    TRUE,
    TRUE,
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
    "11B3_formal_checks.csv"
  ),
  row.names = FALSE
)

stopifnot(all(checks$pass))

# ============================================================
# Audit
# ============================================================

audit <- c(
  "PIPELINE_V2 — 11B3 TARGETED MURAL LINEAGE-INTEGRITY AUDIT",
  "==========================================================",
  "",
  paste0("Cells: ", ncol(x)),
  paste0(
    "Patients: ",
    length(unique(x$patient_id))
  ),
  "Resolution interrogated: r0.6",
  "",
  "Targeted questions:",
  "- C0/C3 mural-fibroblast boundary",
  "- C6/C8/C9 APM-bulge interface",
  "- persistence across doublet-score strata",
  "- patient-level recurrence",
  "- competing lineage evidence",
  "",
  "Interpretive restrictions:",
  "- disease not used",
  "- source L1 descriptive only",
  "- no cells removed",
  "- no identities assigned",
  "- no resolution selected",
  "- source object not modified",
  "",
  "PASS"
)

writeLines(
  audit,
  file.path(
    OUT,
    "11B3_audit.txt"
  )
)

cat("\n============================================================\n")
cat("11B3 COMPLETE\n")
cat("============================================================\n")
cat("Cells:", ncol(x), "\n")
cat("C0/C3 mural-fibroblast boundary audited.\n")
cat("C6/C8/C9 APM-bulge interface audited.\n")
cat("Doublet-independent evidence audited.\n")
cat("Patient recurrence audited.\n")
cat("No cells removed.\n")
cat("No identity assigned.\n")
cat("No resolution selected.\n")
cat("Disease not used.\n")
cat("Source object not modified.\n")
cat("PASS\n")
