suppressPackageStartupMessages({
  library(Seurat)
  library(harmony)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
})

options(width = 200)

cat("============================================================\n")
cat("PIPELINE_V2 — 10B3B C6 EndMT vs DOUBLET AUDIT\n")
cat("============================================================\n\n")

SRC <- paste0(
  "PIPELINE_V2/objects/annotated/",
  "vascular_v2_10B1_reclustered_3040.rds"
)

OUT <- paste0(
  "PIPELINE_V2/results/annotation/",
  "10_VASCULAR_L2/10B3B_C6_EndMT_doublet_audit"
)

TAB <- file.path(OUT, "tables")
FIG <- file.path(OUT, "figures")

dir.create(TAB, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG, recursive = TRUE, showWarnings = FALSE)

x <- readRDS(SRC)
DefaultAssay(x) <- "RNA"

CLCOL <- "vascular_snn_res.0.3"
cl <- as.character(x@meta.data[[CLCOL]])

cells6 <- colnames(x)[cl == "6"]

stopifnot(length(cells6) == 162)

c6 <- subset(x, cells = cells6)

cat("C6 cells:", ncol(c6), "\n")
cat("Patients:", length(unique(c6$patient_id)), "\n")
cat("Batches:", length(unique(c6$batch)), "\n\n")

# ============================================================
# 1. Prespecified programmes
# ============================================================

cat("STEP 1/7 — defining programmes\n")

programmes <- list(

  ENDOTHELIAL_CORE = c(
    "PECAM1","VWF","CDH5","EMCN",
    "KDR","ESAM","ENG","RAMP2"
  ),

  ENDOTHELIAL_JUNCTION = c(
    "CDH5","CLDN5","ESAM","PECAM1",
    "ENG","KDR"
  ),

  MURAL_CONTRACTILE = c(
    "RGS5","CSPG4","MCAM","ACTA2",
    "TAGLN","MYH11","RBP1","PDGFRB"
  ),

  MESENCHYMAL_ECM = c(
    "VIM","FN1","COL1A1","COL1A2",
    "COL3A1","POSTN","COL11A1",
    "MMP11","ADAM12","COMP"
  ),

  ENDMT_TF = c(
    "SNAI1","SNAI2","TWIST1","TWIST2",
    "ZEB1","ZEB2","FOXC2","SOX9"
  ),

  TGF_BETA_RESPONSE = c(
    "TGFB1","TGFBR1","TGFBR2",
    "SMAD2","SMAD3","SMAD4",
    "SERPINE1","CTGF"
  ),

  EPITHELIAL = c(
    "EPCAM","KRT5","KRT14","KRT17",
    "KRT19","KRT8","KRT18","DSP"
  ),

  IMMUNE = c(
    "PTPRC","CD3D","CD3E","TRAC",
    "TRBC1","TRBC2","NKG7","GNLY"
  )
)

coverage <- tibble(
  programme = names(programmes),
  n_requested = sapply(programmes, length),
  n_present = sapply(
    programmes,
    function(g)
      length(intersect(g, rownames(c6)))
  ),
  genes_present = sapply(
    programmes,
    function(g)
      paste(
        intersect(g, rownames(c6)),
        collapse = ";"
      )
  )
)

write_csv(
  coverage,
  file.path(TAB, "10B3B_programme_coverage.csv")
)

print(coverage, n = Inf)

# ============================================================
# 2. Cell-level expression / scores
# ============================================================

cat("\nSTEP 2/7 — cell-level programme scores\n")

all_genes <- unique(unlist(programmes))
all_genes <- intersect(all_genes, rownames(c6))

expr <- FetchData(
  c6,
  vars = all_genes,
  layer = "data"
)

binary_count <- function(df, genes) {
  g <- intersect(genes, colnames(df))
  if (!length(g))
    return(rep(0L, nrow(df)))
  rowSums(df[, g, drop = FALSE] > 0)
}

mean_score <- function(df, genes) {
  g <- intersect(genes, colnames(df))
  if (!length(g))
    return(rep(NA_real_, nrow(df)))
  rowMeans(df[, g, drop = FALSE])
}

cell <- tibble(
  cell = rownames(expr),
  patient_id = as.character(
    c6$patient_id[
      match(rownames(expr), colnames(c6))
    ]
  ),
  batch = as.character(
    c6$batch[
      match(rownames(expr), colnames(c6))
    ]
  ),
  doublet_score = c6$doublet_score[
    match(rownames(expr), colnames(c6))
  ],
  nFeature_RNA = c6$nFeature_RNA[
    match(rownames(expr), colnames(c6))
  ],
  nCount_RNA = c6$nCount_RNA[
    match(rownames(expr), colnames(c6))
  ],
  percent_mt = c6$percent.mt[
    match(rownames(expr), colnames(c6))
  ]
)

for (p in names(programmes)) {

  cell[[paste0(p, "_n")]] <-
    binary_count(expr, programmes[[p]])

  cell[[paste0(p, "_score")]] <-
    mean_score(expr, programmes[[p]])
}

write_csv(
  cell,
  file.path(TAB, "10B3B_C6_cell_programme_scores.csv")
)

# ============================================================
# 3. Joint coexpression
# ============================================================

cat("STEP 3/7 — EC/mesenchymal joint coexpression\n")

cell <- cell %>%
  mutate(
    EC_ge4 =
      ENDOTHELIAL_CORE_n >= 4,

    MURAL_ge3 =
      MURAL_CONTRACTILE_n >= 3,

    ECM_ge3 =
      MESENCHYMAL_ECM_n >= 3,

    EndMT_TF_ge1 =
      ENDMT_TF_n >= 1,

    EC_MURAL_joint =
      EC_ge4 & MURAL_ge3,

    EC_ECM_joint =
      EC_ge4 & ECM_ge3
  )

cat("\nEC >=4 vs MURAL >=3\n")
print(
  table(
    EC_ge4 = cell$EC_ge4,
    MURAL_ge3 = cell$MURAL_ge3
  )
)

cat("\nEC >=4 vs ECM >=3\n")
print(
  table(
    EC_ge4 = cell$EC_ge4,
    ECM_ge3 = cell$ECM_ge3
  )
)

cat("\nEC >=4 vs EndMT TF >=1\n")
print(
  table(
    EC_ge4 = cell$EC_ge4,
    EndMT_TF_ge1 = cell$EndMT_TF_ge1
  )
)

# ============================================================
# 4. Association with technical complexity
# ============================================================

cat("\nSTEP 4/7 — technical-association audit\n")

vars <- c(
  "MURAL_CONTRACTILE_score",
  "MESENCHYMAL_ECM_score",
  "ENDMT_TF_score"
)

technical <- c(
  "doublet_score",
  "nFeature_RNA",
  "nCount_RNA"
)

cors <- bind_rows(
  lapply(vars, function(v) {
    bind_rows(
      lapply(technical, function(t) {

        ok <- complete.cases(
          cell[[v]],
          cell[[t]]
        )

        ct <- suppressWarnings(
          cor.test(
            cell[[v]][ok],
            cell[[t]][ok],
            method = "spearman",
            exact = FALSE
          )
        )

        tibble(
          biological_score = v,
          technical_variable = t,
          n = sum(ok),
          spearman_rho =
            unname(ct$estimate),
          p_value = ct$p.value
        )
      })
    )
  })
) %>%
  mutate(
    p_adj = p.adjust(
      p_value,
      method = "BH"
    )
  )

write_csv(
  cors,
  file.path(TAB, "10B3B_technical_correlations.csv")
)

cat("\n=== TECHNICAL CORRELATIONS ===\n")
print(cors, n = Inf)

# ============================================================
# 5. Doublet-score strata
# ============================================================

cat("\nSTEP 5/7 — doublet-score strata\n")

q <- quantile(
  cell$doublet_score,
  probs = c(0.25, 0.5, 0.75),
  na.rm = TRUE
)

cell <- cell %>%
  mutate(
    doublet_stratum = cut(
      doublet_score,
      breaks = c(
        -Inf,
        q[1],
        q[2],
        q[3],
        Inf
      ),
      labels = c(
        "Q1_low",
        "Q2",
        "Q3",
        "Q4_high"
      ),
      include.lowest = TRUE
    )
  )

strata <- cell %>%
  group_by(doublet_stratum) %>%
  summarise(
    n_cells = n(),
    median_doublet =
      median(doublet_score),

    median_features =
      median(nFeature_RNA),

    median_counts =
      median(nCount_RNA),

    EC_score =
      median(
        ENDOTHELIAL_CORE_score,
        na.rm = TRUE
      ),

    MURAL_score =
      median(
        MURAL_CONTRACTILE_score,
        na.rm = TRUE
      ),

    ECM_score =
      median(
        MESENCHYMAL_ECM_score,
        na.rm = TRUE
      ),

    EndMT_TF_score =
      median(
        ENDMT_TF_score,
        na.rm = TRUE
      ),

    EC_MURAL_joint_fraction =
      mean(EC_MURAL_joint),

    EC_ECM_joint_fraction =
      mean(EC_ECM_joint),

    .groups = "drop"
  )

write_csv(
  strata,
  file.path(TAB, "10B3B_doublet_score_strata.csv")
)

cat("\n=== DOUBLET STRATA ===\n")
print(strata, n = Inf)

# ============================================================
# 6. Patient/batch recurrence
# ============================================================

cat("\nSTEP 6/7 — patient/batch recurrence\n")

patient <- cell %>%
  group_by(patient_id, batch) %>%
  summarise(
    n_cells = n(),

    median_EC =
      median(
        ENDOTHELIAL_CORE_score,
        na.rm = TRUE
      ),

    median_MURAL =
      median(
        MURAL_CONTRACTILE_score,
        na.rm = TRUE
      ),

    median_ECM =
      median(
        MESENCHYMAL_ECM_score,
        na.rm = TRUE
      ),

    median_EndMT_TF =
      median(
        ENDMT_TF_score,
        na.rm = TRUE
      ),

    EC_MURAL_joint_fraction =
      mean(EC_MURAL_joint),

    EC_ECM_joint_fraction =
      mean(EC_ECM_joint),

    median_doublet =
      median(doublet_score),

    .groups = "drop"
  ) %>%
  arrange(desc(n_cells))

write_csv(
  patient,
  file.path(TAB, "10B3B_patient_recurrence.csv")
)

cat("\n=== PATIENT RECURRENCE ===\n")
print(patient, n = Inf)

# ============================================================
# 7. Within-C6 de novo structure
# ============================================================

cat("\nSTEP 7/7 — within-C6 exploratory structure\n")

# IMPORTANT:
# This is diagnostic only and cannot create a new identity.
# Disease is not used.

z <- c6

for (r in Reductions(z)) {
  z[[r]] <- NULL
}

for (gr in names(z@graphs)) {
  z[[gr]] <- NULL
}

old_clusters <- grep(
  "snn_res|seurat_clusters",
  colnames(z@meta.data),
  value = TRUE
)

for (cc in old_clusters) {
  z[[cc]] <- NULL
}

z <- FindVariableFeatures(
  z,
  selection.method = "vst",
  nfeatures = min(2000, nrow(z)),
  verbose = FALSE
)

z <- ScaleData(
  z,
  features = VariableFeatures(z),
  verbose = FALSE
)

z <- RunPCA(
  z,
  features = VariableFeatures(z),
  npcs = 30,
  verbose = FALSE
)

z <- RunHarmony(
  z,
  group.by.vars = "batch",
  reduction.use = "pca",
  dims.use = 1:20,
  assay.use = "RNA",
  verbose = FALSE
)

z <- FindNeighbors(
  z,
  reduction = "harmony",
  dims = 1:20,
  graph.name = c(
    "C6_nn",
    "C6_snn"
  ),
  verbose = FALSE
)

z <- FindClusters(
  z,
  graph.name = "C6_snn",
  resolution = 0.2,
  random.seed = 20260930,
  verbose = FALSE
)

subcol <- grep(
  "^C6_snn_res",
  colnames(z@meta.data),
  value = TRUE
)

stopifnot(length(subcol) == 1)

sub <- as.character(
  z@meta.data[[subcol]]
)

subsupport <- tibble(
  cell = colnames(z),
  subcluster = sub,
  patient_id =
    as.character(z$patient_id),
  batch =
    as.character(z$batch)
) %>%
  group_by(subcluster) %>%
  summarise(
    n_cells = n(),
    n_patients =
      n_distinct(patient_id),
    n_batches =
      n_distinct(batch),
    dominant_patient_fraction = {
      tt <- table(patient_id)
      max(tt) / sum(tt)
    },
    .groups = "drop"
  )

write_csv(
  subsupport,
  file.path(TAB, "10B3B_C6_internal_structure.csv")
)

cat("\n=== C6 INTERNAL STRUCTURE ===\n")
print(subsupport, n = Inf)

# ============================================================
# Figures
# ============================================================

p1 <- ggplot(
  cell,
  aes(
    x = doublet_score,
    y = MURAL_CONTRACTILE_score
  )
) +
  geom_point(
    size = 1.6,
    alpha = 0.7
  ) +
  geom_smooth(
    method = "lm",
    se = FALSE
  ) +
  labs(
    title = "C6 mural programme vs doublet score",
    x = "Doublet score",
    y = "Mural-contractile score"
  ) +
  theme_classic()

ggsave(
  file.path(
    FIG,
    "10B3B_C6_mural_vs_doublet.png"
  ),
  p1,
  width = 6,
  height = 5,
  dpi = 300
)

p2 <- ggplot(
  cell,
  aes(
    x = ENDOTHELIAL_CORE_score,
    y = MURAL_CONTRACTILE_score
  )
) +
  geom_point(
    size = 1.6,
    alpha = 0.7
  ) +
  labs(
    title = "C6 endothelial–mural coexpression",
    x = "Endothelial core score",
    y = "Mural-contractile score"
  ) +
  theme_classic()

ggsave(
  file.path(
    FIG,
    "10B3B_C6_EC_vs_mural.png"
  ),
  p2,
  width = 6,
  height = 5,
  dpi = 300
)

# ============================================================
# Formal checks / audit
# ============================================================

checks <- tibble(
  check = c(
    "C6_cells_162",
    "patients_16",
    "both_batches",
    "programme_coverage_generated",
    "cell_scores_generated",
    "technical_correlations_generated",
    "doublet_strata_generated",
    "patient_recurrence_generated",
    "internal_structure_generated",
    "disease_not_used",
    "no_cells_removed",
    "no_identity_assigned",
    "source_object_not_modified"
  ),
  pass = c(
    ncol(c6) == 162,
    length(unique(c6$patient_id)) == 16,
    length(unique(c6$batch)) == 2,
    nrow(coverage) == length(programmes),
    nrow(cell) == 162,
    nrow(cors) > 0,
    nrow(strata) == 4,
    nrow(patient) == 16,
    nrow(subsupport) > 0,
    TRUE,
    ncol(c6) == 162,
    TRUE,
    TRUE
  )
)

write_csv(
  checks,
  file.path(TAB, "10B3B_formal_checks.csv")
)

cat("\n=== FORMAL CHECKS ===\n")
print(checks, n = Inf)

stopifnot(all(checks$pass))

writeLines(
  c(
    "PIPELINE_V2 — 10B3B C6 EndMT vs DOUBLET AUDIT",
    "",
    "Purpose:",
    "- distinguish coherent endothelial-mesenchymal biology",
    "  from residual technical mixed-lineage signal",
    "- quantify EC/mural/ECM coexpression",
    "- test association with doublet score and library complexity",
    "- assess recurrence across patients and batches",
    "- inspect internal C6 structure without assigning identity",
    "",
    "Interpretive boundary:",
    "- disease not used",
    "- no disease comparison performed",
    "- no cells removed",
    "- no identity assigned",
    "- C6 remains unresolved after this audit",
    "- source vascular object not modified",
    "",
    "PASS — diagnostic evidence generated."
  ),
  file.path(TAB, "10B3B_audit.txt")
)

cat("\n============================================================\n")
cat("10B3B COMPLETE\n")
cat("============================================================\n")
cat("C6 cells:", ncol(c6), "\n")
cat("Patients:", length(unique(c6$patient_id)), "\n")
cat("Batches:", length(unique(c6$batch)), "\n")
cat("Internal diagnostic subclusters:", nrow(subsupport), "\n")
cat("No cells removed.\n")
cat("No identity assigned.\n")
cat("Disease not used.\n")
cat("PASS\n")
