suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(readr)
})

cat("\n============================================================\n")
cat("PIPELINE_V2 — 08F IMMUNE UNIVERSE REFINEMENT\n")
cat("============================================================\n\n")

# ==========================================================================
# PURPOSE
#
# 1. Correct the cluster-level programme standardisation used in 08E.
# 2. Adjudicate obvious non-immune populations conservatively.
# 3. Resolve mixed cluster C9 at single-cell level.
# 4. Flag possible epithelial-immune doublets separately.
#
# IMPORTANT:
# - No cells are removed.
# - No L1 labels are modified.
# - No immune L2 labels are assigned.
# - No clustering/reduction is recalculated.
# - Output is an adjudication table only.
# ==========================================================================

IN <- "PIPELINE_V2/objects/annotated/immune_v2_08B_reclustered.rds"

OUT <- "PIPELINE_V2/results/annotation/08_IMMUNE_L2/08F_universe_refinement"
TDIR <- file.path(OUT,"tables")

dir.create(TDIR,recursive=TRUE,showWarnings=FALSE)

CC <- "immune_snn_res.0.6"

EXPECTED_CELLS <- 4762
EXPECTED_GENES <- 23412
EXPECTED_CLUSTERS <- 17

# ==========================================================================
# 1. LOAD
# ==========================================================================

cat("STEP 1/9 — loading 08B immune universe\n")

x <- readRDS(IN)
DefaultAssay(x) <- "RNA"

stopifnot(
  ncol(x)==EXPECTED_CELLS,
  nrow(x)==EXPECTED_GENES,
  CC %in% colnames(x@meta.data)
)

x$immune_cluster_06 <- as.character(
  x@meta.data[[CC]]
)

cluster_levels <- as.character(0:16)

stopifnot(
  all(cluster_levels %in% unique(x$immune_cluster_06)),
  length(unique(x$immune_cluster_06))==EXPECTED_CLUSTERS
)

cat("Cells:",ncol(x),"\n")
cat("Genes:",nrow(x),"\n")
cat("Patients:",length(unique(x$patient_id)),"\n")
cat("Clusters:",length(unique(x$immune_cluster_06)),"\n")

# ==========================================================================
# 2. DEFINE STRICT LINEAGE PROGRAMMES
#
# Deliberately restricted to relatively lineage-informative genes.
# State genes are not used to decide lineage.
# ==========================================================================

cat("\nSTEP 2/9 — defining strict lineage programmes\n")

programmes <- list(

  Leukocyte = c(
    "PTPRC","LST1","CD3D","CD3E",
    "TRAC","LYZ"
  ),

  T_cell = c(
    "CD3D","CD3E","TRAC","CD247"
  ),

  NK = c(
    "NKG7","GNLY","KLRD1",
    "NCR1","FGFBP2","PRF1"
  ),

  B_cell = c(
    "CD79A","CD79B","MS4A1",
    "PAX5","CD37"
  ),

  Plasma = c(
    "MZB1","JCHAIN","TNFRSF17",
    "SDC1","DERL3"
  ),

  Myeloid = c(
    "LST1","LYZ","CTSS",
    "FCER1G","TYROBP"
  ),

  Langerhans = c(
    "CD207","CD1A","FCER1A",
    "CD1E"
  ),

  Mast = c(
    "TPSAB1","TPSB2","CPA3",
    "KIT","MS4A2","HDC"
  ),

  Epithelial = c(
    "KRT5","KRT14","KRT15",
    "KRT17","KRT19",
    "COL17A1","EPCAM","KRTDAP"
  ),

  Neural_Schwann = c(
    "PMP2","MPZ","PLP1",
    "SOX10","S100B","CDH19"
  ),

  Fibroblast = c(
    "COL1A1","COL1A2",
    "DCN","LUM","PDGFRA"
  ),

  Endothelial = c(
    "PECAM1","VWF","EMCN",
    "KDR","CLDN5"
  )
)

coverage <- bind_rows(
  lapply(names(programmes),function(nm){

    requested <- programmes[[nm]]
    present <- intersect(requested,rownames(x))

    data.frame(
      programme=nm,
      n_present=length(present),
      n_requested=length(requested),
      coverage=length(present)/length(requested),
      missing=paste(
        setdiff(requested,present),
        collapse=";"
      ),
      stringsAsFactors=FALSE
    )
  })
)

cat("\n=== PROGRAMME COVERAGE ===\n")
print(
  as.data.frame(coverage),
  row.names=FALSE
)

write_csv(
  coverage,
  file.path(
    TDIR,
    "08F_programme_coverage.csv"
  )
)

stopifnot(
  all(coverage$n_present>0)
)

# ==========================================================================
# 3. EXPRESSION MATRIX
# ==========================================================================

cat("\nSTEP 3/9 — extracting expression matrix\n")

genes <- intersect(
  unique(unlist(programmes)),
  rownames(x)
)

dat <- LayerData(
  x,
  assay="RNA",
  layer="data"
)[genes,,drop=FALSE]

stopifnot(
  ncol(dat)==ncol(x),
  identical(
    colnames(dat),
    colnames(x)
  )
)

cat("Programme genes used:",nrow(dat),"\n")

# ==========================================================================
# 4. CORRECTED CLUSTER-LEVEL PROGRAMME SCORES
# ==========================================================================

cat("\nSTEP 4/9 — corrected cluster-level programme scoring\n")

cluster_score_rows <- list()

for(cl in cluster_levels){

  idx <- which(x$immune_cluster_06 == cl)
  stopifnot(length(idx) > 0)

  for(nm in names(programmes)){

    g <- intersect(programmes[[nm]], rownames(dat))

    vals <- Matrix::colMeans(
      dat[g,idx,drop=FALSE]
    )

    cluster_score_rows[[length(cluster_score_rows)+1]] <- data.frame(
      cluster=cl,
      programme=nm,
      n_cells=length(idx),
      mean_programme=mean(vals),
      median_programme=median(vals),
      pct_positive=100*mean(vals > 0),
      stringsAsFactors=FALSE
    )
  }
}

cluster_long <- bind_rows(cluster_score_rows)

# --------------------------------------------------------------------------
# Build explicit cluster x programme matrix
# --------------------------------------------------------------------------

mean_matrix <- xtabs(
  mean_programme ~ cluster + programme,
  data=cluster_long
)

mean_matrix <- as.matrix(mean_matrix)

mean_matrix <- mean_matrix[
  cluster_levels,
  names(programmes),
  drop=FALSE
]

stopifnot(
  nrow(mean_matrix) == EXPECTED_CLUSTERS,
  ncol(mean_matrix) == length(programmes)
)

# --------------------------------------------------------------------------
# CORRECT z-score:
# standardise EACH programme across the 17 clusters.
# --------------------------------------------------------------------------

z_matrix <- apply(
  mean_matrix,
  2,
  function(v){

    s <- sd(v)

    if(is.na(s) || s == 0){
      rep(0,length(v))
    } else {
      (v - mean(v)) / s
    }
  }
)

z_matrix <- as.matrix(z_matrix)

rownames(z_matrix) <- rownames(mean_matrix)
colnames(z_matrix) <- colnames(mean_matrix)

z_long <- as.data.frame(
  as.table(z_matrix),
  stringsAsFactors=FALSE
)

colnames(z_long) <- c(
  "cluster",
  "programme",
  "programme_z"
)

z_long$cluster <- as.character(z_long$cluster)
z_long$programme <- as.character(z_long$programme)
z_long$programme_z <- as.numeric(z_long$programme_z)

cluster_corrected <- cluster_long %>%
  left_join(
    z_long,
    by=c("cluster","programme")
  )

stopifnot(
  nrow(cluster_corrected) ==
    EXPECTED_CLUSTERS * length(programmes),
  !anyNA(cluster_corrected$programme_z)
)

write_csv(
  cluster_corrected,
  file.path(
    TDIR,
    "08F_corrected_cluster_programme_scores.csv"
  )
)

write_csv(
  data.frame(
    cluster=rownames(mean_matrix),
    mean_matrix,
    check.names=FALSE
  ),
  file.path(
    TDIR,
    "08F_cluster_programme_mean_matrix.csv"
  )
)

write_csv(
  data.frame(
    cluster=rownames(z_matrix),
    z_matrix,
    check.names=FALSE
  ),
  file.path(
    TDIR,
    "08F_cluster_programme_z_matrix.csv"
  )
)

# ==========================================================================
# 5. CORRECTED PROGRAMME RANKING
# ==========================================================================

cat("\nSTEP 5/9 — corrected programme ranking\n")

immune_programmes <- c(
  "Leukocyte",
  "T_cell",
  "NK",
  "B_cell",
  "Plasma",
  "Myeloid",
  "Langerhans",
  "Mast"
)

nonimmune_programmes <- c(
  "Epithelial",
  "Neural_Schwann",
  "Fibroblast",
  "Endothelial"
)

immune_rank <- cluster_corrected %>%
  filter(programme %in% immune_programmes) %>%
  group_by(cluster) %>%
  arrange(desc(programme_z), .by_group=TRUE) %>%
  mutate(immune_rank=row_number()) %>%
  ungroup()

nonimmune_rank <- cluster_corrected %>%
  filter(programme %in% nonimmune_programmes) %>%
  group_by(cluster) %>%
  arrange(desc(programme_z), .by_group=TRUE) %>%
  mutate(nonimmune_rank=row_number()) %>%
  ungroup()

top_immune <- immune_rank %>%
  filter(immune_rank == 1) %>%
  transmute(
    cluster,
    top_immune_programme=programme,
    top_immune_z=programme_z,
    top_immune_mean=mean_programme,
    top_immune_pct=pct_positive
  )

top_nonimmune <- nonimmune_rank %>%
  filter(nonimmune_rank == 1) %>%
  transmute(
    cluster,
    top_nonimmune_programme=programme,
    top_nonimmune_z=programme_z,
    top_nonimmune_mean=mean_programme,
    top_nonimmune_pct=pct_positive
  )

leukocyte <- cluster_corrected %>%
  filter(programme == "Leukocyte") %>%
  transmute(
    cluster,
    leukocyte_z=programme_z,
    leukocyte_mean=mean_programme,
    leukocyte_pct=pct_positive
  )

corrected_summary <- top_immune %>%
  left_join(top_nonimmune,by="cluster") %>%
  left_join(leukocyte,by="cluster") %>%
  mutate(
    immune_minus_nonimmune_z =
      top_immune_z - top_nonimmune_z,
    cluster_numeric=as.numeric(cluster)
  ) %>%
  arrange(cluster_numeric) %>%
  select(-cluster_numeric)

cat("\n=== CORRECTED CLUSTER PROGRAMME SUMMARY ===\n")

print(
  as.data.frame(corrected_summary),
  row.names=FALSE
)

write_csv(
  corrected_summary,
  file.path(
    TDIR,
    "08F_corrected_cluster_summary.csv"
  )
)

# ==========================================================================
# 6. Z-SCORE SANITY CHECK
# ==========================================================================

cat("\nSTEP 6/9 — z-score sanity check\n")

z_variation <- cluster_corrected %>%
  group_by(programme) %>%
  summarise(
    n_unique_z=n_distinct(round(programme_z,8)),
    mean_z=mean(programme_z),
    sd_z=sd(programme_z),
    min_z=min(programme_z),
    max_z=max(programme_z),
    .groups="drop"
  )

cat("\n=== Z-SCORE SANITY CHECK ===\n")

print(
  as.data.frame(z_variation),
  row.names=FALSE
)

write_csv(
  z_variation,
  file.path(
    TDIR,
    "08F_zscore_sanity_check.csv"
  )
)

stopifnot(
  all(z_variation$n_unique_z > 1),
  all(abs(z_variation$mean_z) < 1e-10),
  all(abs(z_variation$sd_z - 1) < 1e-8)
)

# --------------------------------------------------------------------------
# Focused corrected report
# --------------------------------------------------------------------------

focus_clusters <- c(
  "3","7","8","9",
  "10","11","12",
  "13","14","15","16"
)

cat("\n============================================================\n")
cat("FOCUSED CORRECTED PROGRAMME REPORT\n")
cat("============================================================\n")

for(cl in focus_clusters){

  cat("\n------------------------------------------------------------\n")
  cat("CLUSTER",cl,"\n")

  z <- cluster_corrected %>%
    filter(cluster == cl) %>%
    arrange(desc(programme_z)) %>%
    select(
      programme,
      programme_z,
      mean_programme,
      pct_positive
    )

  print(
    as.data.frame(z),
    row.names=FALSE
  )
}


# ==========================================================================
# 7. CELL-LEVEL LINEAGE EVIDENCE
# ==========================================================================

cat("\nSTEP 7/9 — cell-level lineage evidence for C9/C12\n")

audit_clusters <- c("9","12")

audit_cells <- colnames(x)[
  x$immune_cluster_06 %in% audit_clusters
]

audit_dat <- dat[,audit_cells,drop=FALSE]

# --------------------------------------------------------------------------
# Helper: mean expression of a gene set per cell
# --------------------------------------------------------------------------

cell_score <- function(mat, genes){

  g <- intersect(genes,rownames(mat))

  stopifnot(length(g)>0)

  Matrix::colMeans(
    mat[g,,drop=FALSE]
  )
}

# Strict cell-level scores

score_leukocyte <- cell_score(
  audit_dat,
  c("PTPRC","LST1","CD3D","CD3E","TRAC","LYZ")
)

score_T <- cell_score(
  audit_dat,
  c("CD3D","CD3E","TRAC","CD247")
)

score_myeloid <- cell_score(
  audit_dat,
  c("LST1","LYZ","CTSS","FCER1G","TYROBP")
)

score_epithelial <- cell_score(
  audit_dat,
  c(
    "KRT5","KRT14","KRT15",
    "KRT17","KRT19",
    "COL17A1","EPCAM","KRTDAP"
  )
)

score_schwann <- cell_score(
  audit_dat,
  c(
    "PMP2","MPZ","PLP1",
    "SOX10","S100B","CDH19"
  )
)

# --------------------------------------------------------------------------
# Gene-level binary support
# --------------------------------------------------------------------------

expr_gene <- function(g){

  if(!g %in% rownames(audit_dat)){
    return(
      rep(0,ncol(audit_dat))
    )
  }

  as.numeric(
    audit_dat[g,] > 0
  )
}

immune_gene_count <-
  expr_gene("PTPRC") +
  expr_gene("LST1") +
  expr_gene("CD3D") +
  expr_gene("CD3E") +
  expr_gene("TRAC") +
  expr_gene("LYZ")

epithelial_gene_count <-
  expr_gene("KRT5") +
  expr_gene("KRT14") +
  expr_gene("KRT15") +
  expr_gene("KRT17") +
  expr_gene("KRT19") +
  expr_gene("COL17A1") +
  expr_gene("EPCAM") +
  expr_gene("KRTDAP")

schwann_gene_count <-
  expr_gene("PMP2") +
  expr_gene("MPZ") +
  expr_gene("PLP1") +
  expr_gene("SOX10") +
  expr_gene("S100B") +
  expr_gene("CDH19")

# ==========================================================================
# 8. CONSERVATIVE CELL-LEVEL ADJUDICATION
#
# Rules intentionally require coherent multi-gene evidence.
#
# REVIEW_MIXED_DOUBLET is NOT equivalent to removal.
# It means that immune and epithelial evidence coexist strongly enough that
# the cell requires explicit downstream doublet review.
# ==========================================================================

md <- x@meta.data[audit_cells,,drop=FALSE]

cell_audit <- data.frame(
  cell=audit_cells,
  cluster=as.character(md$immune_cluster_06),
  source_L1=as.character(md$L1),
  patient_id=as.character(md$patient_id),
  disease=as.character(md$disease),
  batch=as.character(md$batch),
  doublet_score=as.numeric(md$doublet_score),

  leukocyte_score=as.numeric(score_leukocyte),
  T_score=as.numeric(score_T),
  myeloid_score=as.numeric(score_myeloid),
  epithelial_score=as.numeric(score_epithelial),
  schwann_score=as.numeric(score_schwann),

  immune_gene_count=as.integer(immune_gene_count),
  epithelial_gene_count=as.integer(epithelial_gene_count),
  schwann_gene_count=as.integer(schwann_gene_count),

  stringsAsFactors=FALSE
)

# --------------------------------------------------------------------------
# Evidence flags
# --------------------------------------------------------------------------

cell_audit <- cell_audit %>%
  mutate(

    immune_supported =
      immune_gene_count >= 2 &
      leukocyte_score > 0,

    epithelial_supported =
      epithelial_gene_count >= 2 &
      epithelial_score > 0,

    schwann_supported =
      schwann_gene_count >= 3 &
      schwann_score > 0,

    mixed_immune_epithelial =
      immune_supported &
      epithelial_supported,

    lineage_winner = case_when(

      schwann_supported &
        schwann_score > leukocyte_score &
        schwann_score > epithelial_score
        ~ "Neural_Schwann",

      immune_supported &
        !epithelial_supported
        ~ "Immune",

      epithelial_supported &
        !immune_supported
        ~ "Epithelial",

      mixed_immune_epithelial
        ~ "Mixed_immune_epithelial",

      TRUE
        ~ "Uncertain"
    )
  )

# --------------------------------------------------------------------------
# Conservative disposition
#
# C12:
#   coherent Schwann cells return to Neural/Schwann.
#   anything else remains reviewable.
#
# C9:
#   immune-only cells stay immune.
#   epithelial-only cells return epithelial.
#   mixed cells are NOT automatically deleted.
# --------------------------------------------------------------------------

cell_audit <- cell_audit %>%
  mutate(

    refined_disposition = case_when(

      cluster=="12" &
        lineage_winner=="Neural_Schwann"
        ~ "RETURN_NEURAL_SCHWANN",

      cluster=="12"
        ~ "REVIEW_UNCERTAIN",

      cluster=="9" &
        lineage_winner=="Immune"
        ~ "KEEP_IMMUNE",

      cluster=="9" &
        lineage_winner=="Epithelial"
        ~ "RETURN_EPITHELIAL",

      cluster=="9" &
        lineage_winner=="Mixed_immune_epithelial"
        ~ "REVIEW_MIXED_DOUBLET",

      cluster=="9" &
        lineage_winner=="Neural_Schwann"
        ~ "RETURN_NEURAL_SCHWANN",

      TRUE
        ~ "REVIEW_UNCERTAIN"
    )
  )

write_csv(
  cell_audit,
  file.path(
    TDIR,
    "08F_C9_C12_cell_level_adjudication.csv"
  )
)

# ==========================================================================
# CELL-LEVEL SUMMARY
# ==========================================================================

cell_summary <- cell_audit %>%
  count(
    cluster,
    refined_disposition,
    name="n_cells"
  ) %>%
  group_by(cluster) %>%
  mutate(
    fraction=n_cells/sum(n_cells)
  ) %>%
  ungroup()

cat("\n=== C9/C12 REFINED DISPOSITION ===\n")

print(
  as.data.frame(cell_summary),
  row.names=FALSE
)

write_csv(
  cell_summary,
  file.path(
    TDIR,
    "08F_C9_C12_disposition_summary.csv"
  )
)

# ==========================================================================
# PATIENT SUPPORT OF EACH DISPOSITION
# ==========================================================================

patient_disposition <- cell_audit %>%
  count(
    cluster,
    refined_disposition,
    patient_id,
    disease,
    batch,
    name="n_cells"
  ) %>%
  arrange(
    cluster,
    refined_disposition,
    disease,
    desc(n_cells)
  )

write_csv(
  patient_disposition,
  file.path(
    TDIR,
    "08F_C9_C12_disposition_by_patient.csv"
  )
)

# ==========================================================================
# MIXED-CELL DOUBLET AUDIT
# ==========================================================================

doublet_audit <- cell_audit %>%
  group_by(
    cluster,
    refined_disposition
  ) %>%
  summarise(
    n_cells=n(),
    median_doublet_score=median(
      doublet_score,
      na.rm=TRUE
    ),
    q25_doublet_score=quantile(
      doublet_score,
      0.25,
      na.rm=TRUE
    ),
    q75_doublet_score=quantile(
      doublet_score,
      0.75,
      na.rm=TRUE
    ),
    median_immune_genes=median(
      immune_gene_count
    ),
    median_epithelial_genes=median(
      epithelial_gene_count
    ),
    median_schwann_genes=median(
      schwann_gene_count
    ),
    .groups="drop"
  )

cat("\n=== C9/C12 DOUBLET / LINEAGE QC ===\n")

print(
  as.data.frame(doublet_audit),
  row.names=FALSE
)

write_csv(
  doublet_audit,
  file.path(
    TDIR,
    "08F_C9_C12_doublet_lineage_QC.csv"
  )
)

# ==========================================================================
# C9 DETAILED REPORT
# ==========================================================================

cat("\n============================================================\n")
cat("CLUSTER 9 — CELL-LEVEL ADJUDICATION\n")
cat("============================================================\n")

c9 <- cell_audit %>%
  filter(cluster=="9")

cat("Cells:",nrow(c9),"\n\n")

print(
  as.data.frame(
    c9 %>%
      count(
        refined_disposition,
        name="n_cells"
      ) %>%
      mutate(
        fraction=n_cells/sum(n_cells)
      )
  ),
  row.names=FALSE
)

cat("\nC9 patient support by disposition:\n")

print(
  as.data.frame(
    c9 %>%
      count(
        refined_disposition,
        patient_id,
        disease,
        name="n_cells"
      ) %>%
      arrange(
        refined_disposition,
        desc(n_cells)
      )
  ),
  row.names=FALSE
)

# ==========================================================================
# C12 DETAILED REPORT
# ==========================================================================

cat("\n============================================================\n")
cat("CLUSTER 12 — CELL-LEVEL ADJUDICATION\n")
cat("============================================================\n")

c12 <- cell_audit %>%
  filter(cluster=="12")

cat("Cells:",nrow(c12),"\n\n")

print(
  as.data.frame(
    c12 %>%
      count(
        refined_disposition,
        name="n_cells"
      ) %>%
      mutate(
        fraction=n_cells/sum(n_cells)
      )
  ),
  row.names=FALSE
)


# ==========================================================================
# 9. MASTER REFINEMENT TABLE + FORMAL INVARIANTS
# ==========================================================================

cat("\nSTEP 8/9 — master refinement table\n")

all_cells <- colnames(x)

master <- data.frame(
  cell=all_cells,
  immune_cluster_06=as.character(x$immune_cluster_06),
  source_L1=as.character(x$L1),
  patient_id=as.character(x$patient_id),
  disease=as.character(x$disease),
  batch=as.character(x$batch),
  stringsAsFactors=FALSE
)

# Default: retain cells in immune universe.
master$proposed_disposition <- "KEEP_IMMUNE"

# C9/C12 receive the cell-level adjudication generated above.
disp_lookup <- setNames(
  cell_audit$refined_disposition,
  cell_audit$cell
)

idx_audit <- master$cell %in% names(disp_lookup)

master$proposed_disposition[idx_audit] <-
  unname(
    disp_lookup[
      master$cell[idx_audit]
    ]
  )

stopifnot(
  nrow(master)==EXPECTED_CELLS,
  !anyNA(master$proposed_disposition),
  length(unique(master$cell))==EXPECTED_CELLS
)

write_csv(
  master,
  file.path(
    TDIR,
    "08F_master_cell_disposition.csv"
  )
)

# --------------------------------------------------------------------------
# Proposed disposition summary
# --------------------------------------------------------------------------

master_summary <- master %>%
  count(
    proposed_disposition,
    name="n_cells"
  ) %>%
  mutate(
    pct=100*n_cells/sum(n_cells)
  ) %>%
  arrange(desc(n_cells))

cat("\n=== PROPOSED IMMUNE-UNIVERSE DISPOSITION ===\n")

print(
  as.data.frame(master_summary),
  row.names=FALSE
)

write_csv(
  master_summary,
  file.path(
    TDIR,
    "08F_master_disposition_summary.csv"
  )
)

# --------------------------------------------------------------------------
# Proposed immune universe only
#
# REVIEW categories are NOT automatically retained in the final universe.
# They remain unresolved until explicit adjudication after inspection.
# --------------------------------------------------------------------------

n_keep <- sum(
  master$proposed_disposition=="KEEP_IMMUNE"
)

n_return_epi <- sum(
  master$proposed_disposition=="RETURN_EPITHELIAL"
)

n_return_schwann <- sum(
  master$proposed_disposition=="RETURN_NEURAL_SCHWANN"
)

n_review_mixed <- sum(
  master$proposed_disposition=="REVIEW_MIXED_DOUBLET"
)

n_review_uncertain <- sum(
  master$proposed_disposition=="REVIEW_UNCERTAIN"
)

stopifnot(
  n_keep +
    n_return_epi +
    n_return_schwann +
    n_review_mixed +
    n_review_uncertain ==
    EXPECTED_CELLS
)

# --------------------------------------------------------------------------
# Patient support among cells provisionally retained as immune
# --------------------------------------------------------------------------

keep_support <- master %>%
  filter(
    proposed_disposition=="KEEP_IMMUNE"
  ) %>%
  count(
    patient_id,
    disease,
    batch,
    name="n_cells"
  ) %>%
  arrange(
    disease,
    desc(n_cells)
  )

write_csv(
  keep_support,
  file.path(
    TDIR,
    "08F_KEEP_IMMUNE_by_patient.csv"
  )
)

cat("\n=== PROVISIONALLY RETAINED IMMUNE CELLS BY PATIENT ===\n")

print(
  as.data.frame(keep_support),
  row.names=FALSE
)

# ==========================================================================
# FORMAL CHECKS
# ==========================================================================

cat("\nSTEP 9/9 — formal checks and audit trail\n")

checks <- data.frame(
  check=c(
    "input_cells_4762",
    "input_genes_23412",
    "working_clusters_17",
    "master_rows_4762",
    "master_unique_cells_4762",
    "all_cells_have_disposition",
    "corrected_zscores_vary_by_cluster",
    "corrected_zscores_mean_zero",
    "corrected_zscores_sd_one",
    "C9_cell_audit_complete",
    "C12_cell_audit_complete",
    "C9_plus_C12_audit_256",
    "disposition_counts_sum_4762",
    "no_global_L1_modified",
    "no_cells_removed_in_08F",
    "no_frozen_object_written"
  ),
  pass=c(
    ncol(x)==EXPECTED_CELLS,
    nrow(x)==EXPECTED_GENES,
    length(cluster_levels)==EXPECTED_CLUSTERS,
    nrow(master)==EXPECTED_CELLS,
    length(unique(master$cell))==EXPECTED_CELLS,
    !anyNA(master$proposed_disposition),
    all(z_variation$n_unique_z>1),
    all(abs(z_variation$mean_z)<1e-10),
    all(abs(z_variation$sd_z-1)<1e-8),
    sum(cell_audit$cluster=="9")==183,
    sum(cell_audit$cluster=="12")==73,
    nrow(cell_audit)==256,
    sum(master_summary$n_cells)==EXPECTED_CELLS,
    TRUE,
    TRUE,
    TRUE
  ),
  stringsAsFactors=FALSE
)

cat("\n=== FORMAL INVARIANTS ===\n")

print(
  checks,
  row.names=FALSE
)

write_csv(
  checks,
  file.path(
    TDIR,
    "08F_formal_checks.csv"
  )
)

stopifnot(
  all(checks$pass)
)

# ==========================================================================
# DECISION REPORT
# ==========================================================================

decision_report <- c(
  "PIPELINE_V2 — 08F immune universe refinement",
  paste0("Generated: ",Sys.time()),
  "Input: immune_v2_08B_reclustered.rds",
  "Working resolution: 0.6",
  paste0("Input cells: ",ncol(x)),
  paste0("Input genes: ",nrow(x)),
  paste0("Input clusters: ",length(cluster_levels)),
  "",
  "PURPOSE",
  "Correct the programme-standardisation defect identified in 08E and perform conservative cell-level adjudication of clusters 9 and 12.",
  "",
  "CORRECTED PROGRAMME SCORING",
  "Programme means were calculated independently for each cluster.",
  "Each programme was then standardized across the 17 immune working clusters.",
  "Disease was not used for lineage assignment.",
  "",
  "C9/C12 REFINEMENT",
  paste0("KEEP_IMMUNE: ",n_keep),
  paste0("RETURN_EPITHELIAL: ",n_return_epi),
  paste0("RETURN_NEURAL_SCHWANN: ",n_return_schwann),
  paste0("REVIEW_MIXED_DOUBLET: ",n_review_mixed),
  paste0("REVIEW_UNCERTAIN: ",n_review_uncertain),
  "",
  "IMPORTANT",
  "RETURN categories are proposed lineage corrections, not deletions.",
  "REVIEW_MIXED_DOUBLET is a review flag and does not itself imply cell removal.",
  "REVIEW_UNCERTAIN is unresolved and requires explicit adjudication.",
  "No cell was removed by this script.",
  "No cell was reassigned in the frozen global atlas.",
  "No frozen immune object was written.",
  "No L1 metadata were modified.",
  "",
  "NEXT STEP",
  "Inspect 08F outputs and explicitly adjudicate C9/C12 before rebuilding and freezing the refined immune universe."
)

writeLines(
  decision_report,
  file.path(
    TDIR,
    "08F_decision_report.txt"
  )
)

cat("\n============================================================\n")
cat("08F COMPLETE\n")
cat("============================================================\n")
cat("Input cells:",ncol(x),"\n")
cat("KEEP_IMMUNE:",n_keep,"\n")
cat("RETURN_EPITHELIAL:",n_return_epi,"\n")
cat("RETURN_NEURAL_SCHWANN:",n_return_schwann,"\n")
cat("REVIEW_MIXED_DOUBLET:",n_review_mixed,"\n")
cat("REVIEW_UNCERTAIN:",n_review_uncertain,"\n")
cat("No cells removed.\n")
cat("No global L1 labels modified.\n")
cat("No frozen immune object written.\n")
cat("PASS — refinement evidence complete; explicit adjudication required.\n")

