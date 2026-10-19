# =============================================================================
# PIPELINE_V2 — 05B6_cluster25_doublet_score_calibration.R
# Calibrate cluster-25 residual KRT-high cells against the ACTUAL scDblFinder
# score distributions from the 89,668 post-QC cells, including the 6,775 cells
# historically removed as doublets.
#
# AUDIT ONLY:
# - no cells removed
# - no identities changed
# - no PCA/Harmony/UMAP/clustering recalculated
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(ggplot2)
})

options(stringsAsFactors = FALSE)
set.seed(20260924)
source("PIPELINE_V2/config/paths.R")

cat("\n============================================================\n")
cat("PIPELINE_V2 — 05B6 DOUBLET SCORE CALIBRATION\n")
cat("============================================================\n\n")

# Historical A02 contains all post-QC cells with scDblFinder score/class attached
HIST_A02 <- "04_ANALYSIS/12_scrnaseq/DE_NOVO_16/01_QC/A02_16_singlet_objects.rds"
# If the singlet-only object cannot contain removed doublets, use post-QC input and
# reconstruct score/class from the current 04 script output is impossible. We
# therefore search first for an all-cell object carrying the score metadata.
POSTQC <- file.path(V2_ROOT,"objects","checkpoints","03_postQC_objects.rds")
SINGLETS <- file.path(V2_ROOT,"objects","checkpoints","04_singlet_objects.rds")
ATLAS <- file.path(V2_FROZEN,"global_atlas_v2_frozen.rds")
SOUPX <- file.path(V2_ROOT,"objects","ambient_corrected","global_atlas_v2_SoupX_expression.rds")

OUT <- file.path(V2_ROOT,"results","QC","ambient_RNA","05B6_doublet_calibration")
TD <- file.path(OUT,"tables")
FD <- file.path(OUT,"figures")
dir.create(TD,recursive=TRUE,showWarnings=FALSE)
dir.create(FD,recursive=TRUE,showWarnings=FALSE)

atlas <- readRDS(ATLAS)
soupx <- readRDS(SOUPX)
stopifnot(ncol(atlas)==82893, identical(colnames(atlas),colnames(soupx)))

# ---------------------------------------------------------------------------
# Recover score distributions.
# The V2 04 singlet checkpoint contains scores only for retained cells.
# To calibrate against removed doublets, attempt to recover scores/classes from
# historical artifacts. Fail safely rather than inventing them.
# ---------------------------------------------------------------------------
cat("STEP 1/5 — recovering scDblFinder score distributions\n")

sing <- readRDS(SINGLETS)
sing_df <- do.call(rbind,lapply(names(sing),function(pt){
  z <- sing[[pt]]
  stopifnot("doublet_score" %in% colnames(z@meta.data))
  data.frame(
    cell_id=colnames(z),
    patient_id=pt,
    doublet_score=as.numeric(z$doublet_score),
    class="retained_singlet",
    stringsAsFactors=FALSE
  )
}))
rownames(sing_df) <- NULL
stopifnot(nrow(sing_df)==82893)

# Try candidate all-cell files carrying historical doublet metadata.
candidates <- c(
  "04_ANALYSIS/12_scrnaseq/DE_NOVO_16/01_QC/A01_16_postQC_objects.rds",
  POSTQC
)

all_df <- NULL
source_used <- NA_character_

for(f in candidates){
  if(!file.exists(f)) next
  zz <- readRDS(f)
  has_score <- all(vapply(zz,function(z) "doublet_score" %in% colnames(z@meta.data),logical(1)))
  has_class <- all(vapply(zz,function(z) "doublet_class" %in% colnames(z@meta.data),logical(1)))
  if(has_score && has_class){
    all_df <- do.call(rbind,lapply(names(zz),function(pt){
      z <- zz[[pt]]
      data.frame(
        cell_id=colnames(z),patient_id=pt,
        doublet_score=as.numeric(z$doublet_score),
        class=as.character(z$doublet_class),
        stringsAsFactors=FALSE
      )
    }))
    source_used <- f
    break
  }
}

# If all-cell scores were not persisted, we can still calibrate KRT-high against
# all retained singlets. We explicitly report that removed-doublet calibration
# is unavailable rather than silently approximating it.
removed_available <- !is.null(all_df)

if(removed_available){
  cat("Recovered all-cell score/class metadata from:\n ",source_used,"\n")
  cat("Rows:",nrow(all_df),"\n")
  print(table(all_df$class))
} else {
  cat("WARNING: no persisted all-cell object with doublet_score + doublet_class found.\n")
  cat("Calibration against removed doublets will be marked unavailable.\n")
}

# ---------------------------------------------------------------------------
# Define cluster-25 KRT-high/low from corrected expression exactly as 05B5.
# ---------------------------------------------------------------------------
cat("\nSTEP 2/5 — defining cluster-25 KRT states\n")

DefaultAssay(soupx) <- "RNA"
cells25 <- colnames(soupx)[as.character(soupx$RNA_snn_res.0.6)=="25"]
KRT <- intersect(c("KRT5","KRT14","KRT15","KRT17","KRT1","KRT10"),rownames(soupx))
cc <- LayerData(soupx,assay="RNA",layer="counts")[KRT,cells25,drop=FALSE]
nk <- Matrix::colSums(cc>0)
khi <- nk>=3
stopifnot(length(cells25)==330,sum(khi)==147)

state <- data.frame(
  cell_id=cells25,
  state=ifelse(khi,"cluster25_KRT_high","cluster25_KRT_low"),
  stringsAsFactors=FALSE
)
state$doublet_score <- atlas$doublet_score[match(state$cell_id,colnames(atlas))]
stopifnot(!anyNA(state$doublet_score))

# ---------------------------------------------------------------------------
# Percentile calibration.
# ---------------------------------------------------------------------------
cat("\nSTEP 3/5 — percentile calibration\n")

sing_scores <- sing_df$doublet_score
sing_ecdf <- ecdf(sing_scores)

state$percentile_among_retained_singlets <- 100*sing_ecdf(state$doublet_score)

if(removed_available){
  # Normalize class spelling defensively.
  is_doublet <- tolower(all_df$class)=="doublet"
  removed_scores <- all_df$doublet_score[is_doublet]
  retained_scores_all <- all_df$doublet_score[!is_doublet]
  removed_ecdf <- ecdf(removed_scores)
  state$percentile_among_removed_doublets <- 100*removed_ecdf(state$doublet_score)
} else {
  removed_scores <- numeric(0)
  state$percentile_among_removed_doublets <- NA_real_
}

write.csv(state,file.path(TD,"05B6_cluster25_cell_score_calibration.csv"),row.names=FALSE)

qfun <- function(v){
  qs <- quantile(v,c(.01,.05,.10,.25,.50,.75,.90,.95,.99),na.rm=TRUE,names=FALSE)
  data.frame(q01=qs[1],q05=qs[2],q10=qs[3],q25=qs[4],median=qs[5],
             q75=qs[6],q90=qs[7],q95=qs[8],q99=qs[9])
}

dist_list <- list(
  all_retained_singlets=sing_scores,
  cluster25_KRT_low=state$doublet_score[state$state=="cluster25_KRT_low"],
  cluster25_KRT_high=state$doublet_score[state$state=="cluster25_KRT_high"]
)
if(removed_available) dist_list$removed_doublets <- removed_scores

quantiles <- do.call(rbind,lapply(names(dist_list),function(nm){
  cbind(group=nm,n=length(dist_list[[nm]]),qfun(dist_list[[nm]]))
}))
rownames(quantiles) <- NULL
write.csv(quantiles,file.path(TD,"05B6_doublet_score_quantiles.csv"),row.names=FALSE)

# Threshold-free overlap summaries.
summary_state <- do.call(rbind,lapply(split(state,state$state),function(z){
  data.frame(
    state=z$state[1],n=nrow(z),
    median_score=median(z$doublet_score),
    median_percentile_singlets=median(z$percentile_among_retained_singlets),
    pct_above_singlet_p90=100*mean(z$doublet_score>=quantile(sing_scores,.90)),
    pct_above_singlet_p95=100*mean(z$doublet_score>=quantile(sing_scores,.95)),
    pct_above_singlet_p99=100*mean(z$doublet_score>=quantile(sing_scores,.99)),
    median_percentile_removed=if(removed_available) median(z$percentile_among_removed_doublets) else NA
  )
}))
write.csv(summary_state,file.path(TD,"05B6_cluster25_state_calibration_summary.csv"),row.names=FALSE)

# ---------------------------------------------------------------------------
# Patient-level calibration.
# ---------------------------------------------------------------------------
cat("\nSTEP 4/5 — patient-level calibration\n")

state$patient_id <- atlas$patient_id[match(state$cell_id,colnames(atlas))]
pt <- do.call(rbind,lapply(split(state,state$patient_id),function(z){
  hi <- z$state=="cluster25_KRT_high"; lo <- !hi
  data.frame(
    patient_id=z$patient_id[1],n=nrow(z),n_high=sum(hi),n_low=sum(lo),
    median_score_high=if(any(hi)) median(z$doublet_score[hi]) else NA,
    median_score_low=if(any(lo)) median(z$doublet_score[lo]) else NA,
    pct_high_above_global_singlet_p95=if(any(hi)) 100*mean(z$doublet_score[hi]>=quantile(sing_scores,.95)) else NA
  )
}))
write.csv(pt,file.path(TD,"05B6_patient_calibration.csv"),row.names=FALSE)

# ---------------------------------------------------------------------------
# Figures.
# ---------------------------------------------------------------------------
cat("\nSTEP 5/5 — figures\n")

plot_df <- data.frame(
  score=c(sing_scores,
          state$doublet_score[state$state=="cluster25_KRT_low"],
          state$doublet_score[state$state=="cluster25_KRT_high"],
          if(removed_available) removed_scores else numeric(0)),
  group=c(rep("All retained singlets",length(sing_scores)),
          rep("Cluster 25 KRT-low",sum(state$state=="cluster25_KRT_low")),
          rep("Cluster 25 KRT-high",sum(state$state=="cluster25_KRT_high")),
          if(removed_available) rep("Removed doublets",length(removed_scores)) else character(0))
)

p1 <- ggplot(plot_df,aes(score,group))+
  geom_violin(scale="width",trim=TRUE,fill="grey90",color="grey30",linewidth=.3)+
  geom_boxplot(width=.12,outlier.shape=NA,fill="white",linewidth=.3)+
  theme_classic(base_size=10)+
  labs(x="scDblFinder score",y=NULL,title="Calibration of cluster 25 against scDblFinder distributions")
ggsave(file.path(FD,"05B6_doublet_score_distributions.pdf"),p1,width=8,height=4.5,bg="white")
ggsave(file.path(FD,"05B6_doublet_score_distributions.png"),p1,width=8,height=4.5,dpi=600,bg="white")

# ECDF gives a less density-dependent comparison.
p2 <- ggplot(plot_df,aes(score,color=group))+
  stat_ecdf(linewidth=.8)+theme_classic(base_size=10)+
  labs(x="scDblFinder score",y="Empirical cumulative fraction",color=NULL,
       title="scDblFinder score ECDF")
ggsave(file.path(FD,"05B6_doublet_score_ECDF.pdf"),p2,width=7,height=5,bg="white")
ggsave(file.path(FD,"05B6_doublet_score_ECDF.png"),p2,width=7,height=5,dpi=600,bg="white")

writeLines(c(
  "PIPELINE_V2 — 05B6 doublet score calibration",
  paste0("Generated: ",format(Sys.time(),"%Y-%m-%d %H:%M:%S")),
  paste0("Retained singlets: ",length(sing_scores)),
  paste0("Cluster25 KRT-low: ",sum(state$state=="cluster25_KRT_low")),
  paste0("Cluster25 KRT-high: ",sum(state$state=="cluster25_KRT_high")),
  paste0("Removed-doublet scores available: ",removed_available),
  if(removed_available) paste0("All-cell score source: ",source_used) else
    "Removed-doublet scores were not persisted in an accessible all-cell checkpoint; no approximation was made.",
  "No score threshold was invented.",
  "No cells removed and canonical atlas unchanged."
),file.path(TD,"05B6_audit.txt"))

cat("\n=== SCORE QUANTILES ===\n")
print(quantiles,row.names=FALSE)
cat("\n=== CLUSTER 25 CALIBRATION ===\n")
print(summary_state,row.names=FALSE)
cat("\nRemoved-doublet calibration available:",removed_available,"\n")
if(removed_available) cat("Removed doublets recovered:",length(removed_scores),"\n")
cat("\nPASS — calibration audit only; canonical atlas unchanged.\n")
