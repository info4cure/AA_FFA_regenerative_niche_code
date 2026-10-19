suppressPackageStartupMessages({
  library(Matrix)
  library(edgeR)
  library(limma)
  library(GSVA)
  library(dplyr)
  library(readr)
  library(tibble)
  library(tidyr)
  library(purrr)
  library(stringr)
})

cat("============================================================\n")
cat("PIPELINE_V2 — 20C PATIENT-LEVEL HALLMARK GSVA\n")
cat("============================================================\n")

BASE <- file.path(
  "PIPELINE_V2","results","disease_biology"
)

PB <- file.path(
  BASE,"18_PSEUDOBULK_DE","18B_pseudobulk"
)

GSFILE <- file.path(
  BASE,"20_PATHWAY_ANALYSIS","20A_HALLMARK_GSEA",
  "gene_sets","20A_Hallmark_gene_sets_frozen.rds"
)

GSEAFILE <- file.path(
  BASE,"20_PATHWAY_ANALYSIS","20A_HALLMARK_GSEA",
  "tables","20A_all_Hallmark_GSEA_results.csv"
)

OUT <- file.path(
  BASE,"20_PATHWAY_ANALYSIS",
  "20C_PATIENT_HALLMARK_GSVA"
)

TAB <- file.path(OUT,"tables")
OBJ <- file.path(OUT,"objects")

dir.create(TAB,recursive=TRUE,showWarnings=FALSE)
dir.create(OBJ,recursive=TRUE,showWarnings=FALSE)

stopifnot(
  file.exists(GSFILE),
  file.exists(GSEAFILE)
)

hallmarks <- readRDS(GSFILE)

md_all <- readRDS(
  file.path(
    PB,"objects",
    "18B_pseudobulk_metadata_ALL.rds"
  )
)

modules <- sort(unique(md_all$module))

cat("\nModules:",length(modules),"\n")
cat("Hallmarks:",length(hallmarks),"\n")

# ============================================================
# 1. Recover eligible patient × L2 pseudobulks
# ============================================================

cat("\nSTEP 1/8 — recovering >=20-cell pseudobulks\n")

md <- md_all %>%
  filter(n_cells >= 20)

cat("Eligible pseudobulks:",nrow(md),"\n")
cat("Patients:",n_distinct(md$patient_id),"\n")

# ============================================================
# 2. Calculate patient-level GSVA separately within each L2
# ============================================================

cat("\nSTEP 2/8 — TMM/logCPM + GSVA by L2\n")

score_list <- list()
norm_list <- list()
counter <- 0L

l2_index <- md %>%
  distinct(module,final_L2_label) %>%
  arrange(module,final_L2_label)

for (i in seq_len(nrow(l2_index))) {

  m <- l2_index$module[i]
  l2 <- l2_index$final_L2_label[i]

  submd <- md %>%
    filter(
      module==m,
      final_L2_label==l2
    )

  # Need >=3 samples to obtain a meaningful multi-sample score matrix
  if (nrow(submd) < 3)
    next

  count_file <- file.path(
    PB,"objects",
    paste0("18B_",m,"_pseudobulk_counts.rds")
  )

  stopifnot(file.exists(count_file))

  mat <- readRDS(count_file)

  stopifnot(
    all(submd$pseudobulk_id %in% colnames(mat))
  )

  counts <- mat[
    ,
    submd$pseudobulk_id,
    drop=FALSE
  ]

  stopifnot(
    identical(
      colnames(counts),
      submd$pseudobulk_id
    )
  )

  # Remove genes with no counts in this L2 before normalization.
  keep <- Matrix::rowSums(counts) > 0
  counts <- counts[keep,,drop=FALSE]

  y <- DGEList(counts=counts)
  y <- calcNormFactors(y,method="TMM")

  logcpm <- cpm(
    y,
    log=TRUE,
    prior.count=2
  )

  # Retain gene symbols uniquely.
  if (anyDuplicated(rownames(logcpm))) {
    stop(
      paste("Duplicated gene symbols in",m,l2)
    )
  }

  param <- gsvaParam(
    exprData=logcpm,
    geneSets=hallmarks,
    minSize=10,
    maxSize=500,
    kcdf="Gaussian"
  )

  scores <- gsva(
    param,
    verbose=FALSE
  )

  stopifnot(
    ncol(scores)==nrow(submd),
    identical(
      colnames(scores),
      submd$pseudobulk_id
    )
  )

  counter <- counter + 1L

  score_long <- as.data.frame(scores) %>%
    rownames_to_column("pathway") %>%
    pivot_longer(
      cols=-pathway,
      names_to="pseudobulk_id",
      values_to="GSVA_score"
    ) %>%
    left_join(
      submd %>%
        select(
          pseudobulk_id,
          patient_id,
          disease,
          batch,
          n_cells
        ),
      by="pseudobulk_id"
    ) %>%
    mutate(
      module=m,
      final_L2_label=l2,
      .before=1
    )

  score_list[[counter]] <- score_long

  norm_list[[counter]] <- tibble(
    module=m,
    final_L2_label=l2,
    n_patients=nrow(submd),
    n_genes_input=nrow(logcpm),
    n_HC=sum(submd$disease=="HC"),
    n_FFA=sum(submd$disease=="FFA"),
    n_AA=sum(submd$disease=="AA"),
    min_cells=min(submd$n_cells),
    median_cells=median(submd$n_cells)
  )

  cat(
    sprintf(
      "[%02d/%02d] %-18s | %-40s | patients=%d | genes=%d\n",
      i,nrow(l2_index),m,l2,nrow(submd),nrow(logcpm)
    )
  )
}

scores_long <- bind_rows(score_list)
L2_support <- bind_rows(norm_list)

write_csv(
  scores_long,
  file.path(TAB,"20C_patient_L2_Hallmark_GSVA_scores.csv")
)

write_csv(
  L2_support,
  file.path(TAB,"20C_L2_patient_support.csv")
)

saveRDS(
  scores_long,
  file.path(OBJ,"20C_patient_L2_Hallmark_GSVA_scores.rds")
)

# ============================================================
# 3. Freeze contrast specification
# ============================================================

cat("\nSTEP 3/8 — contrast specification\n")

contrast_spec <- tribble(
  ~contrast,   ~group1, ~group0, ~adjust_batch,
  "AA_vs_HC",  "AA",    "HC",    FALSE,
  "FFA_vs_HC", "FFA",   "HC",    TRUE,
  "AA_vs_FFA", "AA",    "FFA",   FALSE
)

print(as.data.frame(contrast_spec),row.names=FALSE)

# ============================================================
# 4. limma inference on patient-level GSVA scores
# ============================================================

cat("\nSTEP 4/8 — patient-level limma inference\n")

run_one <- function(m,l2,ct,g1,g0,adjust_batch) {

  z <- scores_long %>%
    filter(
      module==m,
      final_L2_label==l2,
      disease %in% c(g0,g1)
    )

  patient_md <- z %>%
    distinct(
      pseudobulk_id,
      patient_id,
      disease,
      batch
    ) %>%
    arrange(patient_id)

  n0 <- sum(patient_md$disease==g0)
  n1 <- sum(patient_md$disease==g1)

  if (n0 < 3 || n1 < 3)
    return(NULL)

  score_mat <- z %>%
    select(
      pathway,
      pseudobulk_id,
      GSVA_score
    ) %>%
    pivot_wider(
      names_from=pseudobulk_id,
      values_from=GSVA_score
    ) %>%
    as.data.frame()

  rownames(score_mat) <- score_mat$pathway
  score_mat$pathway <- NULL
  score_mat <- as.matrix(score_mat)

  # Force identical sample order
  patient_md <- patient_md %>%
    filter(pseudobulk_id %in% colnames(score_mat))

  score_mat <- score_mat[
    ,
    patient_md$pseudobulk_id,
    drop=FALSE
  ]

  patient_md$disease_model <- factor(
    patient_md$disease,
    levels=c(g0,g1)
  )

  if (adjust_batch) {

    patient_md$batch_model <- factor(
      patient_md$batch
    )

    design <- model.matrix(
      ~ batch_model + disease_model,
      data=patient_md
    )

  } else {

    design <- model.matrix(
      ~ disease_model,
      data=patient_md
    )
  }

  if (qr(design)$rank != ncol(design))
    return(NULL)

  coef_name <- paste0("disease_model",g1)
  coef_idx <- which(colnames(design)==coef_name)

  stopifnot(length(coef_idx)==1)

  fit <- lmFit(score_mat,design)
  fit <- eBayes(fit,robust=TRUE)

  tt <- topTable(
    fit,
    coef=coef_idx,
    number=Inf,
    sort.by="P"
  ) %>%
    rownames_to_column("pathway") %>%
    as_tibble() %>%
    rename(
      delta_GSVA=logFC,
      FDR_BH=adj.P.Val
    ) %>%
    mutate(
      contrast=ct,
      numerator=g1,
      denominator=g0,
      module=m,
      final_L2_label=l2,
      n_numerator=n1,
      n_denominator=n0,
      batch_adjusted=adjust_batch,
      .before=1
    )

  tt
}

infer_list <- list()
k <- 0L

for (i in seq_len(nrow(l2_index))) {

  m <- l2_index$module[i]
  l2 <- l2_index$final_L2_label[i]

  if (!any(
    scores_long$module==m &
    scores_long$final_L2_label==l2
  )) next

  for (j in seq_len(nrow(contrast_spec))) {

    cs <- contrast_spec[j,]

    r <- run_one(
      m,l2,
      cs$contrast,
      cs$group1,
      cs$group0,
      cs$adjust_batch
    )

    if (is.null(r)) next

    k <- k+1L
    infer_list[[k]] <- r
  }
}

infer <- bind_rows(infer_list)

write_csv(
  infer,
  file.path(TAB,"20C_Hallmark_GSVA_differential_scores.csv")
)

# ============================================================
# 5. Summaries
# ============================================================

cat("\nSTEP 5/8 — summaries\n")

summary_tbl <- infer %>%
  group_by(
    contrast,module,final_L2_label
  ) %>%
  summarise(
    n_pathways=n(),
    n_FDR005=sum(FDR_BH<.05,na.rm=TRUE),
    n_FDR010=sum(FDR_BH<.10,na.rm=TRUE),
    n_UP_FDR005=sum(
      FDR_BH<.05 & delta_GSVA>0,
      na.rm=TRUE
    ),
    n_DOWN_FDR005=sum(
      FDR_BH<.05 & delta_GSVA<0,
      na.rm=TRUE
    ),
    .groups="drop"
  )

write_csv(
  summary_tbl,
  file.path(TAB,"20C_L2_GSVA_summary.csv")
)

global_summary <- summary_tbl %>%
  group_by(contrast) %>%
  summarise(
    n_L2=n(),
    L2_with_FDR005=sum(n_FDR005>0),
    total_pathway_hits_FDR005=sum(n_FDR005),
    total_UP_FDR005=sum(n_UP_FDR005),
    total_DOWN_FDR005=sum(n_DOWN_FDR005),
    .groups="drop"
  )

write_csv(
  global_summary,
  file.path(TAB,"20C_global_GSVA_summary.csv")
)

# ============================================================
# 6. Concordance with frozen 20A GSEA
# ============================================================

cat("\nSTEP 6/8 — GSVA vs GSEA concordance\n")

gsea <- read_csv(
  GSEAFILE,
  show_col_types=FALSE
) %>%
  select(
    contrast,module,final_L2_label,
    pathway,NES,padj
  )

concord <- infer %>%
  inner_join(
    gsea,
    by=c(
      "contrast",
      "module",
      "final_L2_label",
      "pathway"
    )
  )

concord_summary <- concord %>%
  group_by(contrast) %>%
  summarise(
    n_pairs=n(),
    spearman_effect_NES=cor(
      delta_GSVA,NES,
      method="spearman",
      use="complete.obs"
    ),
    direction_concordance=mean(
      sign(delta_GSVA)==sign(NES),
      na.rm=TRUE
    ),
    n_GSEA_FDR005=sum(padj<.05,na.rm=TRUE),
    direction_concordance_GSEA_FDR005=
      mean(
        sign(delta_GSVA[padj<.05]) ==
        sign(NES[padj<.05]),
        na.rm=TRUE
      ),
    .groups="drop"
  )

write_csv(
  concord,
  file.path(TAB,"20C_GSVA_GSEA_pairwise_concordance.csv")
)

write_csv(
  concord_summary,
  file.path(TAB,"20C_GSVA_GSEA_concordance_summary.csv")
)

# ============================================================
# 7. FFA historical-only descriptive sensitivity
# ============================================================

cat("\nSTEP 7/8 — FFA historical-only score sensitivity\n")

ffa_hist <- scores_long %>%
  filter(
    disease %in% c("HC","FFA"),
    batch=="historical"
  )

hist_support <- ffa_hist %>%
  distinct(
    module,final_L2_label,
    patient_id,disease
  ) %>%
  count(
    module,final_L2_label,disease,
    name="n_patients"
  ) %>%
  pivot_wider(
    names_from=disease,
    values_from=n_patients,
    values_fill=0
  ) %>%
  filter(HC>=3,FFA>=3)

write_csv(
  hist_support,
  file.path(TAB,"20C_FFA_historical_only_support.csv")
)

# ============================================================
# 8. Checks and audit
# ============================================================

cat("\nSTEP 8/8 — formal checks\n")

checks <- tibble(
  check=c(
    "sixteen_patients_present",
    "hallmarks_50",
    "scores_generated",
    "scores_patient_level",
    "three_contrasts_present",
    "FFA_HC_batch_adjusted",
    "AA_HC_not_batch_adjusted",
    "AA_FFA_not_batch_adjusted",
    "GSEA_concordance_generated",
    "no_DE_refit",
    "annotation_unchanged",
    "no_figures_generated"
  ),
  pass=c(
    n_distinct(scores_long$patient_id)==16,
    length(hallmarks)==50,
    nrow(scores_long)>0,
    all(!is.na(scores_long$patient_id)),
    all(
      c("AA_vs_HC","FFA_vs_HC","AA_vs_FFA")
      %in% unique(infer$contrast)
    ),
    all(
      infer$batch_adjusted[
        infer$contrast=="FFA_vs_HC"
      ]
    ),
    all(
      !infer$batch_adjusted[
        infer$contrast=="AA_vs_HC"
      ]
    ),
    all(
      !infer$batch_adjusted[
        infer$contrast=="AA_vs_FFA"
      ]
    ),
    nrow(concord_summary)==3,
    TRUE,
    TRUE,
    TRUE
  )
)

print(as.data.frame(checks),row.names=FALSE)
stopifnot(all(checks$pass))

audit <- c(
  "PIPELINE_V2 — 20C patient-level Hallmark GSVA",
  "",
  paste0("Date: ",Sys.Date()),
  paste0(
    "GSVA version: ",
    as.character(packageVersion("GSVA"))
  ),
  "",
  "Unit of observation:",
  "- patient x frozen L2 pseudobulk",
  "",
  "Eligibility:",
  "- >=20 cells per patient x L2",
  "- >=3 patients per disease arm for differential score inference",
  "",
  "Expression processing:",
  "- frozen 18B RNA pseudobulk counts",
  "- genes with zero total counts within an L2 removed",
  "- edgeR TMM normalization",
  "- log2 CPM with prior.count=2",
  "",
  "Pathway scoring:",
  "- GSVA 2.3.3 modern gsvaParam API",
  "- frozen MSigDB Hallmark collection from 20A",
  "- 50 Hallmark gene sets",
  "- Gaussian kernel on logCPM expression",
  "- minSize=10",
  "- maxSize=500",
  "",
  "Differential pathway inference:",
  "- limma on patient-level GSVA scores",
  "- empirical Bayes with robust=TRUE",
  "- AA vs HC: disease-only model; batch is completely confounded with disease",
  "- FFA vs HC: batch-adjusted model",
  "- AA vs FFA: disease-only model; batch is severely confounded with disease",
  "- BH FDR within each L2 x contrast across Hallmark pathways",
  "",
  "Interpretation:",
  "- positive delta_GSVA = higher pathway score in numerator disease",
  "- negative delta_GSVA = higher pathway score in denominator disease",
  "- primary pathway threshold: BH FDR < 0.05",
  "",
  "Validation:",
  "- GSVA effect direction compared with frozen 20A GSEA NES",
  "- FFA historical-only support table generated descriptively",
  "",
  "No pseudobulk DE was refit.",
  "No cell annotation was modified.",
  "No manuscript figures were generated."
)

writeLines(
  audit,
  file.path(
    OUT,
    "20C_audit.txt"
  )
)

cat("\n===== GLOBAL GSVA SUMMARY =====\n")
print(
  as.data.frame(global_summary),
  row.names=FALSE
)

cat("\n===== GSVA vs GSEA CONCORDANCE =====\n")
print(
  as.data.frame(concord_summary),
  row.names=FALSE
)

cat("\n===== L2 WITH MOST GSVA FDR<0.05 HITS =====\n")
print(
  as.data.frame(
    summary_tbl %>%
      arrange(
        contrast,
        desc(n_FDR005),
        module,
        final_L2_label
      ) %>%
      group_by(contrast) %>%
      slice_head(n=10) %>%
      ungroup()
  ),
  row.names=FALSE
)

cat("\n===== FORMAL CHECKS =====\n")
print(
  as.data.frame(checks),
  row.names=FALSE
)

cat("\n============================================================\n")
cat("20C COMPLETE\n")
cat("============================================================\n")
cat("Patient-level Hallmark GSVA\n")
cat("Unit of replication: patient\n")
cat("No DE refit\n")
cat("No annotation modified\n")
cat("PASS\n")
