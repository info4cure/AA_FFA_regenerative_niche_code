suppressPackageStartupMessages({
  library(edgeR)
  library(dplyr)
  library(readr)
  library(tibble)
  library(purrr)
  library(stringr)
})

cat("============================================================\n")
cat("PIPELINE_V2 — 18D3 FFA vs HC BATCH MODEL SENSITIVITY\n")
cat("============================================================\n")

BASE <- file.path(
  "PIPELINE_V2", "results", "disease_biology", "18_PSEUDOBULK_DE"
)

PB_BASE <- file.path(BASE, "18B_pseudobulk")
PRIMARY <- file.path(BASE, "18C_DE")

OUT <- file.path(BASE, "18D_interpretation", "18D3_FFA_batch_sensitivity")
TAB <- file.path(OUT, "tables")
FULL <- file.path(OUT, "full_results")

dir.create(TAB, recursive=TRUE, showWarnings=FALSE)
dir.create(FULL, recursive=TRUE, showWarnings=FALSE)

# ============================================================
# 1. Load frozen 18B pseudobulk
# ============================================================

cat("\nSTEP 1/7 — loading frozen pseudobulk\n")

pb_meta <- readRDS(
  file.path(PB_BASE, "objects", "18B_pseudobulk_metadata_ALL.rds")
)

modules <- sort(unique(pb_meta$module))

pb <- list()

for (m in modules) {
  f <- file.path(
    PB_BASE, "objects",
    paste0("18B_", m, "_pseudobulk_counts.rds")
  )
  stopifnot(file.exists(f))
  pb[[m]] <- readRDS(f)
}

# Eligibility is frozen from 18A
elig <- read_csv(
  file.path(
    BASE, "18A_preflight", "tables",
    "18A_pseudobulk_contrast_eligibility.csv"
  ),
  show_col_types=FALSE
)

cat("Modules:", length(modules), "\n")

# ============================================================
# 2. Identify FFA-vs-HC eligible L2
# ============================================================

cat("\nSTEP 2/7 — recovering FFA vs HC eligibility\n")

cat("Eligibility columns:\n")
print(names(elig))

# Avoid assuming exact eligibility-column implementation.
# Recover analyses directly from the frozen 18C FFA-vs-HC result files.

primary_files <- list.files(
  PRIMARY,
  pattern="^18C_FFA_vs_HC__.*\\.csv$",
  recursive=TRUE,
  full.names=TRUE
)

stopifnot(length(primary_files) == 18)

primary_index <- map_dfr(primary_files, function(f) {

  x <- read_csv(f, show_col_types=FALSE)

  stopifnot(
    "module" %in% names(x),
    "final_L2_label" %in% names(x)
  )

  tibble(
    primary_file=f,
    module=unique(x$module),
    final_L2_label=unique(x$final_L2_label)
  )
})

stopifnot(
  nrow(primary_index) == 18,
  !anyDuplicated(
    paste(primary_index$module, primary_index$final_L2_label)
  )
)

print(primary_index %>% select(module, final_L2_label), n=Inf)

# ============================================================
# 3. edgeR function
# ============================================================

run_model <- function(module_name, l2_label,
                      model_type=c("adjusted","unadjusted","historical")) {

  model_type <- match.arg(model_type)

  mat <- pb[[module_name]]

  md <- pb_meta %>%
    filter(
      module == module_name,
      final_L2_label == l2_label,
      disease %in% c("HC","FFA"),
      n_cells >= 20
    )

  if (model_type == "historical") {
    md <- md %>%
      filter(batch == "historical")
  }

  md <- md %>%
    mutate(
      disease_model=factor(disease, levels=c("HC","FFA"))
    )

  n_HC  <- sum(md$disease_model == "HC")
  n_FFA <- sum(md$disease_model == "FFA")

  if (n_HC < 3 || n_FFA < 3)
    return(NULL)

  stopifnot(all(md$pseudobulk_id %in% colnames(mat)))

  counts <- mat[, md$pseudobulk_id, drop=FALSE]

  stopifnot(identical(colnames(counts), md$pseudobulk_id))

  if (model_type == "adjusted") {

    md$batch_model <- factor(md$batch)

    design <- model.matrix(
      ~ batch_model + disease_model,
      data=md
    )

  } else {

    design <- model.matrix(
      ~ disease_model,
      data=md
    )
  }

  stopifnot(qr(design)$rank == ncol(design))

  coef_idx <- which(
    colnames(design) == "disease_modelFFA"
  )

  stopifnot(length(coef_idx) == 1)

  y <- DGEList(counts=counts)

  keep <- filterByExpr(
    y,
    design=design
  )

  stopifnot(sum(keep) > 0)

  y <- y[
    keep,,
    keep.lib.sizes=FALSE
  ]

  y <- calcNormFactors(
    y,
    method="TMM"
  )

  y <- estimateDisp(
    y,
    design=design,
    robust=TRUE
  )

  fit <- glmQLFit(
    y,
    design=design,
    robust=TRUE
  )

  qlf <- glmQLFTest(
    fit,
    coef=coef_idx
  )

  res <- topTags(
    qlf,
    n=Inf,
    sort.by="PValue"
  )$table %>%
    rownames_to_column("gene") %>%
    as_tibble() %>%
    rename(FDR_BH=FDR) %>%
    mutate(
      model=model_type,
      module=module_name,
      final_L2_label=l2_label,
      n_HC=n_HC,
      n_FFA=n_FFA,
      n_genes_tested=sum(keep),
      .before=1
    )

  res
}

# ============================================================
# 4. Run three specifications
# ============================================================

cat("\nSTEP 3/7 — fitting three model specifications\n")

all_results <- list()

for (i in seq_len(nrow(primary_index))) {

  m  <- primary_index$module[i]
  l2 <- primary_index$final_L2_label[i]

  cat(
    sprintf(
      "[%02d/18] %s :: %s\n",
      i, m, l2
    )
  )

  for (model in c("adjusted","unadjusted","historical")) {

    r <- run_model(m, l2, model)

    if (is.null(r))
      next

    key <- paste(m, l2, model, sep="__")
    all_results[[key]] <- r

    write_csv(
      r,
      file.path(
        FULL,
        paste0(
          "18D3_",
          model, "__",
          gsub("[^A-Za-z0-9._-]+","_",m), "__",
          gsub("[^A-Za-z0-9._-]+","_",l2),
          ".csv"
        )
      )
    )
  }
}

# ============================================================
# 5. Per-model summary
# ============================================================

cat("\nSTEP 4/7 — summarising DE\n")

summary_tbl <- map_dfr(
  names(all_results),
  function(k) {

    x <- all_results[[k]]

    tibble(
      model=unique(x$model),
      module=unique(x$module),
      final_L2_label=unique(x$final_L2_label),
      n_HC=unique(x$n_HC),
      n_FFA=unique(x$n_FFA),
      n_genes_tested=unique(x$n_genes_tested),

      n_FDR_001=sum(x$FDR_BH < .01, na.rm=TRUE),
      n_FDR_005=sum(x$FDR_BH < .05, na.rm=TRUE),
      n_FDR_010=sum(x$FDR_BH < .10, na.rm=TRUE),

      n_up_FDR_005=sum(
        x$FDR_BH < .05 & x$logFC > 0,
        na.rm=TRUE
      ),

      n_down_FDR_005=sum(
        x$FDR_BH < .05 & x$logFC < 0,
        na.rm=TRUE
      ),

      n_FDR010_LFC05=sum(
        x$FDR_BH < .10 & abs(x$logFC) >= .5,
        na.rm=TRUE
      )
    )
  }
)

write_csv(
  summary_tbl,
  file.path(TAB, "18D3_model_DE_summary.csv")
)

# ============================================================
# 6. Gene-level model concordance
# ============================================================

cat("\nSTEP 5/7 — model concordance\n")

get_res <- function(m,l2,model) {
  key <- paste(m,l2,model,sep="__")
  all_results[[key]]
}

concordance <- map_dfr(
  seq_len(nrow(primary_index)),
  function(i) {

    m  <- primary_index$module[i]
    l2 <- primary_index$final_L2_label[i]

    A <- get_res(m,l2,"adjusted")
    U <- get_res(m,l2,"unadjusted")
    H <- get_res(m,l2,"historical")

    if (is.null(A) || is.null(U))
      return(NULL)

    AU <- inner_join(
      A %>% select(gene, logFC_A=logFC, P_A=PValue, FDR_A=FDR_BH),
      U %>% select(gene, logFC_U=logFC, P_U=PValue, FDR_U=FDR_BH),
      by="gene"
    )

    rho_AU <- cor(
      AU$logFC_A,
      AU$logFC_U,
      method="spearman",
      use="complete.obs"
    )

    pearson_AU <- cor(
      AU$logFC_A,
      AU$logFC_U,
      method="pearson",
      use="complete.obs"
    )

    dir_AU <- mean(
      sign(AU$logFC_A) == sign(AU$logFC_U),
      na.rm=TRUE
    )

    rho_UH <- NA_real_
    rho_AH <- NA_real_
    dir_UH <- NA_real_
    dir_AH <- NA_real_

    if (!is.null(H)) {

      UH <- inner_join(
        U %>% select(gene, logFC_U=logFC),
        H %>% select(gene, logFC_H=logFC),
        by="gene"
      )

      AH <- inner_join(
        A %>% select(gene, logFC_A=logFC),
        H %>% select(gene, logFC_H=logFC),
        by="gene"
      )

      rho_UH <- cor(
        UH$logFC_U,
        UH$logFC_H,
        method="spearman",
        use="complete.obs"
      )

      rho_AH <- cor(
        AH$logFC_A,
        AH$logFC_H,
        method="spearman",
        use="complete.obs"
      )

      dir_UH <- mean(
        sign(UH$logFC_U) == sign(UH$logFC_H),
        na.rm=TRUE
      )

      dir_AH <- mean(
        sign(AH$logFC_A) == sign(AH$logFC_H),
        na.rm=TRUE
      )
    }

    sA <- summary_tbl %>%
      filter(
        model=="adjusted",
        module==m,
        final_L2_label==l2
      )

    sU <- summary_tbl %>%
      filter(
        model=="unadjusted",
        module==m,
        final_L2_label==l2
      )

    sH <- summary_tbl %>%
      filter(
        model=="historical",
        module==m,
        final_L2_label==l2
      )

    tibble(
      module=m,
      final_L2_label=l2,

      genes_common_AU=nrow(AU),

      rho_logFC_AU=rho_AU,
      pearson_logFC_AU=pearson_AU,
      direction_AU=dir_AU,

      rho_logFC_UH=rho_UH,
      rho_logFC_AH=rho_AH,
      direction_UH=dir_UH,
      direction_AH=dir_AH,

      adjusted_FDR005=sA$n_FDR_005,
      adjusted_FDR010=sA$n_FDR_010,

      unadjusted_FDR005=sU$n_FDR_005,
      unadjusted_FDR010=sU$n_FDR_010,

      historical_FDR005=
        if(nrow(sH)) sH$n_FDR_005 else NA_integer_,

      historical_FDR010=
        if(nrow(sH)) sH$n_FDR_010 else NA_integer_
    )
  }
)

write_csv(
  concordance,
  file.path(TAB, "18D3_L2_model_concordance.csv")
)

# ============================================================
# 7. Reproduction check against frozen 18C
# ============================================================

cat("\nSTEP 6/7 — verifying adjusted model reproduces 18C\n")

repro <- map_dfr(
  seq_len(nrow(primary_index)),
  function(i) {

    m <- primary_index$module[i]
    l2 <- primary_index$final_L2_label[i]

    frozen <- read_csv(
      primary_index$primary_file[i],
      show_col_types=FALSE
    )

    fresh <- get_res(m,l2,"adjusted")

    z <- inner_join(
      frozen %>%
        select(
          gene,
          logFC_frozen=logFC,
          FDR_frozen=FDR_BH
        ),
      fresh %>%
        select(
          gene,
          logFC_fresh=logFC,
          FDR_fresh=FDR_BH
        ),
      by="gene"
    )

    tibble(
      module=m,
      final_L2_label=l2,
      frozen_genes=nrow(frozen),
      fresh_genes=nrow(fresh),
      common_genes=nrow(z),

      max_abs_logFC_difference=max(
        abs(z$logFC_frozen-z$logFC_fresh),
        na.rm=TRUE
      ),

      max_abs_FDR_difference=max(
        abs(z$FDR_frozen-z$FDR_fresh),
        na.rm=TRUE
      ),

      exact_gene_count=
        nrow(frozen)==nrow(fresh)
    )
  }
)

write_csv(
  repro,
  file.path(TAB, "18D3_18C_reproduction_check.csv")
)

stopifnot(
  all(repro$exact_gene_count),
  all(repro$max_abs_logFC_difference < 1e-10),
  all(repro$max_abs_FDR_difference < 1e-10)
)

# ============================================================
# Final summaries
# ============================================================

cat("\nSTEP 7/7 — final summary\n")

totals <- summary_tbl %>%
  group_by(model) %>%
  summarise(
    n_L2=n(),
    total_FDR001=sum(n_FDR_001),
    total_FDR005=sum(n_FDR_005),
    total_FDR010=sum(n_FDR_010),
    total_FDR010_LFC05=sum(n_FDR010_LFC05),
    total_up_FDR005=sum(n_up_FDR_005),
    total_down_FDR005=sum(n_down_FDR_005),
    .groups="drop"
  )

write_csv(
  totals,
  file.path(TAB, "18D3_global_model_summary.csv")
)

cat("\n===== GLOBAL MODEL SUMMARY =====\n")
print(as.data.frame(totals), row.names=FALSE)

cat("\n===== L2 MODEL CONCORDANCE =====\n")
print(
  as.data.frame(
    concordance %>%
      select(
        final_L2_label,
        adjusted_FDR005,
        unadjusted_FDR005,
        historical_FDR005,
        rho_logFC_AU,
        rho_logFC_UH,
        direction_AU,
        direction_UH
      )
  ),
  row.names=FALSE
)

cat("\n===== FORMAL CHECKS =====\n")
cat("18 L2 recovered                 :", nrow(primary_index)==18, "\n")
cat("18C reproduced exactly          : TRUE\n")
cat("Primary 18C modified            : FALSE\n")
cat("Annotation modified             : FALSE\n")
cat("Sensitivity only                : TRUE\n")

audit <- c(
  "PIPELINE_V2 — 18D3 FFA vs HC batch-model sensitivity",
  "",
  "Primary 18C was not modified.",
  "Three specifications compared:",
  "A) adjusted: all eligible FFA/HC, ~ batch + disease",
  "B) unadjusted: same patients, ~ disease",
  "C) historical-only: historical patients only, ~ disease",
  "",
  "All models use:",
  "- RNA pseudobulk counts from frozen 18B",
  "- >=20 cells per patient x L2",
  "- edgeR filterByExpr",
  "- TMM normalization",
  "- estimateDisp robust=TRUE",
  "- glmQLFit robust=TRUE",
  "- glmQLFTest",
  "- BH FDR",
  "",
  "Adjusted model was required to reproduce frozen 18C exactly.",
  "Historical-only is a sensitivity analysis, not an independent cohort.",
  "No annotation modified."
)


writeLines(
  audit,
  file.path(OUT, "18D3_audit.txt")
)

cat("\n============================================================\n")
cat("18D3 COMPLETE\n")
cat("============================================================\n")
cat("FFA vs HC batch-model sensitivity only\n")
cat("Primary 18C unchanged\n")
cat("PASS\n")
