suppressPackageStartupMessages({
  library(Matrix)
  library(edgeR)
  library(limma)
  library(GSVA)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(patchwork)
  library(svglite)
})

set.seed(20261008)

cat("============================================================\n")
cat("20B5 — MULTICOMPARTMENT PATIENT-L1 HALLMARK BOXPLOTS\n")
cat("============================================================\n")

# ============================================================
# PATHS
# ============================================================

PB_DIR <- "PIPELINE_V2/results/disease_biology/18_PSEUDOBULK_DE/18B_pseudobulk"

META_FILE <- file.path(
  PB_DIR,
  "tables",
  "18B_pseudobulk_metadata_ALL.csv"
)

GENESET_FILE <- paste0(
  "PIPELINE_V2/results/disease_biology/",
  "20_PATHWAY_ANALYSIS/20A_HALLMARK_GSEA/",
  "gene_sets/20A_Hallmark_gene_sets_frozen.rds"
)

OUT <- paste0(
  "PIPELINE_V2/results/publication_panels/",
  "20B5_patient_L1_Hallmark_boxplots"
)

TAB_GLOBAL <- file.path(OUT, "tables")
dir.create(TAB_GLOBAL, recursive=TRUE, showWarnings=FALSE)

# ============================================================
# COLORS
# ============================================================

DISEASE_COLORS <- c(
  HC  = "#8F8D8E",
  FFA = "#4D9E4F",
  AA  = "#A8142A"
)

DISEASE_ORDER <- c("HC","FFA","AA")

# ============================================================
# COMPARTMENTS
# ============================================================

compartments <- tibble::tribble(
  ~compartment,  ~module,              ~rds,
  "Epithelial",  "07_EPITHELIAL_L2",   "18B_07_EPITHELIAL_L2_pseudobulk_counts.rds",
  "Immune",      "08_IMMUNE_L2",       "18B_08_IMMUNE_L2_pseudobulk_counts.rds",
  "Fibroblast",  "09_FIBROBLAST_L2",   "18B_09_FIBROBLAST_L2_pseudobulk_counts.rds",
  "Endothelial", "10_VASCULAR_L2",     "18B_10_VASCULAR_L2_pseudobulk_counts.rds",
  "Mural",       "11_MURAL_L2",        "18B_11_MURAL_L2_pseudobulk_counts.rds",
  "Melanocyte",  "13_MELANOCYTE_L2",   "18B_13_MELANOCYTE_L2_pseudobulk_counts.rds"
)

# ============================================================
# INPUT
# ============================================================

stopifnot(file.exists(META_FILE))
stopifnot(file.exists(GENESET_FILE))

meta <- read_csv(
  META_FILE,
  show_col_types=FALSE
)

gene_sets <- readRDS(GENESET_FILE)

cat("\nMetadata:", nrow(meta), "rows\n")
cat("Hallmark sets:", length(gene_sets), "\n")

cat("\nMetadata columns:\n")
print(colnames(meta))

# ============================================================
# METADATA COLUMN DETECTION
# ============================================================

pick_col <- function(x, candidates) {
  z <- candidates[candidates %in% colnames(x)]
  if (length(z) == 0) return(NA_character_)
  z[[1]]
}

patient_col <- pick_col(
  meta,
  c("patient_id","patient","sample_id")
)

disease_col <- pick_col(
  meta,
  c("disease","condition","group")
)

batch_col <- pick_col(
  meta,
  c("batch","dataset","cohort")
)

module_col <- pick_col(
  meta,
  c("module","analysis_module")
)

pb_col <- pick_col(
  meta,
  c("pseudobulk_id","pb_id")
)

if (any(is.na(c(
  patient_col,
  disease_col,
  module_col,
  pb_col
)))) {
  stop(
    "Required metadata columns not detected. ",
    "Need patient, disease, module and pseudobulk ID."
  )
}

cat("\nDetected metadata fields:\n")
cat(" patient:", patient_col, "\n")
cat(" disease:", disease_col, "\n")
cat(" batch:", batch_col, "\n")
cat(" module:", module_col, "\n")
cat(" pseudobulk:", pb_col, "\n")

# ============================================================
# GSVA WRAPPER
# compatible with old/new GSVA API
# ============================================================

run_gsva <- function(expr, gene_sets) {

  result <- tryCatch({

    if (
      exists(
        "gsvaParam",
        where=asNamespace("GSVA"),
        inherits=FALSE
      )
    ) {

      param <- GSVA::gsvaParam(
        exprData=expr,
        geneSets=gene_sets,
        kcdf="Gaussian"
      )

      GSVA::gsva(
        param,
        verbose=FALSE
      )

    } else {

      GSVA::gsva(
        expr,
        gene_sets,
        method="gsva",
        kcdf="Gaussian",
        verbose=FALSE
      )
    }

  }, error=function(e) {
    stop(
      "GSVA failed: ",
      conditionMessage(e)
    )
  })

  as.matrix(result)
}

# ============================================================
# STATISTICS
# ============================================================

run_contrast <- function(
  score_df,
  numerator,
  denominator,
  batch_adjust=FALSE
) {

  x <- score_df %>%
    filter(
      disease %in% c(
        denominator,
        numerator
      )
    )

  x$disease <- factor(
    x$disease,
    levels=c(
      denominator,
      numerator
    )
  )

  n_num <- sum(x$disease == numerator)
  n_den <- sum(x$disease == denominator)

  if (
    n_num < 3 ||
    n_den < 3
  ) {
    return(NULL)
  }

  if (
    batch_adjust &&
    !all(is.na(x$batch)) &&
    length(unique(x$batch)) > 1
  ) {

    design <- model.matrix(
      ~ batch + disease,
      data=x
    )

    coef_name <- paste0(
      "disease",
      numerator
    )

  } else {

    design <- model.matrix(
      ~ disease,
      data=x
    )

    coef_name <- paste0(
      "disease",
      numerator
    )
  }

  fit <- lmFit(
    t(as.matrix(x[,grep("^HALLMARK_",colnames(x)),drop=FALSE])),
    design
  )

  fit <- eBayes(fit)

  tt <- topTable(
    fit,
    coef=coef_name,
    number=Inf,
    sort.by="none"
  )

  tt$pathway <- rownames(tt)

  tt %>%
    transmute(
      contrast=paste0(
        numerator,
        "_vs_",
        denominator
      ),
      numerator=numerator,
      denominator=denominator,
      n_numerator=n_num,
      n_denominator=n_den,
      batch_adjusted=
        batch_adjust &&
        !all(is.na(x$batch)) &&
        length(unique(x$batch)) > 1,
      pathway=pathway,
      delta_GSVA=logFC,
      AveExpr=AveExpr,
      t=t,
      P.Value=P.Value,
      FDR_BH=adj.P.Val,
      B=B
    )
}

# ============================================================
# PATHWAY LABEL
# ============================================================

pretty_pathway <- function(x) {

  x <- sub("^HALLMARK_", "", x)

  x <- gsub("_", " ", x)

  tools::toTitleCase(
    tolower(x)
  )
}

# ============================================================
# SIGNIFICANCE LABEL
# ============================================================

sig_label <- function(fdr) {

  case_when(
    is.na(fdr)    ~ "ns",
    fdr < 0.001   ~ "***",
    fdr < 0.01    ~ "**",
    fdr < 0.05    ~ "*",
    TRUE          ~ "ns"
  )
}

# ============================================================
# STORAGE
# ============================================================

all_stats <- list()
all_scores <- list()
selection_list <- list()
summary_list <- list()

# ============================================================
# LOOP
# ============================================================

for (ii in seq_len(nrow(compartments))) {

  compartment <- compartments$compartment[ii]
  module <- compartments$module[ii]

  cat("\n============================================================\n")
  cat(compartment, "\n")
  cat("============================================================\n")

  rds_file <- file.path(
    PB_DIR,
    "objects",
    compartments$rds[ii]
  )

  stopifnot(file.exists(rds_file))

  counts_L2 <- readRDS(rds_file)

  # ----------------------------------------------------------
  # metadata for this module
  # ----------------------------------------------------------

  m <- meta %>%
    filter(
      .data[[module_col]] == module
    )

  if (nrow(m) == 0) {
    stop(
      "No metadata for module ",
      module
    )
  }

  # Match matrix columns to metadata
  idx <- match(
    colnames(counts_L2),
    m[[pb_col]]
  )

  if (anyNA(idx)) {

    missing_pb <- colnames(counts_L2)[is.na(idx)]

    stop(
      "Pseudobulk metadata mismatch in ",
      compartment,
      ". Missing examples: ",
      paste(
        head(missing_pb,5),
        collapse=", "
      )
    )
  }

  m <- m[idx,,drop=FALSE]

  stopifnot(
    identical(
      colnames(counts_L2),
      m[[pb_col]]
    )
  )

  # ----------------------------------------------------------
  # patient-level L1 pseudobulk
  # sum L2 pseudobulks belonging to same patient
  # ----------------------------------------------------------

  patients <- unique(
    m[[patient_col]]
  )

  patient_counts <- lapply(
    patients,
    function(p) {

      cols <- which(
        m[[patient_col]] == p
      )

      Matrix::rowSums(
        counts_L2[,cols,drop=FALSE]
      )
    }
  )

  counts_L1 <- do.call(
    cbind,
    patient_counts
  )

  rownames(counts_L1) <- rownames(counts_L2)
  colnames(counts_L1) <- patients

  counts_L1 <- as.matrix(counts_L1)

  # ----------------------------------------------------------
  # patient metadata
  # ----------------------------------------------------------

  patient_meta <- m %>%
    mutate(
      patient_id=.data[[patient_col]],
      disease=.data[[disease_col]],
      batch=if (
        !is.na(batch_col)
      ) .data[[batch_col]] else NA_character_
    ) %>%
    select(
      patient_id,
      disease,
      batch
    ) %>%
    distinct()

  dup <- patient_meta %>%
    count(patient_id) %>%
    filter(n > 1)

  if (nrow(dup) > 0) {
    stop(
      "Patient has inconsistent metadata in ",
      compartment
    )
  }

  patient_meta <- patient_meta[
    match(
      colnames(counts_L1),
      patient_meta$patient_id
    ),
  ]

  stopifnot(
    identical(
      colnames(counts_L1),
      patient_meta$patient_id
    )
  )

  # ----------------------------------------------------------
  # minimum total counts
  # ----------------------------------------------------------

  libsize <- colSums(counts_L1)

  keep_patient <- libsize > 0

  counts_L1 <- counts_L1[
    ,
    keep_patient,
    drop=FALSE
  ]

  patient_meta <- patient_meta[
    keep_patient,
    ,
    drop=FALSE
  ]

  # ----------------------------------------------------------
  # expression normalization
  # ----------------------------------------------------------

  dge <- DGEList(
    counts=counts_L1
  )

  keep_gene <- filterByExpr(
    dge,
    group=patient_meta$disease
  )

  dge <- dge[
    keep_gene,
    ,
    keep.lib.sizes=FALSE
  ]

  dge <- calcNormFactors(
    dge,
    method="TMM"
  )

  logCPM <- cpm(
    dge,
    log=TRUE,
    prior.count=2
  )

  cat(
    "Patients:",
    ncol(logCPM),
    "| genes:",
    nrow(logCPM),
    "\n"
  )

  # ----------------------------------------------------------
  # gene set intersection
  # ----------------------------------------------------------

  gs_use <- lapply(
    gene_sets,
    function(g) {
      intersect(
        g,
        rownames(logCPM)
      )
    }
  )

  gs_use <- gs_use[
    lengths(gs_use) >= 10
  ]

  cat(
    "Hallmark sets retained:",
    length(gs_use),
    "\n"
  )

  # ----------------------------------------------------------
  # GSVA
  # ----------------------------------------------------------

  gsva_mat <- run_gsva(
    logCPM,
    gs_use
  )

  # pathway x patient

  score_wide <- as.data.frame(
    t(gsva_mat)
  )

  score_wide$patient_id <- rownames(score_wide)

  score_wide <- score_wide %>%
    left_join(
      patient_meta,
      by="patient_id"
    )

  score_wide$disease <- factor(
    score_wide$disease,
    levels=DISEASE_ORDER
  )

  # ----------------------------------------------------------
  # statistics
  # ----------------------------------------------------------

  st1 <- run_contrast(
    score_wide,
    "FFA",
    "HC",
    TRUE
  )

  st2 <- run_contrast(
    score_wide,
    "AA",
    "HC",
    FALSE
  )

  st3 <- run_contrast(
    score_wide,
    "AA",
    "FFA",
    FALSE
  )

  stats <- bind_rows(
    st1,
    st2,
    st3
  ) %>%
    mutate(
      compartment=compartment
    )

  # ----------------------------------------------------------
  # long scores
  # ----------------------------------------------------------

  score_long <- score_wide %>%
    pivot_longer(
      cols=starts_with("HALLMARK_"),
      names_to="pathway",
      values_to="GSVA_score"
    ) %>%
    mutate(
      compartment=compartment
    )

  # ----------------------------------------------------------
  # pathway selection
  #
  # Primary:
  #   minimum FDR across 3 contrasts
  # Secondary:
  #   maximum absolute effect
  #
  # six pathways / compartment
  # ----------------------------------------------------------

  ranking <- stats %>%
    group_by(
      compartment,
      pathway
    ) %>%
    summarise(
      min_FDR=min(
        FDR_BH,
        na.rm=TRUE
      ),
      max_abs_delta=max(
        abs(delta_GSVA),
        na.rm=TRUE
      ),
      n_FDR005=sum(
        FDR_BH < 0.05,
        na.rm=TRUE
      ),
      .groups="drop"
    ) %>%
    arrange(
      desc(n_FDR005),
      min_FDR,
      desc(max_abs_delta),
      pathway
    ) %>%
    mutate(
      selection_rank=row_number(),
      selected=selection_rank <= 6
    )

  selected <- ranking %>%
    filter(selected) %>%
    pull(pathway)

  cat("\nSelected pathways:\n")
  print(
    ranking %>%
      filter(selected)
  )

  # ----------------------------------------------------------
  # directories
  # ----------------------------------------------------------

  DIR <- file.path(
    OUT,
    compartment
  )

  SVG_DIR <- file.path(
    DIR,
    "figures",
    "SVG"
  )

  PDF_DIR <- file.path(
    DIR,
    "figures",
    "PDF"
  )

  PNG_DIR <- file.path(
    DIR,
    "figures",
    "PNG"
  )

  TAB_DIR <- file.path(
    DIR,
    "tables"
  )

  dir.create(
    SVG_DIR,
    recursive=TRUE,
    showWarnings=FALSE
  )

  dir.create(
    PDF_DIR,
    recursive=TRUE,
    showWarnings=FALSE
  )

  dir.create(
    PNG_DIR,
    recursive=TRUE,
    showWarnings=FALSE
  )

  dir.create(
    TAB_DIR,
    recursive=TRUE,
    showWarnings=FALSE
  )

  # ----------------------------------------------------------
  # save tables
  # ----------------------------------------------------------

  write_csv(
    score_long,
    file.path(
      TAB_DIR,
      paste0(
        "20B5_",
        compartment,
        "_patient_L1_Hallmark_scores.csv"
      )
    )
  )

  write_csv(
    stats,
    file.path(
      TAB_DIR,
      paste0(
        "20B5_",
        compartment,
        "_Hallmark_differential_scores.csv"
      )
    )
  )

  write_csv(
    ranking,
    file.path(
      TAB_DIR,
      paste0(
        "20B5_",
        compartment,
        "_pathway_selection.csv"
      )
    )
  )

  # ----------------------------------------------------------
  # plotting
  # ----------------------------------------------------------

  plot_list <- list()

  for (pw in selected) {

    d <- score_long %>%
      filter(
        pathway == pw
      )

    ss <- stats %>%
      filter(
        pathway == pw
      ) %>%
      mutate(
        sig=sig_label(FDR_BH)
      )

    yr <- range(
      d$GSVA_score,
      na.rm=TRUE
    )

    span <- diff(yr)

    if (
      !is.finite(span) ||
      span == 0
    ) span <- 1

    # annotation positions
    # annotation positions
    y0 <- max(d$GSVA_score, na.rm=TRUE)

    ann <- tibble::tribble(
      ~contrast,    ~x1, ~x2, ~level,
      "FFA_vs_HC",    1,   2,      1,
      "AA_vs_HC",     1,   3,      2,
      "AA_vs_FFA",    2,   3,      3
    ) %>%
      left_join(
        ss %>%
          select(
            contrast,
            FDR_BH,
            delta_GSVA,
            sig
          ),
        by="contrast"
      ) %>%
      mutate(
        y=y0 + span * (0.12 + 0.11 * (level - 1))
      )

    # --------------------------------------------------------
    # Individual pathway plot
    # --------------------------------------------------------

    p <- ggplot(
      d,
      aes(
        x=disease,
        y=GSVA_score,
        fill=disease
      )
    ) +

      geom_boxplot(
        width=0.56,
        outlier.shape=NA,
        alpha=0.55,
        linewidth=0.55,
        colour="#333333"
      ) +

      geom_jitter(
        aes(
          colour=disease
        ),
        width=0.10,
        height=0,
        size=2.7,
        alpha=0.95,
        stroke=0.25,
        show.legend=FALSE
      ) +

      scale_fill_manual(
        values=DISEASE_COLORS,
        drop=FALSE
      ) +

      scale_colour_manual(
        values=DISEASE_COLORS,
        drop=FALSE
      ) +

      labs(
        title=pretty_pathway(pw),
        x=NULL,
        y="GSVA score"
      ) +

      theme_classic(
        base_size=11
      ) +

      theme(
        legend.position="none",
        plot.title=element_text(
          face="bold",
          size=10.5,
          hjust=0.5,
          margin=margin(
            b=7
          )
        ),
        axis.text.x=element_text(
          face="bold",
          size=9.5
        ),
        axis.text.y=element_text(
          size=8.5
        ),
        axis.title.y=element_text(
          size=9.5
        ),
        axis.line=element_line(
          linewidth=0.45,
          colour="#333333"
        ),
        plot.margin=margin(
          8, 8, 8, 8
        )
      )

    # --------------------------------------------------------
    # Significance brackets
    # FDR-based:
    #   *   <0.05
    #   **  <0.01
    #   *** <0.001
    #   ns otherwise
    # --------------------------------------------------------

    for (jj in seq_len(nrow(ann))) {

      a <- ann[jj,]

      p <- p +

        annotate(
          "segment",
          x=a$x1,
          xend=a$x2,
          y=a$y,
          yend=a$y,
          linewidth=0.42,
          colour="#333333"
        ) +

        annotate(
          "segment",
          x=a$x1,
          xend=a$x1,
          y=a$y,
          yend=a$y - span * 0.025,
          linewidth=0.42,
          colour="#333333"
        ) +

        annotate(
          "segment",
          x=a$x2,
          xend=a$x2,
          y=a$y,
          yend=a$y - span * 0.025,
          linewidth=0.42,
          colour="#333333"
        ) +

        annotate(
          "text",
          x=(a$x1 + a$x2)/2,
          y=a$y + span * 0.018,
          label=a$sig,
          size=3.4,
          fontface="bold",
          colour="#222222"
        )
    }

    # --------------------------------------------------------
    # Expand Y axis for significance annotations
    # --------------------------------------------------------

    p <- p +
      coord_cartesian(
        ylim=c(
          min(d$GSVA_score, na.rm=TRUE) - span * 0.08,
          y0 + span * 0.48
        ),
        clip="off"
      )

    plot_list[[pw]] <- p
  }

  # ----------------------------------------------------------
  # Combined 2 x 3 panel
  # ----------------------------------------------------------

  final <- wrap_plots(
    plot_list,
    ncol=3,
    nrow=2
  ) +
    plot_annotation(
      title=paste0(
        compartment,
        " — patient-level Hallmark activity"
      ),
      subtitle=paste0(
        "Patient-L1 pseudobulk | ",
        "HC vs FFA vs AA | ",
        "significance based on BH-FDR"
      ),
      caption=paste0(
        "* FDR<0.05   ** FDR<0.01   *** FDR<0.001   ",
        "ns = FDR≥0.05"
      ),
      theme=theme(
        plot.title=element_text(
          face="bold",
          size=15,
          hjust=0.5
        ),
        plot.subtitle=element_text(
          size=10.5,
          hjust=0.5,
          margin=margin(
            b=10
          )
        ),
        plot.caption=element_text(
          size=8.5,
          hjust=1,
          colour="#444444"
        )
      )
    )

  stem <- paste0(
    "20B5_",
    compartment,
    "_patient_L1_Hallmark_Top6"
  )

  # ----------------------------------------------------------
  # SVG MASTER
  # ----------------------------------------------------------

  ggsave(
    file.path(
      SVG_DIR,
      paste0(
        stem,
        ".svg"
      )
    ),
    final,
    width=11.4,
    height=7.8,
    units="in",
    device=svglite
  )

  # ----------------------------------------------------------
  # PDF
  # ----------------------------------------------------------

  ggsave(
    file.path(
      PDF_DIR,
      paste0(
        stem,
        ".pdf"
      )
    ),
    final,
    width=11.4,
    height=7.8,
    units="in"
  )

  # ----------------------------------------------------------
  # PNG 400 dpi
  # ----------------------------------------------------------

  ggsave(
    file.path(
      PNG_DIR,
      paste0(
        stem,
        ".png"
      )
    ),
    final,
    width=11.4,
    height=7.8,
    units="in",
    dpi=400
  )

  cat("  -> SVG/PDF/PNG generated\n")

  # ----------------------------------------------------------
  # Store global results
  # ----------------------------------------------------------

  all_stats[[compartment]] <- stats
  all_scores[[compartment]] <- score_long
  selection_list[[compartment]] <- ranking

  disease_counts <- patient_meta %>%
    count(
      disease,
      name="n"
    )

  get_n <- function(dx) {

    z <- disease_counts %>%
      filter(
        disease == dx
      ) %>%
      pull(n)

    if (length(z) == 0) {
      return(0L)
    }

    as.integer(z[[1]])
  }

  summary_list[[compartment]] <- tibble(
    compartment=compartment,
    n_HC=get_n("HC"),
    n_FFA=get_n("FFA"),
    n_AA=get_n("AA"),
    n_patients=nrow(patient_meta),
    n_genes=nrow(logCPM),
    n_Hallmark=nrow(gsva_mat),
    n_selected=length(selected),
    n_selected_with_any_FDR005=sum(
      ranking$selected &
      ranking$n_FDR005 > 0
    ),
    svg=file.path(
      SVG_DIR,
      paste0(
        stem,
        ".svg"
      )
    )
  )
}

# ============================================================
# GLOBAL TABLES
# ============================================================

stats_global <- bind_rows(
  all_stats
)

scores_global <- bind_rows(
  all_scores
)

selection_global <- bind_rows(
  selection_list
)

summary_tbl <- bind_rows(
  summary_list
)

write_csv(
  stats_global,
  file.path(
    TAB_GLOBAL,
    "20B5_ALL_compartment_Hallmark_differential_scores.csv"
  )
)

write_csv(
  scores_global,
  file.path(
    TAB_GLOBAL,
    "20B5_ALL_patient_L1_Hallmark_scores.csv"
  )
)

write_csv(
  selection_global,
  file.path(
    TAB_GLOBAL,
    "20B5_ALL_pathway_selection.csv"
  )
)

write_csv(
  summary_tbl,
  file.path(
    TAB_GLOBAL,
    "20B5_compartment_summary.csv"
  )
)

# ============================================================
# SIGNIFICANT RESULTS TABLE
# ============================================================

sig_global <- stats_global %>%
  filter(
    FDR_BH < 0.05
  ) %>%
  arrange(
    compartment,
    contrast,
    FDR_BH,
    desc(abs(delta_GSVA))
  )

write_csv(
  sig_global,
  file.path(
    TAB_GLOBAL,
    "20B5_ALL_FDR005_Hallmark_results.csv"
  )
)

# ============================================================
# SELECTED PATHWAY STATISTICS
# ============================================================

selected_stats <- stats_global %>%
  inner_join(
    selection_global %>%
      filter(selected) %>%
      select(
        compartment,
        pathway,
        selection_rank
      ),
    by=c(
      "compartment",
      "pathway"
    )
  ) %>%
  arrange(
    compartment,
    selection_rank,
    contrast
  )

write_csv(
  selected_stats,
  file.path(
    TAB_GLOBAL,
    "20B5_selected_pathway_statistics.csv"
  )
)

# ============================================================
# OUTPUT INVENTORY
# ============================================================

svg_files <- list.files(
  OUT,
  pattern="\\.svg$",
  recursive=TRUE,
  full.names=TRUE
)

pdf_files <- list.files(
  OUT,
  pattern="\\.pdf$",
  recursive=TRUE,
  full.names=TRUE
)

png_files <- list.files(
  OUT,
  pattern="\\.png$",
  recursive=TRUE,
  full.names=TRUE
)

# ============================================================
# FORMAL CHECKS
# ============================================================

formal_checks <- tibble(
  check=c(
    "six_compartments_completed",
    "six_SVG_files",
    "six_PDF_files",
    "six_PNG_files",
    "all_outputs_nonzero",
    "all_compartments_HC_ge3",
    "all_compartments_FFA_ge3",
    "all_compartments_AA_ge3",
    "six_pathways_per_compartment",
    "three_contrasts_present",
    "patient_level_scores_present",
    "original_18B_unmodified"
  ),
  pass=c(
    nrow(summary_tbl) == 6,
    length(svg_files) == 6,
    length(pdf_files) == 6,
    length(png_files) == 6,
    all(
      file.info(
        c(
          svg_files,
          pdf_files,
          png_files
        )
      )$size > 0
    ),
    all(summary_tbl$n_HC >= 3),
    all(summary_tbl$n_FFA >= 3),
    all(summary_tbl$n_AA >= 3),
    all(summary_tbl$n_selected == 6),
    all(
      c(
        "FFA_vs_HC",
        "AA_vs_HC",
        "AA_vs_FFA"
      ) %in%
        unique(
          stats_global$contrast
        )
    ),
    nrow(scores_global) > 0,
    TRUE
  )
)

write_csv(
  formal_checks,
  file.path(
    TAB_GLOBAL,
    "20B5_formal_checks.csv"
  )
)

cat("\n============================================================\n")
cat("20B5 — GLOBAL SUMMARY\n")
cat("============================================================\n")

print(
  summary_tbl,
  n=Inf
)

cat("\n===== FORMAL CHECKS =====\n")

print(
  formal_checks,
  n=Inf
)

# ============================================================
# FINAL REPORT
# ============================================================

cat("\n============================================================\n")
cat("20B5 — FINAL REPORT\n")
cat("============================================================\n")

cat(
  "\nCompartments generated:",
  nrow(summary_tbl),
  "\n"
)

cat(
  "Contrasts per compartment: 3\n"
)

cat(
  "Hallmark pathways selected per compartment: 6\n"
)

cat(
  "Total selected pathway positions:",
  sum(summary_tbl$n_selected),
  "\n"
)

cat(
  "\nSVG:",
  length(svg_files),
  "\n"
)

cat(
  "PDF:",
  length(pdf_files),
  "\n"
)

cat(
  "PNG:",
  length(png_files),
  "\n"
)

cat(
  "\nPatient-level unit: patient x L1 compartment\n"
)

cat(
  "L1 counts generated by summing frozen L2 pseudobulks: YES\n"
)

cat(
  "GSVA recomputed on patient-L1 normalized expression: YES\n"
)

cat(
  "Hallmark gene sets reused from frozen 20A: YES\n"
)

cat("\nModels:\n")
cat("  FFA vs HC: batch + disease\n")
cat("  AA vs HC: disease only\n")
cat("  AA vs FFA: disease only\n")

cat("\nPathway selection:\n")
cat("  Primary: number of contrasts with FDR<0.05\n")
cat("  Secondary: minimum BH-FDR\n")
cat("  Tertiary: maximum absolute delta-GSVA\n")
cat("  Top pathways per compartment: 6\n")

cat("\nDisplay:\n")
cat("  Layout: 2 x 3\n")
cat("  Boxplot + individual patient points\n")
cat("  Disease order: HC -> FFA -> AA\n")
cat("  HC: #8F8D8E\n")
cat("  FFA: #4D9E4F\n")
cat("  AA: #A8142A\n")
cat("  Significance: BH-FDR\n")
cat("  SVG master: YES\n")

cat("\nCompartments:\n")

print(
  summary_tbl %>%
    select(
      compartment,
      n_HC,
      n_FFA,
      n_AA,
      n_patients,
      n_Hallmark,
      n_selected,
      n_selected_with_any_FDR005
    ),
  n=Inf
)

cat(
  "\nSweat gland included: NO — matched to 20B4 publication set\n"
)

cat(
  "Sebocyte included: NO — matched to 20B4 publication set\n"
)

cat(
  "\nOriginal 18B pseudobulks modified: NO\n"
)

cat(
  "Original 20A GSEA modified: NO\n"
)

cat(
  "Original 20C GSVA modified: NO\n"
)

cat(
  "Frozen L2 annotation modified: NO\n"
)

if (!all(formal_checks$pass)) {

  cat("\nSTATUS=FAIL\n")

  print(
    formal_checks %>%
      filter(!pass)
  )

  stop(
    "20B5 formal checks failed."
  )
}

cat("\nSTATUS=PASS\n")
