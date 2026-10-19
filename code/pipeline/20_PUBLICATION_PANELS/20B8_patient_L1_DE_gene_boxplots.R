suppressPackageStartupMessages({
  library(Matrix)
  library(edgeR)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(patchwork)
  library(svglite)
  library(stringr)
})

set.seed(20261008)

cat("============================================================\n")
cat("20B8 — PATIENT-L1 PSEUDOBULK DE GENE BOXPLOTS\n")
cat("============================================================\n")

# ============================================================
# PATHS
# ============================================================

PB_DIR <- paste0(
  "PIPELINE_V2/results/disease_biology/",
  "18_PSEUDOBULK_DE/18B_pseudobulk"
)

META_FILE <- file.path(
  PB_DIR,
  "tables",
  "18B_pseudobulk_metadata_ALL.csv"
)

DE_ROOT <- paste0(
  "PIPELINE_V2/results/publication_panels/",
  "20B4_all_L1_DE_heatmaps"
)

OUT <- paste0(
  "PIPELINE_V2/results/publication_panels/",
  "20B8_patient_L1_DE_gene_boxplots"
)

GLOBAL_TAB <- file.path(
  OUT,
  "tables"
)

dir.create(
  GLOBAL_TAB,
  recursive=TRUE,
  showWarnings=FALSE
)

# ============================================================
# PARAMETERS — FROZEN FOR 20B8
# ============================================================

FDR_THRESHOLD <- 0.10
ABS_LOGFC_THRESHOLD <- 1
LOGCPM_THRESHOLD <- 1
MIN_L1_CELLS <- 20
MAX_GENES <- 6

contrast_order <- c(
  "FFA_vs_HC",
  "AA_vs_HC",
  "AA_vs_FFA"
)

disease_order <- c(
  "HC",
  "FFA",
  "AA"
)

DISEASE_COLORS <- c(
  HC  = "#8F8D8E",
  FFA = "#4D9E4F",
  AA  = "#A8142A"
)

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
# INPUT METADATA
# ============================================================

stopifnot(
  file.exists(META_FILE)
)

meta_all <- read_csv(
  META_FILE,
  show_col_types=FALSE
)

required_meta <- c(
  "module",
  "pseudobulk_id",
  "patient_id",
  "disease",
  "batch",
  "n_cells"
)

stopifnot(
  all(
    required_meta %in%
      colnames(meta_all)
  )
)

# ============================================================
# PUBLICATION GENE FILTER
# ============================================================

publication_gene <- function(gene) {

  !grepl(
    "^MIR",
    gene,
    ignore.case=TRUE
  ) &
  !grepl(
    "^LINC",
    gene,
    ignore.case=TRUE
  ) &
  !grepl(
    "-AS[0-9]*$",
    gene,
    ignore.case=TRUE
  )
}

# ============================================================
# FDR SYMBOL
# ============================================================

fdr_symbol <- function(x) {

  case_when(
    is.na(x)   ~ "",
    x < 0.001  ~ "***",
    x < 0.01   ~ "**",
    x < 0.05   ~ "*",
    x < 0.10   ~ "\u2020",
    TRUE       ~ ""
  )
}

# ============================================================
# LOAD COMPLETE 20B4 DE
# ============================================================

load_compartment_DE <- function(compartment) {

  out <- list()

  for (cc in contrast_order) {

    f <- file.path(
      DE_ROOT,
      compartment,
      "DE",
      paste0(
        "20B4_",
        compartment,
        "_",
        cc,
        "_full_DE.csv"
      )
    )

    stopifnot(
      file.exists(f)
    )

    z <- read_csv(
      f,
      show_col_types=FALSE
    ) %>%
      filter(
        publication_gene(gene)
      ) %>%
      mutate(
        abs_logFC=abs(logFC)
      )

    out[[cc]] <- z
  }

  bind_rows(out)
}

# ============================================================
# GENE SELECTION
#
# Eligibility:
# FDR <0.10
# |logFC| >=1
# logCPM >=1
#
# Selection:
# first pass: up to 2 genes from each contrast
# second pass: fill remaining positions from global ranking
#
# A gene appears only once per compartment.
# ============================================================

select_genes <- function(de) {

  eligible <- de %>%
    filter(
      FDR_BH < FDR_THRESHOLD,
      abs_logFC >= ABS_LOGFC_THRESHOLD,
      logCPM >= LOGCPM_THRESHOLD
    ) %>%
    arrange(
      FDR_BH,
      desc(abs_logFC),
      desc(logCPM),
      gene
    )

  if (nrow(eligible) == 0) {

    return(
      list(
        eligible=eligible,
        selected=tibble()
      )
    )
  }

  selected <- tibble(
    gene=character()
  )

  # ----------------------------------------------------------
  # Pass 1 — contrast diversity
  # max 2 unique genes per contrast
  # ----------------------------------------------------------

  for (cc in contrast_order) {

    z <- eligible %>%
      filter(
        contrast == cc,
        !gene %in% selected$gene
      ) %>%
      slice_head(
        n=2
      )

    if (nrow(z) > 0) {

      selected <- bind_rows(
        selected,
        z
      )
    }

    if (n_distinct(selected$gene) >= MAX_GENES) {
      break
    }
  }

  # ----------------------------------------------------------
  # Pass 2 — fill remaining positions globally
  # ----------------------------------------------------------

  if (n_distinct(selected$gene) < MAX_GENES) {

    remaining_n <-
      MAX_GENES -
      n_distinct(
        selected$gene
      )

    fill <- eligible %>%
      filter(
        !gene %in% selected$gene
      ) %>%
      group_by(
        gene
      ) %>%
      arrange(
        FDR_BH,
        desc(abs_logFC),
        .by_group=TRUE
      ) %>%
      slice_head(
        n=1
      ) %>%
      ungroup() %>%
      arrange(
        FDR_BH,
        desc(abs_logFC),
        desc(logCPM),
        gene
      ) %>%
      slice_head(
        n=remaining_n
      )

    selected <- bind_rows(
      selected,
      fill
    )
  }

  selected <- selected %>%
    distinct(
      gene,
      .keep_all=TRUE
    ) %>%
    slice_head(
      n=MAX_GENES
    ) %>%
    mutate(
      selection_rank=row_number()
    )

  list(
    eligible=eligible,
    selected=selected
  )
}

# ============================================================
# STORAGE
# ============================================================

selection_all <- list()
eligible_all <- list()
scores_all <- list()
summary_all <- list()
audit_all <- list()

# ============================================================
# MAIN LOOP
# ============================================================

for (ii in seq_len(nrow(compartments))) {

  compartment <- compartments$compartment[ii]
  module <- compartments$module[ii]

  cat("\n============================================================\n")
  cat(compartment, "\n")
  cat("============================================================\n")

  # ----------------------------------------------------------
  # Complete DE + gene selection
  # ----------------------------------------------------------

  de <- load_compartment_DE(
    compartment
  )

  sel <- select_genes(
    de
  )

  eligible <- sel$eligible
  selected <- sel$selected

  eligible_all[[compartment]] <- eligible

  cat(
    "Eligible DE rows:",
    nrow(eligible),
    "\n"
  )

  cat(
    "Unique eligible genes:",
    n_distinct(
      eligible$gene
    ),
    "\n"
  )

  cat(
    "Selected genes:",
    ifelse(
      nrow(selected) == 0,
      "NONE",
      paste(
        selected$gene,
        collapse=" | "
      )
    ),
    "\n"
  )

  # ----------------------------------------------------------
  # Record zero-signal compartments
  # ----------------------------------------------------------

  if (nrow(selected) == 0) {

    summary_all[[compartment]] <- tibble(
      compartment=compartment,
      n_HC=NA_integer_,
      n_FFA=NA_integer_,
      n_AA=NA_integer_,
      n_eligible_rows=nrow(eligible),
      n_eligible_genes=0L,
      n_selected=0L,
      figure_generated=FALSE
    )

    audit_all[[compartment]] <- tibble(
      compartment=compartment,
      status="NO_ELIGIBLE_GENES_FDR010"
    )

    next
  }

  # ----------------------------------------------------------
  # Read frozen L2 pseudobulk counts
  # ----------------------------------------------------------

  rds_file <- file.path(
    PB_DIR,
    "objects",
    compartments$rds[ii]
  )

  stopifnot(
    file.exists(rds_file)
  )

  mat <- readRDS(
    rds_file
  )

  # ----------------------------------------------------------
  # Module metadata
  # ----------------------------------------------------------

  meta <- meta_all %>%
    filter(
      module == .env$module
    )

  idx <- match(
    colnames(mat),
    meta$pseudobulk_id
  )

  if (anyNA(idx)) {

    stop(
      "Metadata mismatch in ",
      compartment
    )
  }

  meta <- meta[
    idx,
    ,
    drop=FALSE
  ]

  stopifnot(
    identical(
      colnames(mat),
      meta$pseudobulk_id
    )
  )

  # ----------------------------------------------------------
  # Reproduce 20B4 patient L1 support
  # ----------------------------------------------------------

  patient_info <- meta %>%
    group_by(
      patient_id,
      disease,
      batch
    ) %>%
    summarise(
      total_cells=sum(
        n_cells,
        na.rm=TRUE
      ),
      n_L2_present=n_distinct(
        final_L2_label
      ),
      .groups="drop"
    ) %>%
    filter(
      total_cells >= MIN_L1_CELLS
    ) %>%
    mutate(
      disease=factor(
        disease,
        levels=disease_order
      )
    ) %>%
    arrange(
      disease,
      patient_id
    )

  patient_order <- patient_info$patient_id

  # ----------------------------------------------------------
  # Reproduce 20B4 patient-count aggregation
  # ----------------------------------------------------------

  patient_counts <- do.call(
    cbind,
    lapply(
      patient_order,
      function(pid) {

        idx <- which(
          meta$patient_id == pid
        )

        Matrix::rowSums(
          mat[
            ,
            idx,
            drop=FALSE
          ]
        )
      }
    )
  )

  rownames(patient_counts) <- rownames(mat)
  colnames(patient_counts) <- patient_order

  patient_counts <- as.matrix(
    patient_counts
  )

  # ----------------------------------------------------------
  # EXACT 20B4 display normalization
  # ----------------------------------------------------------

  y_all <- DGEList(
    counts=patient_counts
  )

  y_all <- calcNormFactors(
    y_all,
    method="TMM"
  )

  logCPM <- cpm(
    y_all,
    log=TRUE,
    prior.count=2
  )

  # ----------------------------------------------------------
  # Check selected genes exist
  # ----------------------------------------------------------

  selected_genes_before <- as.character(
    selected$gene
  )

  genes_present <- selected_genes_before[
    selected_genes_before %in%
      rownames(logCPM)
  ]

  genes_missing <- setdiff(
    selected_genes_before,
    rownames(logCPM)
  )

  cat(
    "Selected genes present in logCPM:",
    length(genes_present),
    "/",
    length(selected_genes_before),
    "\n"
  )

  if (length(genes_missing) > 0) {

    cat(
      "Missing selected genes:",
      paste(
        genes_missing,
        collapse=" | "
      ),
      "\n"
    )
  }

  if (length(genes_present) == 0) {

    stop(
      "Selected genes absent from logCPM in ",
      compartment
    )
  }

  selected <- selected %>%
    filter(
      .data$gene %in%
        genes_present
    )

  stopifnot(
    nrow(selected) > 0,
    all(
      selected$gene %in%
        rownames(logCPM)
    )
  )

  # ----------------------------------------------------------
  # Expression long table
  # ----------------------------------------------------------

  expr <- as.data.frame(
    t(
      logCPM[
        selected$gene,
        ,
        drop=FALSE
      ]
    )
  )

  expr$patient_id <- rownames(expr)

  expr_long <- expr %>%
    pivot_longer(
      cols=-patient_id,
      names_to="gene",
      values_to="logCPM"
    ) %>%
    left_join(
      patient_info %>%
        mutate(
          disease=as.character(
            disease
          )
        ),
      by="patient_id"
    ) %>%
    mutate(
      disease=factor(
        disease,
        levels=disease_order
      ),
      compartment=compartment
    )

  scores_all[[compartment]] <- expr_long

  # ----------------------------------------------------------
  # All significant DE comparisons for selected genes
  # ----------------------------------------------------------

  selected_stats <- de %>%
    filter(
      gene %in% selected$gene,
      FDR_BH < FDR_THRESHOLD
    ) %>%
    mutate(
      symbol=fdr_symbol(
        FDR_BH
      )
    )

  # ----------------------------------------------------------
  # Directories
  # ----------------------------------------------------------

  COMP_DIR <- file.path(
    OUT,
    compartment
  )

  SVG_DIR <- file.path(
    COMP_DIR,
    "figures",
    "SVG"
  )

  PDF_DIR <- file.path(
    COMP_DIR,
    "figures",
    "PDF"
  )

  PNG_DIR <- file.path(
    COMP_DIR,
    "figures",
    "PNG"
  )

  TAB_DIR <- file.path(
    COMP_DIR,
    "tables"
  )

  for (d in c(
    SVG_DIR,
    PDF_DIR,
    PNG_DIR,
    TAB_DIR
  )) {

    dir.create(
      d,
      recursive=TRUE,
      showWarnings=FALSE
    )
  }

  # ----------------------------------------------------------
  # Save tables
  # ----------------------------------------------------------

  write_csv(
    selected,
    file.path(
      TAB_DIR,
      paste0(
        "20B8_",
        compartment,
        "_selected_genes.csv"
      )
    )
  )

  write_csv(
    selected_stats,
    file.path(
      TAB_DIR,
      paste0(
        "20B8_",
        compartment,
        "_selected_gene_DE_statistics.csv"
      )
    )
  )

  write_csv(
    expr_long,
    file.path(
      TAB_DIR,
      paste0(
        "20B8_",
        compartment,
        "_patient_logCPM.csv"
      )
    )
  )

  # ----------------------------------------------------------
  # Individual gene plots
  # ----------------------------------------------------------

  plot_list <- list()

  for (gene_i in selected$gene) {

    dd <- expr_long %>%
      filter(
        gene == gene_i
      )

    ss <- selected_stats %>%
      filter(
        gene == gene_i
      )

    yr <- range(
      dd$logCPM,
      na.rm=TRUE
    )

    span <- diff(yr)

    if (
      !is.finite(span) ||
      span == 0
    ) {
      span <- 1
    }

    # --------------------------------------------------------
    # Bracket coordinates
    # --------------------------------------------------------

    bracket_map <- tibble::tribble(
      ~contrast,    ~x1, ~x2, ~level,
      "FFA_vs_HC",    1,   2,      1,
      "AA_vs_HC",     1,   3,      2,
      "AA_vs_FFA",    2,   3,      3
    )

    ann <- bracket_map %>%
      inner_join(
        ss %>%
          select(
            contrast,
            FDR_BH,
            logFC,
            symbol
          ),
        by="contrast"
      ) %>%
      arrange(
        level
      )

    ymax <- max(
      dd$logCPM,
      na.rm=TRUE
    )

    if (nrow(ann) > 0) {

      ann <- ann %>%
        mutate(
          display_level=row_number(),
          y=ymax +
            span *
            (
              0.12 +
              0.12 *
              (
                display_level - 1
              )
            )
        )
    }

    # --------------------------------------------------------
    # Boxplot
    # --------------------------------------------------------

    p <- ggplot(
      dd,
      aes(
        x=disease,
        y=logCPM,
        fill=disease
      )
    ) +

      geom_boxplot(
        width=0.56,
        outlier.shape=NA,
        alpha=0.48,
        linewidth=0.55,
        colour="#333333"
      ) +

      geom_jitter(
        aes(
          colour=disease
        ),
        width=0.09,
        height=0,
        size=2.8,
        alpha=0.95,
        stroke=0.30,
        show.legend=FALSE
      ) +

      scale_fill_manual(
        values=DISEASE_COLORS,
        drop=FALSE
      ) +

      scale_colour_manual(
        values=DISEASE_COLORS
        ,
        drop=FALSE
      ) +

      labs(
        title=gene_i,
        x=NULL,
        y="TMM-normalized logCPM"
      ) +

      theme_classic(
        base_size=11
      ) +

      theme(
        legend.position="none",

        plot.title=element_text(
          face="bold",
          size=11.5,
          hjust=0.5
        ),

        axis.text.x=element_text(
          face="bold",
          size=9.5
        ),

        axis.text.y=element_text(
          size=8.5
        ),

        axis.title.y=element_text(
          size=9
        ),

        axis.line=element_line(
          linewidth=0.45,
          colour="#333333"
        ),

        plot.margin=margin(
          8,
          10,
          8,
          10
        )
      )

    # --------------------------------------------------------
    # Add only FDR <0.10 brackets
    # --------------------------------------------------------

    if (nrow(ann) > 0) {

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
            y=a$y + span * 0.022,
            label=a$symbol,
            size=3.6,
            fontface="bold",
            colour="#222222"
          )
      }
    }

    # --------------------------------------------------------
    # Y-axis expansion
    # --------------------------------------------------------

    if (nrow(ann) > 0) {

      upper <- max(
        ann$y,
        na.rm=TRUE
      ) + span * 0.12

    } else {

      upper <- ymax + span * 0.10
    }

    lower <- min(
      dd$logCPM,
      na.rm=TRUE
    ) - span * 0.08

    p <- p +
      coord_cartesian(
        ylim=c(
          lower,
          upper
        ),
        clip="off"
      )

    plot_list[[gene_i]] <- p
  }

  # ----------------------------------------------------------
  # Combined compartment panel
  # ----------------------------------------------------------

  n_selected <- length(
    plot_list
  )

  ncol_panel <- ifelse(
    n_selected <= 2,
    n_selected,
    3
  )

  nrow_panel <- ceiling(
    n_selected /
      ncol_panel
  )

  final <- wrap_plots(
    plot_list,
    ncol=ncol_panel,
    nrow=nrow_panel
  ) +

    plot_annotation(
      title=paste0(
        compartment,
        " — patient-level L1 pseudobulk expression"
      ),

      subtitle=paste0(
        "Genes selected from frozen 20B4 DE: ",
        "FDR<0.10, |logFC|>=1, logCPM>=1"
      ),

      caption=paste0(
        "Each point represents one patient. ",
        "Brackets reproduce frozen 20B4 pseudobulk comparisons; ",
        "*** FDR<0.001, ** FDR<0.01, * FDR<0.05, ",
        "\u2020 FDR 0.05-0.10."
      ),

      theme=theme(
        plot.title=element_text(
          face="bold",
          size=15,
          hjust=0.5
        ),

        plot.subtitle=element_text(
          size=9.5,
          hjust=0.5
        ),

        plot.caption=element_text(
          size=8,
          hjust=1
        )
      )
    )

  # ----------------------------------------------------------
  # Adaptive figure dimensions
  # ----------------------------------------------------------

  fig_width <- ifelse(
    ncol_panel == 1,
    4.2,
    ifelse(
      ncol_panel == 2,
      7.6,
      10.8
    )
  )

  fig_height <- ifelse(
    nrow_panel == 1,
    4.8,
    7.4
  )

  stem <- paste0(
    "20B8_",
    compartment,
    "_patient_L1_DE_gene_boxplots"
  )

  # ----------------------------------------------------------
  # SAVE SVG MASTER
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
    width=fig_width,
    height=fig_height,
    units="in",
    device=svglite
  )

  # ----------------------------------------------------------
  # SAVE PDF
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
    width=fig_width,
    height=fig_height,
    units="in"
  )

  # ----------------------------------------------------------
  # SAVE PNG
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
    width=fig_width,
    height=fig_height,
    units="in",
    dpi=400
  )

  cat(
    "  -> SVG/PDF/PNG generated\n"
  )

  # ----------------------------------------------------------
  # Store selected genes
  # ----------------------------------------------------------

  selection_all[[compartment]] <- selected %>%
    mutate(
      compartment=compartment
    )

  # ----------------------------------------------------------
  # Summary
  # ----------------------------------------------------------

  disease_counts <- patient_info %>%
    mutate(
      disease=as.character(
        disease
      )
    ) %>%
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

    as.integer(
      z[[1]]
    )
  }

  summary_all[[compartment]] <- tibble(
    compartment=compartment,

    n_HC=get_n("HC"),
    n_FFA=get_n("FFA"),
    n_AA=get_n("AA"),

    n_eligible_rows=nrow(
      eligible
    ),

    n_eligible_genes=n_distinct(
      eligible$gene
    ),

    n_selected=nrow(
      selected
    ),

    figure_generated=TRUE,

    svg=file.path(
      SVG_DIR,
      paste0(
        stem,
        ".svg"
      )
    )
  )

  audit_all[[compartment]] <- tibble(
    compartment=compartment,
    status="FIGURE_GENERATED"
  )
}

# ============================================================
# GLOBAL TABLES
# ============================================================

selection_global <- bind_rows(
  selection_all
)

eligible_global <- bind_rows(
  lapply(
    names(eligible_all),
    function(comp) {

      eligible_all[[comp]] %>%
        mutate(
          compartment_source=comp
        )
    }
  )
)

scores_global <- bind_rows(
  scores_all
)

summary_tbl <- bind_rows(
  summary_all
)

audit_tbl <- bind_rows(
  audit_all
)

write_csv(
  selection_global,
  file.path(
    GLOBAL_TAB,
    "20B8_ALL_selected_genes.csv"
  )
)

write_csv(
  eligible_global,
  file.path(
    GLOBAL_TAB,
    "20B8_ALL_eligible_DE_candidates.csv"
  )
)

write_csv(
  scores_global,
  file.path(
    GLOBAL_TAB,
    "20B8_ALL_patient_logCPM_selected_genes.csv"
  )
)

write_csv(
  summary_tbl,
  file.path(
    GLOBAL_TAB,
    "20B8_compartment_summary.csv"
  )
)

write_csv(
  audit_tbl,
  file.path(
    GLOBAL_TAB,
    "20B8_compartment_audit.csv"
  )
)

# ============================================================
# SELECTION LEDGER
#
# For every selected gene, retain ALL three 20B4 contrasts,
# not only significant ones. This makes the final biological
# pattern fully auditable.
# ============================================================

selection_ledger <- list()

if (nrow(selection_global) > 0) {

  for (comp in unique(
    selection_global$compartment
  )) {

    genes_comp <- selection_global %>%
      filter(
        compartment == comp
      ) %>%
      pull(gene)

    de_comp <- load_compartment_DE(
      comp
    )

    selection_ledger[[comp]] <- de_comp %>%
      filter(
        gene %in% genes_comp
      ) %>%
      select(
        compartment,
        contrast,
        numerator,
        denominator,
        batch_adjusted,
        gene,
        logFC,
        logCPM,
        PValue,
        FDR_BH
      ) %>%
      arrange(
        match(
          gene,
          genes_comp
        ),
        match(
          contrast,
          contrast_order
        )
      )
  }
}

selection_ledger_global <- bind_rows(
  selection_ledger
)

write_csv(
  selection_ledger_global,
  file.path(
    GLOBAL_TAB,
    "20B8_selected_gene_all_contrasts_ledger.csv"
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
# EXPECTED FIGURE COUNT
#
# Mural is expected to have zero eligible genes under the
# current frozen criteria. Figure count is therefore determined
# from summary_tbl rather than hard-coded.
# ============================================================

expected_figures <- sum(
  summary_tbl$figure_generated,
  na.rm=TRUE
)

# ============================================================
# FORMAL CHECKS
# ============================================================

checks <- tibble(
  check=c(
    "six_compartments_audited",
    "figure_count_matches_eligible_compartments",
    "SVG_count_correct",
    "PDF_count_correct",
    "PNG_count_correct",
    "all_selected_FDR_lt_010",
    "all_selected_abs_logFC_ge1",
    "all_selected_logCPM_ge1",
    "max_six_genes_per_compartment",
    "all_patient_support_ge20",
    "Mural_not_artificially_filled",
    "all_outputs_nonzero"
  ),

  pass=c(
    nrow(summary_tbl) == 6,

    length(svg_files) ==
      expected_figures,

    length(svg_files) ==
      expected_figures,

    length(pdf_files) ==
      expected_figures,

    length(png_files) ==
      expected_figures,

    if (
      nrow(selection_global) > 0
    ) {
      all(
        selection_global$FDR_BH <
          FDR_THRESHOLD
      )
    } else {
      TRUE
    },

    if (
      nrow(selection_global) > 0
    ) {
      all(
        abs(
          selection_global$logFC
        ) >=
          ABS_LOGFC_THRESHOLD
      )
    } else {
      TRUE
    },

    if (
      nrow(selection_global) > 0
    ) {
      all(
        selection_global$logCPM >=
          LOGCPM_THRESHOLD
      )
    } else {
      TRUE
    },

    all(
      summary_tbl$n_selected <=
        MAX_GENES
    ),

    if (
      nrow(scores_global) > 0
    ) {
      all(
        scores_global$total_cells >=
          MIN_L1_CELLS
      )
    } else {
      TRUE
    },

    summary_tbl %>%
      filter(
        compartment == "Mural"
      ) %>%
      summarise(
        ok=
          n_selected == 0 &&
          !figure_generated
      ) %>%
      pull(ok),

    if (
      length(
        c(
          svg_files,
          pdf_files,
          png_files
        )
      ) > 0
    ) {
      all(
        file.info(
          c(
            svg_files,
            pdf_files,
            png_files
          )
        )$size > 0
      )
    } else {
      FALSE
    }
  )
)

write_csv(
  checks,
  file.path(
    GLOBAL_TAB,
    "20B8_formal_checks.csv"
  )
)

# ============================================================
# FINAL REPORT
# ============================================================

cat("\n============================================================\n")
cat("20B8 — FINAL REPORT\n")
cat("============================================================\n")

cat(
  "\nSource DE: frozen 20B4 patient-L1 pseudobulk DE\n"
)

cat(
  "Display expression: TMM-normalized patient-L1 logCPM\n"
)

cat(
  "New DE testing performed: NO\n"
)

cat(
  "\nEligibility:\n"
)

cat(
  "  FDR <",
  FDR_THRESHOLD,
  "\n"
)

cat(
  "  |logFC| >=",
  ABS_LOGFC_THRESHOLD,
  "\n"
)

cat(
  "  logCPM >=",
  LOGCPM_THRESHOLD,
  "\n"
)

cat(
  "  minimum L1 cells per patient:",
  MIN_L1_CELLS,
  "\n"
)

cat(
  "  MIR/LINC/antisense excluded: YES\n"
)

cat(
  "\nSelection:\n"
)

cat(
  "  Maximum genes per compartment:",
  MAX_GENES,
  "\n"
)

cat(
  "  Contrast-diversity first pass: YES\n"
)

cat(
  "  Artificial filling with non-eligible genes: NO\n"
)

cat(
  "\nStatistics displayed:\n"
)

cat(
  "  Frozen 20B4 FDR only: YES\n"
)

cat(
  "  New Wilcoxon/t-test on plotted values: NO\n"
)

cat(
  "  *** FDR<0.001\n"
)

cat(
  "  **  FDR<0.01\n"
)

cat(
  "  *   FDR<0.05\n"
)

cat(
  "  \u2020   FDR 0.05-0.10\n"
)

cat(
  "\nFigures generated:",
  expected_figures,
  "\n"
)

cat(
  "SVG:",
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
  "\nCompartments:\n"
)

print(
  summary_tbl,
  n=Inf
)

cat(
  "\nSelected genes:\n"
)

if (nrow(selection_global) > 0) {

  print(
    selection_global %>%
      select(
        compartment,
        selection_rank,
        gene,
        contrast,
        logFC,
        logCPM,
        FDR_BH
      ) %>%
      arrange(
        compartment,
        selection_rank
      ),
    n=Inf
  )

} else {

  cat(
    "NONE\n"
  )
}

cat(
  "\n===== FORMAL CHECKS =====\n"
)

print(
  checks,
  n=Inf
)

cat(
  "\nOriginal 18B pseudobulks modified: NO\n"
)

cat(
  "Frozen 20B4 DE modified: NO\n"
)

cat(
  "Frozen L2 annotation modified: NO\n"
)

if (!all(checks$pass)) {

  cat(
    "\nSTATUS=FAIL\n"
  )

  print(
    checks %>%
      filter(
        !pass
      )
  )

  stop(
    "20B8 formal checks failed."
  )
}

cat(
  "\nSTATUS=PASS\n"
)
