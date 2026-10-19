#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Matrix)
  library(edgeR)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(readr)
  library(ggplot2)
  library(patchwork)
  library(svglite)
})

cat("============================================================\n")
cat("20B4 — MULTICOMPARTMENT PATIENT-L1 DE HEATMAPS\n")
cat("============================================================\n")

ROOT <- "PIPELINE_V2"

PB_ROOT <- file.path(
  ROOT,
  "results/disease_biology/18_PSEUDOBULK_DE/18B_pseudobulk/objects"
)

META_FILE <- file.path(
  PB_ROOT,
  "18B_pseudobulk_metadata_ALL.rds"
)

OUT <- file.path(
  ROOT,
  "results/publication_panels/20B4_all_L1_DE_heatmaps"
)

dir.create(OUT, recursive=TRUE, showWarnings=FALSE)

spec <- tribble(
  ~compartment,  ~module,               ~counts_file,
  "Epithelial",  "07_EPITHELIAL_L2",    "18B_07_EPITHELIAL_L2_pseudobulk_counts.rds",
  "Immune",      "08_IMMUNE_L2",        "18B_08_IMMUNE_L2_pseudobulk_counts.rds",
  "Fibroblast",  "09_FIBROBLAST_L2",    "18B_09_FIBROBLAST_L2_pseudobulk_counts.rds",
  "Endothelial", "10_VASCULAR_L2",      "18B_10_VASCULAR_L2_pseudobulk_counts.rds",
  "Mural",       "11_MURAL_L2",         "18B_11_MURAL_L2_pseudobulk_counts.rds",
  "Melanocyte",  "13_MELANOCYTE_L2",    "18B_13_MELANOCYTE_L2_pseudobulk_counts.rds"
)

contrast_spec <- tribble(
  ~contrast,    ~group1, ~group0, ~adjust_batch,
  "FFA_vs_HC",  "FFA",   "HC",    TRUE,
  "AA_vs_HC",   "AA",    "HC",    FALSE,
  "AA_vs_FFA",  "AA",    "FFA",   FALSE
)

DISEASE_COLS <- c(
  HC="#8F8D8E",
  FFA="#4D9E4F",
  AA="#A8142A"
)

pb_meta_all <- readRDS(META_FILE)

# ============================================================
# HELPERS
# ============================================================

select_top20 <- function(res, cc) {

  eligible <- res %>%
    filter(
      !grepl("^MIR",gene,ignore.case=TRUE),
      !grepl("^LINC",gene,ignore.case=TRUE),
      !grepl("-AS[0-9]+$",gene,ignore.case=TRUE)
    )

  up <- eligible %>%
    filter(logFC > 0) %>%
    arrange(PValue,desc(abs(logFC))) %>%
    slice_head(n=20) %>%
    mutate(
      direction="UP",
      rank_direction=row_number(),
      contrast=cc
    )

  down <- eligible %>%
    filter(logFC < 0) %>%
    arrange(PValue,desc(abs(logFC))) %>%
    slice_head(n=20) %>%
    mutate(
      direction="DOWN",
      rank_direction=row_number(),
      contrast=cc
    )

  stopifnot(
    nrow(up)==20,
    nrow(down)==20
  )

  bind_rows(up,down) %>%
    mutate(
      row_id=paste0(
        contrast,"__",
        direction,"__",
        sprintf("%02d",rank_direction),"__",
        gene
      )
    )
}

run_DE <- function(
  counts,
  md,
  cc,
  g1,
  g0,
  adjust_batch,
  compartment
) {

  md2 <- md %>%
    filter(
      disease %in% c(g0,g1),
      total_cells >= 20
    ) %>%
    mutate(
      disease_model=factor(
        disease,
        levels=c(g0,g1)
      )
    )

  n0 <- sum(md2$disease_model==g0)
  n1 <- sum(md2$disease_model==g1)

  if (n0 < 3 || n1 < 3) {
    return(NULL)
  }

  cts <- counts[
    ,
    md2$patient_id,
    drop=FALSE
  ]

  if (adjust_batch) {

    md2$batch_model <- factor(md2$batch)

    design <- model.matrix(
      ~ batch_model + disease_model,
      data=md2
    )

  } else {

    design <- model.matrix(
      ~ disease_model,
      data=md2
    )
  }

  stopifnot(
    qr(design)$rank==ncol(design)
  )

  coef_name <- paste0(
    "disease_model",
    g1
  )

  coef_idx <- which(
    colnames(design)==coef_name
  )

  stopifnot(length(coef_idx)==1)

  y <- DGEList(counts=cts)

  keep <- filterByExpr(
    y,
    design=design
  )

  stopifnot(sum(keep)>0)

  y <- y[
    keep,
    ,
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

  topTags(
    qlf,
    n=Inf,
    sort.by="PValue"
  )$table %>%
    rownames_to_column("gene") %>%
    as_tibble() %>%
    rename(FDR_BH=FDR) %>%
    mutate(
      contrast=cc,
      numerator=g1,
      denominator=g0,
      compartment=compartment,
      n_numerator=n1,
      n_denominator=n0,
      batch_adjusted=adjust_batch,
      n_genes_tested=sum(keep),
      .before=1
    )
}

# ============================================================
# COMPARTMENT LOOP
# ============================================================

global_summary <- list()

for (ii in seq_len(nrow(spec))) {

  compartment <- spec$compartment[ii]
  module <- spec$module[ii]

  cat("\n============================================================\n")
  cat(compartment,"\n")
  cat("============================================================\n")

  COUT <- file.path(
    OUT,
    compartment
  )

  DE_DIR  <- file.path(COUT,"DE")
  SEL_DIR <- file.path(COUT,"selected_genes")
  SRC_DIR <- file.path(COUT,"source_data")
  TAB_DIR <- file.path(COUT,"tables")
  SVG_DIR <- file.path(COUT,"figures/SVG")
  PDF_DIR <- file.path(COUT,"figures/PDF")
  PNG_DIR <- file.path(COUT,"figures/PNG")

  for (d in c(
    DE_DIR,SEL_DIR,SRC_DIR,TAB_DIR,
    SVG_DIR,PDF_DIR,PNG_DIR
  )) {
    dir.create(
      d,
      recursive=TRUE,
      showWarnings=FALSE
    )
  }

  mat <- readRDS(
    file.path(
      PB_ROOT,
      spec$counts_file[ii]
    )
  )

  meta <- pb_meta_all %>%
    filter(module==!!module)

  meta <- meta[
    match(
      colnames(mat),
      meta$pseudobulk_id
    ),
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
  # L2 -> patient L1
  # ----------------------------------------------------------

  patient_info <- meta %>%
    group_by(
      patient_id
    ) %>%
    summarise(
      disease=first(disease),
      batch=first(batch),
      total_cells=sum(n_cells),
      n_L2_present=n(),
      .groups="drop"
    ) %>%
    filter(
      total_cells >= 20
    ) %>%
    mutate(
      disease=factor(
        disease,
        levels=c("HC","FFA","AA")
      )
    ) %>%
    arrange(
      disease,
      patient_id
    )

  patient_order <- patient_info$patient_id

  patient_counts <- do.call(
    cbind,
    lapply(
      patient_order,
      function(pid) {

        idx <- which(
          meta$patient_id==pid
        )

        Matrix::rowSums(
          mat[,idx,drop=FALSE]
        )
      }
    )
  )

  rownames(patient_counts) <- rownames(mat)
  colnames(patient_counts) <- patient_order
  patient_counts <- as.matrix(patient_counts)

  # ----------------------------------------------------------
  # Display normalization across eligible patients
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

  write_csv(
    patient_info,
    file.path(
      TAB_DIR,
      paste0(
        "20B4_",
        compartment,
        "_patient_metadata.csv"
      )
    )
  )

  # ----------------------------------------------------------
  # DE
  # ----------------------------------------------------------

  DE <- list()

  for (kk in seq_len(nrow(contrast_spec))) {

    cc <- contrast_spec[kk,]

    cat(
      "  ",
      cc$contrast,
      " | batch adjusted=",
      cc$adjust_batch,
      "\n",
      sep=""
    )

    res <- run_DE(
      patient_counts,
      patient_info,
      cc$contrast,
      cc$group1,
      cc$group0,
      cc$adjust_batch,
      compartment
    )

    stopifnot(!is.null(res))

    DE[[cc$contrast]] <- res

    write_csv(
      res,
      file.path(
        DE_DIR,
        paste0(
          "20B4_",
          compartment,"_",
          cc$contrast,
          "_full_DE.csv"
        )
      )
    )
  }

  # ----------------------------------------------------------
  # TOP20
  # ----------------------------------------------------------

  selected <- setNames(
    lapply(
      names(DE),
      function(cc) {
        select_top20(
          DE[[cc]],
          cc
        )
      }
    ),
    names(DE)
  )

  selected_all <- bind_rows(selected)

  stopifnot(
    nrow(selected_all)==120,
    !any(
      grepl(
        "^MIR|^LINC|-AS[0-9]+$",
        selected_all$gene,
        ignore.case=TRUE
      )
    )
  )

  write_csv(
    selected_all,
    file.path(
      SEL_DIR,
      paste0(
        "20B4_",
        compartment,
        "_horizontal_Top20_selected_genes.csv"
      )
    )
  )

  # ----------------------------------------------------------
  # MATRIX
  # ----------------------------------------------------------

  build_matrix <- function(sel) {

    z <- do.call(
      rbind,
      lapply(
        sel$gene,
        function(g) {

          v <- as.numeric(
            logCPM[
              g,
              patient_order
            ]
          )

          zz <- as.numeric(
            scale(v)
          )

          zz[!is.finite(zz)] <- 0

          pmax(
            pmin(zz,2.5),
            -2.5
          )
        }
      )
    )

    rownames(z) <- sel$row_id
    colnames(z) <- patient_order

    z
  }

  matrices <- lapply(
    selected,
    build_matrix
  )

  # ----------------------------------------------------------
  # GROUP BOUNDARIES
  # ----------------------------------------------------------

  n_HC <- sum(
    patient_info$disease=="HC"
  )

  n_FFA <- sum(
    patient_info$disease=="FFA"
  )

  n_AA <- sum(
    patient_info$disease=="AA"
  )

  boundaries <- c(
    n_HC + 0.5,
    n_HC + n_FFA + 0.5
  )

  # ----------------------------------------------------------
  # PANEL FUNCTION
  # ----------------------------------------------------------

  make_panel <- function(
    cc,
    title_text
  ) {

    sel <- selected[[cc]]
    z <- matrices[[cc]]

    dat <- as.data.frame(z) %>%
      rownames_to_column("row_id") %>%
      pivot_longer(
        -row_id,
        names_to="patient_id",
        values_to="z"
      ) %>%
      mutate(
        patient_id=factor(
          patient_id,
          levels=patient_order
        ),
        row_id=factor(
          row_id,
          levels=rev(sel$row_id)
        )
      )

    labels <- setNames(
      sel$gene,
      sel$row_id
    )

    ann <- patient_info %>%
      mutate(
        patient_id=factor(
          patient_id,
          levels=patient_order
        ),
        y=1
      )

    p_ann <- ggplot(
      ann,
      aes(
        patient_id,
        y,
        fill=disease
      )
    ) +
      geom_tile() +
      scale_fill_manual(
        values=DISEASE_COLS
      ) +
      scale_x_discrete(
        expand=c(0,0)
      ) +
      theme_void() +
      theme(
        legend.position="none",
        plot.margin=margin(0,0,0,0)
      )

    p <- ggplot(
      dat,
      aes(
        patient_id,
        row_id,
        fill=z
      )
    ) +

      geom_tile() +

      geom_vline(
        xintercept=boundaries,
        colour="white",
        linewidth=.45
      ) +

      geom_hline(
        yintercept=20.5,
        colour="black",
        linewidth=.65
      ) +

      scale_fill_gradient2(
        low="#2166AC",
        mid="#F7F7F7",
        high="#B2182B",
        midpoint=0,
        limits=c(-2.5,2.5),
        name="Row z-score"
      ) +

      scale_x_discrete(
        expand=c(0,0)
      ) +

      scale_y_discrete(
        labels=labels,
        expand=c(0,0)
      ) +

      coord_fixed(
        ratio=1
      ) +

      labs(
        title=title_text,
        x=NULL,
        y=NULL
      ) +

      theme_minimal(
        base_size=7
      ) +

      theme(
        panel.grid=element_blank(),

        axis.text.x=element_text(
          angle=45,
          hjust=1,
          size=5.2
        ),

        axis.text.y=element_text(
          size=5.1,
          colour="black"
        ),

        plot.title=element_text(
          face="bold",
          hjust=.5,
          size=9
        ),

        legend.position="bottom",

        plot.margin=margin(
          0,0,0,0
        )
      )

    p_ann / p +
      plot_layout(
        heights=c(.025,.975)
      )
  }

  p1 <- make_panel(
    "FFA_vs_HC",
    "FFA vs HC"
  )

  p2 <- make_panel(
    "AA_vs_HC",
    "AA vs HC"
  )

  p3 <- make_panel(
    "AA_vs_FFA",
    "AA vs FFA"
  )

  final <- (
    p1 | p2 | p3
  ) +
    plot_layout(
      widths=c(1,1,1),
      guides="collect"
    ) &
    theme(
      legend.position="bottom"
    )

  final <- final +
    plot_annotation(
      title=paste0(
        compartment,
        " - differential expression"
      ),
      subtitle=paste0(
        "Top 20 UP and Top 20 DOWN per contrast; ",
        "MIR*, LINC* and -AS# genes excluded"
      ),
      theme=theme(
        plot.title=element_text(
          face="bold",
          hjust=.5,
          size=11
        ),
        plot.subtitle=element_text(
          hjust=.5,
          size=7
        )
      )
    )

  stem <- paste0(
    "20B4_",
    compartment,
    "_horizontal_three_contrasts_Top20_FILTERED"
  )

  ggsave(
    file.path(
      SVG_DIR,
      paste0(stem,".svg")
    ),
    final,
    width=8.7,
    height=6.3,
    units="in",
    device=svglite
  )

  ggsave(
    file.path(
      PDF_DIR,
      paste0(stem,".pdf")
    ),
    final,
    width=8.7,
    height=6.3,
    units="in"
  )

  ggsave(
    file.path(
      PNG_DIR,
      paste0(stem,".png")
    ),
    final,
    width=8.7,
    height=6.3,
    units="in",
    dpi=400
  )

  global_summary[[compartment]] <- tibble(
    compartment=compartment,
    n_HC=n_HC,
    n_FFA=n_FFA,
    n_AA=n_AA,
    n_patients=nrow(patient_info),
    n_selected_genes=120,
    svg=file.path(
      SVG_DIR,
      paste0(stem,".svg")
    )
  )

  cat(
    "  -> SVG generated\n"
  )
}

#
# ============================================================
# GLOBAL SUMMARY
# ============================================================

cat("\n============================================================\n")
cat("20B4 — GLOBAL SUMMARY\n")
cat("============================================================\n")

summary_tbl <- bind_rows(
  global_summary
)

stopifnot(
  nrow(summary_tbl) == 6,
  n_distinct(summary_tbl$compartment) == 6
)

write_csv(
  summary_tbl,
  file.path(
    OUT,
    "20B4_multicompartment_summary.csv"
  )
)

print(
  summary_tbl,
  n=Inf
)

# ============================================================
# FORMAL OUTPUT CHECKS
# ============================================================

cat("\n===== FORMAL CHECKS =====\n")

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

de_files <- list.files(
  OUT,
  pattern="_full_DE\\.csv$",
  recursive=TRUE,
  full.names=TRUE
)

selected_files <- list.files(
  OUT,
  pattern="_horizontal_Top20_selected_genes\\.csv$",
  recursive=TRUE,
  full.names=TRUE
)

checks <- tibble(
  check=c(
    "six_compartments_completed",
    "six_SVG_files",
    "six_PDF_files",
    "six_PNG_files",
    "eighteen_DE_results",
    "six_selected_gene_tables",
    "all_outputs_nonzero",
    "all_compartments_HC_ge3",
    "all_compartments_FFA_ge3",
    "all_compartments_AA_ge3",
    "all_gene_panels_120_rows"
  ),
  pass=c(
    nrow(summary_tbl)==6,
    length(svg_files)==6,
    length(pdf_files)==6,
    length(png_files)==6,
    length(de_files)==18,
    length(selected_files)==6,

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
    all(summary_tbl$n_selected_genes == 120)
  )
)

print(
  checks,
  n=Inf
)

write_csv(
  checks,
  file.path(
    OUT,
    "20B4_formal_checks.csv"
  )
)

if (!all(checks$pass)) {

  stop(
    "20B4 formal checks failed: ",
    paste(
      checks$check[
        !checks$pass
      ],
      collapse=", "
    )
  )
}

# ============================================================
# ANALYSIS SPECIFICATION
# ============================================================

analysis_specification <- tibble(
  item=c(
    "aggregation",
    "minimum_L1_cells_per_patient",
    "contrasts",
    "FFA_vs_HC_model",
    "AA_vs_HC_model",
    "AA_vs_FFA_model",
    "gene_filtering_for_DE",
    "normalization",
    "dispersion",
    "DE_test",
    "multiple_testing",
    "publication_gene_filter",
    "publication_ranking",
    "genes_per_direction",
    "heatmap_expression",
    "heatmap_scaling",
    "heatmap_clipping",
    "column_clustering",
    "column_order",
    "heatmap_layout",
    "master_vector_format"
  ),
  value=c(
    "L2 pseudobulk counts summed within patient to L1 compartment",
    "20",
    "FFA_vs_HC; AA_vs_HC; AA_vs_FFA",
    "~ batch_model + disease_model",
    "~ disease_model; batch not identifiable",
    "~ disease_model; batch severely confounded",
    "edgeR filterByExpr within L1 x contrast",
    "TMM",
    "estimateDisp robust=TRUE",
    "glmQLFit robust=TRUE; glmQLFTest",
    "BH within L1 x contrast",
    "exclude MIR*, LINC*, and genes ending -AS# before display ranking",
    "PValue ascending within positive/negative logFC; abs(logFC) tie-break",
    "Top20 UP + Top20 DOWN per contrast",
    "TMM-normalized patient-level L1 logCPM",
    "row z-score across eligible patients",
    "[-2.5, 2.5]",
    "NO",
    "HC then FFA then AA; patient ID within disease",
    "three horizontal contrast heatmaps",
    "SVG"
  )
)

write_csv(
  analysis_specification,
  file.path(
    OUT,
    "20B4_analysis_specification.csv"
  )
)

# ============================================================
# OUTPUT INVENTORY
# ============================================================

inventory <- tibble(
  file=c(
    svg_files,
    pdf_files,
    png_files,
    de_files,
    selected_files
  )
) %>%
  mutate(
    bytes=file.info(file)$size,
    extension=tolower(
      tools::file_ext(file)
    )
  )

write_csv(
  inventory,
  file.path(
    OUT,
    "20B4_output_inventory.csv"
  )
)

# ============================================================
# HASHES
# ============================================================

hash_targets <- inventory$file[
  file.exists(
    inventory$file
  )
]

hashes <- tools::md5sum(
  hash_targets
)

hash_tbl <- tibble(
  file=names(hashes),
  md5=unname(hashes)
)

write_csv(
  hash_tbl,
  file.path(
    OUT,
    "20B4_output_hashes.csv"
  )
)

# ============================================================
# FINAL REPORT
# ============================================================

cat("\n============================================================\n")
cat("20B4 — FINAL REPORT\n")
cat("============================================================\n")

cat("\nCompartments generated:", nrow(summary_tbl), "\n")
cat("Contrasts per compartment: 3\n")
cat("Top-ranked genes per contrast: 20 UP + 20 DOWN\n")
cat("Displayed gene positions per compartment: 120\n")

cat("\nSVG:", length(svg_files), "\n")
cat("PDF:", length(pdf_files), "\n")
cat("PNG:", length(png_files), "\n")
cat("DE tables:", length(de_files), "\n")

cat("\nPublication filters:\n")
cat("  MIR*: excluded\n")
cat("  LINC*: excluded\n")
cat("  -AS#: excluded\n")

cat("\nModels:\n")
cat("  FFA vs HC: batch + disease\n")
cat("  AA vs HC: disease only\n")
cat("  AA vs FFA: disease only\n")

cat("\nDisplay:\n")
cat("  Three heatmaps horizontally\n")
cat("  Cells square\n")
cat("  Minimal inter-heatmap spacing\n")
cat("  Columns clustered: NO\n")
cat("  Patient order: HC -> FFA -> AA\n")
cat("  SVG master: YES\n")

cat("\nCompartments:\n")
print(
  summary_tbl %>%
    select(
      compartment,
      n_HC,
      n_FFA,
      n_AA,
      n_patients
    ),
  n=Inf
)

cat("\nSweat gland included: NO — deferred low support\n")
cat("Sebocyte included: NO — FFA n=2 after L1 >=20 filter\n")

cat("\nOriginal 18B pseudobulks modified: NO\n")
cat("Frozen L2 annotation modified: NO\n")
cat("L2-level 18C DE modified: NO\n")
cat("New L1-level DE performed: YES\n")

cat("\nSTATUS=PASS\n")
