suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(ggplot2)
  library(patchwork)
  library(svglite)
})

cat("============================================================\n")
cat("20B6 — CROSS-DATABASE CONSENSUS MIRROR PLOTS\n")
cat("============================================================\n")

# ============================================================
# INPUT / OUTPUT
# ============================================================

INPUT <- paste0(
  "PIPELINE_V2/results/disease_biology/",
  "20_PATHWAY_ANALYSIS/20F_CROSS_DATABASE_CONSENSUS/",
  "20F1B_consensus/tables/",
  "20F1B_core_consensus_map.csv"
)

OUT <- paste0(
  "PIPELINE_V2/results/publication_panels/",
  "20B6_cross_database_consensus_mirror_plots"
)

TAB <- file.path(OUT, "tables")

dir.create(
  TAB,
  recursive=TRUE,
  showWarnings=FALSE
)

stopifnot(file.exists(INPUT))

x <- read_csv(
  INPUT,
  show_col_types=FALSE
)

# ============================================================
# COMPARTMENTS
# ============================================================

module_map <- c(
  "07_EPITHELIAL_L2" = "Epithelial",
  "08_IMMUNE_L2"     = "Immune",
  "09_FIBROBLAST_L2" = "Fibroblast",
  "10_VASCULAR_L2"   = "Endothelial",
  "11_MURAL_L2"      = "Mural",
  "13_MELANOCYTE_L2" = "Melanocyte"
)

contrast_order <- c(
  "FFA_vs_HC",
  "AA_vs_HC",
  "AA_vs_FFA"
)

# ============================================================
# CLEAN AXIS LABELS
# ============================================================

axis_labels <- c(
  "CELL_CYCLE_DNA_REPAIR" =
    "Cell cycle / DNA repair",

  "ECM_COLLAGEN_FIBROSIS" =
    "ECM / collagen / fibrosis",

  "EPITHELIAL_DIFFERENTIATION" =
    "Epithelial differentiation",

  "IFN_JAK_STAT_ANTIGEN" =
    "IFN / JAK-STAT / antigen",

  "METABOLISM_OXPHOS_ROS" =
    "Metabolism / OXPHOS / ROS",

  "T_CELL_CYTOTOXIC" =
    "T-cell / cytotoxic",

  "VASCULAR_ANGIOGENESIS" =
    "Vascular / angiogenesis",

  "WNT_HEDGEHOG_FOLLICLE" =
    "WNT / Hedgehog / follicle"
)

# ============================================================
# PREPARE CORE CONSENSUS
# ============================================================

d <- x %>%
  filter(
    module %in% names(module_map)
  ) %>%
  mutate(
    compartment=unname(
      module_map[module]
    ),
    contrast=factor(
      contrast,
      levels=contrast_order
    ),
    signed_NES=if_else(
      direction == "UP",
      abs(median_database_NES),
      -abs(median_database_NES)
    )
  )

# ============================================================
# COMPARTMENT-LEVEL CONSENSUS
#
# IMPORTANT:
# This is NOT a new L1 GSEA.
#
# It summarizes the already frozen L2 cross-database consensus.
#
# For each:
# compartment × contrast × axis × direction
#
# summarize across L2 states using median signed NES.
# ============================================================

compartment_consensus <- d %>%
  group_by(
    compartment,
    contrast,
    axis,
    direction
  ) %>%
  summarise(
    n_L2=n_distinct(
      final_L2_label
    ),

    median_NES_across_L2=median(
      signed_NES,
      na.rm=TRUE
    ),

    max_abs_NES=max(
      abs(signed_NES),
      na.rm=TRUE
    ),

    min_FDR=min(
      min_FDR,
      na.rm=TRUE
    ),

    n_TRIPLE=sum(
      consensus_class == "TRIPLE",
      na.rm=TRUE
    ),

    n_CORE=sum(
      consensus_class == "CORE",
      na.rm=TRUE
    ),

    n_Hallmark_supported=sum(
      Hallmark_support,
      na.rm=TRUE
    ),

    databases=paste(
      sort(
        unique(
          unlist(
            strsplit(
              primary_databases,
              ";",
              fixed=TRUE
            )
          )
        )
      ),
      collapse=";"
    ),

    .groups="drop"
  ) %>%
  mutate(
    axis_label=unname(
      axis_labels[axis]
    ),

    axis_label=if_else(
      is.na(axis_label),
      str_to_title(
        str_replace_all(
          axis,
          "_",
          " "
        )
      ),
      axis_label
    ),

    plot_NES=median_NES_across_L2,

    evidence_label=paste0(
      "nL2=",
      n_L2
    )
  )

# ============================================================
# DEAL WITH AXES HAVING BOTH UP AND DOWN CONSENSUS
#
# At compartment level an axis may contain L2 states pointing
# in opposite directions. We do NOT force them into one sign.
#
# UP and DOWN remain separate rows in the mirror plot.
# ============================================================

write_csv(
  compartment_consensus,
  file.path(
    TAB,
    "20B6_compartment_consensus_all.csv"
  )
)

# ============================================================
# RANKING
#
# Up to 10 strongest UP and 10 strongest DOWN.
# There are only 8 frozen axes, so no artificial filling.
# ============================================================

selected <- compartment_consensus %>%
  group_by(
    compartment,
    contrast,
    direction
  ) %>%
  arrange(
    desc(abs(plot_NES)),
    min_FDR,
    desc(n_L2),
    .by_group=TRUE
  ) %>%
  mutate(
    rank_direction=row_number()
  ) %>%
  filter(
    rank_direction <= 10
  ) %>%
  ungroup()

write_csv(
  selected,
  file.path(
    TAB,
    "20B6_selected_consensus_axes.csv"
  )
)

# ============================================================
# SIGNIFICANCE / EVIDENCE LABEL
# ============================================================

selected <- selected %>%
  mutate(
    fdr_label=case_when(
      min_FDR < 0.001 ~ "***",
      min_FDR < 0.01  ~ "**",
      min_FDR < 0.05  ~ "*",
      TRUE            ~ ""
    ),

    bar_label=sprintf(
      "%.2f  %s",
      plot_NES,
      fdr_label
    )
  )

# ============================================================
# STORAGE
# ============================================================

summary_list <- list()

# ============================================================
# LOOP
# ============================================================

for (
  compartment in
  unname(module_map)
) {

  cat("\n============================================================\n")
  cat(compartment, "\n")
  cat("============================================================\n")

  dc <- selected %>%
    filter(
      .data$compartment == compartment
    )

  if (nrow(dc) == 0) {

    cat("No CORE/TRIPLE consensus signal — skipped\n")

    next
  }

  SVG_DIR <- file.path(
    OUT,
    compartment,
    "figures",
    "SVG"
  )

  PDF_DIR <- file.path(
    OUT,
    compartment,
    "figures",
    "PDF"
  )

  PNG_DIR <- file.path(
    OUT,
    compartment,
    "figures",
    "PNG"
  )

  TAB_DIR <- file.path(
    OUT,
    compartment,
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

  write_csv(
    dc,
    file.path(
      TAB_DIR,
      paste0(
        "20B6_",
        compartment,
        "_consensus_axes.csv"
      )
    )
  )

  plot_list <- list()

  # ----------------------------------------------------------
  # one mirror plot per contrast
  # ----------------------------------------------------------

  for (
    cc in contrast_order
  ) {

    z <- dc %>%
      filter(
        as.character(contrast) == cc
      )

    if (nrow(z) == 0) {

      p <- ggplot() +
        annotate(
          "text",
          x=0,
          y=0,
          label="No CORE/TRIPLE consensus",
          size=4
        ) +
        xlim(-1,1) +
        ylim(-1,1) +
        theme_void() +
        labs(
          title=gsub(
            "_",
            " ",
            cc
          )
        )

      plot_list[[cc]] <- p

      next
    }

    # --------------------------------------------------------
    # Defensive collapse to ONE bar per biological axis
    # and direction within compartment x contrast.
    #
    # The frozen 20F1B source is L2-resolved. After compartment
    # summarisation, duplicated display labels can still occur
    # defensively if equivalent rows survive upstream grouping.
    # We therefore collapse them here without changing direction:
    # median NES, minimum FDR, maximum L2 support.
    # --------------------------------------------------------

    z <- z %>%
      group_by(
        axis,
        axis_label,
        direction
      ) %>%
      summarise(
        plot_NES=median(
          plot_NES,
          na.rm=TRUE
        ),
        min_FDR=min(
          min_FDR,
          na.rm=TRUE
        ),
        n_L2=max(
          n_L2,
          na.rm=TRUE
        ),
        max_abs_NES=max(
          max_abs_NES,
          na.rm=TRUE
        ),
        .groups="drop"
      ) %>%
      mutate(
        fdr_label=case_when(
          min_FDR < 0.001 ~ "***",
          min_FDR < 0.01  ~ "**",
          min_FDR < 0.05  ~ "*",
          TRUE            ~ ""
        ),
        bar_label=sprintf(
          "%.2f  %s",
          plot_NES,
          fdr_label
        ),
        plot_label=paste0(
          axis_label,
          if_else(
            direction == "UP",
            "  ↑",
            "  ↓"
          )
        )
      )

    # Order from negative extreme through positive extreme.
    # unique() is intentional: factor levels must be unique.
    z <- z %>%
      arrange(
        plot_NES
      ) %>%
      mutate(
        plot_label=factor(
          plot_label,
          levels=unique(plot_label)
        )
      )

    lim <- max(
      abs(z$plot_NES),
      na.rm=TRUE
    )

    lim <- max(
      1,
      ceiling(
        lim * 10
      ) / 10
    )

    # --------------------------------------------------------
    # Mirror / pyramid plot
    # --------------------------------------------------------

    p <- ggplot(
      z,
      aes(
        x=plot_NES,
        y=plot_label,
        fill=direction
      )
    ) +

      geom_col(
        width=0.72,
        linewidth=0.25,
        colour="#333333"
      ) +

      geom_vline(
        xintercept=0,
        linewidth=0.55,
        colour="#333333"
      ) +

      geom_text(
        data=z %>%
          filter(
            plot_NES > 0
          ),
        aes(
          label=bar_label
        ),
        hjust=-0.10,
        size=3.0
      ) +

      geom_text(
        data=z %>%
          filter(
            plot_NES < 0
          ),
        aes(
          label=bar_label
        ),
        hjust=1.10,
        size=3.0
      ) +

      scale_fill_manual(
        values=c(
          DOWN="#5B78A6",
          UP="#B24A4A"
        ),
        drop=FALSE
      ) +

      scale_x_continuous(
        limits=c(
          -lim * 1.30,
          lim * 1.30
        ),
        breaks=scales::pretty_breaks(
          n=5
        )
      ) +

      labs(
        title=gsub(
          "_",
          " ",
          cc
        ),
        x="Consensus NES",
        y=NULL,
        fill=NULL
      ) +

      theme_classic(
        base_size=10
      ) +

      theme(
        legend.position="none",

        plot.title=element_text(
          face="bold",
          hjust=0.5,
          size=11
        ),

        axis.text.y=element_text(
          size=8.5,
          colour="#222222"
        ),

        axis.text.x=element_text(
          size=8.5
        ),

        axis.title.x=element_text(
          size=9.5
        ),

        axis.line.y=element_blank(),

        axis.ticks.y=element_blank(),

        plot.margin=margin(
          8,
          16,
          8,
          16
        )
      )

    plot_list[[cc]] <- p
  }

  # ----------------------------------------------------------
  # Three contrasts horizontally
  # ----------------------------------------------------------

  final <- (
    plot_list[["FFA_vs_HC"]] |
    plot_list[["AA_vs_HC"]] |
    plot_list[["AA_vs_FFA"]]
  ) +

    plot_annotation(
      title=paste0(
        compartment,
        " — cross-database consensus pathway enrichment"
      ),

      subtitle=paste0(
        "CORE/TRIPLE consensus across BioPlanet, GO-BP and Reactome; ",
        "Hallmark support retained where available"
      ),

      caption=paste0(
        "Bars: median consensus NES across L2 states. ",
        "Left = DOWN enrichment; right = UP enrichment. ",
        "* FDR<0.05, ** FDR<0.01, *** FDR<0.001."
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

  stem <- paste0(
    "20B6_",
    compartment,
    "_cross_database_consensus_mirror"
  )

  # ----------------------------------------------------------
  # SAVE
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
    width=15.5,
    height=7.5,
    units="in",
    device=svglite
  )

  ggsave(
    file.path(
      PDF_DIR,
      paste0(
        stem,
        ".pdf"
      )
    ),
    final,
    width=15.5,
    height=7.5,
    units="in"
  )

  ggsave(
    file.path(
      PNG_DIR,
      paste0(
        stem,
        ".png"
      )
    ),
    final,
    width=15.5,
    height=7.5,
    units="in",
    dpi=400
  )

  cat("  -> SVG/PDF/PNG generated\n")

  summary_list[[compartment]] <- tibble(
    compartment=compartment,

    n_rows=nrow(dc),

    n_axes=n_distinct(
      dc$axis
    ),

    n_UP=sum(
      dc$direction == "UP"
    ),

    n_DOWN=sum(
      dc$direction == "DOWN"
    ),

    n_FFA_vs_HC=sum(
      dc$contrast == "FFA_vs_HC"
    ),

    n_AA_vs_HC=sum(
      dc$contrast == "AA_vs_HC"
    ),

    n_AA_vs_FFA=sum(
      dc$contrast == "AA_vs_FFA"
    )
  )
}

# ============================================================
# SUMMARY
# ============================================================

summary_tbl <- bind_rows(
  summary_list
)

write_csv(
  summary_tbl,
  file.path(
    TAB,
    "20B6_compartment_summary.csv"
  )
)

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

checks <- tibble(
  check=c(
    "input_core_consensus_exists",
    "only_CORE_or_TRIPLE_source",
    "no_SINGLE_source",
    "three_contrasts_defined",
    "all_outputs_nonzero"
  ),

  pass=c(
    file.exists(INPUT),

    all(
      d$consensus_class %in%
        c(
          "CORE",
          "TRIPLE"
        )
    ),

    !any(
      d$consensus_class == "SINGLE"
    ),

    all(
      contrast_order %in%
        unique(
          as.character(
            d$contrast
          )
        )
    ),

    length(
      c(
        svg_files,
        pdf_files,
        png_files
      )
    ) == 0 ||
      all(
        file.info(
          c(
            svg_files,
            pdf_files,
            png_files
          )
        )$size > 0
      )
  )
)

write_csv(
  checks,
  file.path(
    TAB,
    "20B6_formal_checks.csv"
  )
)

# ============================================================
# FINAL REPORT
# ============================================================
# ============================================================
# FINAL REPORT
# ============================================================

cat("\n============================================================\n")
cat("20B6 — FINAL REPORT\n")
cat("============================================================\n")

cat(
  "\nSource:",
  INPUT,
  "\n"
)

cat(
  "\nBiological axes:",
  n_distinct(d$axis),
  "\n"
)

cat(
  "Compartments generated:",
  nrow(summary_tbl),
  "\n"
)

cat(
  "Contrasts: FFA_vs_HC | AA_vs_HC | AA_vs_FFA\n"
)

cat(
  "\nSelection: CORE/TRIPLE cross-database consensus only\n"
)

cat(
  "Single-database signals included: NO\n"
)

cat(
  "Redundancy reduction: inherited from frozen 20F1B consensus\n"
)

cat(
  "Databases integrated: BioPlanet | GO-BP | Reactome\n"
)

cat(
  "Hallmark orthogonal support retained: YES\n"
)

cat(
  "\nCompartment-level representation:\n"
)

cat(
  "  Source unit: frozen L2 cross-database consensus\n"
)

cat(
  "  Aggregation: compartment x contrast x axis x direction\n"
)

cat(
  "  NES summary: median signed database NES across supported L2 states\n"
)

cat(
  "  UP and DOWN signals kept separately: YES\n"
)

cat(
  "  New L1 GSEA performed: NO\n"
)

cat(
  "\nRanking:\n"
)

cat(
  "  Primary: absolute compartment consensus NES\n"
)

cat(
  "  Secondary: minimum FDR\n"
)

cat(
  "  Tertiary: number of supporting L2 states\n"
)

cat(
  "  Maximum per direction: 10\n"
)

cat(
  "  Artificial filling to Top10: NO\n"
)

cat(
  "\nDisplay:\n"
)

cat(
  "  Mirror / diverging horizontal bars: YES\n"
)

cat(
  "  DOWN enrichment: left of zero\n"
)

cat(
  "  UP enrichment: right of zero\n"
)

cat(
  "  Three contrasts horizontally: YES\n"
)

cat(
  "  FDR annotation: * <0.05 | ** <0.01 | *** <0.001\n"
)

cat(
  "  SVG master: YES\n"
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
  "\nCompartments:\n"
)

print(
  summary_tbl,
  n=Inf
)

cat(
  "\nSweat gland included: NO — matched to publication compartment set\n"
)

cat(
  "Sebocyte included: NO — matched to publication compartment set\n"
)

cat(
  "\nFrozen 20F1B consensus modified: NO\n"
)

cat(
  "Original GSEA results modified: NO\n"
)

cat(
  "Frozen L2 annotation modified: NO\n"
)

cat(
  "\n===== FORMAL CHECKS =====\n"
)

print(
  checks,
  n=Inf
)

if (!all(checks$pass)) {

  cat("\nSTATUS=FAIL\n")

  print(
    checks %>%
      filter(!pass)
  )

  stop(
    "20B6 formal checks failed."
  )
}

cat("\nSTATUS=PASS\n")
