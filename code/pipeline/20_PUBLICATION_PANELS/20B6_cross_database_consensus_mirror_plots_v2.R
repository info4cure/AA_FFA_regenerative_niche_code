suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(ggplot2)
  library(patchwork)
  library(svglite)
})

cat("============================================================\n")
cat("20B6 v2 — DISEASE-vs-HC CONSENSUS ENRICHMENT PYRAMIDS\n")
cat("============================================================\n")

INPUT <- paste0(
  "PIPELINE_V2/results/disease_biology/",
  "20_PATHWAY_ANALYSIS/20F_CROSS_DATABASE_CONSENSUS/",
  "20F1B_consensus/tables/",
  "20F1B_core_consensus_map.csv"
)

OUT <- paste0(
  "PIPELINE_V2/results/publication_panels/",
  "20B6_cross_database_consensus_mirror_plots_v2"
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
# PUBLICATION COMPARTMENTS
# ============================================================

module_map <- c(
  "07_EPITHELIAL_L2" = "Epithelial",
  "08_IMMUNE_L2"     = "Immune",
  "09_FIBROBLAST_L2" = "Fibroblast",
  "10_VASCULAR_L2"   = "Endothelial",
  "11_MURAL_L2"      = "Mural",
  "13_MELANOCYTE_L2" = "Melanocyte"
)

# Disease-vs-control plots requested
contrast_order <- c(
  "FFA_vs_HC",
  "AA_vs_HC"
)

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
# INPUT FILTER
# CORE/TRIPLE only — already guaranteed by frozen source
# ============================================================

d <- x %>%
  filter(
    module %in% names(module_map),
    contrast %in% contrast_order
  ) %>%
  mutate(
    compartment=unname(
      module_map[module]
    ),
    signed_NES=median_database_NES
  )

# ============================================================
# DIRECTION SANITY CHECK
# ============================================================

direction_check <- d %>%
  mutate(
    direction_ok=
      (direction == "UP"   & signed_NES > 0) |
      (direction == "DOWN" & signed_NES < 0)
  )

if (!all(direction_check$direction_ok)) {

  bad <- direction_check %>%
    filter(!direction_ok)

  print(bad)

  stop(
    "Direction/NES inconsistency detected."
  )
}

# ============================================================
# ONE CONSENSUS VALUE PER:
# compartment x contrast x biological axis
#
# IMPORTANT:
# each frozen L2 CORE/TRIPLE signal contributes its signed NES.
# Median across L2 states gives one compartment-level direction.
# ============================================================

consensus <- d %>%
  group_by(
    compartment,
    contrast,
    axis
  ) %>%
  summarise(
    consensus_NES=median(
      signed_NES,
      na.rm=TRUE
    ),

    n_L2=n_distinct(
      final_L2_label
    ),

    n_signals=n(),

    min_FDR=min(
      min_FDR,
      na.rm=TRUE
    ),

    n_CORE=sum(
      consensus_class == "CORE",
      na.rm=TRUE
    ),

    n_TRIPLE=sum(
      consensus_class == "TRIPLE",
      na.rm=TRUE
    ),

    n_Hallmark_supported=sum(
      Hallmark_support,
      na.rm=TRUE
    ),

    .groups="drop"
  ) %>%
  mutate(
    abs_NES=abs(
      consensus_NES
    ),

    direction=case_when(
      consensus_NES > 0 ~ "UP",
      consensus_NES < 0 ~ "DOWN",
      TRUE              ~ "ZERO"
    ),

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
    )
  )

write_csv(
  consensus,
  file.path(
    TAB,
    "20B6v2_all_compartment_consensus_axes.csv"
  )
)

# ============================================================
# PUBLICATION FILTER
#
# CORE/TRIPLE source already implies significant component
# evidence.
#
# Require a meaningful compartment-level directional effect:
# |median consensus NES| >= 1.
#
# This removes axes whose opposing L2 signals cancel near zero.
# ============================================================

eligible <- consensus %>%
  filter(
    is.finite(consensus_NES),
    abs_NES >= 1
  )

write_csv(
  eligible,
  file.path(
    TAB,
    "20B6v2_eligible_consensus_axes_absNES_ge1.csv"
  )
)

# ============================================================
# TOP 10 PER DIRECTION
# no artificial filling
# ============================================================

selected <- eligible %>%
  group_by(
    compartment,
    contrast,
    direction
  ) %>%
  arrange(
    desc(abs_NES),
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
  ungroup() %>%
  mutate(
    fdr_symbol=case_when(
      min_FDR < 0.001 ~ "***",
      min_FDR < 0.01  ~ "**",
      min_FDR < 0.05  ~ "*",
      TRUE            ~ ""
    ),

    value_label=paste0(
      sprintf(
        "%.2f",
        abs_NES
      ),
      " ",
      fdr_symbol
    )
  )

write_csv(
  selected,
  file.path(
    TAB,
    "20B6v2_selected_axes.csv"
  )
)

# ============================================================
# COLORS REPRESENT DIRECTION, NOT DISEASE
# ============================================================

direction_colors <- c(
  DOWN="#5F86BD",
  UP="#E96B72"
)

# ============================================================
# PLOT FUNCTION
# ============================================================

make_pyramid <- function(
  df,
  contrast,
  disease
) {

  z <- df %>%
    filter(
      .data$contrast == .env$contrast
    )

  if (nrow(z) == 0) {

    return(
      ggplot() +
        annotate(
          "text",
          x=0,
          y=0,
          label=paste0(
            "No robust CORE/TRIPLE consensus\n",
            "(|consensus NES| ≥ 1)"
          ),
          size=4
        ) +
        xlim(-1,1) +
        ylim(-1,1) +
        theme_void() +
        labs(
          title=paste0(
            disease,
            " vs HC"
          )
        )
    )
  }

  # ----------------------------------------------------------
  # Separate labels on the left and right, as in a true
  # enrichment pyramid.
  #
  # Each biological axis occurs only once because consensus was
  # calculated before direction assignment.
  # ----------------------------------------------------------

  z <- z %>%
    arrange(
      consensus_NES
    ) %>%
    mutate(
      row_id=row_number()
    )

  lim <- max(
    abs(z$consensus_NES),
    na.rm=TRUE
  )

  lim <- max(
    1,
    ceiling(
      lim * 10
    ) / 10
  )

  label_pad <- lim * 0.08

  p <- ggplot(
    z,
    aes(
      x=consensus_NES,
      y=reorder(
        axis_label,
        row_id
      ),
      fill=direction
    )
  ) +

    geom_col(
      width=0.70,
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
          direction == "UP"
        ),
      aes(
        label=value_label
      ),
      hjust=-0.12,
      size=3.2,
      fontface="bold"
    ) +

    geom_text(
      data=z %>%
        filter(
          direction == "DOWN"
        ),
      aes(
        label=value_label
      ),
      hjust=1.12,
      size=3.2,
      fontface="bold"
    ) +

    scale_fill_manual(
      values=direction_colors,
      breaks=c(
        "DOWN",
        "UP"
      ),
      labels=c(
        paste0(
          "Down in ",
          disease
        ),
        paste0(
          "Up in ",
          disease
        )
      )
    ) +

    scale_x_continuous(
      limits=c(
        -lim * 1.32,
        lim * 1.32
      ),
      breaks=scales::pretty_breaks(
        n=5
      ),
      labels=function(v) {
        abs(v)
      }
    ) +

    labs(
      title=paste0(
        disease,
        " vs HC"
      ),
      subtitle=paste0(
        "Down in ",
        disease,
        "                       Up in ",
        disease
      ),
      x="Consensus NES",
      y=NULL,
      fill=NULL
    ) +

    theme_classic(
      base_size=11
    ) +

    theme(
      legend.position="top",

      legend.justification="center",

      legend.text=element_text(
        size=10
      ),

      plot.title=element_text(
        face="bold",
        size=13,
        hjust=0.5
      ),

      plot.subtitle=element_text(
        face="bold",
        size=9.5,
        hjust=0.5
      ),

      axis.text.y=element_text(
        size=9.5,
        colour="#222222"
      ),

      axis.text.x=element_text(
        size=9
      ),

      axis.title.x=element_text(
        size=10
      ),

      axis.line.y=element_blank(),

      axis.ticks.y=element_blank(),

      plot.margin=margin(
        8,
        18,
        8,
        18
      )
    )

  p
}

# ============================================================
# GENERATE ONE FIGURE PER COMPARTMENT
# ============================================================

summary_list <- list()

for (
  comp in unname(
    module_map
  )
) {

  cat("\n============================================================\n")
  cat(comp, "\n")
  cat("============================================================\n")

  dc <- selected %>%
    filter(
      .data$compartment == .env$comp
    )

  SVG_DIR <- file.path(
    OUT,
    comp,
    "figures",
    "SVG"
  )

  PDF_DIR <- file.path(
    OUT,
    comp,
    "figures",
    "PDF"
  )

  PNG_DIR <- file.path(
    OUT,
    comp,
    "figures",
    "PNG"
  )

  TAB_DIR <- file.path(
    OUT,
    comp,
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
        "20B6v2_",
        comp,
        "_selected_axes.csv"
      )
    )
  )

  p_ffa <- make_pyramid(
    dc,
    "FFA_vs_HC",
    "FFA"
  )

  p_aa <- make_pyramid(
    dc,
    "AA_vs_HC",
    "AA"
  )

  final <- (
    p_ffa |
    p_aa
  ) +

    plot_annotation(
      title=paste0(
        comp,
        " — consensus pathway enrichment"
      ),

      subtitle=paste0(
        "Cross-database CORE/TRIPLE consensus: ",
        "BioPlanet, GO-BP and Reactome; ",
        "Hallmark support retained where available"
      ),

      caption=paste0(
        "Compartment summary of frozen L2 consensus. ",
        "Bars show median signed NES across supported L2 states. ",
        "Only |consensus NES| ≥ 1 shown. ",
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
    "20B6v2_",
    comp,
    "_disease_vs_HC_consensus_pyramid"
  )

  ggsave(
    file.path(
      SVG_DIR,
      paste0(
        stem,
        ".svg"
      )
    ),
    final,
    width=13.5,
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
    width=13.5,
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
    width=13.5,
    height=7.5,
    units="in",
    dpi=400
  )

  cat("  -> SVG/PDF/PNG generated\n")

  summary_list[[comp]] <- tibble(
    compartment=comp,

    FFA_DOWN=sum(
      dc$contrast == "FFA_vs_HC" &
      dc$direction == "DOWN"
    ),

    FFA_UP=sum(
      dc$contrast == "FFA_vs_HC" &
      dc$direction == "UP"
    ),

    AA_DOWN=sum(
      dc$contrast == "AA_vs_HC" &
      dc$direction == "DOWN"
    ),

    AA_UP=sum(
      dc$contrast == "AA_vs_HC" &
      dc$direction == "UP"
    ),

    n_total=nrow(dc)
  )
}

# ============================================================
# GLOBAL SUMMARY
# ============================================================

summary_tbl <- bind_rows(
  summary_list
)

write_csv(
  summary_tbl,
  file.path(
    TAB,
    "20B6v2_compartment_summary.csv"
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
    "six_compartments_generated",
    "six_SVG",
    "six_PDF",
    "six_PNG",
    "no_duplicate_axis_per_compartment_contrast",
    "all_displayed_absNES_ge1",
    "only_disease_vs_HC",
    "all_outputs_nonzero"
  ),

  pass=c(
    nrow(summary_tbl) == 6,
    length(svg_files) == 6,
    length(pdf_files) == 6,
    length(png_files) == 6,

    selected %>%
      count(
        compartment,
        contrast,
        axis
      ) %>%
      summarise(
        ok=all(n == 1)
      ) %>%
      pull(ok),

    all(
      abs(
        selected$consensus_NES
      ) >= 1
    ),

    all(
      selected
      $contrast %in% c(
        "FFA_vs_HC",
        "AA_vs_HC"
      )
    ),

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
    "20B6v2_formal_checks.csv"
  )
)

# ============================================================
# FINAL REPORT
# ============================================================

cat("\n============================================================\n")
cat("20B6 v2 — FINAL REPORT\n")
cat("============================================================\n")

cat(
  "\nSource:",
  INPUT,
  "\n"
)

cat(
  "\nCompartments generated:",
  nrow(summary_tbl),
  "\n"
)

cat(
  "Reference group: HC\n"
)

cat(
  "Contrasts displayed: FFA vs HC | AA vs HC\n"
)

cat(
  "\nDirection:\n"
)

cat(
  "  NES < 0: DOWN in disease vs HC\n"
)

cat(
  "  NES > 0: UP in disease vs HC\n"
)

cat(
  "\nConsensus definition:\n"
)

cat(
  "  Source: frozen 20F1B CORE/TRIPLE cross-database consensus\n"
)

cat(
  "  Databases: BioPlanet | GO-BP | Reactome\n"
)

cat(
  "  Hallmark support retained where available: YES\n"
)

cat(
  "  Single-database signals included: NO\n"
)

cat(
  "  One biological axis per compartment x contrast: YES\n"
)

cat(
  "  L2 aggregation: median signed NES\n"
)

cat(
  "  Opposing L2 signals allowed to cancel: YES\n"
)

cat(
  "  Minimum displayed |consensus NES|: 1\n"
)

cat(
  "\nSelection:\n"
)

cat(
  "  Maximum UP axes per contrast: 10\n"
)

cat(
  "  Maximum DOWN axes per contrast: 10\n"
)

cat(
  "  Artificial filling: NO\n"
)

cat(
  "\nDisplay:\n"
)

cat(
  "  Left: DOWN in disease\n"
)

cat(
  "  Right: UP in disease\n"
)

cat(
  "  Blue: DOWN\n"
)

cat(
  "  Red: UP\n"
)

cat(
  "  Disease-vs-HC panels: 2\n"
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
  "\nSweat gland included: NO\n"
)

cat(
  "Sebocyte included: NO\n"
)

cat(
  "\nFrozen 20F1B consensus modified: NO\n"
)

cat(
  "Original GSEA modified: NO\n"
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
      filter(
        !pass
      )
  )

  stop(
    "20B6 v2 formal checks failed."
  )
}

cat("\nSTATUS=PASS\n")
