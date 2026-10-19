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
cat("20B7 — MULTICOMPARTMENT CONSENSUS NES STAR PLOTS\n")
cat("============================================================\n")

# ============================================================
# PATHS
# ============================================================

INPUT <- paste0(
  "PIPELINE_V2/results/publication_panels/",
  "20B6_cross_database_consensus_mirror_plots_v2/",
  "tables/20B6v2_all_compartment_consensus_axes.csv"
)

OUT <- paste0(
  "PIPELINE_V2/results/publication_panels/",
  "20B7_multicompartment_consensus_starplots"
)

SVG_DIR <- file.path(OUT, "figures", "SVG")
PDF_DIR <- file.path(OUT, "figures", "PDF")
PNG_DIR <- file.path(OUT, "figures", "PNG")
TAB_DIR <- file.path(OUT, "tables")

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

stopifnot(
  file.exists(INPUT)
)

x <- read_csv(
  INPUT,
  show_col_types=FALSE
)

cat(
  "\nInput rows:",
  nrow(x),
  "\n"
)

cat(
  "Compartments:",
  paste(
    sort(unique(x$compartment)),
    collapse=" | "
  ),
  "\n"
)

cat(
  "Axes:",
  paste(
    sort(unique(x$axis)),
    collapse=" | "
  ),
  "\n"
)

# ============================================================
# CONTRASTS
# ============================================================

contrast_map <- c(
  "AA_vs_HC"  = "AA vs HC",
  "FFA_vs_HC" = "FFA vs HC"
)

x <- x %>%
  filter(
    contrast %in% names(contrast_map)
  )

# ============================================================
# CANONICAL AXIS ORDER
# ============================================================

axis_order <- c(
  "IFN_JAK_STAT_ANTIGEN",
  "T_CELL_CYTOTOXIC",
  "CELL_CYCLE_DNA_REPAIR",
  "METABOLISM_OXPHOS_ROS",
  "ECM_COLLAGEN_FIBROSIS",
  "VASCULAR_ANGIOGENESIS",
  "WNT_HEDGEHOG_FOLLICLE",
  "EPITHELIAL_DIFFERENTIATION"
)

axis_labels <- c(
  "IFN_JAK_STAT_ANTIGEN" =
    "IFN / JAK-STAT /\nantigen",

  "T_CELL_CYTOTOXIC" =
    "T-cell /\ncytotoxic",

  "CELL_CYCLE_DNA_REPAIR" =
    "Cell cycle /\nDNA repair",

  "METABOLISM_OXPHOS_ROS" =
    "Metabolism /\nOXPHOS / ROS",

  "ECM_COLLAGEN_FIBROSIS" =
    "ECM / collagen /\nfibrosis",

  "VASCULAR_ANGIOGENESIS" =
    "Vascular /\nangiogenesis",

  "WNT_HEDGEHOG_FOLLICLE" =
    "WNT / Hedgehog /\nfollicle",

  "EPITHELIAL_DIFFERENTIATION" =
    "Epithelial\ndifferentiation"
)

# ============================================================
# COMPARTMENT ORDER
# ============================================================

compartment_order <- c(
  "Epithelial",
  "Immune",
  "Fibroblast",
  "Endothelial",
  "Mural",
  "Melanocyte"
)

# ============================================================
# COMPARTMENT COLORS
#
# Epithelial retains canonical project L1 color.
# Remaining colors are compartment-identity colors,
# deliberately distinct from disease colors.
# ============================================================

compartment_colors <- c(
  "Epithelial"  = "#8A034D",
  "Immune"      = "#3C5488",
  "Fibroblast"  = "#D07C2F",
  "Endothelial" = "#2C8C89",
  "Mural"       = "#7A5195",
  "Melanocyte"  = "#6B6B6B"
)

# ============================================================
# COMPLETE GRID
#
# Missing compartment × axis combinations are NOT interpreted
# as NES=0. They remain NA.
# ============================================================

plot_data <- x %>%
  select(
    compartment,
    contrast,
    axis,
    consensus_NES,
    n_L2,
    min_FDR
  ) %>%
  complete(
    compartment=compartment_order,
    contrast=names(contrast_map),
    axis=axis_order
  ) %>%
  mutate(
    compartment=factor(
      compartment,
      levels=compartment_order
    ),

    axis=factor(
      axis,
      levels=axis_order
    ),

    axis_label=unname(
      axis_labels[
        as.character(axis)
      ]
    )
  )

write_csv(
  plot_data,
  file.path(
    TAB_DIR,
    "20B7_starplot_input_complete.csv"
  )
)

# ============================================================
# GLOBAL RADIAL SCALE
#
# Same NES scale for AA and FFA so geometries are comparable.
# We transform signed NES to positive radial coordinates:
#
# radial_value = NES + offset
#
# The circle at radial_value = offset corresponds to NES=0.
# Therefore:
#   outside zero circle = positive enrichment
#   inside zero circle  = negative enrichment
# ============================================================

max_abs_NES <- max(
  abs(plot_data$consensus_NES),
  na.rm=TRUE
)

NES_LIMIT <- ceiling(
  max_abs_NES * 2
) / 2

NES_LIMIT <- max(
  NES_LIMIT,
  2
)

RADIAL_OFFSET <- NES_LIMIT

plot_data <- plot_data %>%
  mutate(
    radial_value=consensus_NES + RADIAL_OFFSET
  )

cat(
  "\nGlobal |NES| limit:",
  NES_LIMIT,
  "\n"
)

cat(
  "Zero-NES radial circle:",
  RADIAL_OFFSET,
  "\n"
)

# ============================================================
# RADIAL GRID
# ============================================================

nes_breaks <- seq(
  -NES_LIMIT,
  NES_LIMIT,
  by=0.5
)

nes_breaks <- nes_breaks[
  nes_breaks >= -NES_LIMIT &
  nes_breaks <= NES_LIMIT
]

radial_breaks <- nes_breaks + RADIAL_OFFSET

# ============================================================
# STAR PLOT FUNCTION
# ============================================================

make_star <- function(
  dat,
  contrast_id
) {

  z <- dat %>%
    filter(
      contrast == contrast_id
    )

  # ----------------------------------------------------------
  # close each polygon/line by repeating first axis
  # ----------------------------------------------------------

  z_line <- z %>%
    filter(
      !is.na(consensus_NES)
    ) %>%
    arrange(
      compartment,
      axis
    ) %>%
    group_by(
      compartment
    ) %>%
    group_modify(
      ~{
        if (nrow(.x) < 2) {
          return(.x)
        }

        bind_rows(
          .x,
          .x[1,,drop=FALSE]
        )
      }
    ) %>%
    ungroup()

  # ----------------------------------------------------------
  # Background axis spokes
  # ----------------------------------------------------------

  spokes <- tibble(
    axis=factor(
      axis_order,
      levels=axis_order
    ),
    r0=0,
    r1=2 * RADIAL_OFFSET
  )

  # ----------------------------------------------------------
  # Main plot
  # ----------------------------------------------------------

  p <- ggplot() +

    # radial spokes
    geom_segment(
      data=spokes,
      aes(
        x=axis,
        xend=axis,
        y=r0,
        yend=r1
      ),
      linewidth=0.30,
      colour="#D7D7D7"
    ) +

    # compartment lines
    geom_path(
      data=z_line,
      aes(
        x=axis,
        y=radial_value,
        group=compartment,
        colour=compartment
      ),
      linewidth=0.95,
      alpha=0.82,
      na.rm=TRUE
    ) +

    # patient-compartment consensus points
    geom_point(
      data=z,
      aes(
        x=axis,
        y=radial_value,
        colour=compartment
      ),
      size=2.6,
      alpha=0.95,
      na.rm=TRUE
    ) +

    # emphasize zero-NES circle
    geom_hline(
      yintercept=RADIAL_OFFSET,
      linewidth=0.80,
      colour="#333333"
    ) +

    coord_polar(
      start=0
    ) +

    scale_colour_manual(
      values=compartment_colors,
      drop=FALSE
    ) +

    scale_x_discrete(
      labels=function(v) {
        unname(
          axis_labels[v]
        )
      },
      drop=FALSE
    ) +

    scale_y_continuous(
      limits=c(
        0,
        2 * RADIAL_OFFSET
      ),

      breaks=radial_breaks,

      labels=function(r) {

        nes <- r - RADIAL_OFFSET

        ifelse(
          abs(nes) < 1e-8,
          "0",
          sprintf(
            "%+.1f",
            nes
          )
        )
      }
    ) +

    labs(
      title=unname(
        contrast_map[
          contrast_id
        ]
      ),

      subtitle=paste0(
        "Signed consensus NES — ",
        "inside zero ring = down; outside = up"
      ),

      x=NULL,
      y=NULL,
      colour="Compartment"
    ) +

    theme_minimal(
      base_size=11
    ) +

    theme(
      panel.grid.major.x=element_blank(),

      panel.grid.minor=element_blank(),

      panel.grid.major.y=element_line(
        colour="#E5E5E5",
        linewidth=0.35
      ),

      axis.text.x=element_text(
        size=8.6,
        face="bold",
        colour="#222222"
      ),

      axis.text.y=element_text(
        size=7.5,
        colour="#666666"
      ),

      plot.title=element_text(
        face="bold",
        size=14,
        hjust=0.5
      ),

      plot.subtitle=element_text(
        size=9,
        hjust=0.5
      ),

      legend.position="bottom",

      legend.title=element_text(
        face="bold"
      ),

      legend.text=element_text(
        size=9
      ),

      plot.margin=margin(
        15,
        25,
        15,
        25
      )
    )

  p
}

# ============================================================
# GENERATE
# ============================================================

p_AA <- make_star(
  plot_data,
  "AA_vs_HC"
)

p_FFA <- make_star(
  plot_data,
  "FFA_vs_HC"
)

# ============================================================
# INDIVIDUAL FIGURES
# ============================================================

ggsave(
  file.path(
    SVG_DIR,
    "20B7_AA_vs_HC_consensus_star.svg"
  ),
  p_AA,
  width=9,
  height=8,
  units="in",
  device=svglite
)

ggsave(
  file.path(
    PDF_DIR,
    "20B7_AA_vs_HC_consensus_star.pdf"
  ),
  p_AA,
  width=9,
  height=8,
  units="in"
)

ggsave(
  file.path(
    PNG_DIR,
    "20B7_AA_vs_HC_consensus_star.png"
  ),
  p_AA,
  width=9,
  height=8,
  units="in",
  dpi=400
)

ggsave(
  file.path(
    SVG_DIR,
    "20B7_FFA_vs_HC_consensus_star.svg"
  ),
  p_FFA,
  width=9,
  height=8,
  units="in",
  device=svglite
)

ggsave(
  file.path(
    PDF_DIR,
    "20B7_FFA_vs_HC_consensus_star.pdf"
  ),
  p_FFA,
  width=9,
  height=8,
  units="in"
)

ggsave(
  file.path(
    PNG_DIR,
    "20B7_FFA_vs_HC_consensus_star.png"
  ),
  p_FFA,
  width=9,
  height=8,
  units="in",
  dpi=400
)

# ============================================================
# COMBINED AA + FFA
# ============================================================

combined <- (
  p_AA |
  p_FFA
) +
  plot_layout(
    guides="collect"
  ) +
  plot_annotation(
    title="Multicompartment consensus pathway architecture",

    subtitle=paste0(
      "Cross-database CORE/TRIPLE consensus across ",
      "BioPlanet, GO-BP and Reactome"
    ),

    caption=paste0(
      "Signed NES summarized across supported L2 states. ",
      "The bold circular reference denotes NES=0; ",
      "outward displacement indicates positive enrichment and ",
      "inward displacement negative enrichment."
    ),

    theme=theme(
      plot.title=element_text(
        face="bold",
        size=16,
        hjust=0.5
      ),

      plot.subtitle=element_text(
        size=10,
        hjust=0.5
      ),

      plot.caption=element_text(
        size=8.5,
        hjust=1
      ),

      legend.position="bottom"
    )
  )

ggsave(
  file.path(
    SVG_DIR,
    "20B7_AA_FFA_vs_HC_multicompartment_consensus_star.svg"
  ),
  combined,
  width=16.5,
  height=8.5,
  units="in",
  device=svglite
)

ggsave(
  file.path(
    PDF_DIR,
    "20B7_AA_FFA_vs_HC_multicompartment_consensus_star.pdf"
  ),
  combined,
  width=16.5,
  height=8.5,
  units="in"
)

ggsave(
  file.path(
    PNG_DIR,
    "20B7_AA_FFA_vs_HC_multicompartment_consensus_star.png"
  ),
  combined,
  width=16.5,
  height=8.5,
  units="in",
  dpi=400
)

# ============================================================
# FORMAL CHECKS
# ============================================================

svg_files <- list.files(
  SVG_DIR,
  pattern="\\.svg$",
  full.names=TRUE
)

pdf_files <- list.files(
  PDF_DIR,
  pattern="\\.pdf$",
  full.names=TRUE
)

png_files <- list.files(
  PNG_DIR,
  pattern="\\.png$",
  full.names=TRUE
)

checks <- tibble(
  check=c(
    "two_disease_vs_HC_contrasts",
    "eight_frozen_axes",
    "six_compartments_defined",
    "signed_NES_retained",
    "missing_values_not_imputed_zero",
    "three_SVG",
    "three_PDF",
    "three_PNG",
    "all_outputs_nonzero"
  ),

  pass=c(
    all(
      c(
        "AA_vs_HC",
        "FFA_vs_HC"
      ) %in%
        unique(
          plot_data$contrast
        )
    ),

    length(
      axis_order
    ) == 8,

    length(
      compartment_order
    ) == 6,

    any(
      plot_data$consensus_NES < 0,
      na.rm=TRUE
    ) &&
      any(
        plot_data$consensus_NES > 0,
        na.rm=TRUE
      ),

    any(
      is.na(
        plot_data$consensus_NES
      )
    ),

    length(svg_files) == 3,
    length(pdf_files) == 3,
    length(png_files) == 3,

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
    TAB_DIR,
    "20B7_formal_checks.csv"
  )
)

# ============================================================
# SUMMARY
# ============================================================

summary_tbl <- plot_data %>%
  group_by(
    contrast,
    compartment
  ) %>%
  summarise(
    n_axes_available=sum(
      !is.na(
        consensus_NES
      )
    ),

    n_positive=sum(
      consensus_NES > 0,
      na.rm=TRUE
    ),

    n_negative=sum(
      consensus_NES < 0,
      na.rm=TRUE
    ),

    median_NES=median(
      consensus_NES,
      na.rm=TRUE
    ),

    .groups="drop"
  )

write_csv(
  summary_tbl,
  file.path(
    TAB_DIR,
    "20B7_compartment_summary.csv"
  )
)

# ============================================================
# FINAL REPORT
# ============================================================
# ============================================================
# FINAL REPORT
# ============================================================

cat("\n============================================================\n")
cat("20B7 - FINAL REPORT\n")
cat("============================================================\n")

cat(
  "\nInput:",
  INPUT,
  "\n"
)

cat(
  "Contrasts: AA vs HC | FFA vs HC\n"
)

cat(
  "Frozen biological axes:",
  length(axis_order),
  "\n"
)

cat(
  "Compartments defined:",
  length(compartment_order),
  "\n"
)

cat(
  "Compartments with consensus data:",
  n_distinct(
    x$compartment
  ),
  "\n"
)

cat(
  "\nSigned NES retained: YES\n"
)

cat(
  "NES=0 represented by reference ring: YES\n"
)

cat(
  "Positive NES outside reference ring: YES\n"
)

cat(
  "Negative NES inside reference ring: YES\n"
)

cat(
  "Missing compartment-axis combinations imputed as zero: NO\n"
)

cat(
  "\nRepresentation:\n"
)

cat(
  "  Each radial axis = one frozen biological consensus axis\n"
)

cat(
  "  Each line = one cellular compartment\n"
)

cat(
  "  Color = cellular compartment\n"
)

cat(
  "  Common radial NES scale across AA and FFA: YES\n"
)

cat(
  "  Direct AA-vs-FFA comparison displayed: NO\n"
)

cat(
  "\nFigures:\n"
)

cat(
  "  AA vs HC individual star plot: YES\n"
)

cat(
  "  FFA vs HC individual star plot: YES\n"
)

cat(
  "  Combined AA + FFA figure: YES\n"
)

cat(
  "  SVG:",
  length(svg_files),
  "\n"
)

cat(
  "  PDF:",
  length(pdf_files),
  "\n"
)

cat(
  "  PNG:",
  length(png_files),
  "\n"
)

cat(
  "\nSummary:\n"
)

print(
  summary_tbl,
  n=Inf
)

cat(
  "\n===== FORMAL CHECKS =====\n"
)

print(
  checks,
  n=Inf
)

cat(
  "\nFrozen 20F1B consensus modified: NO\n"
)

cat(
  "20B6 v2 modified: NO\n"
)

cat(
  "Original GSEA modified: NO\n"
)

cat(
  "Frozen L2 annotation modified: NO\n"
)

if (!all(checks$pass)) {

  cat("\nSTATUS=FAIL\n")

  print(
    checks %>%
      filter(!pass)
  )

  stop(
    "20B7 formal checks failed."
  )
}

cat("\nSTATUS=PASS\n")
