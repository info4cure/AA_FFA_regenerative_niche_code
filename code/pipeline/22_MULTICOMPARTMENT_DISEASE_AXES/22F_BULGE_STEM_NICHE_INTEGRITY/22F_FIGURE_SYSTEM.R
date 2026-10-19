
# ============================================================
# 22F FIGURE SYSTEM
# Canonical visual system for Bulge/Isthmus publication figures
# ============================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tibble)
})

# ============================================================
# 1. CANONICAL PALETTES
# ============================================================

K3_COL <- c(
  T3="#2F7189",
  T2="#4B4C7A",
  T1="#B36B00"
)

DISEASE_COL <- c(
  HC="#8F8D8E",
  FFA="#4D9E4F",
  AA="#A8142A"
)

MODULE_COL <- c(
  M1="#4E79A7",
  M2="#59A14F",
  M3="#E15759",
  M4="#F28E2B",
  M5="#B07AA1"
)

DIRECTION_COL <- c(
  "Late-T1 higher"="#2166AC",
  "Y37 higher"="#B2182B"
)

# ============================================================
# 2. CANONICAL TRUNK STATES
# ============================================================

BULGE_STATE_REFERENCE <- tibble(
  state=factor(
    c("T3","T2","T1"),
    levels=c("T3","T2","T1")
  ),
  xmin=c(
    0,
    0.438,
    0.844
  ),
  xmax=c(
    0.438,
    0.844,
    1
  )
)

BULGE_T2_BOUNDARY <- 0.438
BULGE_T1_BOUNDARY <- 0.844

# ============================================================
# 3. CANONICAL BIOLOGICAL LABELS
# ============================================================

BULGE_LABELS <- c(
  T3="T3",
  T2="T2",
  T1="T1",
  late_T1="Proliferative T1",
  Y12="Y12 branch point",
  Y8="Minor post-Y12 arm (Y8)",
  Y37="Post-proliferative T1 (Y37)"
)

# Important:
# T1 remains the frozen K=3 macroterritory.
# "Proliferative T1" and "Post-proliferative T1 (Y37)"
# are contextual labels for topological/molecular substates.
# They DO NOT redefine the frozen K=3 partition.

# ============================================================
# 4. TYPOGRAPHY
# ============================================================

FS <- list(
  figure_title=15,
  figure_subtitle=9.5,

  panel_title=10.5,
  panel_subtitle=8.2,

  axis_title=8.5,
  axis_text=7.5,

  legend_title=8,
  legend_text=7,

  annotation=7.5,
  small=6.8,
  caption=7
)

# ============================================================
# 5. LINE / POINT GEOMETRY
# ============================================================

GEOM <- list(
  umap_point_size=0.48,
  umap_point_alpha=0.95,
  umap_point_stroke=0.16,

  graph_linewidth=0.75,

  boundary_linewidth=0.45,

  heatmap_separator=0.28,

  patient_point_size=2.0
)

# ============================================================
# 6. BASE PUBLICATION THEME
# ============================================================

theme_bulge <- function(
  base_size=8.5
) {

  theme_classic(
    base_size=base_size
  ) +

    theme(
      plot.title=
        element_text(
          face="bold",
          size=FS$panel_title,
          hjust=0,
          margin=margin(
            b=2
          )
        ),

      plot.subtitle=
        element_text(
          size=FS$panel_subtitle,
          hjust=0,
          margin=margin(
            b=5
          )
        ),

      axis.title=
        element_text(
          size=FS$axis_title
        ),

      axis.text=
        element_text(
          size=FS$axis_text,
          colour="#252525"
        ),

      legend.title=
        element_text(
          size=FS$legend_title,
          face="bold"
        ),

      legend.text=
        element_text(
          size=FS$legend_text
        ),

      strip.text=
        element_text(
          size=FS$axis_text,
          face="bold"
        ),

      plot.margin=
        margin(
          5,
          6,
          5,
          5
        )
    )
}

# ============================================================
# 7. FIGURE-LEVEL ANNOTATION THEME
# ============================================================

theme_bulge_annotation <- function() {

  theme(
    plot.title=
      element_text(
        face="bold",
        size=FS$figure_title,
        hjust=0,
        margin=margin(
          b=3
        )
      ),

    plot.subtitle=
      element_text(
        size=FS$figure_subtitle,
        hjust=0,
        margin=margin(
          b=8
        )
      ),

    plot.caption=
      element_text(
        size=FS$caption,
        colour="#555555",
        hjust=0,
        margin=margin(
          t=7
        )
      )
  )
}

# ============================================================
# 8. STANDARD T3/T2/T1 BAND
# ============================================================

bulge_state_band <- function(
  show_labels=TRUE,
  show_axis=TRUE,
  title=NULL,
  subtitle=NULL
) {

  p <- ggplot(
    BULGE_STATE_REFERENCE
  ) +

    geom_rect(
      aes(
        xmin=xmin,
        xmax=xmax,
        ymin=0,
        ymax=1,
        fill=state
      ),
      colour="white",
      linewidth=0.5
    ) +

    scale_fill_manual(
      values=K3_COL,
      drop=FALSE
    ) +

    scale_x_continuous(
      limits=c(0,1),
      breaks=seq(
        0,
        1,
        0.2
      ),
      expand=c(0,0)
    ) +

    coord_cartesian(
      ylim=c(0,1),
      clip="off"
    ) +

    labs(
      title=title,
      subtitle=subtitle,
      x=if(show_axis)
        "Normalized trunk position"
      else
        NULL,
      y=NULL
    ) +

    theme_bulge() +

    theme(
      axis.text.y=
        element_blank(),

      axis.ticks.y=
        element_blank(),

      axis.line.y=
        element_blank(),

      legend.position="none"
    )

  if(show_labels) {

    p <- p +

      geom_text(
        aes(
          x=(xmin+xmax)/2,
          y=0.5,
          label=state
        ),
        colour="white",
        fontface="bold",
        size=4
      )
  }

  p
}

# ============================================================
# 9. ADD STATE BOUNDARIES TO ANY 0-1 TRAJECTORY PLOT
# ============================================================

add_bulge_state_boundaries <- function(
  p,
  include_T3_T2=TRUE,
  include_T2_T1=TRUE,
  colour="#353535",
  alpha=0.75
) {

  if(include_T3_T2) {

    p <- p +

      geom_vline(
        xintercept=BULGE_T2_BOUNDARY,
        linetype=2,
        linewidth=GEOM$boundary_linewidth,
        colour=colour,
        alpha=alpha
      )
  }

  if(include_T2_T1) {

    p <- p +

      geom_vline(
        xintercept=BULGE_T1_BOUNDARY,
        linetype=2,
        linewidth=GEOM$boundary_linewidth,
        colour=colour,
        alpha=alpha
      )
  }

  p
}

# ============================================================
# 10. ADD T3/T2/T1 LABELS ABOVE A TRAJECTORY
# ============================================================

add_bulge_state_labels <- function(
  p,
  y,
  colour="#252525",
  size=3
) {

  labs_df <- BULGE_STATE_REFERENCE |>
    mutate(
      x=(xmin+xmax)/2
    )

  p +

    geom_text(
      data=labs_df,
      aes(
        x=x,
        y=y,
        label=state
      ),
      inherit.aes=FALSE,
      colour=colour,
      fontface="bold",
      size=size
    )
}

# ============================================================
# 11. DISEASE SCALE
# ============================================================

scale_fill_disease <- function(...) {

  scale_fill_manual(
    values=DISEASE_COL,
    drop=FALSE,
    ...
  )
}

scale_colour_disease <- function(...) {

  scale_colour_manual(
    values=DISEASE_COL,
    drop=FALSE,
    ...
  )
}

# ============================================================
# 12. K3 SCALE
# ============================================================

scale_fill_K3 <- function(...) {

  scale_fill_manual(
    values=K3_COL,
    drop=FALSE,
    ...
  )
}

scale_colour_K3 <- function(...) {

  scale_colour_manual(
    values=K3_COL,
    drop=FALSE,
    ...
  )
}

# ============================================================
# 13. STANDARD HEATMAP SCALE
# ============================================================

scale_fill_bulge_z <- function(
  limits=c(-2.5,2.5),
  name="Row z"
) {

  scale_fill_gradient2(
    low="#2166AC",
    mid="white",
    high="#B2182B",
    midpoint=0,
    limits=limits,
    oob=scales::squish,
    name=name
  )
}

# ============================================================
# 14. STANDARD DELTA SCALE
# ============================================================

scale_fill_bulge_delta <- function(
  limits,
  name="Delta"
) {

  scale_fill_gradient2(
    low="#2166AC",
    mid="white",
    high="#B2182B",
    midpoint=0,
    limits=c(
      -abs(limits),
      abs(limits)
    ),
    oob=scales::squish,
    name=name
  )
}

# ============================================================
# 15. STANDARD FIGURE DIMENSIONS
# ============================================================

FIGURE_SIZE <- list(

  full_width=c(
    width=16,
    height=10
  ),

  integrated=c(
    width=16,
    height=13
  ),

  trajectory=c(
    width=15.5,
    height=9
  ),

  landscape=c(
    width=16,
    height=12
  )
)

# ============================================================
# 16. CANONICAL EXPORT
# ============================================================

save_bulge_figure <- function(
  plot,
  stem,
  svg_dir,
  pdf_dir,
  png_dir,
  width,
  height,
  dpi=400
) {

  dir.create(
    svg_dir,
    recursive=TRUE,
    showWarnings=FALSE
  )

  dir.create(
    pdf_dir,
    recursive=TRUE,
    showWarnings=FALSE
  )

  dir.create(
    png_dir,
    recursive=TRUE,
    showWarnings=FALSE
  )

  ggsave(
    file.path(
      svg_dir,
      paste0(
        stem,
        ".svg"
      )
    ),
    plot,
    width=width,
    height=height,
    units="in",
    limitsize=FALSE
  )

  ggsave(
    file.path(
      pdf_dir,
      paste0(
        stem,
        ".pdf"
      )
    ),
    plot,
    width=width,
    height=height,
    units="in",
    limitsize=FALSE
  )

  ggsave(
    file.path(
      png_dir,
      paste0(
        stem,
        ".png"
      )
    ),
    plot,
    width=width,
    height=height,
    units="in",
    dpi=dpi,
    limitsize=FALSE
  )

  invisible(
    TRUE
  )
}

# ============================================================
# 17. SYSTEM AUDIT
# ============================================================

stopifnot(
  identical(
    K3_COL,
    c(
      T3="#2F7189",
      T2="#4B4C7A",
      T1="#B36B00"
    )
  ),

  identical(
    DISEASE_COL,
    c(
      HC="#8F8D8E",
      FFA="#4D9E4F",
      AA="#A8142A"
    )
  ),

  nrow(
    BULGE_STATE_REFERENCE
  )==3,

  BULGE_T2_BOUNDARY==0.438,

  BULGE_T1_BOUNDARY==0.844
)

cat(
  "22F_FIGURE_SYSTEM loaded successfully\n"
)

cat(
  "K3: T3 -> T2 -> T1\n"
)

cat(
  "Trunk boundaries: 0.438 / 0.844\n"
)

cat(
  "Disease palette: HC / FFA / AA frozen\n"
)
