#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(readr)
  library(ggplot2)
  library(patchwork)
  library(svglite)
})

cat("============================================================\n")
cat("20B3D — PLOT ONLY: HORIZONTAL TOP20 HEATMAPS\n")
cat("============================================================\n")

ROOT <- paste0(
  "PIPELINE_V2/results/publication_panels/",
  "20B3_patient_L1_DE_heatmaps/Fibroblast"
)

DE_DIR  <- file.path(ROOT,"DE")
SRC_DIR <- file.path(ROOT,"source_data")
SEL_DIR <- file.path(ROOT,"selected_genes")
TAB_DIR <- file.path(ROOT,"tables")

SVG_DIR <- file.path(ROOT,"figures/SVG")
PDF_DIR <- file.path(ROOT,"figures/PDF")
PNG_DIR <- file.path(ROOT,"figures/PNG")

for (d in c(
  SEL_DIR,TAB_DIR,
  SVG_DIR,PDF_DIR,PNG_DIR
)) {
  dir.create(d,recursive=TRUE,showWarnings=FALSE)
}

# ============================================================
# 1. INPUTS
# ============================================================

cat("STEP 1/5 — loading existing DE and expression\n")

logCPM <- read.csv(
  file.path(
    SRC_DIR,
    "20B3_Fibroblast_patient_logCPM_TMM.csv"
  ),
  row.names=1,
  check.names=FALSE
)

meta <- read_csv(
  file.path(
    TAB_DIR,
    "20B3_Fibroblast_patient_metadata.csv"
  ),
  show_col_types=FALSE
)

contrasts <- c(
  "FFA_vs_HC",
  "AA_vs_HC",
  "AA_vs_FFA"
)

DE <- setNames(
  lapply(
    contrasts,
    function(cc) {
      read_csv(
        file.path(
          DE_DIR,
          paste0(
            "20B3_Fibroblast_",
            cc,
            "_full_DE.csv"
          )
        ),
        show_col_types=FALSE
      )
    }
  ),
  contrasts
)

stopifnot(
  ncol(logCPM)==16,
  all(contrasts %in% names(DE))
)

meta <- meta %>%
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

patient_order <- meta$patient_id

stopifnot(
  setequal(
    patient_order,
    colnames(logCPM)
  )
)

logCPM <- logCPM[
  ,
  patient_order,
  drop=FALSE
]

# ============================================================
# 2. FILTER + TOP20
# ============================================================

cat("STEP 2/5 — filtering genes and selecting Top20\n")

select_top20 <- function(res,cc) {

  eligible <- res %>%
    filter(
      !grepl("^MIR",gene,ignore.case=TRUE),
      !grepl("^LINC",gene,ignore.case=TRUE),
      !grepl("-AS[0-9]+$",gene,ignore.case=TRUE)
    )

  up <- eligible %>%
    filter(logFC>0) %>%
    arrange(PValue,desc(abs(logFC))) %>%
    slice_head(n=20) %>%
    mutate(
      direction="UP",
      rank=row_number(),
      contrast=cc
    )

  down <- eligible %>%
    filter(logFC<0) %>%
    arrange(PValue,desc(abs(logFC))) %>%
    slice_head(n=20) %>%
    mutate(
      direction="DOWN",
      rank=row_number(),
      contrast=cc
    )

  stopifnot(
    nrow(up)==20,
    nrow(down)==20
  )

  bind_rows(up,down) %>%
    mutate(
      row_id=paste(
        contrast,
        direction,
        sprintf("%02d",rank),
        gene,
        sep="__"
      )
    )
}

selected <- setNames(
  lapply(
    contrasts,
    function(cc) {
      select_top20(
        DE[[cc]],
        cc
      )
    }
  ),
  contrasts
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
    "20B3D_Fibroblast_horizontal_Top20_selected_genes.csv"
  )
)

# ============================================================
# 3. MATRICES
# ============================================================

cat("STEP 3/5 — building z-score matrices\n")

build_matrix <- function(sel) {

  out <- do.call(
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

        z <- as.numeric(scale(v))

        z[!is.finite(z)] <- 0

        pmax(
          pmin(z,2.5),
          -2.5
        )
      }
    )
  )

  rownames(out) <- sel$row_id
  colnames(out) <- patient_order

  out
}

matrices <- setNames(
  lapply(
    contrasts,
    function(cc) {
      build_matrix(
        selected[[cc]]
      )
    }
  ),
  contrasts
)

# ============================================================
# 4. PLOTS
# ============================================================

cat("STEP 4/5 — drawing three horizontal heatmaps\n")

DISEASE_COLS <- c(
  HC="#8F8D8E",
  FFA="#4D9E4F",
  AA="#A8142A"
)

make_panel <- function(cc,title_text) {

  sel <- selected[[cc]]
  mat <- matrices[[cc]]

  dat <- as.data.frame(mat) %>%
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

  ann <- meta %>%
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
    scale_x_discrete(expand=c(0,0)) +
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
      xintercept=c(5.5,10.5),
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

    theme_minimal(base_size=7) +

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

      legend.direction="horizontal",

      legend.title.position="top",

      legend.key.width=
        grid::unit(14,"mm"),

      legend.key.height=
        grid::unit(2.5,"mm"),

      plot.margin=margin(0,0,0,0)
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

# collect_guides gives ONE common colour scale.
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
    title="Fibroblast - differential expression",
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
  "20B3D_Fibroblast_horizontal_",
  "three_contrasts_Top20_FILTERED"
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

# ============================================================
# 5. AUDIT
# ============================================================

cat("STEP 5/5 — checks\n")

out_svg <- file.path(
  SVG_DIR,
  paste0(stem,".svg")
)

stopifnot(
  file.exists(out_svg),
  file.info(out_svg)$size > 0,
  nrow(selected_all)==120,
  all(
    table(
      selected_all$contrast,
      selected_all$direction
    )==20
  )
)

cat("\nHeatmaps: 3\n")
cat("Genes/heatmap: 40 (20 UP + 20 DOWN)\n")
cat("Patients/heatmap: 16\n")
cat("MIR*: excluded\n")
cat("LINC*: excluded\n")
cat("-AS#: excluded\n")
cat("Cells square: YES\n")
cat("Column clustering: NO\n")
cat("Patient order: HC -> FFA -> AA\n")
cat("SVG:",out_svg,"\n")
cat("\nSTATUS=PASS\n")
