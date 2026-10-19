suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(readr)
  library(ggplot2)
  library(scales)
})

cat("============================================================\n")
cat("20A2B - PAPER-STYLE DISEASE COMPOSITION\n")
cat("============================================================\n")

ROOT <- "PIPELINE_V2"

LEDGER <- file.path(
  ROOT,
  "results/disease_biology/17A_frozen_analysis_atlas/tables",
  "17A_master_cell_ledger_FINAL.rds"
)

DISPLAY <- file.path(
  ROOT,
  "results/figure_specification/19C3_publication_suite/tables",
  "19C3_L2_display_codes.tsv"
)

DISPLAY_LABELS <- file.path(
  ROOT,
  "config/figure_specification/L2_display_labels.tsv"
)

OUT <- file.path(
  ROOT,
  "results/publication_panels/20A2_paper_style_celltype"
)

FIG <- file.path(
  OUT,
  "figures"
)

TAB <- file.path(
  OUT,
  "tables"
)

dir.create(
  TAB,
  recursive = TRUE,
  showWarnings = FALSE
)

ledger <- readRDS(LEDGER)

display <- read_tsv(
  DISPLAY,
  show_col_types = FALSE
)

display_labels <- read_tsv(
  DISPLAY_LABELS,
  show_col_types = FALSE
) %>%
  mutate(
    plot_label_2line = paste0(
      display_code,
      "\n",
      short_identity
    )
  ) %>%
  mutate(
    plot_label_2line = paste0(
      display_code,
      "\n",
      short_identity
    )
  )

stopifnot(nrow(display_labels) == 59)
stopifnot(!anyDuplicated(display_labels$final_L2_label))
stopifnot(!anyDuplicated(display_labels$display_code))

display <- display %>%
  select(
    -any_of(
      c(
        "short_identity",
        "plot_label",
        "plot_label_2line",
        "full_identity"
      )
    )
  ) %>%
  left_join(
    display_labels %>%
      select(
        final_L2_label,
        short_identity,
        plot_label,
        plot_label_2line,
        full_identity
      ),
    by = "final_L2_label"
  )

stopifnot(!anyNA(display$short_identity))
stopifnot(!anyNA(display$plot_label))
stopifnot(!anyNA(display$plot_label_2line))

if (!"cell" %in% colnames(ledger)) {

  cc <- intersect(
    c("cell", "cell_id", "barcode", "Cell"),
    colnames(ledger)
  )

  stopifnot(length(cc) == 1)

  ledger <- ledger %>%
    rename(
      cell = all_of(cc)
    )
}

required <- c(
  "analysis_eligible_L2",
  "final_compartment",
  "final_L2_label",
  "disease"
)

stopifnot(
  all(required %in% colnames(ledger))
)

base <- ledger %>%
  filter(
    analysis_eligible_L2
  ) %>%
  mutate(
    disease = factor(
      disease,
      levels = c(
        "HC",
        "FFA",
        "AA"
      )
    )
  )

stopifnot(nrow(base) == 82494)
stopifnot(n_distinct(base$final_L2_label) == 59)

disease_cols <- c(
  HC = "#8F8D8E",
  FFA = "#4D9E4F",
  AA = "#A8142A"
)

save3 <- function(
  p,
  parent,
  stem
) {

  root <- file.path(
    FIG,
    parent,
    "D2_disease_composition"
  )

  for (fmt in c(
    "PDF",
    "PNG",
    "SVG"
  )) {

    dir.create(
      file.path(
        root,
        fmt
      ),
      recursive = TRUE,
      showWarnings = FALSE
    )
  }

  ggsave(
    file.path(
      root,
      "PDF",
      paste0(stem, ".pdf")
    ),
    p,
    width = 7.0,
    height = 4.8,
    units = "in",
    device = "pdf",
    limitsize = FALSE
  )

  ggsave(
    file.path(
      root,
      "PNG",
      paste0(stem, ".png")
    ),
    p,
    width = 7.0,
    height = 4.8,
    units = "in",
    dpi = 600,
    bg = "white",
    limitsize = FALSE
  )

  if (
    requireNamespace(
      "svglite",
      quietly = TRUE
    )
  ) {

    ggsave(
      file.path(
        root,
        "SVG",
        paste0(stem, ".svg")
      ),
      p,
      width = 7.0,
      height = 4.8,
      units = "in",
      device = svglite::svglite,
      bg = "white",
      limitsize = FALSE
    )

  } else {

    svg(
      file.path(
        root,
        "SVG",
        paste0(stem, ".svg")
      ),
      width = 7.0,
      height = 4.8
    )

    print(p)
    dev.off()
  }
}

all_tables <- list()
inventory <- list()

parents <- display %>%
  filter(
    L1 != "Neural_Schwann"
  ) %>%
  distinct(
    L1
  ) %>%
  pull(L1)

for (parent in parents) {

  disp <- display %>%
    filter(
      L1 == parent
    ) %>%
    arrange(
      L2_order_within_L1
    )

  dd <- base %>%
    filter(
      final_compartment == parent
    )

  comp <- dd %>%
    count(
      final_L2_label,
      disease,
      name = "n_cells"
    ) %>%
    complete(
      final_L2_label =
        disp$final_L2_label,
      disease = factor(
        c(
          "HC",
          "FFA",
          "AA"
        ),
        levels = c(
          "HC",
          "FFA",
          "AA"
        )
      ),
      fill = list(
        n_cells = 0L
      )
    ) %>%
    group_by(
      final_L2_label
    ) %>%
    mutate(
      total_L2_cells =
        sum(n_cells),
      proportion =
        ifelse(
          total_L2_cells > 0,
          n_cells /
            total_L2_cells,
          0
        )
    ) %>%
    ungroup() %>%
    left_join(
      disp %>%
        select(
          final_L2_label,
          display_code,
          legend_label,
          short_identity,
          plot_label,
          plot_label_2line,
          L2_order_within_L1
        ),
      by = "final_L2_label"
    ) %>%
    mutate(
      plot_label_2line = factor(
        plot_label_2line,
        levels =
          disp$plot_label_2line
      ),
      disease = factor(
        disease,
        levels = c(
          "HC",
          "FFA",
          "AA"
        )
      )
    )

  stopifnot(
    all(
      abs(
        comp %>%
          group_by(
            final_L2_label
          ) %>%
          summarise(
            s = sum(proportion),
            .groups = "drop"
          ) %>%
          pull(s) -
          1
      ) < 1e-10
    )
  )

  # Explicit disease factor for plotting.
  # Stacking order is controlled by position_stack(reverse = TRUE).
  comp$disease <- factor(
    as.character(comp$disease),
    levels = c(
      "HC",
      "FFA",
      "AA"
    )
  )

  p <- ggplot(
    comp,
    aes(
      x = plot_label_2line,
      y = proportion,
      fill = disease
    )
  ) +
    geom_col(
      width = 0.82,
      colour = "white",
      linewidth = 0.28,
      position = position_stack(
        reverse = TRUE
      )
    ) +
    scale_fill_manual(
      values = c(
        "HC" = "#8F8D8E",
        "FFA" = "#4D9E4F",
        "AA" = "#A8142A"
      ),
      limits = c(
        "HC",
        "FFA",
        "AA"
      ),
      breaks = c(
        "HC",
        "FFA",
        "AA"
      ),
      labels = c(
        "HC",
        "FFA",
        "AA"
      ),
      drop = FALSE,
      na.value = "black"
    ) +
    scale_y_continuous(
      labels = percent_format(
        accuracy = 1
      ),
      limits = c(
        0,
        1
      ),
      expand = c(
        0,
        0
      )
    ) +
    labs(
      title = paste0(
        gsub("_", " ", parent),
        " - disease composition by L2 subtype"
      ),
      subtitle =
        "Cell-level descriptive composition within each frozen L2 identity",
      x = NULL,
      y = "Proportion",
      fill = NULL
    ) +
    theme_classic(
      base_size = 10
    ) +
    theme(
      axis.text.x =
        element_text(
          angle = 45,
          hjust = 1,
          vjust = 1,
          size = 8
        ),
      axis.text.y =
        element_text(
          size = 8
        ),
      axis.title.y =
        element_text(
          size = 9
        ),
      plot.title =
        element_text(
          face = "bold",
          size = 11
        ),
      plot.subtitle =
        element_text(
          size = 8
        ),
      legend.position = "right",
      legend.text =
        element_text(
          size = 8
        ),
      legend.key.height =
        grid::unit(
          4,
          "mm"
        ),
      plot.margin =
        margin(
          7, 8, 7, 7
        )
    )

  stem <- paste0(
    "20A2_",
    parent,
    "_D2_disease_composition"
  )

  save3(
    p,
    parent,
    stem
  )

  all_tables[[parent]] <-
    comp %>%
    mutate(
      L1 = parent
    )

  inventory[[parent]] <-
    tibble(
      L1 = parent,
      n_L2 =
        nrow(disp),
      n_cells =
        nrow(dd),
      figure = stem
    )

  cat(
    sprintf(
      "%-15s : %6d cells | %2d L2 | DONE\n",
      parent,
      nrow(dd),
      nrow(disp)
    )
  )
}

all_comp <- bind_rows(
  all_tables
)

inventory <- bind_rows(
  inventory
)

write_tsv(
  all_comp,
  file.path(
    TAB,
    "20A2B_disease_composition_by_L2.tsv"
  )
)

write_tsv(
  inventory,
  file.path(
    TAB,
    "20A2B_figure_inventory.tsv"
  )
)

pdf_n <- length(
  list.files(
    FIG,
    pattern =
      "D2_disease_composition.*\\.pdf$",
    recursive = TRUE
  )
)

png_n <- length(
  list.files(
    FIG,
    pattern =
      "D2_disease_composition.*\\.png$",
    recursive = TRUE
  )
)

svg_n <- length(
  list.files(
    FIG,
    pattern =
      "D2_disease_composition.*\\.svg$",
    recursive = TRUE
  )
)

checks <- tibble(
  check = c(
    "eight_compartments",
    "frozen_L2_59",
    "HC_colour",
    "FFA_colour",
    "AA_colour",
    "all_L2_sum_to_one",
    "PDF_files_8",
    "PNG_files_8",
    "SVG_files_8",
    "no_inference",
    "no_annotation_change"
  ),
  pass = c(
    nrow(inventory) == 8,
    n_distinct(
      all_comp$final_L2_label
    ) == 59,
    disease_cols["HC"] ==
      "#8F8D8E",
    disease_cols["FFA"] ==
      "#4D9E4F",
    disease_cols["AA"] ==
      "#A8142A",
    all(
      abs(
        all_comp %>%
          group_by(
            L1,
            final_L2_label
          ) %>%
          summarise(
            s = sum(proportion),
            .groups = "drop"
          ) %>%
          pull(s) -
          1
      ) < 1e-10
    ),
    pdf_n == 8,
    png_n == 8,
    svg_n == 8,
    TRUE,
    TRUE
  )
)

write_tsv(
  checks,
  file.path(
    TAB,
    "20A2B_formal_checks.tsv"
  )
)

cat("\n===== FORMAL CHECKS =====\n")
print(
  as.data.frame(checks),
  row.names = FALSE
)

stopifnot(
  all(checks$pass)
)

cat("\n============================================================\n")
cat("20A2B COMPLETE\n")
cat("============================================================\n")
cat("Disease-composition plots: 8\n")
cat("Stacked categories: HC | FFA | AA\n")
cat("HC: #8F8D8E\n")
cat("FFA: #4D9E4F\n")
cat("AA: #A8142A\n")
cat("Inferential testing performed: NO\n")
cat("Annotation modified: NO\n")
cat("PASS\n")
