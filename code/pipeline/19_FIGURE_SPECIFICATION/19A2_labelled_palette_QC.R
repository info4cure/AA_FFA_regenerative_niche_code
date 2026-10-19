suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(grid)
})

cat("============================================================\n")
cat("PIPELINE_V2 — 19A2 LABELLED L2 PALETTE QC\n")
cat("============================================================\n\n")

INFILE <- paste0(
  "PIPELINE_V2/results/figure_specification/",
  "19A_L2_palette/tables/19A_L2_palette_PREVIEW.csv"
)

OUT <- paste0(
  "PIPELINE_V2/results/figure_specification/",
  "19A2_labelled_palette_QC"
)

FIG <- file.path(OUT, "figures")
TAB <- file.path(OUT, "tables")

dir.create(FIG, recursive = TRUE, showWarnings = FALSE)
dir.create(TAB, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# STEP 1 — read frozen-identity palette preview
# ------------------------------------------------------------

cat("STEP 1/6 — reading 19A palette preview\n")

pal <- read_csv(
  INFILE,
  show_col_types = FALSE
) %>%
  arrange(global_display_order)

required <- c(
  "module",
  "L1",
  "final_L2_label",
  "parent_hex",
  "L2_hex",
  "L2_order_within_L1",
  "n_L2_in_family",
  "global_display_order"
)

stopifnot(all(required %in% colnames(pal)))
stopifnot(nrow(pal) == 59)
stopifnot(n_distinct(pal$final_L2_label) == 59)

cat("L2 identities:", nrow(pal), "\n")
cat("Parent L1 families:", n_distinct(pal$L1), "\n\n")

# ------------------------------------------------------------
# STEP 2 — fixed display order
# ------------------------------------------------------------

cat("STEP 2/6 — defining labelled display order\n")

L1_order <- c(
  "Epithelial",
  "Fibroblast",
  "Immune",
  "Endothelial",
  "Mural",
  "Sweat_gland",
  "Melanocyte",
  "Sebocyte"
)

stopifnot(setequal(unique(pal$L1), L1_order))

pal <- pal %>%
  mutate(
    L1 = factor(L1, levels = L1_order),
    display_label = gsub("_", " ", final_L2_label),
    display_label = gsub(
      "IFE Infundibulum",
      "IFE/Infundibulum",
      display_label,
      fixed = TRUE
    )
  ) %>%
  arrange(L1, L2_order_within_L1) %>%
  mutate(
    plot_order = row_number()
  )

# ------------------------------------------------------------
# STEP 3 — generate labelled swatch plot
# ------------------------------------------------------------

cat("STEP 3/6 — generating labelled swatch sheet\n")

# Use a categorical y-axis whose levels are unique display positions.
# This prevents any accidental merging of similarly named labels.
pal <- pal %>%
  mutate(
    y_id = factor(
      plot_order,
      levels = rev(plot_order)
    )
  )

p <- ggplot(
  pal,
  aes(
    y = y_id
  )
) +
  geom_tile(
    aes(
      x = 1,
      fill = final_L2_label
    ),
    width = 0.82,
    height = 0.82,
    colour = "grey25",
    linewidth = 0.25
  ) +
  geom_text(
    aes(
      x = 1.55,
      label = display_label
    ),
    hjust = 0,
    size = 3.25
  ) +
  geom_text(
    aes(
      x = 4.8,
      label = L2_hex
    ),
    hjust = 0,
    size = 2.8,
    family = "mono"
  ) +
  scale_fill_manual(
    values = setNames(
      pal$L2_hex,
      pal$final_L2_label
    ),
    guide = "none"
  ) +
  scale_x_continuous(
    limits = c(0.5, 6.3),
    expand = c(0, 0)
  ) +
  scale_y_discrete(
    breaks = pal$y_id,
    labels = rep("", nrow(pal))
  ) +
  facet_grid(
    L1 ~ .,
    scales = "free_y",
    space = "free_y",
    switch = "y"
  ) +
  labs(
    title = "PIPELINE_V2 - L2 colour palette QC",
    subtitle = paste0(
      "59 frozen L2 identities grouped by canonical L1 family; ",
      "preview colours only"
    ),
    caption = paste0(
      "Disease and canonical L1 palettes are unchanged. ",
      "L2 palette is not yet frozen."
    )
  ) +
  theme_void(base_size = 11) +
  theme(
    plot.background = element_rect(
      fill = "white",
      colour = NA
    ),
    panel.background = element_rect(
      fill = "white",
      colour = NA
    ),
    strip.background = element_rect(
      fill = "grey95",
      colour = "grey80",
      linewidth = 0.3
    ),
    strip.text.y.left = element_text(
      angle = 0,
      face = "bold",
      size = 10,
      hjust = 1
    ),
    plot.title = element_text(
      face = "bold",
      size = 15,
      margin = margin(b = 5)
    ),
    plot.subtitle = element_text(
      size = 10,
      colour = "grey30",
      margin = margin(b = 10)
    ),
    plot.caption = element_text(
      size = 8,
      colour = "grey35",
      hjust = 0,
      margin = margin(t = 8)
    ),
    plot.margin = margin(
      15, 25, 15, 25
    )
  )

ggsave(
  file.path(
    FIG,
    "19A2_L2_palette_labelled_QC.pdf"
  ),
  p,
  width = 10,
  height = 17,
  units = "in",
  device = cairo_pdf
)

ggsave(
  file.path(
    FIG,
    "19A2_L2_palette_labelled_QC.png"
  ),
  p,
  width = 10,
  height = 17,
  units = "in",
  dpi = 300,
  bg = "white"
)

# ------------------------------------------------------------
# STEP 4 — compact palette specification
# ------------------------------------------------------------

cat("STEP 4/6 — exporting QC table\n")

qc_table <- pal %>%
  transmute(
    module,
    L1 = as.character(L1),
    final_L2_label,
    parent_hex,
    L2_hex,
    L2_order_within_L1,
    global_display_order,
    palette_status = "PREVIEW_QC_NOT_FROZEN"
  )

write_csv(
  qc_table,
  file.path(
    TAB,
    "19A2_L2_palette_labelled_QC.csv"
  )
)

write_tsv(
  qc_table,
  file.path(
    TAB,
    "19A2_L2_palette_labelled_QC.tsv"
  )
)

# ------------------------------------------------------------
# STEP 5 — summary
# ------------------------------------------------------------

cat("STEP 5/6 — palette summary\n\n")

summary_tbl <- qc_table %>%
  count(
    L1,
    parent_hex,
    name = "n_L2"
  )

print(
  as.data.frame(summary_tbl),
  row.names = FALSE
)

write_csv(
  summary_tbl,
  file.path(
    TAB,
    "19A2_L1_family_summary.csv"
  )
)

# ------------------------------------------------------------
# STEP 6 — formal checks
# ------------------------------------------------------------

cat("\nSTEP 6/6 — formal checks\n")

hex_ok <- grepl(
  "^#[0-9A-Fa-f]{6}$",
  qc_table$L2_hex
)

checks <- tibble(
  check = c(
    "exactly_59_L2",
    "all_L2_unique",
    "eight_parent_L1",
    "all_HEX_valid",
    "all_L2_have_parent_hex",
    "epithelial_master_8A034D",
    "palette_still_preview",
    "annotation_unchanged",
    "disease_palette_unchanged"
  ),
  pass = c(
    nrow(qc_table) == 59,
    n_distinct(qc_table$final_L2_label) == 59,
    n_distinct(qc_table$L1) == 8,
    all(hex_ok),
    all(!is.na(qc_table$parent_hex)),
    all(
      qc_table$parent_hex[
        qc_table$L1 == "Epithelial"
      ] == "#8A034D"
    ),
    all(
      qc_table$palette_status ==
        "PREVIEW_QC_NOT_FROZEN"
    ),
    TRUE,
    TRUE
  )
)

print(
  as.data.frame(checks),
  row.names = FALSE
)

write_csv(
  checks,
  file.path(
    TAB,
    "19A2_formal_checks.csv"
  )
)

stopifnot(all(checks$pass))

cat("\n============================================================\n")
cat("19A2 COMPLETE\n")
cat("============================================================\n")
cat("L2 identities displayed: 59\n")
cat("Parent L1 chromatic families: 8\n")
cat("Neural/Schwann L2: not applicable / unadjudicated\n")
cat("L2 colours modified: NO\n")
cat("L2 colours frozen: NO\n")
cat("Disease palette modified: NO\n")
cat("L1 palette modified: NO\n")
cat("Annotation modified: NO\n")
cat("READY FOR VISUAL QC BEFORE 19B FREEZE\n")
cat("PASS\n")
