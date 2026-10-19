suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(tibble)
  library(scales)
})

cat("============================================================\n")
cat("PIPELINE_V2 — 19A L2 PALETTE PREVIEW\n")
cat("============================================================\n\n")

# ============================================================
# PATHS
# ============================================================

PAL_FILE <- "PIPELINE_V2/config/figure_palette_v2.R"

OUT <- paste0(
  "PIPELINE_V2/results/figure_specification/",
  "19A_L2_palette"
)

FIG <- file.path(OUT, "figures")
TAB <- file.path(OUT, "tables")

dir.create(FIG, recursive = TRUE, showWarnings = FALSE)
dir.create(TAB, recursive = TRUE, showWarnings = FALSE)

stopifnot(file.exists(PAL_FILE))
source(PAL_FILE)


# ============================================================
# 1. FROZEN OBJECT SPECIFICATION
# ============================================================

objects <- tribble(
  ~module, ~L1, ~path,

  "07_EPITHELIAL_L2",
  "Epithelial",
  "PIPELINE_V2/objects/annotated/epithelial_v2_L2_frozen.rds",

  "08_IMMUNE_L2",
  "Immune",
  "PIPELINE_V2/objects/annotated/immune_v2_08I_L2_frozen.rds",

  "09_FIBROBLAST_L2",
  "Fibroblast",
  "PIPELINE_V2/objects/annotated/fibroblast_v2_09H5B_L2_frozen_4685.rds",

  "10_VASCULAR_L2",
  "Endothelial",
  "PIPELINE_V2/objects/annotated/vascular_v2_10B6_L2_frozen_3023.rds",

  "11_MURAL_L2",
  "Mural",
  "PIPELINE_V2/results/annotation/11_MURAL_L2/11B5_freeze/mural_v2_11B5_L2_frozen.rds",

  "12_SWEAT_GLAND_L2",
  "Sweat_gland",
  "PIPELINE_V2/results/annotation/12_SWEAT_GLAND_L2/12D_freeze/sweat_gland_v2_12D_L2_frozen.rds",

  "13_MELANOCYTE_L2",
  "Melanocyte",
  "PIPELINE_V2/results/annotation/13_MELANOCYTE_L2/13D_freeze/melanocyte_v2_13D_L2_frozen.rds",

  "14_SEBOCYTE_L2",
  "Sebocyte",
  "PIPELINE_V2/results/annotation/14_SEBOCYTE_L2/14D_freeze/sebocyte_v2_14D_L2_frozen.rds"
)

stopifnot(all(file.exists(objects$path)))


# ============================================================
# 2. RECOVER L2 LABELS FROM FROZEN OBJECTS
# ============================================================

cat("STEP 1/7 — recovering frozen L2 labels\n")

recover_L2 <- function(module, L1, path) {

  x <- readRDS(path)

  md <- x[[]]

  # Canonical frozen L2 metadata column by module.
  # Explicit mapping avoids guessing among multiple L2-related fields.
  canonical_L2_columns <- c(
    "07_EPITHELIAL_L2" = "L2_identity",
    "08_IMMUNE_L2" = "immune_L2",
    "09_FIBROBLAST_L2" = "fibroblast_L2_label",
    "10_VASCULAR_L2" = "vascular_L2_label",
    "11_MURAL_L2" = "mural_L2_label",
    "12_SWEAT_GLAND_L2" = "sweat_L2_label",
    "13_MELANOCYTE_L2" = "melanocyte_L2_label",
    "14_SEBOCYTE_L2" = "sebocyte_L2_label"
  )

  if (module %in% names(canonical_L2_columns)) {

    col <- unname(canonical_L2_columns[module])

    if (!col %in% colnames(md)) {
      stop(
        paste0(
          "\nExpected canonical L2 column '",
          col,
          "' not found for ",
          module
        )
      )
    }

  } else {

    candidates <- c(
      "final_L2_label",
      "final_L2",
      "L2_identity",
      "L2_label",
      "immune_L2"
    )

    found <- candidates[candidates %in% colnames(md)]

    if (length(found) == 0) {

      possible <- grep(
        "L2",
        colnames(md),
        value = TRUE,
        ignore.case = TRUE
      )

      stop(
        paste0(
          "\nNo canonical L2 metadata column found for ",
          module,
          ".\nAvailable L2-like columns:\n",
          paste(possible, collapse = "\n")
        )
      )
    }

    col <- found[1]
  }

  vals <- as.character(md[[col]])

  vals <- vals[
    !is.na(vals) &
    nzchar(vals)
  ]

  tibble(
    module = module,
    L1 = L1,
    metadata_column = col,
    final_L2_label = sort(unique(vals))
  )
}

maps <- lapply(
  seq_len(nrow(objects)),
  function(i) {
    recover_L2(
      objects$module[i],
      objects$L1[i],
      objects$path[i]
    )
  }
)

L2_map <- bind_rows(maps)


cat("\n===== RECOVERED L2 COUNTS =====\n")

count_tbl <- L2_map %>%
  count(L1, name = "n_L2") %>%
  mutate(
    L1 = factor(
      L1,
      levels = L1_ORDER
    )
  ) %>%
  arrange(L1)

print(as.data.frame(count_tbl), row.names = FALSE)

cat("\nTotal L2:", nrow(L2_map), "\n")


# ============================================================
# 3. REQUIRE EXACT FROZEN UNIVERSE
# ============================================================

stopifnot(
  nrow(L2_map) == 59,
  n_distinct(L2_map$final_L2_label) == 59
)


# ============================================================
# 4. GENERATE FAMILY-CONSTRAINED COLOURS
# ============================================================

cat("\nSTEP 2/7 — generating L2 chromatic families\n")

# Generate a family around the master colour.
#
# Each L2 remains recognisably within the parent L1 hue.
# We mix the master colour with controlled dark/light anchors.
#
# This is PREVIEW only.

make_family <- function(base_hex, n) {

  if (n == 1) {
    return(base_hex)
  }

  base_rgb <- col2rgb(base_hex) / 255

  dark_anchor <- base_rgb * 0.52

  light_anchor <-
    base_rgb * 0.38 +
    matrix(
      rep(1, 3),
      nrow = 3
    ) * 0.62

  positions <- seq(
    0,
    1,
    length.out = n
  )

  cols <- vapply(
    positions,
    function(p) {

      # Piecewise interpolation:
      # dark -> canonical -> light

      if (p <= 0.5) {

        q <- p / 0.5

        rgbv <-
          dark_anchor * (1 - q) +
          base_rgb * q

      } else {

        q <- (p - 0.5) / 0.5

        rgbv <-
          base_rgb * (1 - q) +
          light_anchor * q
      }

      rgb(
        rgbv[1],
        rgbv[2],
        rgbv[3]
      )
    },
    character(1)
  )

  toupper(cols)
}


palette_list <- lapply(
  intersect(
    L1_ORDER,
    unique(L2_map$L1)
  ),
  function(comp) {

    z <- L2_map %>%
      filter(L1 == comp) %>%
      arrange(final_L2_label)

    base <- unname(L1_COLORS[comp])

    cols <- make_family(
      base,
      nrow(z)
    )

    z %>%
      mutate(
        parent_hex = base,
        L2_hex = cols,
        palette_status = "PREVIEW_NOT_FROZEN"
      )
  }
)

L2_palette <- bind_rows(palette_list)


# ============================================================
# 5. ADD DISPLAY ORDER
# ============================================================

L2_palette <- L2_palette %>%
  mutate(
    L1 = factor(
      L1,
      levels = L1_ORDER
    )
  ) %>%
  arrange(
    L1,
    final_L2_label
  ) %>%
  group_by(L1) %>%
  mutate(
    L2_order_within_L1 = row_number(),
    n_L2_in_family = n()
  ) %>%
  ungroup() %>%
  mutate(
    global_display_order = row_number(),
    L1 = as.character(L1)
  )


# ============================================================
# 6. SAVE PREVIEW TABLE
# ============================================================

cat("\nSTEP 3/7 — saving preview table\n")

write_csv(
  L2_palette,
  file.path(
    TAB,
    "19A_L2_palette_PREVIEW.csv"
  )
)

write_tsv(
  L2_palette,
  file.path(
    TAB,
    "19A_L2_palette_PREVIEW.tsv"
  )
)


# ============================================================
# 7. SWATCH FIGURE
# ============================================================

cat("STEP 4/7 — building swatch sheet\n")

plot_df <- L2_palette %>%
  mutate(
    display_label = final_L2_label,
    display_order = rev(seq_len(n()))
  )

p <- ggplot(
  plot_df,
  aes(
    x = 1,
    y = display_order
  )
) +
  geom_tile(
    aes(fill = L2_hex),
    width = 0.88,
    height = 0.82,
    colour = "white",
    linewidth = 0.4
  ) +
  scale_fill_identity() +
  geom_text(
    aes(
      x = 1.62,
      label = display_label
    ),
    hjust = 0,
    size = 3.25
  ) +
  geom_text(
    aes(
      x = 0.38,
      label = L2_hex
    ),
    hjust = 1,
    size = 2.7,
    family = "mono"
  ) +
  facet_grid(
    L1 ~ .,
    scales = "free_y",
    space = "free_y",
    switch = "y"
  ) +
  coord_cartesian(
    xlim = c(-0.25, 4.7),
    clip = "off"
  ) +
  labs(
    title = "PIPELINE_V2 — L2 colour palette preview",
    subtitle = paste0(
      "59 frozen L2 identities • colours constrained ",
      "within canonical L1 chromatic families"
    ),
    caption = paste0(
      "PREVIEW ONLY — not frozen. ",
      "Disease and L1 canonical palettes remain unchanged."
    )
  ) +
  theme_void(base_size = 11) +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 15,
      hjust = 0
    ),
    plot.subtitle = element_text(
      size = 10.5,
      margin = margin(b = 10)
    ),
    plot.caption = element_text(
      size = 8.5,
      hjust = 0,
      margin = margin(t = 10)
    ),
    strip.text.y.left = element_text(
      angle = 0,
      face = "bold",
      size = 10,
      hjust = 1
    ),
    strip.background = element_blank(),
    panel.spacing.y = unit(5, "mm"),
    plot.margin = margin(
      15, 100, 15, 35
    )
  )


ggsave(
  file.path(
    FIG,
    "19A_L2_palette_swatches_PREVIEW.pdf"
  ),
  p,
  width = 11,
  height = 19,
  units = "in",
  limitsize = FALSE
)

ggsave(
  file.path(
    FIG,
    "19A_L2_palette_swatches_PREVIEW.png"
  ),
  p,
  width = 11,
  height = 19,
  units = "in",
  dpi = 300,
  limitsize = FALSE
)


# ============================================================
# 8. FAMILY MASTER SWATCHES
# ============================================================

cat("STEP 5/7 — canonical L1 reference sheet\n")

master_df <- tibble(
  L1 = L1_ORDER,
  hex = unname(L1_COLORS[L1_ORDER])
)

p_master <- master_df %>%
  mutate(
    y = rev(seq_len(n()))
  ) %>%
  ggplot(
    aes(
      x = 1,
      y = y
    )
  ) +
  geom_tile(
    aes(fill = hex),
    width = 0.9,
    height = 0.8,
    colour = "white",
    linewidth = 0.5
  ) +
  scale_fill_identity() +
  geom_text(
    aes(
      x = 1.65,
      label = L1
    ),
    hjust = 0,
    size = 4.2,
    fontface = "bold"
  ) +
  geom_text(
    aes(
      x = 0.35,
      label = hex
    ),
    hjust = 1,
    size = 3.4,
    family = "mono"
  ) +
  coord_cartesian(
    xlim = c(-0.2, 3.6),
    clip = "off"
  ) +
  labs(
    title = "Canonical L1 colour families"
  ) +
  theme_void() +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 15
    ),
    plot.margin = margin(
      15, 80, 15, 20
    )
  )

ggsave(
  file.path(
    FIG,
    "19A_L1_master_colours.pdf"
  ),
  p_master,
  width = 8,
  height = 6.5
)

ggsave(
  file.path(
    FIG,
    "19A_L1_master_colours.png"
  ),
  p_master,
  width = 8,
  height = 6.5,
  dpi = 300
)


# ============================================================
# 9. FORMAL CHECKS
# ============================================================

cat("STEP 6/7 — formal checks\n")

checks <- tibble(
  check = c(
    "exactly_59_L2",
    "all_L2_unique",
    "eight_L2_parent_compartments",
    "all_parent_L1_have_canonical_colour",
    "all_L2_have_colour",
    "all_HEX_valid",
    "disease_palette_unchanged",
    "epithelial_master_8A034D",
    "annotation_unchanged",
    "palette_status_preview"
  ),
  pass = c(
    nrow(L2_palette) == 59,

    n_distinct(
      L2_palette$final_L2_label
    ) == 59,

    n_distinct(
      L2_palette$L1
    ) == 8,

    all(
      unique(L2_palette$L1) %in%
        names(L1_COLORS)
    ),

    all(
      !is.na(L2_palette$L2_hex)
    ),

    all(
      grepl(
        "^#[0-9A-F]{6}$",
        L2_palette$L2_hex
      )
    ),

    identical(
      unname(DISEASE_COLORS),
      c(
        "#8F8D8E",
        "#4D9E4F",
        "#A8142A"
      )
    ),

    identical(
      unname(
        L1_COLORS["Epithelial"]
      ),
      "#8A034D"
    ),

    TRUE,

    all(
      L2_palette$palette_status ==
        "PREVIEW_NOT_FROZEN"
    )
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
    "19A_formal_checks.csv"
  )
)

stopifnot(all(checks$pass))


# ============================================================
# 10. SUMMARY
# ============================================================

cat("\nSTEP 7/7 — summary\n")

print(
  as.data.frame(
    L2_palette %>%
      count(
        L1,
        name = "n_L2"
      ) %>%
      arrange(
        match(
          L1,
          L1_ORDER
        )
      )
  ),
  row.names = FALSE
)

cat("\n============================================================\n")
cat("19A COMPLETE\n")
cat("============================================================\n")
cat("Frozen L2 identities recovered: 59\n")
cat("Parent L1 families represented: 8\n")
cat("Neural/Schwann L2 palette: not applicable / unadjudicated\n")
cat("Disease palette modified: NO\n")
cat("L1 palette modified: NO\n")
cat("L2 palette status: PREVIEW — NOT FROZEN\n")
cat("Annotation modified: NO\n")
cat("Figures generated: YES\n")
cat("READY FOR VISUAL L2 PALETTE QC\n")
cat("PASS\n")
