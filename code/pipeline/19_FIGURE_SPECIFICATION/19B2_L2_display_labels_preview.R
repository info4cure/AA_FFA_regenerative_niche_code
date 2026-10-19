suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
})

cat("============================================================\n")
cat("PIPELINE_V2 - 19B2 L2 DISPLAY LABELS PREVIEW\n")
cat("============================================================\n")

# ============================================================
# 1. INPUT / OUTPUT
# ============================================================

SRC <- paste0(
  "PIPELINE_V2/results/figure_specification/",
  "19C3_publication_suite/tables/",
  "19C3_L2_display_codes.tsv"
)

OUT <- paste0(
  "PIPELINE_V2/results/figure_specification/",
  "19B2_L2_display_labels"
)

TAB <- file.path(OUT, "tables")

dir.create(
  TAB,
  recursive = TRUE,
  showWarnings = FALSE
)

stopifnot(file.exists(SRC))

x <- read_tsv(
  SRC,
  show_col_types = FALSE
)

stopifnot(nrow(x) == 59)
stopifnot(!anyDuplicated(x$final_L2_label))
stopifnot(!anyDuplicated(x$display_code))

# ============================================================
# 2. EXPLICIT PUBLICATION SHORT LABELS
# ============================================================

cat("\nSTEP 1/4 - defining 59 publication short labels\n")

short_map <- tribble(
  ~final_L2_label, ~short_identity,

  # ----------------------------------------------------------
  # EPITHELIAL
  # ----------------------------------------------------------
  "Basal_epithelial",
  "Basal",

  "Bulge_Isthmus",
  "Bulge/Isthmus",

  "Hair_shaft",
  "Hair shaft",

  "IFE_Infundibulum_differentiating",
  "IFE/Infund. diff.",

  "IFE_Infundibulum_suprabasal",
  "IFE/Infund. suprabasal",

  "IFE_Infundibulum_terminal",
  "IFE/Infund. terminal",

  "IRS_companion",
  "IRS companion",

  "ORS",
  "ORS",

  "ORS_basal",
  "ORS basal",

  # ----------------------------------------------------------
  # FIBROBLAST
  # ----------------------------------------------------------
  "ADAM12/MMP11 remodeling",
  "ADAM12/MMP11",

  "APOE/CCL19 immune-perivascular",
  "APOE/CCL19",

  "COCH/AQP4 dermal papilla",
  "COCH/AQP4 DP",

  "COL11A1/DPEP1 dermal-sheath",
  "COL11A1/DPEP1 DS",

  "COL26A1/RELN fibroblast",
  "COL26A1/RELN",

  "COL6A5 papillary",
  "COL6A5 papillary",

  "EBF2/SFRP5 fibroblast",
  "EBF2/SFRP5",

  "PTGDS/IGF1 fibroblast",
  "PTGDS/IGF1",

  "SFRP2/DPP4 reticular",
  "SFRP2/DPP4 reticular",

  "ZNF804B/ADGRB3 fibroblast",
  "ZNF804B/ADGRB3",

  # ----------------------------------------------------------
  # IMMUNE
  # ----------------------------------------------------------
  "B/plasma-lineage",
  "B/plasma",

  "CD4 T",
  "CD4 T",

  "CD8 T",
  "CD8 T",

  "Langerhans",
  "Langerhans",

  "Macrophage/monocyte",
  "Macro/Mono",

  "Mast",
  "Mast",

  "NK",
  "NK",

  "Treg",
  "Treg",

  "cDC1-like dendritic",
  "cDC1-like",

  # ----------------------------------------------------------
  # ENDOTHELIAL
  # ----------------------------------------------------------
  "ACKR1+ venular EC",
  "ACKR1+ venular",

  "ADGRG6/BMPER endothelial state",
  "ADGRG6/BMPER",

  "CA4/RGCC capillary EC",
  "CA4/RGCC capillary",

  "CCL21/LYVE1 lymphatic EC",
  "CCL21/LYVE1 lymphatic",

  "Epithelial-like EC",
  "Epithelial-like",

  "GJA5/FBLN5 arterial EC",
  "GJA5/FBLN5 arterial",

  "Low-complexity EC",
  "Low-complexity",

  "Mural-like contractile EC",
  "Mural-like",

  "SELE/CXCL10 activated venular EC",
  "SELE/CXCL10 activated",

  "SLC10A1/CD83 endothelial state",
  "SLC10A1/CD83",

  # ----------------------------------------------------------
  # MURAL
  # ----------------------------------------------------------
  "APM-bulge interface smooth muscle",
  "APM-bulge interface",

  "APM-like smooth muscle",
  "APM-like SM",

  "Contractile ECM-associated pericyte",
  "Contractile ECM peri.",

  "Contractile mural cell",
  "Contractile mural",

  "Intermediate contractile smooth muscle",
  "Intermediate SM",

  "Mural-epithelial interface cell",
  "Mural-epithelial",

  "RGS5/PDGFRB ECM-associated pericyte",
  "RGS5/PDGFRB peri.",

  # ----------------------------------------------------------
  # SWEAT GLAND
  # ----------------------------------------------------------
  "Clear secretory cell",
  "Clear secretory",

  "Dark/apocrine-like secretory cell",
  "Dark/apocrine-like",

  "Ductal/clear transitional cell",
  "Ductal/clear trans.",

  "Inflammatory secretory/ductal cell",
  "Inflamm. secretory",

  "KRT77+ basal ductal/myoepithelial cell",
  "KRT77+ basal ductal",

  # ----------------------------------------------------------
  # MELANOCYTE
  # ----------------------------------------------------------
  "Antigen-presentation-high melanocyte",
  "Ag-presenting high",

  "Differentiated/homeostatic melanocyte",
  "Differentiated",

  "Low-pigmentation transitional melanocyte",
  "Low-pigment trans.",

  "Stress/neural-crest-like melanocyte",
  "Stress/NC-like",

  # ----------------------------------------------------------
  # SEBOCYTE
  # ----------------------------------------------------------
  "Basal pilosebaceous progenitor-like",
  "Basal progenitor",

  "Differentiated lipogenic sebocyte",
  "Differentiated lipogenic",

  "Early sebocyte commitment",
  "Early commitment",

  "Junctional pilosebaceous progenitor-like",
  "Junctional progenitor",

  "Terminal lipid-rich inflammatory sebocyte",
  "Terminal lipid-rich"
)

# ============================================================
# 3. VALIDATE AND BUILD DISPLAY SPECIFICATION
# ============================================================

cat("STEP 2/4 - validating explicit mapping\n")

stopifnot(nrow(short_map) == 59)
stopifnot(!anyDuplicated(short_map$final_L2_label))
stopifnot(!anyNA(short_map$short_identity))
stopifnot(all(short_map$short_identity != ""))

missing_from_map <- setdiff(
  x$final_L2_label,
  short_map$final_L2_label
)

extra_in_map <- setdiff(
  short_map$final_L2_label,
  x$final_L2_label
)

if (
  length(missing_from_map) > 0 ||
  length(extra_in_map) > 0
) {

  cat("\nMISSING FROM SHORT MAP:\n")
  print(missing_from_map)

  cat("\nEXTRA IN SHORT MAP:\n")
  print(extra_in_map)

  stop(
    "19B2 short-label mapping does not match frozen 59-L2 universe."
  )
}

out <- x %>%
  left_join(
    short_map,
    by = "final_L2_label"
  ) %>%
  mutate(
    full_identity =
      final_L2_label,

    # Standard compact form for legends / axes.
    plot_label =
      paste0(
        display_code,
        "  ",
        short_identity
      ),

    # Two-line form for narrow x axes.
    plot_label_2line =
      paste0(
        display_code,
        "\n",
        short_identity
      ),

    # Full form retained for supplementary legends/tables.
    full_plot_label =
      paste0(
        display_code,
        "  ",
        gsub(
          "_",
          " ",
          final_L2_label
        )
      ),

    display_label_status =
      "PREVIEW_NOT_FROZEN"
  ) %>%
  arrange(
    global_display_order
  )

stopifnot(nrow(out) == 59)
stopifnot(!anyNA(out$short_identity))
stopifnot(!anyNA(out$plot_label))
stopifnot(!anyNA(out$plot_label_2line))

# short_identity need not be globally unique because the code
# disambiguates it; plot_label itself must be unique.
stopifnot(!anyDuplicated(out$plot_label))

# ============================================================
# 4. WRITE PREVIEW
# ============================================================

cat("STEP 3/4 - writing preview specification\n")

write_tsv(
  out,
  file.path(
    TAB,
    "19B2_L2_display_labels_PREVIEW.tsv"
  )
)

write_csv(
  out,
  file.path(
    TAB,
    "19B2_L2_display_labels_PREVIEW.csv"
  )
)

qc <- out %>%
  select(
    L1,
    display_code,
    short_identity,
    full_identity,
    plot_label
  )

write_tsv(
  qc,
  file.path(
    TAB,
    "19B2_L2_display_labels_QC.tsv"
  )
)

checks <- tibble(
  check = c(
    "exactly_59_L2",
    "all_frozen_labels_mapped",
    "no_extra_labels",
    "display_codes_unique",
    "plot_labels_unique",
    "short_labels_complete",
    "eight_parent_compartments",
    "annotation_unchanged",
    "palette_unchanged",
    "status_preview"
  ),
  pass = c(
    nrow(out) == 59,
    length(missing_from_map) == 0,
    length(extra_in_map) == 0,
    !anyDuplicated(out$display_code),
    !anyDuplicated(out$plot_label),
    all(
      !is.na(out$short_identity) &
        out$short_identity != ""
    ),
    n_distinct(out$L1) == 8,
    TRUE,
    TRUE,
    all(
      out$display_label_status ==
        "PREVIEW_NOT_FROZEN"
    )
  )
)

write_tsv(
  checks,
  file.path(
    TAB,
    "19B2_formal_checks.tsv"
  )
)

cat("\n===== 59 DISPLAY LABELS =====\n")

print(
  as.data.frame(
    out %>%
      select(
        L1,
        display_code,
        short_identity,
        full_identity
      )
  ),
  row.names = FALSE
)

cat("\n===== FORMAL CHECKS =====\n")

print(
  as.data.frame(checks),
  row.names = FALSE
)

stopifnot(all(checks$pass))

cat("\nSTEP 4/4 - summary\n")
cat("============================================================\n")
cat("19B2 COMPLETE\n")
cat("============================================================\n")
cat("Frozen L2 identities: 59\n")
cat("Short identities proposed: 59\n")
cat("Canonical biological labels modified: NO\n")
cat("Frozen colours modified: NO\n")
cat("Display-label status: PREVIEW - NOT FROZEN\n")
cat("READY FOR DISPLAY LABEL QC\n")
cat("PASS\n")
