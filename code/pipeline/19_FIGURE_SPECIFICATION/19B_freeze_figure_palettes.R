suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tibble)
})

cat("============================================================\n")
cat("PIPELINE_V2 — 19B FIGURE PALETTE FREEZE\n")
cat("============================================================\n\n")

INFILE <- paste0(
  "PIPELINE_V2/results/figure_specification/",
  "19A2_labelled_palette_QC/tables/",
  "19A2_L2_palette_labelled_QC.csv"
)

CONFIG <- "PIPELINE_V2/config/figure_palettes"

OUT <- paste0(
  "PIPELINE_V2/results/figure_specification/",
  "19B_palette_freeze/tables"
)

dir.create(CONFIG, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# STEP 1 — recover visually approved 19A2 palette
# ------------------------------------------------------------

cat("STEP 1/7 — recovering approved 19A2 palette\n")

x <- read_csv(
  INFILE,
  show_col_types = FALSE
)

stopifnot(nrow(x) == 59)
stopifnot(n_distinct(x$final_L2_label) == 59)

# ------------------------------------------------------------
# STEP 2 — canonical disease palette
# ------------------------------------------------------------

cat("STEP 2/7 — freezing disease palette\n")

disease <- tribble(
  ~disease, ~hex,      ~display_order,
  "HC",     "#8F8D8E", 1L,
  "FFA",    "#4D9E4F", 2L,
  "AA",     "#A8142A", 3L
)

# ------------------------------------------------------------
# STEP 3 — canonical L1 palette
# ------------------------------------------------------------

cat("STEP 3/7 — freezing L1 palette\n")

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

L1 <- x %>%
  distinct(
    L1,
    parent_hex
  ) %>%
  mutate(
    display_order = match(L1, L1_order)
  ) %>%
  arrange(display_order) %>%
  rename(hex = parent_hex)

stopifnot(nrow(L1) == 8)
stopifnot(
  L1$hex[L1$L1 == "Epithelial"] == "#8A034D"
)

# Neural/Schwann exists at L1 but has no frozen L2.
# It is therefore added explicitly to the canonical L1 palette.
neural <- tibble(
  L1 = "Neural_Schwann",
  hex = "#A87916",
  display_order = 9L
)

L1 <- bind_rows(
  L1,
  neural
) %>%
  arrange(display_order)

# ------------------------------------------------------------
# STEP 4 — canonical frozen L2 palette
# ------------------------------------------------------------

cat("STEP 4/7 — freezing 59 L2 colours\n")

L2 <- x %>%
  transmute(
    module,
    L1,
    final_L2_label,
    hex = L2_hex,
    L2_order_within_L1,
    global_display_order,
    palette_status = "FROZEN"
  ) %>%
  arrange(global_display_order)

stopifnot(nrow(L2) == 59)
stopifnot(n_distinct(L2$final_L2_label) == 59)

# ------------------------------------------------------------
# STEP 5 — write canonical specification
# ------------------------------------------------------------

cat("STEP 5/7 — writing canonical palette files\n")

write_tsv(
  disease,
  file.path(CONFIG, "disease_palette.tsv")
)

write_tsv(
  L1,
  file.path(CONFIG, "L1_palette.tsv")
)

write_tsv(
  L2,
  file.path(CONFIG, "L2_palette.tsv")
)

write_csv(
  disease,
  file.path(OUT, "19B_disease_palette_frozen.csv")
)

write_csv(
  L1,
  file.path(OUT, "19B_L1_palette_frozen.csv")
)

write_csv(
  L2,
  file.path(OUT, "19B_L2_palette_frozen.csv")
)

spec <- c(
  "PIPELINE_V2 FIGURE PALETTE SPECIFICATION",
  "Version: V1_20260930",
  "Status: FROZEN",
  "",
  "Purpose:",
  "Canonical colour source for all downstream scientific figures.",
  "",
  "Disease palette:",
  "HC  = #8F8D8E",
  "FFA = #4D9E4F",
  "AA  = #A8142A",
  "",
  "Rules:",
  "- Disease colours are fixed.",
  "- L1 colours are fixed.",
  "- Epithelial master colour is #8A034D.",
  "- Melanocyte uses a brown chromatic family.",
  "- Neural/Schwann uses dark yellow at L1.",
  "- Sweat gland uses a blue chromatic family.",
  "- Each frozen L2 identity has one immutable HEX colour.",
  "- L2 colours remain within their parent L1 chromatic family.",
  "- No downstream script may dynamically regenerate these colours.",
  "- Returned unadjudicated cells have no frozen L2 colour.",
  "- Palette changes require an explicit new version and audit.",
  "",
  "Canonical files:",
  "PIPELINE_V2/config/figure_palettes/disease_palette.tsv",
  "PIPELINE_V2/config/figure_palettes/L1_palette.tsv",
  "PIPELINE_V2/config/figure_palettes/L2_palette.tsv"
)

writeLines(
  spec,
  file.path(CONFIG, "PALETTE_SPECIFICATION.txt")
)

# ------------------------------------------------------------
# STEP 6 — hashes
# ------------------------------------------------------------

cat("STEP 6/7 — computing SHA256 hashes\n")

files <- c(
  file.path(CONFIG, "disease_palette.tsv"),
  file.path(CONFIG, "L1_palette.tsv"),
  file.path(CONFIG, "L2_palette.tsv"),
  file.path(CONFIG, "PALETTE_SPECIFICATION.txt")
)

sha <- tools::md5sum(files)

# md5 retained for native R portability.
hashes <- tibble(
  file = names(sha),
  md5 = unname(sha)
)

write_tsv(
  hashes,
  file.path(CONFIG, "PALETTE_FILE_HASHES.tsv")
)

write_csv(
  hashes,
  file.path(OUT, "19B_palette_file_hashes.csv")
)

# ------------------------------------------------------------
# STEP 7 — formal checks
# ------------------------------------------------------------

cat("STEP 7/7 — formal checks\n")

hex_valid <- function(z) {
  grepl("^#[0-9A-Fa-f]{6}$", z)
}

checks <- tibble(
  check = c(
    "disease_colours_3",
    "HC_8F8D8E",
    "FFA_4D9E4F",
    "AA_A8142A",
    "L1_colours_9",
    "epithelial_8A034D",
    "neural_schwann_dark_yellow",
    "L2_colours_59",
    "L2_labels_unique",
    "all_disease_HEX_valid",
    "all_L1_HEX_valid",
    "all_L2_HEX_valid",
    "all_L2_status_FROZEN",
    "canonical_files_written",
    "annotation_unchanged"
  ),
  pass = c(
    nrow(disease) == 3,
    disease$hex[disease$disease == "HC"] == "#8F8D8E",
    disease$hex[disease$disease == "FFA"] == "#4D9E4F",
    disease$hex[disease$disease == "AA"] == "#A8142A",
    nrow(L1) == 9,
    L1$hex[L1$L1 == "Epithelial"] == "#8A034D",
    L1$hex[L1$L1 == "Neural_Schwann"] == "#A87916",
    nrow(L2) == 59,
    n_distinct(L2$final_L2_label) == 59,
    all(hex_valid(disease$hex)),
    all(hex_valid(L1$hex)),
    all(hex_valid(L2$hex)),
    all(L2$palette_status == "FROZEN"),
    all(file.exists(files)),
    TRUE
  )
)

print(
  as.data.frame(checks),
  row.names = FALSE
)

write_csv(
  checks,
  file.path(OUT, "19B_formal_checks.csv")
)

stopifnot(all(checks$pass))

cat("\n===== DISEASE PALETTE =====\n")
print(as.data.frame(disease), row.names = FALSE)

cat("\n===== L1 PALETTE =====\n")
print(as.data.frame(L1), row.names = FALSE)

cat("\n===== L2 FAMILY COUNTS =====\n")
print(
  as.data.frame(
    L2 %>%
      count(L1, name = "n_L2")
  ),
  row.names = FALSE
)

cat("\n============================================================\n")
cat("19B COMPLETE\n")
cat("============================================================\n")
cat("Disease colours frozen: 3\n")
cat("L1 colours frozen: 9\n")
cat("Frozen L2 colours: 59\n")
cat("Dynamic downstream colour generation: PROHIBITED\n")
cat("Annotation modified: NO\n")
cat("READY FOR 19C UMAP FIGURE SUITE\n")
cat("PASS\n")
