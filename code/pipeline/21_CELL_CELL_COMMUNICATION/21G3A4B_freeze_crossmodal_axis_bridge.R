suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
})

cat("============================================================\n")
cat("21G3A4B — FREEZE COMMUNICATION ↔ PATHWAY AXIS BRIDGE\n")
cat("============================================================\n")

BASE <- paste0(
  "PIPELINE_V2/results/disease_biology/",
  "21_CELL_CELL_COMMUNICATION/",
  "21G_MULTICOMPARTMENT_COMMUNICATION/"
)

OUT <- file.path(
  BASE,
  "21G4_CROSSMODAL_INTEGRATION"
)

TAB <- file.path(OUT, "tables")

dir.create(
  TAB,
  recursive=TRUE,
  showWarnings=FALSE
)

COMM_FILE <- file.path(
  BASE,
  "21G3_BIOLOGICAL_READOUT",
  "tables",
  "21G2C_primary_A_all.csv"
)

PATH_FILE <- paste0(
  "PIPELINE_V2/results/disease_biology/",
  "20_PATHWAY_ANALYSIS/",
  "20F_CROSS_DATABASE_CONSENSUS/",
  "20F2B_integrated_prioritization/tables/",
  "20F2B_integrated_state_prioritization.csv"
)

stopifnot(
  file.exists(COMM_FILE),
  file.exists(PATH_FILE)
)

comm <- read_csv(
  COMM_FILE,
  show_col_types=FALSE
)

pathway <- read_csv(
  PATH_FILE,
  show_col_types=FALSE
)

cat("\nSTEP 1/6 — recovering frozen axis universes\n")

comm_axes <- sort(unique(comm$axis))
path_axes <- sort(unique(pathway$axis))

cat("Communication axes:", length(comm_axes), "\n")
cat("20F2B pathway axes:", length(path_axes), "\n")

stopifnot(
  length(comm_axes) == 12,
  length(path_axes) == 8
)

cat("\nSTEP 2/6 — defining prespecified biological bridge\n")

bridge <- tribble(
  ~communication_axis,
  ~pathway_axis,
  ~bridge_class,
  ~bridge_interpretation,

  "INFLAMMATORY_IFN_CHEMOKINE",
  "IFN_JAK_STAT_ANTIGEN",
  "DIRECT",
  "Interferon/chemokine inflammatory communication linked to IFN-JAK-STAT and antigen-response transcriptional programmes",

  "CYTOTOXIC_TCELL",
  "T_CELL_CYTOTOXIC",
  "DIRECT",
  "Cytotoxic communication linked to T-cell cytotoxic transcriptional programmes",

  "ECM_INTEGRIN_REMODELING",
  "ECM_COLLAGEN_FIBROSIS",
  "RELATED",
  "ECM-integrin adhesive communication related to extracellular-matrix, collagen and fibrotic transcriptional remodelling",

  "TGF_BETA_FIBROTIC",
  "ECM_COLLAGEN_FIBROSIS",
  "DIRECT",
  "TGF-beta-associated fibrotic communication linked to extracellular-matrix, collagen and fibrotic transcriptional remodelling",

  "VASCULAR_ANGIOGENIC",
  "VASCULAR_ANGIOGENESIS",
  "DIRECT",
  "Angiogenic communication linked to vascular and angiogenic transcriptional programmes",

  "WNT_FOLLICULAR_REGENERATIVE",
  "WNT_HEDGEHOG_FOLLICLE",
  "DIRECT",
  "WNT-associated regenerative communication linked to WNT-Hedgehog and hair-follicle transcriptional programmes",

  "BMP_HOMEOSTATIC_DIFFERENTIATION",
  "EPITHELIAL_DIFFERENTIATION",
  "RELATED",
  "BMP-associated homeostatic differentiation communication related to epithelial differentiation programmes",

  "GROWTH_SURVIVAL",
  NA_character_,
  "NO_PRESPECIFIED_BRIDGE",
  "Broad growth-survival communication has no sufficiently specific frozen 20F2B pathway counterpart",

  "IGF_TROPHIC",
  NA_character_,
  "NO_PRESPECIFIED_BRIDGE",
  "IGF trophic communication has no sufficiently specific frozen 20F2B pathway counterpart",

  "FGF_TROPHIC",
  NA_character_,
  "NO_PRESPECIFIED_BRIDGE",
  "FGF trophic communication has no sufficiently specific frozen 20F2B pathway counterpart",

  "IMMUNE_REGULATORY",
  NA_character_,
  "NO_PRESPECIFIED_BRIDGE",
  "Immune-regulatory communication is not equivalent to the frozen IFN or cytotoxic pathway axes",

  "TNF_IL1_INFLAMMATORY",
  NA_character_,
  "NO_PRESPECIFIED_BRIDGE",
  "TNF/IL1 inflammatory communication is not equivalent to the frozen IFN-JAK-STAT pathway axis"
)

cat("\n===== FROZEN AXIS BRIDGE =====\n")
print(bridge, n=Inf, width=Inf)

cat("\nSTEP 3/6 — checking complete communication-axis coverage\n")

coverage <- tibble(
  communication_axis=comm_axes
) %>%
  left_join(
    bridge,
    by="communication_axis"
  )

print(
  coverage,
  n=Inf,
  width=Inf
)

stopifnot(
  nrow(coverage) == 12,
  all(!is.na(coverage$bridge_class))
)

cat("\nSTEP 4/6 — checking mapped pathway axes\n")

mapped <- bridge %>%
  filter(
    bridge_class !=
      "NO_PRESPECIFIED_BRIDGE"
  )

stopifnot(
  all(
    mapped$pathway_axis %in%
      path_axes
  )
)

cat("\nMapped communication axes:",
    nrow(mapped), "\n")

cat("Unmapped communication axes:",
    sum(
      bridge$bridge_class ==
        "NO_PRESPECIFIED_BRIDGE"
    ),
    "\n"
)

cat("\n===== BRIDGE CLASS COUNTS =====\n")

print(
  bridge %>%
    count(
      bridge_class,
      name="n_axes"
    ),
  n=Inf
)

cat("\nSTEP 5/6 — freezing bridge\n")

write_csv(
  bridge,
  file.path(
    TAB,
    "21G3A4B_frozen_communication_pathway_axis_bridge.csv"
  )
)

write_csv(
  coverage,
  file.path(
    TAB,
    "21G3A4B_axis_bridge_coverage.csv"
  )
)

cat("\nSTEP 6/6 — formal checks\n")

formal_checks <- tribble(
  ~check, ~pass,

  "twelve_communication_axes_retained",
  length(comm_axes) == 12,

  "eight_pathway_axes_retained",
  length(path_axes) == 8,

  "communication_axis_unique_in_bridge",
  !anyDuplicated(bridge$communication_axis),

  "all_communication_axes_classified",
  setequal(
    comm_axes,
    bridge$communication_axis
  ),

  "all_mapped_pathway_axes_exist",
  all(
    mapped$pathway_axis %in%
      path_axes
  ),

  "direct_and_related_distinguished",
  all(
    c("DIRECT","RELATED") %in%
      bridge$bridge_class
  ),

  "unmapped_axes_explicit",
  any(
    bridge$bridge_class ==
      "NO_PRESPECIFIED_BRIDGE"
  ),

  "no_disease_direction_used",
  TRUE,

  "no_effect_used_for_mapping",
  TRUE,

  "no_p_value_used_for_mapping",
  TRUE,

  "no_q_value_used_for_mapping",
  TRUE,

  "no_crossmodal_score",
  TRUE,

  "no_statistical_inference",
  TRUE,

  "no_analytical_refit",
  TRUE
)

write_csv(
  formal_checks,
  file.path(
    TAB,
    "21G3A4B_formal_checks.csv"
  )
)

cat("\n===== FORMAL CHECKS =====\n")
print(formal_checks, n=Inf, width=Inf)

stopifnot(
  all(formal_checks$pass)
)

audit <- c(
  "21G3A4B — FROZEN COMMUNICATION ↔ PATHWAY AXIS BRIDGE",
  "",
  paste0(
    "Communication axes = ",
    length(comm_axes)
  ),
  paste0(
    "20F2B pathway axes = ",
    length(path_axes)
  ),
  paste0(
    "Mapped communication axes = ",
    nrow(mapped)
  ),
  paste0(
    "No-prespecified-bridge axes = ",
    sum(
      bridge$bridge_class ==
        "NO_PRESPECIFIED_BRIDGE"
    )
  ),
  "",
  "BRIDGE PRINCIPLES",
  "-----------------",
  "Bridge defined before cross-modal concordance analysis.",
  "Disease direction was not used to define mappings.",
  "Communication effect sizes were not used to define mappings.",
  "P-values and q-values were not used to define mappings.",
  "DIRECT and RELATED mappings remain analytically distinct.",
  "Absence of a pathway bridge is not interpreted as biological absence.",
  "No cross-modal score is defined.",
  "",
  "ANALYTICAL INTEGRITY",
  "--------------------",
  "No new statistical inference.",
  "No communication refit.",
  "No pathway refit.",
  "No DE refit.",
  "No DA refit.",
  "No annotation change.",
  "",
  "PASS"
)

writeLines(
  audit,
  file.path(
    OUT,
    "21G3A4B_audit.txt"
  )
)

cat("\n============================================================\n")
cat("21G3A4B COMPLETE\n")
cat("============================================================\n")
cat("Communication-pathway axis bridge frozen before concordance analysis\n")
cat("12 communication axes retained\n")
cat(nrow(mapped), "axes mapped to frozen 20F2B biology\n")
cat(sum(
  bridge$bridge_class ==
    "NO_PRESPECIFIED_BRIDGE"
), "axes explicitly retained without pathway bridge\n")
cat("DIRECT and RELATED mappings distinguished\n")
cat("No disease direction used for mapping\n")
cat("No cross-modal score\n")
cat("No statistical inference\n")
cat("No analytical refit\n")
cat("PASS\n")
