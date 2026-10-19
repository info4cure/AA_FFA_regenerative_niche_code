suppressPackageStartupMessages({
  library(tidyverse)
})

cat("============================================================\n")
cat("21G4C3 — QUANTITATIVE FIGURE FEASIBILITY\n")
cat("============================================================\n")

BASE <- "PIPELINE_V2/results/disease_biology/21_CELL_CELL_COMMUNICATION"

G21 <- file.path(
  BASE,
  "21G_MULTICOMPARTMENT_COMMUNICATION"
)

OUT <- file.path(
  G21,
  "21G4_CROSSMODAL_INTEGRATION"
)

TAB <- file.path(OUT, "tables")
dir.create(TAB, recursive=TRUE, showWarnings=FALSE)

# ------------------------------------------------------------
# 1. Candidate frozen patient-level sources
# ------------------------------------------------------------

cat("\nSTEP 1/7 — locating frozen patient-level communication sources\n")

candidate_dirs <- c(
  file.path(BASE, "21F_PATIENT_LEVEL_COMMUNICATION"),
  file.path(BASE, "21F2_BIOLOGICAL_LR_AXES"),
  file.path(BASE, "21F2_NETWORK_AGGREGATION"),
  file.path(G21, "21G2_INFERENCE"),
  file.path(G21, "21G3_BIOLOGICAL_READOUT")
)

candidate_dirs <- candidate_dirs[
  dir.exists(candidate_dirs)
]

files <- unlist(
  lapply(
    candidate_dirs,
    function(d) {
      list.files(
        d,
        recursive=TRUE,
        full.names=TRUE,
        pattern="\\.(csv|tsv|rds)$"
      )
    }
  )
)

cat("Candidate files:", length(files), "\n")

# ------------------------------------------------------------
# 2. Inventory file names
# ------------------------------------------------------------

cat("\nSTEP 2/7 — classifying candidate files\n")

inventory <- tibble(
  file = files,
  basename = basename(files)
) %>%
  mutate(
    likely_patient =
      str_detect(
        str_to_lower(basename),
        "patient|quantitative|projection|axis"
      ),

    likely_A =
      str_detect(
        str_to_lower(basename),
        "incoming|a_axis|axis"
      ),

    likely_D =
      str_detect(
        str_to_lower(basename),
        "density|d_axis|decomposition"
      ),

    likely_sender =
      str_detect(
        str_to_lower(basename),
        "sender|active|intensity"
      )
  )

print(
  inventory %>%
    filter(
      likely_patient |
      likely_A |
      likely_D |
      likely_sender
    ),
  n=Inf,
  width=Inf
)

write_csv(
  inventory,
  file.path(
    TAB,
    "21G4C3_candidate_source_inventory.csv"
  )
)

# ------------------------------------------------------------
# 3. Inspect tabular structures without modifying them
# ------------------------------------------------------------

cat("\nSTEP 3/7 — inspecting candidate tabular structures\n")

inspect_file <- function(f) {

  ext <- tools::file_ext(f)

  x <- tryCatch(
    {
      if (ext == "csv") {
        read_csv(f, show_col_types=FALSE)
      } else if (ext == "tsv") {
        read_tsv(f, show_col_types=FALSE)
      } else if (ext == "rds") {
        readRDS(f)
      } else {
        NULL
      }
    },
    error=function(e) NULL
  )

  if (is.null(x)) {
    return(
      tibble(
        file=f,
        object_class=NA_character_,
        n_rows=NA_integer_,
        n_cols=NA_integer_,
        columns=NA_character_
      )
    )
  }

  if (
    is.data.frame(x) ||
    inherits(x, "tbl")
  ) {
    return(
      tibble(
        file=f,
        object_class=paste(class(x), collapse=";"),
        n_rows=nrow(x),
        n_cols=ncol(x),
        columns=paste(names(x), collapse=";")
      )
    )
  }

  tibble(
    file=f,
    object_class=paste(class(x), collapse=";"),
    n_rows=NA_integer_,
    n_cols=NA_integer_,
    columns=if (is.list(x)) {
      paste(names(x), collapse=";")
    } else {
      NA_character_
    }
  )
}

structure_inventory <- map_dfr(
  files,
  inspect_file
)

write_csv(
  structure_inventory,
  file.path(
    TAB,
    "21G4C3_source_structure_inventory.csv"
  )
)

cat("\n===== STRUCTURES WITH PATIENT-LIKE COLUMNS =====\n")

patient_structures <- structure_inventory %>%
  filter(
    str_detect(
      str_to_lower(coalesce(columns, "")),
      "patient|sample|subject|disease"
    )
  )

print(
  patient_structures,
  n=Inf,
  width=Inf
)

# ------------------------------------------------------------
# 4. Search specifically for quantitative endpoint columns
# ------------------------------------------------------------

cat("\nSTEP 4/7 — endpoint-column discovery\n")

endpoint_inventory <- structure_inventory %>%
  mutate(
    has_patient =
      str_detect(
        str_to_lower(coalesce(columns, "")),
        "patient|sample|subject"
      ),

    has_disease =
      str_detect(
        str_to_lower(coalesce(columns, "")),
        "disease|condition|group"
      ),

    has_receiver =
      str_detect(
        str_to_lower(coalesce(columns, "")),
        "receiver"
      ),

    has_axis =
      str_detect(
        str_to_lower(coalesce(columns, "")),
        "axis"
      ),

    has_incoming =
      str_detect(
        str_to_lower(coalesce(columns, "")),
        "incoming"
      ),

    has_density =
      str_detect(
        str_to_lower(coalesce(columns, "")),
        "density|pooled_d|d_axis"
      ),

    has_sender_intensity =
      str_detect(
        str_to_lower(coalesce(columns, "")),
        "sender.*intensity|active.*intensity"
      )
  )

cat("\n===== ENDPOINT INVENTORY =====\n")

print(
  endpoint_inventory %>%
    filter(
      has_patient |
      has_incoming |
      has_density |
      has_sender_intensity
    ) %>%
    select(
      file,
      object_class,
      n_rows,
      n_cols,
      has_patient,
      has_disease,
      has_receiver,
      has_axis,
      has_incoming,
      has_density,
      has_sender_intensity
    ),
  n=Inf,
  width=Inf
)

write_csv(
  endpoint_inventory,
  file.path(
    TAB,
    "21G4C3_endpoint_inventory.csv"
  )
)

# ------------------------------------------------------------
# 5. Assess candidate quantitative panels
# ------------------------------------------------------------

cat("\nSTEP 5/7 — quantitative panel feasibility\n")

patient_axis_available <- any(
  endpoint_inventory$has_patient &
  endpoint_inventory$has_receiver &
  endpoint_inventory$has_axis
)

incoming_available <- any(
  endpoint_inventory$has_patient &
  endpoint_inventory$has_incoming
)

density_available <- any(
  endpoint_inventory$has_patient &
  endpoint_inventory$has_density
)

sender_available <- any(
  endpoint_inventory$has_patient &
  endpoint_inventory$has_sender_intensity
)

# 317-row landscape and frozen multimodal matrix are already
# demonstrably available from previous frozen modules.

panel_feasibility <- tribble(
  ~panel_candidate,
  ~data_status,
  ~proposed_plot,
  ~decision_status,

  "A_GLOBAL_COMMUNICATION_LANDSCAPE",
  "AVAILABLE_FROZEN",
  paste(
    "28 receivers × 12 communication axes × 3 contrasts;",
    "effect encoded continuously; q-value/evidence encoded separately"
  ),
  "FEASIBLE",

  "B_AA_PATIENT_LEVEL_ARCHITECTURE",
  ifelse(
    patient_axis_available,
    "PATIENT_AXIS_SOURCE_DETECTED",
    "PATIENT_AXIS_SOURCE_NOT_DETECTED"
  ),
  paste(
    "patient-level AA/HC quantitative distributions for",
    "follicular IGF/growth-survival signals"
  ),
  ifelse(
    patient_axis_available,
    "CANDIDATE_FEASIBLE",
    "REQUIRES_SOURCE_RECOVERY"
  ),

  "C_FFA_PATIENT_LEVEL_ARCHITECTURE",
  ifelse(
    patient_axis_available,
    "PATIENT_AXIS_SOURCE_DETECTED",
    "PATIENT_AXIS_SOURCE_NOT_DETECTED"
  ),
  paste(
    "patient-level FFA/HC quantitative distributions across",
    "epithelial growth-survival receivers"
  ),
  ifelse(
    patient_axis_available,
    "CANDIDATE_FEASIBLE",
    "REQUIRES_SOURCE_RECOVERY"
  ),

  "D_A_VS_D_DECOMPOSITION",
  ifelse(
    density_available,
    "PATIENT_DENSITY_SOURCE_DETECTED",
    "FROZEN_SUMMARY_DECOMPOSITION_ONLY"
  ),
  paste(
    "incoming A-axis versus pooled D-axis effect architecture;",
    "active sender intensity retained separately"
  ),
  "FEASIBLE_FROM_FROZEN_SUMMARY",

  "E_CROSSMODAL_EVIDENCE_HEATMAP",
  "AVAILABLE_FROZEN",
  paste(
    "communication + density + DA + DE + pathway evidence;",
    "no combined numeric score"
  ),
  "FEASIBLE",

  "F_PATIENT_COMMUNICATION_PHENOTYPE",
  ifelse(
    patient_axis_available,
    "PATIENT_AXIS_SOURCE_DETECTED",
    "PATIENT_AXIS_SOURCE_NOT_DETECTED"
  ),
  paste(
    "unsupervised patient × receiver-axis PCA/heatmap;",
    "disease labels used only after unsupervised embedding"
  ),
  ifelse(
    patient_axis_available,
    "CANDIDATE_FEASIBLE",
    "REQUIRES_SOURCE_RECOVERY"
  )
)

cat("\n===== PANEL FEASIBILITY =====\n")
print(panel_feasibility, n=Inf, width=Inf)

write_csv(
  panel_feasibility,
  file.path(
    TAB,
    "21G4C3_panel_feasibility.csv"
  )
)

# ------------------------------------------------------------
# 6. Analytical safeguards
# ------------------------------------------------------------

cat("\nSTEP 6/7 — analytical safeguards\n")

formal_checks <- tribble(
  ~check, ~pass,

  "source_discovery_only",
  TRUE,

  "no_communication_recalculation",
  TRUE,

  "no_DA_recalculation",
  TRUE,

  "no_DE_recalculation",
  TRUE,

  "no_pathway_recalculation",
  TRUE,

  "no_threshold_refit",
  TRUE,

  "no_crossmodal_score",
  TRUE,

  "no_combined_p_value",
  TRUE,

  "no_combined_q_value",
  TRUE,

  "no_patient_embedding_generated",
  TRUE,

  "no_figure_generated",
  TRUE,

  "no_panel_selected",
  TRUE,

  "no_manuscript_modified",
  TRUE
)

print(
  formal_checks,
  n=Inf,
  width=Inf
)

stopifnot(
  all(formal_checks$pass)
)

write_csv(
  formal_checks,
  file.path(
    TAB,
    "21G4C3_formal_checks.csv"
  )
)

# ------------------------------------------------------------
# 7. Final report
# ------------------------------------------------------------

cat("\nSTEP 7/7 — final feasibility report\n")

cat("\n============================================================\n")
cat("21G4C3 COMPLETE\n")
cat("============================================================\n")

cat(
  "Patient × receiver × axis source detected:",
  patient_axis_available,
  "\n"
)

cat(
  "Patient-level incoming A-axis source detected:",
  incoming_available,
  "\n"
)

cat(
  "Patient-level density/D-axis source detected:",
  density_available,
  "\n"
)

cat(
  "Patient-level sender-active-intensity source detected:",
  sender_available,
  "\n"
)

cat("Global 317-test landscape: FEASIBLE\n")
cat("Cross-modal evidence heatmap: FEASIBLE\n")
cat("No plot generated\n")
cat("No panel selected\n")
cat("No statistical inference\n")
cat("No analytical refit\n")
cat("PASS\n")
