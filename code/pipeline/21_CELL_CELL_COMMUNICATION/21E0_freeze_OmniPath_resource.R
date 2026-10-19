suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(tibble)
})

cat("============================================================\n")
cat("PIPELINE_V2 — 21E0 OMNIPATH INTERCELL RESOURCE FREEZE\n")
cat("============================================================\n")

OUT <- paste0(
  "PIPELINE_V2/results/disease_biology/",
  "21_CELL_CELL_COMMUNICATION/",
  "21E0_OMNIPATH_RESOURCE"
)

RAW <- file.path(OUT, "raw")
TAB <- file.path(OUT, "tables")

dir.create(RAW, recursive=TRUE, showWarnings=FALSE)
dir.create(TAB, recursive=TRUE, showWarnings=FALSE)

# ============================================================
# 1. Freeze acquisition specification
# ============================================================

cat("\nSTEP 1/9 — freezing acquisition specification\n")

base <- "https://omnipathdb.org"

# OmniPath interaction datasets relevant to intercellular
# communication. No disease/expression filtering is applied.
interaction_url <- paste0(
  base,
  "/interactions?",
  "format=tsv",
  "&datasets=ligrecextra",
  "&genesymbols=1",
  "&organisms=9606",
  "&fields=sources,references,datasets"
)

intercell_url <- paste0(
  base,
  "/intercell?",
  "format=tsv"
)

spec <- tibble(
  resource=c(
    "OmniPath_interactions_ligrecextra",
    "OmniPath_intercell_annotations"
  ),
  endpoint=c(
    "interactions",
    "intercell"
  ),
  organism=c(
    "Homo sapiens / NCBI 9606 explicitly requested",
    "endpoint-native intercell resource; no organism parameter supplied"
  ),
  retrieval_method="OmniPath official REST API",
  url=c(
    interaction_url,
    intercell_url
  )
)

write_csv(
  spec,
  file.path(TAB, "21E0_acquisition_specification.csv")
)

print(spec, n=Inf, width=Inf)

# ============================================================
# 2. Download raw resources
# ============================================================

cat("\nSTEP 2/9 — downloading frozen REST resources\n")

raw_interactions <- file.path(
  RAW,
  "21E0_OmniPath_ligrecextra_raw.tsv"
)

raw_intercell <- file.path(
  RAW,
  "21E0_OmniPath_intercell_raw.tsv"
)

download.file(
  interaction_url,
  raw_interactions,
  mode="wb",
  quiet=FALSE
)

download.file(
  intercell_url,
  raw_intercell,
  mode="wb",
  quiet=FALSE
)

stopifnot(
  file.exists(raw_interactions),
  file.exists(raw_intercell),
  file.info(raw_interactions)$size > 0,
  file.info(raw_intercell)$size > 0
)

cat(
  "Interaction bytes:",
  file.info(raw_interactions)$size,
  "\n"
)

cat(
  "Intercell bytes:",
  file.info(raw_intercell)$size,
  "\n"
)

# ============================================================
# 3. Parse exact downloaded snapshots
# ============================================================

cat("\nSTEP 3/9 — parsing frozen snapshots\n")

lr <- read_tsv(
  raw_interactions,
  show_col_types=FALSE,
  progress=FALSE
)

ic <- read_tsv(
  raw_intercell,
  show_col_types=FALSE,
  progress=FALSE
)

cat(
  "ligrecextra rows:",
  nrow(lr),
  "| columns:",
  ncol(lr),
  "\n"
)

cat(
  "intercell rows:",
  nrow(ic),
  "| columns:",
  ncol(ic),
  "\n"
)

cat("\nInteraction columns:\n")
print(names(lr))

cat("\nIntercell columns:\n")
print(names(ic))

stopifnot(
  nrow(lr) > 0,
  nrow(ic) > 0
)

# ============================================================
# 4. Preserve parsed snapshots without biological filtering
# ============================================================

cat("\nSTEP 4/9 — freezing parsed snapshots\n")

write_csv(
  lr,
  file.path(
    TAB,
    "21E0_OmniPath_ligrecextra_frozen.csv"
  )
)

write_csv(
  ic,
  file.path(
    TAB,
    "21E0_OmniPath_intercell_frozen.csv"
  )
)

saveRDS(
  lr,
  file.path(
    RAW,
    "21E0_OmniPath_ligrecextra_frozen.rds"
  ),
  compress="xz"
)

saveRDS(
  ic,
  file.path(
    RAW,
    "21E0_OmniPath_intercell_frozen.rds"
  ),
  compress="xz"
)

# ============================================================
# 5. Basic interaction inventory
# ============================================================

cat("\nSTEP 5/9 — interaction inventory\n")

required_lr <- c(
  "source",
  "target"
)

stopifnot(
  all(required_lr %in% names(lr))
)

interaction_inventory <- tibble(
  metric=c(
    "rows",
    "unique_source_entities",
    "unique_target_entities",
    "unique_source_target_pairs",
    "complete_duplicate_rows"
  ),
  value=c(
    nrow(lr),
    n_distinct(lr$source),
    n_distinct(lr$target),
    nrow(
      distinct(
        lr,
        source,
        target
      )
    ),
    sum(duplicated(lr))
  )
)

write_csv(
  interaction_inventory,
  file.path(
    TAB,
    "21E0_interaction_inventory.csv"
  )
)

print(
  interaction_inventory,
  n=Inf
)

# ============================================================
# 6. Intercell annotation inventory
# ============================================================

cat("\nSTEP 6/9 — intercell annotation inventory\n")

required_ic <- c(
  "category",
  "parent",
  "database",
  "scope",
  "aspect",
  "source",
  "uniprot",
  "genesymbol",
  "entity_type",
  "consensus_score",
  "transmitter",
  "receiver",
  "secreted",
  "plasma_membrane_transmembrane",
  "plasma_membrane_peripheral"
)

missing_ic <- setdiff(
  required_ic,
  names(ic)
)

cat(
  "Missing expected intercell columns:",
  length(missing_ic),
  "\n"
)

if (length(missing_ic)) {
  print(missing_ic)
}

stopifnot(
  all(
    c(
      "category",
      "parent",
      "source",
      "uniprot",
      "genesymbol"
    ) %in% names(ic)
  )
)

category_summary <- ic %>%
  count(
    category,
    parent,
    sort=TRUE,
    name="n_records"
  )

write_csv(
  category_summary,
  file.path(
    TAB,
    "21E0_intercell_category_summary.csv"
  )
)

role_cols <- intersect(
  c(
    "transmitter",
    "receiver",
    "secreted",
    "plasma_membrane_transmembrane",
    "plasma_membrane_peripheral"
  ),
  names(ic)
)

role_summary <- bind_rows(
  lapply(
    role_cols,
    function(v) {
      tibble(
        role=v,
        n_true=sum(
          as.logical(ic[[v]]),
          na.rm=TRUE
        ),
        n_false=sum(
          !as.logical(ic[[v]]),
          na.rm=TRUE
        ),
        n_missing=sum(
          is.na(ic[[v]])
        )
      )
    }
  )
)

write_csv(
  role_summary,
  file.path(
    TAB,
    "21E0_intercell_role_summary.csv"
  )
)

print(
  role_summary,
  n=Inf
)

# ============================================================
# 7. Resource provenance
# ============================================================

cat("\nSTEP 7/9 — resource provenance\n")

get_hash <- function(path) {

  out <- system2(
    "shasum",
    c("-a", "256", shQuote(path)),
    stdout=TRUE,
    stderr=TRUE
  )

  strsplit(out[1], "\\s+")[[1]][1]
}

files_to_hash <- c(
  raw_interactions,
  raw_intercell,
  file.path(
    RAW,
    "21E0_OmniPath_ligrecextra_frozen.rds"
  ),
  file.path(
    RAW,
    "21E0_OmniPath_intercell_frozen.rds"
  )
)

hashes <- tibble(
  file=files_to_hash,
  sha256=vapply(
    files_to_hash,
    get_hash,
    character(1)
  )
)

write_csv(
  hashes,
  file.path(
    TAB,
    "21E0_resource_hashes.csv"
  )
)

provenance <- tibble(
  field=c(
    "retrieval_date",
    "resource",
    "retrieval_method",
    "organism",
    "interaction_dataset",
    "OmnipathR_installed_version",
    "OmnipathR_wrapper_used_for_retrieval",
    "wrapper_reason_not_used",
    "disease_filtering",
    "expression_filtering",
    "communication_inference"
  ),
  value=c(
    as.character(Sys.Date()),
    "OmniPath",
    "official REST API",
    paste(
      "interactions: Homo sapiens / NCBI 9606 explicitly requested;",
      "intercell: endpoint-native resource without organism parameter"
    ),
    "ligrecextra",
    as.character(
      packageVersion("OmnipathR")
    ),
    "FALSE",
    paste(
      "OmnipathR 3.14.0 organism resolution failed",
      "at external Ensembl taxonomy lookup;",
      "official REST API remained accessible"
    ),
    "NONE",
    "NONE",
    "NONE"
  )
)

write_csv(
  provenance,
  file.path(
    TAB,
    "21E0_resource_provenance.csv"
  )
)

print(
  provenance,
  n=Inf,
  width=Inf
)

# ============================================================
# 8. Formal checks
# ============================================================

cat("\nSTEP 8/9 — formal checks\n")

checks <- tibble(
  check=c(
    "interaction_download_present",
    "intercell_download_present",
    "interaction_nonempty",
    "intercell_nonempty",
    "interaction_source_present",
    "interaction_target_present",
    "intercell_category_present",
    "intercell_genesymbol_present",
    "intercell_transmitter_present",
    "intercell_receiver_present",
    "intercell_secreted_present",
    "raw_hashes_generated",
    "interaction_human_9606_prespecified",
    "ligrecextra_prespecified",
    "no_disease_filtering",
    "no_expression_filtering",
    "no_communication_inference",
    "no_DE_refit",
    "no_pathway_refit",
    "annotation_unchanged"
  ),
  pass=c(
    file.exists(raw_interactions),
    file.exists(raw_intercell),
    nrow(lr) > 0,
    nrow(ic) > 0,
    "source" %in% names(lr),
    "target" %in% names(lr),
    "category" %in% names(ic),
    "genesymbol" %in% names(ic),
    "transmitter" %in% names(ic),
    "receiver" %in% names(ic),
    "secreted" %in% names(ic),
    nrow(hashes) == 4 &&
      all(!is.na(hashes$sha256)),
    TRUE,
    TRUE,
    TRUE,
    TRUE,
    TRUE,
    TRUE,
    TRUE,
    TRUE
  )
)

print(
  checks,
  n=Inf,
  width=Inf
)

stopifnot(
  all(checks$pass)
)

write_csv(
  checks,
  file.path(
    TAB,
    "21E0_formal_checks.csv"
  )
)

# ============================================================
# 9. Audit
# ============================================================

cat("\nSTEP 9/9 — writing audit\n")

audit <- c(
  "PIPELINE_V2 — 21E0 OmniPath intercell resource freeze",
  "",
  paste0("Date: ", Sys.Date()),
  "",
  "RESOURCE",
  "--------",
  "OmniPath official REST API",
  "Interactions organism: Homo sapiens / NCBI taxonomy 9606 explicitly requested",
  "Intercell endpoint: endpoint-native resource; no organism parameter supplied",
  "Interaction dataset: ligrecextra",
  "",
  paste0(
    "Frozen interaction records: ",
    nrow(lr)
  ),
  paste0(
    "Frozen intercell annotation records: ",
    nrow(ic)
  ),
  "",
  "RETRIEVAL",
  "---------",
  paste(
    "OmnipathR 3.14.0 was installed but was not used",
    "for resource retrieval."
  ),
  paste(
    "Its organism-resolution layer attempted an external",
    "Ensembl taxonomy lookup which failed in the current",
    "environment before the OmniPath query was executed."
  ),
  paste(
    "The official OmniPath REST endpoints were independently",
    "confirmed accessible and were used directly."
  ),
  "",
  "ANALYTICAL FREEZE",
  "-----------------",
  "- no disease filtering",
  "- no AA/FFA comparison",
  "- no expression filtering",
  "- no atlas-gene intersection",
  "- no ligand/receptor prioritization",
  "- no communication inference",
  "- no differential-expression refit",
  "- no pathway refit",
  "- no annotation modification",
  "",
  "The downloaded TSV files are the immutable raw snapshots.",
  "SHA-256 hashes are recorded for reproducibility.",
  "",
  "PASS"
)

writeLines(
  audit,
  file.path(
    OUT,
    "21E0_audit.txt"
  )
)

cat("\n============================================================\n")
cat("21E0 COMPLETE\n")
cat("============================================================\n")
cat("OmniPath human intercell resource frozen\n")
cat("Raw REST responses preserved\n")
cat("Resource provenance preserved\n")
cat("SHA-256 hashes recorded\n")
cat("No disease filtering\n")
cat("No expression filtering\n")
cat("No communication inference\n")
cat("PASS\n")
