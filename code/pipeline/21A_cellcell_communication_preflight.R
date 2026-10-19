suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(readr)
  library(tibble)
  library(stringr)
  library(purrr)
})

cat("============================================================\n")
cat("PIPELINE_V2 — 21A CELL–CELL COMMUNICATION PREFLIGHT\n")
cat("============================================================\n")

OUT <- paste0(
  "PIPELINE_V2/results/disease_biology/",
  "21_CELL_CELL_COMMUNICATION/21A_PREFLIGHT"
)
dir.create(OUT, recursive=TRUE, showWarnings=FALSE)

# ============================================================
# 1. Locate candidate Seurat objects
# ============================================================

cat("\nSTEP 1/8 — locating candidate Seurat objects\n")

roots <- c(
  "PIPELINE_V2",
  "08_REPRODUCIBILITY"
)

rds_files <- unlist(
  lapply(
    roots[file.exists(roots)],
    function(root) {
      list.files(
        root,
        pattern="\\.rds$",
        recursive=TRUE,
        full.names=TRUE
      )
    }
  )
)

rds_files <- sort(unique(rds_files))

cat("RDS files found:", length(rds_files), "\n")

candidate_files <- rds_files[
  str_detect(
    basename(rds_files),
    regex(
      "final|frozen|annot|atlas|global|integrat|seurat|L2",
      ignore_case=TRUE
    )
  )
]

cat("Filename-based candidates:", length(candidate_files), "\n")

write_lines(
  rds_files,
  file.path(OUT, "21A_all_RDS_files.txt")
)

write_lines(
  candidate_files,
  file.path(OUT, "21A_candidate_RDS_files.txt")
)

# ============================================================
# 2. Inspect candidate objects safely
# ============================================================

cat("\nSTEP 2/8 — inspecting candidate objects\n")

inspect_rds <- function(f) {

  obj <- tryCatch(
    readRDS(f),
    error=function(e) NULL
  )

  if (is.null(obj)) {
    return(
      tibble(
        file=f,
        readable=FALSE,
        is_seurat=FALSE,
        n_cells=NA_real_,
        n_features=NA_real_
      )
    )
  }

  is_s <- inherits(obj, "Seurat")

  if (!is_s) {
    return(
      tibble(
        file=f,
        readable=TRUE,
        is_seurat=FALSE,
        n_cells=NA_real_,
        n_features=NA_real_
      )
    )
  }

  tibble(
    file=f,
    readable=TRUE,
    is_seurat=TRUE,
    n_cells=ncol(obj),
    n_features=nrow(obj)
  )
}

inspection <- map_dfr(
  candidate_files,
  inspect_rds
)

write_csv(
  inspection,
  file.path(OUT, "21A_candidate_object_inventory.csv")
)

print(
  inspection %>%
    filter(is_seurat) %>%
    arrange(desc(n_cells)),
  n=Inf
)

# ============================================================
# 3. Select largest readable Seurat object for metadata audit
# ============================================================

cat("\nSTEP 3/8 — selecting metadata-audit object\n")

seurat_inventory <- inspection %>%
  filter(is_seurat)

if (nrow(seurat_inventory) == 0) {
  stop(
    paste(
      "No readable Seurat object identified automatically.",
      "Inspect 21A_all_RDS_files.txt before proceeding."
    )
  )
}

audit_file <- seurat_inventory %>%
  arrange(desc(n_cells)) %>%
  slice(1) %>%
  pull(file)

cat("Audit object:\n", audit_file, "\n")

obj <- readRDS(audit_file)

cat("Cells   :", ncol(obj), "\n")
cat("Features:", nrow(obj), "\n")

# ============================================================
# 4. Metadata audit
# ============================================================

cat("\nSTEP 4/8 — metadata audit\n")

md <- obj[[]]

cat("Metadata columns:\n")
print(names(md))

expected_exact <- c(
  "patient_id",
  "disease",
  "batch"
)

expected_pattern <- c(
  "L1",
  "L2"
)

exact_present <- expected_exact %in% names(md)

pattern_hits <- map(
  expected_pattern,
  ~ names(md)[
      str_detect(
        names(md),
        regex(.x, ignore_case=TRUE)
      )
    ]
)

metadata_check <- tibble(
  field=expected_exact,
  present=exact_present
)

write_csv(
  metadata_check,
  file.path(OUT, "21A_metadata_required_fields.csv")
)

cat("\nRequired exact metadata:\n")
print(metadata_check)

cat("\nL1-like columns:\n")
print(pattern_hits[[1]])

cat("\nL2-like columns:\n")
print(pattern_hits[[2]])

write_lines(
  pattern_hits[[1]],
  file.path(OUT, "21A_L1_like_columns.txt")
)

write_lines(
  pattern_hits[[2]],
  file.path(OUT, "21A_L2_like_columns.txt")
)

# ============================================================
# 5. Patient / disease / batch structure
# ============================================================

cat("\nSTEP 5/8 — patient/disease/batch structure\n")

if (all(expected_exact %in% names(md))) {

  design <- md %>%
    rownames_to_column("cell") %>%
    distinct(
      patient_id,
      disease,
      batch
    ) %>%
    arrange(
      disease,
      patient_id
    )

  write_csv(
    design,
    file.path(OUT, "21A_patient_design.csv")
  )

  cat("\nPatients by disease:\n")
  print(
    design %>%
      count(disease, name="n_patients")
  )

  cat("\nPatients by disease and batch:\n")
  print(
    design %>%
      count(
        disease,
        batch,
        name="n_patients"
      )
  )

} else {

  design <- tibble()

  cat(
    "Exact patient_id/disease/batch fields not all present;",
    "design summary deferred.\n"
  )
}

# ============================================================
# 6. Assay / layer availability
# ============================================================

cat("\nSTEP 6/8 — assay and layer availability\n")

cat("Assays:\n")
print(Assays(obj))

cat("Default assay:", DefaultAssay(obj), "\n")

assay_info <- map_dfr(
  Assays(obj),
  function(a) {

    aa <- obj[[a]]

    layer_names <- tryCatch(
      Layers(aa),
      error=function(e) character()
    )

    tibble(
      assay=a,
      n_features=nrow(aa),
      layers=paste(layer_names, collapse=";")
    )
  }
)

print(assay_info)

write_csv(
  assay_info,
  file.path(OUT, "21A_assay_layer_inventory.csv")
)

# ============================================================
# 7. Communication-tool environment
# ============================================================

cat("\nSTEP 7/8 — communication-tool environment\n")

packages <- c(
  "CellChat",
  "CellPhoneDB",
  "LIANA",
  "liana",
  "NicheNetR",
  "nichenetr",
  "Seurat",
  "SeuratObject",
  "future",
  "igraph"
)

pkg_info <- tibble(
  package=packages,
  installed=map_lgl(
    packages,
    ~ requireNamespace(.x, quietly=TRUE)
  ),
  version=map_chr(
    packages,
    function(p) {
      if (!requireNamespace(p, quietly=TRUE))
        return(NA_character_)
      as.character(packageVersion(p))
    }
  )
)

print(pkg_info)

write_csv(
  pkg_info,
  file.path(OUT, "21A_package_inventory.csv")
)

# ============================================================
# 8. Formal checks
# ============================================================

cat("\nSTEP 8/8 — formal checks\n")

checks <- tibble(
  check=c(
    "RDS_inventory_generated",
    "Seurat_object_identified",
    "metadata_inventory_generated",
    "assay_inventory_generated",
    "package_inventory_generated",
    "no_communication_inference",
    "no_DE_refit",
    "no_pathway_refit",
    "annotation_unchanged",
    "source_object_not_overwritten"
  ),
  pass=c(
    file.exists(
      file.path(OUT, "21A_all_RDS_files.txt")
    ),
    nrow(seurat_inventory) >= 1,
    file.exists(
      file.path(OUT, "21A_metadata_required_fields.csv")
    ),
    file.exists(
      file.path(OUT, "21A_assay_layer_inventory.csv")
    ),
    file.exists(
      file.path(OUT, "21A_package_inventory.csv")
    ),
    TRUE,
    TRUE,
    TRUE,
    TRUE,
    TRUE
  )
)

print(checks, n=Inf)

stopifnot(all(checks$pass))

audit <- c(
  "PIPELINE_V2 — 21A cell-cell communication preflight",
  "",
  paste0("Date: ", Sys.Date()),
  "",
  paste0("Metadata-audit object: ", audit_file),
  paste0("Cells: ", ncol(obj)),
  paste0("Features: ", nrow(obj)),
  "",
  "Purpose:",
  "- inventory frozen Seurat objects",
  "- identify patient/disease/batch metadata",
  "- identify L1/L2 metadata candidates",
  "- inventory expression assays/layers",
  "- inventory communication-analysis packages",
  "",
  "No cell-cell communication inference performed.",
  "No disease comparison performed.",
  "No differential-expression model refitted.",
  "No pathway analysis refitted.",
  "No cell annotation modified.",
  "No source object overwritten."
)

writeLines(
  audit,
  file.path(OUT, "21A_audit.txt")
)

cat("\n============================================================\n")
cat("21A PREFLIGHT COMPLETE\n")
cat("============================================================\n")
cat("No communication inference performed\n")
cat("No analytical refit\n")
cat("No annotation modified\n")
cat("PASS\n")
