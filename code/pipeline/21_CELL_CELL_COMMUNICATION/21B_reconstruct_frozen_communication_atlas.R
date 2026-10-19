suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(readr)
  library(tibble)
  library(stringr)
  library(purrr)
})

cat("============================================================\n")
cat("PIPELINE_V2 — 21B FROZEN COMMUNICATION ATLAS\n")
cat("============================================================\n")

# ============================================================
# 0. Paths
# ============================================================

GLOBAL <- "PIPELINE_V2/objects/ambient_corrected/global_atlas_v2_SoupX_expression.rds"

L1_OBJ <- "PIPELINE_V2/objects/annotated/global_atlas_v2_L1_frozen.rds"

EPI_OBJ <- "PIPELINE_V2/objects/annotated/epithelial_v2_L2_frozen.rds"

IMM_OBJ <- "PIPELINE_V2/objects/annotated/immune_v2_08I_L2_frozen.rds"

FIB_OBJ <- "PIPELINE_V2/objects/annotated/fibroblast_v2_09H5B_L2_frozen_4685.rds"

VAS_OBJ <- "PIPELINE_V2/objects/annotated/vascular_v2_10B6_L2_frozen_3023.rds"

MUR_OBJ <- paste0(
  "PIPELINE_V2/results/annotation/11_MURAL_L2/",
  "11B5_freeze/mural_v2_11B5_L2_frozen.rds"
)

SWEAT_OBJ <- paste0(
  "PIPELINE_V2/results/annotation/12_SWEAT_GLAND_L2/",
  "12D_freeze/sweat_gland_v2_12D_L2_frozen.rds"
)

MEL_OBJ <- paste0(
  "PIPELINE_V2/results/annotation/13_MELANOCYTE_L2/",
  "13D_freeze/melanocyte_v2_13D_L2_frozen.rds"
)

SEB_OBJ <- paste0(
  "PIPELINE_V2/results/annotation/14_SEBOCYTE_L2/",
  "14D_freeze/sebocyte_v2_14D_L2_frozen.rds"
)

FIB_RETURN <- paste0(
  "PIPELINE_V2/objects/annotated/",
  "epithelial_v2_09H4A_C9_returned_unadjudicated.rds"
)

MEL_RETURN <- paste0(
  "PIPELINE_V2/results/annotation/13_MELANOCYTE_L2/",
  "13C_C4_boundary/",
  "melanocyte_v2_13C_return_to_epithelial_39.rds"
)

OUT <- paste0(
  "PIPELINE_V2/results/disease_biology/",
  "21_CELL_CELL_COMMUNICATION/",
  "21B_FROZEN_COMMUNICATION_ATLAS"
)

TAB <- file.path(OUT, "tables")

OBJ_OUT <- paste0(
  "PIPELINE_V2/objects/frozen/",
  "global_atlas_v2_21B_communication_frozen.rds"
)

dir.create(OUT, recursive=TRUE, showWarnings=FALSE)
dir.create(TAB, recursive=TRUE, showWarnings=FALSE)
dir.create(dirname(OBJ_OUT), recursive=TRUE, showWarnings=FALSE)

required <- c(
  GLOBAL,
  L1_OBJ,
  EPI_OBJ,
  IMM_OBJ,
  FIB_OBJ,
  VAS_OBJ,
  MUR_OBJ,
  SWEAT_OBJ,
  MEL_OBJ,
  SEB_OBJ,
  FIB_RETURN,
  MEL_RETURN
)

stopifnot(all(file.exists(required)))

cat("\nSTEP 1/10 — loading global expression universe\n")

obj <- readRDS(GLOBAL)

stopifnot(
  inherits(obj, "Seurat"),
  ncol(obj) == 82893,
  nrow(obj) == 23412,
  all(c("patient_id","disease","batch") %in% colnames(obj[[]]))
)

DefaultAssay(obj) <- "RNA"

rna_layers <- Layers(obj[["RNA"]])

stopifnot(
  "counts" %in% rna_layers,
  "data" %in% rna_layers
)

cat("Cells   :", ncol(obj), "\n")
cat("Features:", nrow(obj), "\n")
cat("RNA layers:", paste(rna_layers, collapse="; "), "\n")

# ============================================================
# 2. Recover frozen L1
# ============================================================

cat("\nSTEP 2/10 — transferring frozen L1 annotation\n")

l1 <- readRDS(L1_OBJ)

stopifnot(
  inherits(l1, "Seurat"),
  ncol(l1) == 82893,
  "L1" %in% colnames(l1[[]]),
  setequal(colnames(obj), colnames(l1))
)

l1_md <- l1[[]] %>%
  rownames_to_column("barcode") %>%
  select(
    barcode,
    L1,
    any_of(c(
      "L1_status",
      "L1_flag",
      "L1_frozen",
      "annotation_version",
      "annotation_expression_source"
    ))
  )

idx <- match(colnames(obj), l1_md$barcode)

stopifnot(
  !anyNA(idx),
  !anyDuplicated(l1_md$barcode)
)

obj$communication_L1 <- l1_md$L1[idx]

stopifnot(!anyNA(obj$communication_L1))

# ============================================================
# 3. Helper for canonical frozen L2 transfer
# ============================================================

cat("\nSTEP 3/10 — defining exact-barcode L2 transfer\n")

extract_L2 <- function(
  file,
  label_col,
  compartment,
  status_col=NULL
) {

  x <- readRDS(file)

  stopifnot(
    inherits(x, "Seurat"),
    label_col %in% colnames(x[[]])
  )

  md <- x[[]] %>%
    rownames_to_column("barcode")

  if (!is.null(status_col)) {
    stopifnot(status_col %in% names(md))
  }

  out <- tibble(
    barcode=md$barcode,
    communication_compartment=compartment,
    communication_L2=as.character(md[[label_col]]),
    annotation_source=basename(file)
  )

  if (!is.null(status_col)) {
    out$source_L2_status <- as.character(
      md[[status_col]]
    )
  } else {
    out$source_L2_status <- NA_character_
  }

  stopifnot(
    !anyDuplicated(out$barcode),
    !anyNA(out$communication_L2),
    all(out$barcode %in% colnames(obj))
  )

  out
}

# ============================================================
# 4. Canonical frozen L2 sources
# ============================================================

cat("\nSTEP 4/10 — recovering canonical frozen L2 sources\n")

epi <- extract_L2(
  EPI_OBJ,
  "L2_identity",
  "Epithelial",
  "L2_status"
)

imm <- extract_L2(
  IMM_OBJ,
  "immune_L2",
  "Immune",
  "immune_L2_status"
)

fib <- extract_L2(
  FIB_OBJ,
  "fibroblast_L2_label",
  "Fibroblast",
  "fibroblast_L2_status"
)

vas <- extract_L2(
  VAS_OBJ,
  "vascular_L2_label",
  "Vascular",
  "vascular_L2_status"
)

mur <- extract_L2(
  MUR_OBJ,
  "mural_L2_label",
  "Mural",
  "mural_L2_status"
)

sweat <- extract_L2(
  SWEAT_OBJ,
  "sweat_L2_label",
  "Sweat_gland",
  "sweat_L2_status"
)

mel <- extract_L2(
  MEL_OBJ,
  "melanocyte_L2_label",
  "Melanocyte",
  "melanocyte_L2_status"
)

seb <- extract_L2(
  SEB_OBJ,
  "sebocyte_L2_label",
  "Sebocyte",
  "sebocyte_L2_status"
)

source_inventory <- tibble(
  compartment=c(
    "Epithelial",
    "Immune",
    "Fibroblast",
    "Vascular",
    "Mural",
    "Sweat_gland",
    "Melanocyte",
    "Sebocyte"
  ),
  file=c(
    EPI_OBJ,
    IMM_OBJ,
    FIB_OBJ,
    VAS_OBJ,
    MUR_OBJ,
    SWEAT_OBJ,
    MEL_OBJ,
    SEB_OBJ
  ),
  n_cells=c(
    nrow(epi),
    nrow(imm),
    nrow(fib),
    nrow(vas),
    nrow(mur),
    nrow(sweat),
    nrow(mel),
    nrow(seb)
  )
)

print(source_inventory, n=Inf)

# ============================================================
# 5. Recover returned epithelial cells
# ============================================================

cat("\nSTEP 5/10 — recovering returned epithelial boundary cells\n")

fib_ret_obj <- readRDS(FIB_RETURN)
mel_ret_obj <- readRDS(MEL_RETURN)

fib_ret <- colnames(fib_ret_obj)
mel_ret <- colnames(mel_ret_obj)

stopifnot(
  length(fib_ret) == 29,
  length(mel_ret) == 39,
  !anyDuplicated(fib_ret),
  !anyDuplicated(mel_ret),
  length(intersect(fib_ret,mel_ret)) == 0,
  all(fib_ret %in% colnames(obj)),
  all(mel_ret %in% colnames(obj))
)

returned <- bind_rows(
  tibble(
    barcode=fib_ret,
    return_source="FIBROBLAST_REVIEW",
    return_n=29L
  ),
  tibble(
    barcode=mel_ret,
    return_source="MELANOCYTE_REVIEW",
    return_n=39L
  )
)

cat("Fibroblast-returned epithelial:", length(fib_ret), "\n")
cat("Melanocyte-returned epithelial:", length(mel_ret), "\n")
cat("Total returned:", nrow(returned), "\n")

# ============================================================
# 6. Assemble exact-barcode annotation map
# ============================================================

cat("\nSTEP 6/10 — assembling frozen communication annotation map\n")

l2_map <- bind_rows(
  epi,
  imm,
  fib,
  vas,
  mur,
  sweat,
  mel,
  seb
)

# A cell must not inherit two canonical frozen L2 labels.
dup_l2 <- l2_map %>%
  count(barcode, name="n") %>%
  filter(n > 1)

cat("Duplicated canonical L2 barcodes:", nrow(dup_l2), "\n")

stopifnot(nrow(dup_l2) == 0)

# Returned cells are special:
# they must not be forced into an epithelial L2.
returned_overlap <- l2_map %>%
  filter(barcode %in% returned$barcode)

cat(
  "Returned cells overlapping canonical L2 map:",
  nrow(returned_overlap),
  "\n"
)

# Remove returned cells from any inherited L2 assignment.
l2_map_clean <- l2_map %>%
  filter(!barcode %in% returned$barcode)

comm <- tibble(
  barcode=colnames(obj),
  communication_L1=as.character(obj$communication_L1)
) %>%
  left_join(
    l2_map_clean,
    by="barcode"
  ) %>%
  left_join(
    returned %>%
      select(barcode,return_source),
    by="barcode"
  )

# Returned cells remain epithelial at communication L1,
# but explicitly have no adjudicated L2.
comm <- comm %>%
  mutate(
    communication_L1=if_else(
      !is.na(return_source),
      "Epithelial",
      communication_L1
    ),
    communication_annotation_level=case_when(
      !is.na(return_source) ~ "RETURNED_UNADJUDICATED",
      !is.na(communication_L2) ~ "L2",
      TRUE ~ "L1"
    ),
    communication_identity=case_when(
      communication_annotation_level == "L2" ~
        communication_L2,
      communication_annotation_level ==
        "RETURNED_UNADJUDICATED" ~
        paste0(
          "Epithelial_returned_unadjudicated__",
          return_source
        ),
      TRUE ~
        paste0("L1__",communication_L1)
    ),
    communication_L2_frozen=
      communication_annotation_level == "L2"
  )

stopifnot(
  nrow(comm) == 82893,
  !anyDuplicated(comm$barcode),
  !anyNA(comm$communication_identity),
  !anyNA(comm$communication_annotation_level)
)

# ============================================================
# 7. Transfer communication metadata to global object
# ============================================================

cat("\nSTEP 7/10 — transferring communication metadata\n")

idx <- match(colnames(obj),comm$barcode)

stopifnot(!anyNA(idx))

obj$communication_L1 <-
  comm$communication_L1[idx]

obj$communication_compartment <-
  comm$communication_compartment[idx]

obj$communication_L2 <-
  comm$communication_L2[idx]

obj$communication_identity <-
  comm$communication_identity[idx]

obj$communication_annotation_level <-
  comm$communication_annotation_level[idx]

obj$communication_L2_frozen <-
  comm$communication_L2_frozen[idx]

obj$communication_return_source <-
  comm$return_source[idx]

obj$communication_annotation_version <-
  "PIPELINE_V2_21B_FROZEN_COMMUNICATION_ATLAS"

obj$communication_expression_source <-
  "SoupX_corrected"

# ============================================================
# 8. Communication-atlas summaries
# ============================================================

cat("\nSTEP 8/10 — communication atlas summaries\n")

level_summary <- obj[[]] %>%
  count(
    communication_annotation_level,
    name="n_cells"
  ) %>%
  arrange(desc(n_cells))

identity_summary <- obj[[]] %>%
  count(
    communication_L1,
    communication_annotation_level,
    communication_identity,
    name="n_cells"
  ) %>%
  arrange(
    communication_L1,
    desc(n_cells)
  )

patient_identity <- obj[[]] %>%
  count(
    disease,
    patient_id,
    communication_identity,
    name="n_cells"
  )

disease_identity <- patient_identity %>%
  group_by(
    disease,
    communication_identity
  ) %>%
  summarise(
    n_patients=n_distinct(patient_id),
    total_cells=sum(n_cells),
    min_cells=min(n_cells),
    median_cells=median(n_cells),
    max_cells=max(n_cells),
    .groups="drop"
  ) %>%
  arrange(
    disease,
    desc(n_patients),
    desc(total_cells)
  )

patient_support <- obj[[]] %>%
  distinct(
    patient_id,
    disease,
    batch
  ) %>%
  count(
    disease,
    batch,
    name="n_patients"
  )

cat("\n===== ANNOTATION LEVELS =====\n")
print(tibble::as_tibble(level_summary), n=Inf, width=Inf)

cat("\n===== COMMUNICATION IDENTITIES =====\n")
print(tibble::as_tibble(identity_summary), n=Inf, width=Inf)

cat("\n===== PATIENT SUPPORT =====\n")
print(tibble::as_tibble(patient_support), n=Inf, width=Inf)

# ============================================================
# 9. Write frozen object and tables
# ============================================================

cat("\nSTEP 9/10 — writing frozen communication atlas\n")

write_csv(
  source_inventory,
  file.path(
    TAB,
    "21B_annotation_source_inventory.csv"
  )
)

write_csv(
  comm,
  file.path(
    TAB,
    "21B_cell_annotation_map.csv"
  )
)

write_csv(
  returned,
  file.path(
    TAB,
    "21B_returned_epithelial_cells.csv"
  )
)

write_csv(
  level_summary,
  file.path(
    TAB,
    "21B_annotation_level_summary.csv"
  )
)

write_csv(
  identity_summary,
  file.path(
    TAB,
    "21B_identity_summary.csv"
  )
)

write_csv(
  patient_identity,
  file.path(
    TAB,
    "21B_patient_identity_counts.csv"
  )
)

write_csv(
  disease_identity,
  file.path(
    TAB,
    "21B_disease_identity_support.csv"
  )
)

write_csv(
  patient_support,
  file.path(
    TAB,
    "21B_patient_batch_support.csv"
  )
)

saveRDS(
  obj,
  OBJ_OUT,
  compress=FALSE
)

stopifnot(file.exists(OBJ_OUT))

# ============================================================
# 10. Formal checks / audit
# ============================================================

cat("\nSTEP 10/10 — formal checks\n")

md <- obj[[]]

returned_md <- md[
  !is.na(md$communication_return_source),
  ,
  drop=FALSE
]

checks <- tibble(
  check=c(
    "global_expression_universe_82893",
    "RNA_counts_present",
    "RNA_data_present",
    "patient_id_present",
    "disease_present",
    "batch_present",
    "frozen_L1_transferred",
    "canonical_L2_sources_loaded",
    "canonical_L2_barcode_unique",
    "fibroblast_returned_29",
    "melanocyte_returned_39",
    "returned_total_68",
    "returned_cells_no_L2",
    "returned_cells_explicitly_flagged",
    "communication_identity_complete",
    "communication_object_written",
    "no_clustering_performed",
    "no_DE_refit",
    "no_pathway_refit",
    "no_communication_inference",
    "source_objects_not_overwritten"
  ),
  pass=c(
    ncol(obj) == 82893,
    "counts" %in% Layers(obj[["RNA"]]),
    "data" %in% Layers(obj[["RNA"]]),
    "patient_id" %in% names(md),
    "disease" %in% names(md),
    "batch" %in% names(md),
    !anyNA(md$communication_L1),
    nrow(source_inventory) == 8,
    nrow(dup_l2) == 0,
    length(fib_ret) == 29,
    length(mel_ret) == 39,
    nrow(returned) == 68,
    all(
      is.na(
        returned_md$communication_L2
      )
    ),
    all(
      returned_md$communication_annotation_level ==
        "RETURNED_UNADJUDICATED"
    ),
    !anyNA(
      md$communication_identity
    ),
    file.exists(OBJ_OUT),
    TRUE,
    TRUE,
    TRUE,
    TRUE,
    TRUE
  )
)

print(checks,n=Inf)

stopifnot(all(checks$pass))

write_csv(
  checks,
  file.path(
    TAB,
    "21B_formal_checks.csv"
  )
)

audit <- c(
  "PIPELINE_V2 — 21B frozen communication atlas",
  "",
  paste0("Date: ",Sys.Date()),
  "",
  "Global expression universe:",
  "- 82,893 cells",
  "- SoupX-corrected RNA counts/data retained",
  "",
  "Canonical frozen L2 sources:",
  "- epithelial: 07 epithelial L2",
  "- immune: 08I",
  "- fibroblast: 09H5B",
  "- vascular: 10B6",
  "- mural: 11B5",
  "- sweat gland: 12D",
  "- melanocyte: 13D",
  "- sebocyte: 14D",
  "",
  "Returned epithelial boundary cells:",
  "- 29 cells returned from fibroblast review",
  "- 39 cells returned from melanocyte review",
  "- 68 returned cells total",
  "- retained as epithelial at communication L1",
  "- no returned cell assigned to a frozen epithelial L2",
  "- explicitly labelled RETURNED_UNADJUDICATED",
  "",
  "Communication identity construction:",
  "- canonical frozen L2 transferred by exact barcode",
  "- cells without canonical L2 retained at L1",
  "- returned epithelial cells retained as explicit unadjudicated identities",
  "",
  "No clustering performed.",
  "No differential-expression refit.",
  "No pathway refit.",
  "No cell-cell communication inference performed.",
  "No source object overwritten."
)

writeLines(
  audit,
  file.path(
    OUT,
    "21B_audit.txt"
  )
)

cat("\n===== FORMAL CHECKS =====\n")
print(checks,n=Inf)

cat("\n============================================================\n")
cat("21B COMPLETE\n")
cat("============================================================\n")
cat("82,893-cell communication atlas reconstructed\n")
cat("Canonical frozen L2 annotations transferred by exact barcode\n")
cat("29 fibroblast-returned epithelial cells retained unadjudicated\n")
cat("39 melanocyte-returned epithelial cells retained unadjudicated\n")
cat("No returned cell forced into an epithelial L2 identity\n")
cat("No clustering performed\n")
cat("No differential-expression refit\n")
cat("No pathway refit\n")
cat("No communication inference performed\n")
cat("Source objects not overwritten\n")
cat("PASS\n")
