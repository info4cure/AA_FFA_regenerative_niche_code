suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(readr)
})

cat("\n============================================================\n")
cat("PIPELINE_V2 — 07A EPITHELIAL L2 PREFLIGHT\n")
cat("============================================================\n\n")

IN <- "PIPELINE_V2/objects/annotated/global_atlas_v2_L1_frozen.rds"

OUT <- "PIPELINE_V2/results/annotation/07_EPITHELIAL_L2/07A_preflight"
TDIR <- file.path(OUT,"tables")

dir.create(TDIR,recursive=TRUE,showWarnings=FALSE)

x <- readRDS(IN)
DefaultAssay(x) <- "RNA"

# ------------------------------------------------------------
# 1. Canonical invariants
# ------------------------------------------------------------

stopifnot(
  ncol(x)==82893,
  nrow(x)==23412,
  all(x$L1_frozen),
  "L1" %in% colnames(x@meta.data),
  "RNA_snn_res.0.6" %in% colnames(x@meta.data)
)

epi_cells <- colnames(x)[x$L1=="Epithelial"]

stopifnot(length(epi_cells)==63893)

epi <- subset(x,cells=epi_cells)

cat("Global cells:",ncol(x),"\n")
cat("Epithelial cells:",ncol(epi),"\n")
cat("Genes:",nrow(epi),"\n\n")

# ------------------------------------------------------------
# 2. Metadata inventory
# ------------------------------------------------------------

cat("=== METADATA COLUMNS ===\n")
print(colnames(epi@meta.data))

writeLines(
  colnames(epi@meta.data),
  file.path(TDIR,"07A_metadata_columns.txt")
)

# ------------------------------------------------------------
# 3. Reductions inventory
# ------------------------------------------------------------

cat("\n=== STORED REDUCTIONS ===\n")

reds <- names(epi@reductions)

if(length(reds)==0){
  cat("NONE\n")
} else {
  for(r in reds){
    dims <- dim(Embeddings(epi,r))
    cat(
      r,
      ":",
      dims[1],"cells x",
      dims[2],"dimensions\n"
    )
  }
}

writeLines(
  if(length(reds)) reds else "NONE",
  file.path(TDIR,"07A_reductions.txt")
)

# ------------------------------------------------------------
# 4. Disease composition
# ------------------------------------------------------------

needed <- c(
  "patient_id",
  "disease",
  "batch"
)

miss <- setdiff(
  needed,
  colnames(epi@meta.data)
)

if(length(miss)>0){
  stop(
    "Missing metadata: ",
    paste(miss,collapse=", ")
  )
}

disease_tab <- epi@meta.data %>%
  count(disease,name="n_cells") %>%
  mutate(
    pct=100*n_cells/sum(n_cells)
  ) %>%
  arrange(desc(n_cells))

cat("\n=== EPITHELIAL CELLS BY DISEASE ===\n")
print(as.data.frame(disease_tab),row.names=FALSE)

write_csv(
  disease_tab,
  file.path(TDIR,"07A_epithelial_by_disease.csv")
)

# ------------------------------------------------------------
# 5. Patient composition
# ------------------------------------------------------------

patient_tab <- epi@meta.data %>%
  count(
    disease,
    patient_id,
    batch,
    name="n_cells"
  ) %>%
  group_by(disease) %>%
  mutate(
    pct_within_disease=
      100*n_cells/sum(n_cells)
  ) %>%
  ungroup() %>%
  arrange(disease,desc(n_cells))

cat("\n=== EPITHELIAL CELLS BY PATIENT ===\n")
print(as.data.frame(patient_tab),row.names=FALSE)

write_csv(
  patient_tab,
  file.path(TDIR,"07A_epithelial_by_patient.csv")
)

# ------------------------------------------------------------
# 6. Global-cluster composition
# ------------------------------------------------------------

cluster_tab <- epi@meta.data %>%
  mutate(
    global_cluster=
      as.character(.data[["RNA_snn_res.0.6"]])
  ) %>%
  count(
    global_cluster,
    disease,
    name="n_cells"
  ) %>%
  group_by(global_cluster) %>%
  mutate(
    pct_disease=
      100*n_cells/sum(n_cells)
  ) %>%
  ungroup() %>%
  arrange(
    as.integer(global_cluster),
    disease
  )

write_csv(
  cluster_tab,
  file.path(
    TDIR,
    "07A_global_cluster_by_disease.csv"
  )
)

cluster_total <- epi@meta.data %>%
  mutate(
    global_cluster=
      as.character(.data[["RNA_snn_res.0.6"]])
  ) %>%
  count(
    global_cluster,
    name="n_cells"
  ) %>%
  arrange(
    as.integer(global_cluster)
  )

cat("\n=== EPITHELIAL GLOBAL CLUSTERS ===\n")
print(as.data.frame(cluster_total),row.names=FALSE)

write_csv(
  cluster_total,
  file.path(
    TDIR,
    "07A_global_cluster_counts.csv"
  )
)

# ------------------------------------------------------------
# 7. Patient x global-cluster support
# ------------------------------------------------------------

patient_cluster <- epi@meta.data %>%
  mutate(
    global_cluster=
      as.character(.data[["RNA_snn_res.0.6"]])
  ) %>%
  count(
    disease,
    patient_id,
    global_cluster,
    name="n_cells"
  )

write_csv(
  patient_cluster,
  file.path(
    TDIR,
    "07A_patient_global_cluster_counts.csv"
  )
)

support <- patient_cluster %>%
  group_by(global_cluster) %>%
  summarise(
    n_patients=n_distinct(patient_id),
    n_HC=n_distinct(
      patient_id[disease=="HC"]
    ),
    n_AA=n_distinct(
      patient_id[disease=="AA"]
    ),
    n_FFA=n_distinct(
      patient_id[disease=="FFA"]
    ),
    .groups="drop"
  ) %>%
  arrange(
    as.integer(global_cluster)
  )

cat("\n=== PATIENT SUPPORT PER GLOBAL CLUSTER ===\n")
print(as.data.frame(support),row.names=FALSE)

write_csv(
  support,
  file.path(
    TDIR,
    "07A_global_cluster_patient_support.csv"
  )
)

# ------------------------------------------------------------
# 8. QC distributions by disease
# ------------------------------------------------------------

qc_candidates <- c(
  "nCount_RNA",
  "nFeature_RNA",
  "percent.mt",
  "doublet_score"
)

qc_present <- intersect(
  qc_candidates,
  colnames(epi@meta.data)
)

cat("\n=== QC METADATA AVAILABLE ===\n")
cat(
  if(length(qc_present))
    paste(qc_present,collapse=", ")
  else
    "NONE",
  "\n"
)

if(length(qc_present)>0){

  qc_summary <- epi@meta.data %>%
    select(
      disease,
      all_of(qc_present)
    ) %>%
    group_by(disease) %>%
    summarise(
      across(
        all_of(qc_present),
        list(
          median=~median(.x,na.rm=TRUE),
          q25=~quantile(.x,0.25,na.rm=TRUE),
          q75=~quantile(.x,0.75,na.rm=TRUE)
        )
      ),
      .groups="drop"
    )

  write_csv(
    qc_summary,
    file.path(
      TDIR,
      "07A_QC_by_disease.csv"
    )
  )
}

# ------------------------------------------------------------
# 9. Summary for downstream L2 design
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("07A SUMMARY\n")
cat("============================================================\n")
cat("Epithelial cells:",ncol(epi),"\n")
cat(
  "Patients:",
  n_distinct(epi$patient_id),
  "\n"
)
cat(
  "Diseases:",
  paste(
    sort(unique(epi$disease)),
    collapse=", "
  ),
  "\n"
)
cat(
  "Global epithelial clusters:",
  n_distinct(epi$RNA_snn_res.0.6),
  "\n"
)
cat(
  "Stored reductions:",
  if(length(reds))
    paste(reds,collapse=", ")
  else
    "NONE",
  "\n"
)

cat("\nPASS — preflight only; no L2 analysis performed.\n")
