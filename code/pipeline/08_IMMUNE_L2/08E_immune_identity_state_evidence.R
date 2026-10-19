suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(readr)
})

cat("\n============================================================\n")
cat("PIPELINE_V2 — 08E IMMUNE IDENTITY / STATE / CONTAMINATION AUDIT\n")
cat("============================================================\n\n")

IN <- "PIPELINE_V2/objects/annotated/immune_v2_08B_reclustered.rds"
OUT <- "PIPELINE_V2/results/annotation/08_IMMUNE_L2/08E_identity_state"
TDIR <- file.path(OUT,"tables")
dir.create(TDIR,recursive=TRUE,showWarnings=FALSE)

CC <- "immune_snn_res.0.6"
EXPECTED_CELLS <- 4762
EXPECTED_GENES <- 23412

cat("STEP 1/10 — loading immune object\n")

x <- readRDS(IN)
DefaultAssay(x) <- "RNA"

stopifnot(
  ncol(x)==EXPECTED_CELLS,
  nrow(x)==EXPECTED_GENES,
  CC %in% colnames(x@meta.data)
)

x$immune_cluster_06 <- factor(
  as.character(x@meta.data[[CC]]),
  levels=as.character(0:16)
)

stopifnot(!anyNA(x$immune_cluster_06))

cluster_levels <- as.character(0:16)
clusters <- as.character(x$immune_cluster_06)

cat("Cells:",ncol(x),"\n")
cat("Genes:",nrow(x),"\n")
cat("Patients:",length(unique(x$patient_id)),"\n")
cat("Working resolution: 0.6\n")
cat("Clusters:",length(cluster_levels),"\n")

cat("\nSTEP 2/10 — defining programmes\n")

identity_sets <- list(
  Leukocyte_core=c("PTPRC","LST1","CD3D","CD3E"),
  T_core=c("CD3D","CD3E","TRAC","CD247"),
  CD4_T=c("CD3D","CD3E","TRAC","CD4","IL7R","LTB","MAL","TCF7"),
  CD8_T=c("CD3D","CD3E","TRAC","CD8A","CD8B","GZMK","CCL5","GZMH"),
  Treg=c("FOXP3","IL2RA","CTLA4","IKZF2","TNFRSF4","TIGIT","LAYN"),
  NK=c("NKG7","GNLY","KLRD1","NCR1","FGFBP2","PRF1","KLRF1"),
  B_cell=c("CD79A","CD79B","MS4A1","CD37","PAX5","CD74","HLA-DRA"),
  Plasma=c("MZB1","JCHAIN","TNFRSF17","SDC1","DERL3","XBP1"),
  Monocyte=c("LYZ","FCN1","CTSS","S100A8","S100A9","VCAN","LILRB1"),
  Macrophage=c("C1QA","C1QB","C1QC","CD163","MRC1","MERTK","MSR1","STAB1"),
  cDC2=c("CD1C","CLEC10A","FCER1A","CD1E","CD1B","CST3"),
  cDC1=c("CLEC9A","XCR1","WDFY4","BATF3","CADM1","IRF8"),
  pDC=c("CLEC4C","GZMB","IL3RA","TCF4","SERPINF1","IRF7"),
  Langerhans=c("CD207","CD1A","CD1E","FCER1A","HLA-DRA","EPCAM"),
  Mast=c("TPSAB1","TPSB2","CPA3","KIT","MS4A2","HDC","CMA1")
)

nonimmune_sets <- list(
  Epithelial=c(
    "KRT5","KRT14","KRT15","KRT17",
    "KRT19","COL17A1","EPCAM","KRTDAP"
  ),
  Neural_Schwann=c(
    "PMP2","MPZ","PLP1","SOX10",
    "S100B","CDH19","SLC1A3"
  ),
  Fibroblast=c(
    "COL1A1","COL1A2","COL3A1",
    "DCN","LUM","PDGFRA"
  ),
  Endothelial=c(
    "PECAM1","VWF","EMCN","KDR","CLDN5","ERG"
  )
)

state_sets <- list(
  Cycling=c("MKI67","TOP2A","UBE2C","BIRC5","CENPF","RRM2"),
  IFN_response=c("ISG15","IFIT1","IFIT2","IFIT3","MX1","OAS1","RSAD2","IFI6"),
  Cytotoxic=c("NKG7","GNLY","PRF1","GZMB","GZMH","CCL5"),
  Activation_exhaustion=c("PDCD1","LAG3","TIGIT","HAVCR2","TOX","CXCL13","TNFRSF9"),
  Tissue_resident=c("CXCR6","ZNF683","ITGAE","CD69","RGS1")
)

all_sets <- c(identity_sets,nonimmune_sets,state_sets)

class_map <- c(
  setNames(rep("identity",length(identity_sets)),names(identity_sets)),
  setNames(rep("nonimmune",length(nonimmune_sets)),names(nonimmune_sets)),
  setNames(rep("state",length(state_sets)),names(state_sets))
)

# ==========================================================================
# 3. MARKER COVERAGE
# ==========================================================================

cat("\nSTEP 3/10 — marker coverage\n")

coverage <- bind_rows(
  lapply(names(all_sets),function(nm){

    requested <- all_sets[[nm]]
    present <- intersect(requested,rownames(x))

    data.frame(
      programme=nm,
      class=unname(class_map[nm]),
      n_present=length(present),
      n_requested=length(requested),
      coverage=length(present)/length(requested),
      missing=paste(setdiff(requested,present),collapse=";"),
      stringsAsFactors=FALSE
    )
  })
)

cat("\n=== PROGRAMME COVERAGE ===\n")
print(as.data.frame(coverage),row.names=FALSE)

write_csv(
  coverage,
  file.path(TDIR,"08E_programme_coverage.csv")
)

stopifnot(
  all(coverage$n_present > 0)
)

# ==========================================================================
# 4. CLUSTER-LEVEL PROGRAMME EVIDENCE
# ==========================================================================

cat("\nSTEP 4/10 — cluster-level programme evidence\n")

genes <- intersect(
  unique(unlist(all_sets)),
  rownames(x)
)

dat <- LayerData(
  x,
  assay="RNA",
  layer="data"
)[genes,,drop=FALSE]

programme_rows <- list()

for(cl in cluster_levels){

  cells <- which(clusters==cl)
  stopifnot(length(cells)>0)

  for(nm in names(all_sets)){

    g <- intersect(
      all_sets[[nm]],
      rownames(dat)
    )

    stopifnot(length(g)>0)

    vals <- Matrix::colMeans(
      dat[g,cells,drop=FALSE]
    )

    programme_rows[[length(programme_rows)+1]] <- data.frame(
      cluster=cl,
      programme=nm,
      class=unname(class_map[nm]),
      n_cells=length(cells),
      mean_programme=mean(vals),
      median_programme=median(vals),
      pct_positive=100*mean(vals>0),
      stringsAsFactors=FALSE
    )
  }
}

ev <- bind_rows(programme_rows)

ev <- ev %>%
  group_by(programme) %>%
  mutate(
    programme_z=ifelse(
      sd(mean_programme)>0,
      as.numeric(scale(mean_programme)),
      0
    )
  ) %>%
  ungroup()

write_csv(
  ev,
  file.path(TDIR,"08E_cluster_programme_evidence.csv")
)

# ==========================================================================
# 5. SINGLE-GENE EVIDENCE
# ==========================================================================

cat("\nSTEP 5/10 — key single-gene evidence\n")

panel <- unique(c(

  # Leukocyte backbone
  "PTPRC","LST1","CD3D","CD3E","TRAC",

  # T-cell
  "CD4","CD8A","CD8B",
  "IL7R","TCF7","KLRB1","CCR6",

  # Treg
  "FOXP3","IL2RA","CTLA4","IKZF2",

  # NK / cytotoxic
  "NKG7","GNLY","KLRD1","NCR1",
  "FGFBP2","PRF1",

  # B / plasma
  "CD79A","CD79B","MS4A1","PAX5",
  "MZB1","JCHAIN","TNFRSF17","SDC1",

  # Monocyte / macrophage
  "LYZ","FCN1","CTSS","S100A8","S100A9",
  "C1QA","C1QB","C1QC",
  "CD163","MRC1","MERTK","STAB1",

  # DC / Langerhans
  "CD1C","CLEC10A","FCER1A","CD1B",
  "CLEC9A","XCR1","WDFY4","BATF3",
  "CLEC4C","IL3RA",
  "CD207","CD1A",

  # Mast
  "TPSAB1","TPSB2","CPA3",
  "KIT","MS4A2","HDC",

  # States
  "PDCD1","LAG3","CXCL13",
  "CXCR6","ZNF683","ITGAE",
  "MKI67","TOP2A","UBE2C",
  "ISG15","IFIT3","MX1",

  # Epithelial
  "KRT5","KRT14","KRT15","KRT17",
  "KRT19","COL17A1","EPCAM","KRTDAP",

  # Neural / Schwann
  "PMP2","MPZ","PLP1",
  "SOX10","S100B","CDH19",

  # Fibroblast
  "COL1A1","COL1A2","DCN","LUM",

  # Endothelial
  "PECAM1","VWF","EMCN","CLDN5"
))

panel <- intersect(panel,rownames(x))

pd <- LayerData(
  x,
  assay="RNA",
  layer="data"
)[panel,,drop=FALSE]

gene_rows <- list()

for(cl in cluster_levels){

  cells <- which(clusters==cl)

  for(g in panel){

    v <- pd[g,cells]

    gene_rows[[length(gene_rows)+1]] <- data.frame(
      cluster=cl,
      gene=g,
      mean_expression=mean(v),
      pct_detected=100*mean(v>0),
      stringsAsFactors=FALSE
    )
  }
}

gene_ev <- bind_rows(gene_rows)

write_csv(
  gene_ev,
  file.path(TDIR,"08E_key_gene_evidence.csv")
)

# ==========================================================================
# 6. PATIENT / BATCH / QC SUPPORT
# ==========================================================================

cat("\nSTEP 6/10 — patient, batch and QC support\n")

md <- x@meta.data %>%
  mutate(
    cluster=as.character(.data[[CC]])
  )

patient_support <- md %>%
  count(
    cluster,
    patient_id,
    disease,
    batch,
    name="n_cells"
  )

write_csv(
  patient_support,
  file.path(TDIR,"08E_patient_cluster_support.csv")
)

qc <- md %>%
  group_by(cluster) %>%
  summarise(
    n_cells=n(),
    n_patients=n_distinct(patient_id),
    n_HC=n_distinct(patient_id[disease=="HC"]),
    n_AA=n_distinct(patient_id[disease=="AA"]),
    n_FFA=n_distinct(patient_id[disease=="FFA"]),
    n_batches=n_distinct(batch),

    nCount_RNA_median=
      median(nCount_RNA,na.rm=TRUE),

    nFeature_RNA_median=
      median(nFeature_RNA,na.rm=TRUE),

    percent_mt_median=
      median(percent.mt,na.rm=TRUE),

    doublet_score_median=
      median(doublet_score,na.rm=TRUE),

    dominant_patient=
      names(
        sort(
          table(patient_id),
          decreasing=TRUE
        )
      )[1],

    dominant_patient_cells=
      max(table(patient_id)),

    dominant_patient_fraction=
      max(table(patient_id))/n(),

    .groups="drop"
  )

write_csv(
  qc,
  file.path(TDIR,"08E_cluster_QC_patient_support.csv")
)

cat("\n=== CLUSTER SUPPORT ===\n")
print(
  as.data.frame(
    qc %>%
      select(
        cluster,
        n_cells,
        n_patients,
        n_HC,
        n_AA,
        n_FFA,
        n_batches,
        dominant_patient,
        dominant_patient_fraction
      )
  ),
  row.names=FALSE
)

# ==========================================================================
# 7. IMMUNE VS NON-IMMUNE EVIDENCE
#
# Evidence flags only. No reassignment is performed here.
# ==========================================================================

cat("\nSTEP 7/10 — immune versus non-immune evidence\n")

immune_best <- ev %>%
  filter(
    class=="identity",
    programme!="Leukocyte_core"
  ) %>%
  group_by(cluster) %>%
  slice_max(
    programme_z,
    n=1,
    with_ties=FALSE
  ) %>%
  ungroup() %>%
  transmute(
    cluster,
    top_immune_programme=programme,
    top_immune_z=programme_z,
    top_immune_mean=mean_programme,
    top_immune_pct=pct_positive
  )

leukocyte_core <- ev %>%
  filter(programme=="Leukocyte_core") %>%
  transmute(
    cluster,
    leukocyte_core_z=programme_z,
    leukocyte_core_mean=mean_programme,
    leukocyte_core_pct=pct_positive
  )

nonimmune_best <- ev %>%
  filter(class=="nonimmune") %>%
  group_by(cluster) %>%
  slice_max(
    programme_z,
    n=1,
    with_ties=FALSE
  ) %>%
  ungroup() %>%
  transmute(
    cluster,
    top_nonimmune_programme=programme,
    top_nonimmune_z=programme_z,
    top_nonimmune_mean=mean_programme,
    top_nonimmune_pct=pct_positive
  )

disposition <- qc %>%
  left_join(
    immune_best,
    by="cluster"
  ) %>%
  left_join(
    leukocyte_core,
    by="cluster"
  ) %>%
  left_join(
    nonimmune_best,
    by="cluster"
  ) %>%
  mutate(

    immune_minus_nonimmune_z=
      top_immune_z-top_nonimmune_z,

    provisional_disposition=case_when(

      top_nonimmune_programme=="Epithelial" &
        top_nonimmune_z>2 &
        top_nonimmune_z>top_immune_z
        ~ "REVIEW_RETURN_EPITHELIAL",

      top_nonimmune_programme=="Neural_Schwann" &
        top_nonimmune_z>2 &
        top_nonimmune_z>top_immune_z
        ~ "REVIEW_RETURN_NEURAL_SCHWANN",

      top_nonimmune_z>2 &
        top_nonimmune_z>top_immune_z
        ~ "REVIEW_NONIMMUNE",

      TRUE
        ~ "KEEP_IMMUNE_CANDIDATE"
    )
  )

write_csv(
  disposition,
  file.path(
    TDIR,
    "08E_provisional_disposition.csv"
  )
)

cat("\n=== PROVISIONAL CLUSTER DISPOSITION ===\n")

print(
  as.data.frame(
    disposition %>%
      mutate(
        cluster_numeric=as.numeric(cluster)
      ) %>%
      arrange(cluster_numeric) %>%
      select(
        cluster,
        n_cells,
        n_patients,
        leukocyte_core_z,
        leukocyte_core_pct,
        top_immune_programme,
        top_immune_z,
        top_nonimmune_programme,
        top_nonimmune_z,
        immune_minus_nonimmune_z,
        provisional_disposition
      )
  ),
  row.names=FALSE
)

# ==========================================================================
# 8. CELL-LEVEL C9 / C12 LINEAGE AUDIT
#
# Cluster means can conceal mixed populations. C9 and C12 are therefore
# examined at single-cell level using independent immune, epithelial and
# Schwann lineage cores.
# ==========================================================================

cat("\nSTEP 8/10 — cell-level C9/C12 lineage audit\n")

cell_score_sets <- list(

  Immune_core=c(
    "PTPRC",
    "CD3D","CD3E","TRAC",
    "LST1","LYZ"
  ),

  Epithelial_core=c(
    "KRT5","KRT14","KRT15",
    "KRT17","KRT19",
    "EPCAM","COL17A1"
  ),

  Schwann_core=c(
    "PMP2","MPZ","PLP1",
    "SOX10","S100B","CDH19"
  )
)

score_genes <- intersect(
  unique(unlist(cell_score_sets)),
  rownames(x)
)

score_dat <- LayerData(
  x,
  assay="RNA",
  layer="data"
)[score_genes,,drop=FALSE]

cell_scores <- data.frame(
  cell=colnames(x),
  cluster=as.character(x$immune_cluster_06),
  patient_id=x$patient_id,
  disease=x$disease,
  batch=x$batch,
  stringsAsFactors=FALSE
)

for(nm in names(cell_score_sets)){

  g <- intersect(
    cell_score_sets[[nm]],
    rownames(score_dat)
  )

  stopifnot(length(g)>0)

  cell_scores[[nm]] <- Matrix::colMeans(
    score_dat[g,,drop=FALSE]
  )
}

cell_scores <- cell_scores %>%
  mutate(

    immune_minus_epithelial=
      Immune_core-Epithelial_core,

    immune_minus_schwann=
      Immune_core-Schwann_core,

    lineage_winner=case_when(

      Immune_core>=Epithelial_core &
        Immune_core>=Schwann_core
        ~ "Immune",

      Epithelial_core>Immune_core &
        Epithelial_core>=Schwann_core
        ~ "Epithelial",

      Schwann_core>Immune_core &
        Schwann_core>Epithelial_core
        ~ "Neural_Schwann",

      TRUE
        ~ "Uncertain"
    )
  )

c9c12_scores <- cell_scores %>%
  filter(cluster %in% c("9","12"))

stopifnot(
  nrow(c9c12_scores)==
    sum(clusters %in% c("9","12"))
)

write_csv(
  c9c12_scores,
  file.path(
    TDIR,
    "08E_C9_C12_cell_level_lineage_scores.csv"
  )
)

cell_score_summary <- c9c12_scores %>%
  count(
    cluster,
    lineage_winner,
    name="n_cells"
  ) %>%
  group_by(cluster) %>%
  mutate(
    fraction=n_cells/sum(n_cells)
  ) %>%
  ungroup()

cat("\n=== C9/C12 CELL-LEVEL LINEAGE WINNERS ===\n")

print(
  as.data.frame(cell_score_summary),
  row.names=FALSE
)

write_csv(
  cell_score_summary,
  file.path(
    TDIR,
    "08E_C9_C12_cell_level_lineage_summary.csv"
  )
)

# Patient-level reproducibility of the C9/C12 lineage result

cell_score_patient <- c9c12_scores %>%
  count(
    cluster,
    patient_id,
    disease,
    lineage_winner,
    name="n_cells"
  ) %>%
  group_by(
    cluster,
    patient_id,
    disease
  ) %>%
  mutate(
    fraction_within_patient_cluster=
      n_cells/sum(n_cells)
  ) %>%
  ungroup()

write_csv(
  cell_score_patient,
  file.path(
    TDIR,
    "08E_C9_C12_lineage_by_patient.csv"
  )
)

# ==========================================================================
# 9A. C9 EPITHELIAL RETURN AUDIT
# ==========================================================================

cat("\n============================================================\n")
cat("CLUSTER 9 — EPITHELIAL RETURN AUDIT\n")
cat("============================================================\n")

c9_cells <- colnames(x)[
  as.character(x$immune_cluster_06)=="9"
]

c9_md <- x@meta.data[c9_cells,,drop=FALSE]

cat("Cells:",length(c9_cells),"\n")

cat("\nFrozen source L1:\n")
print(
  sort(
    table(c9_md$L1),
    decreasing=TRUE
  )
)

cat("\nPatient support:\n")
print(
  sort(
    table(c9_md$patient_id),
    decreasing=TRUE
  )
)

c9_genes <- gene_ev %>%
  filter(
    cluster=="9",
    gene %in% c(
      "PTPRC","LST1",
      "CD3D","CD3E","TRAC","LYZ",
      "KRT5","KRT14","KRT15",
      "KRT17","KRT19",
      "COL17A1","EPCAM","KRTDAP"
    )
  ) %>%
  arrange(desc(mean_expression))

cat("\nKey C9 lineage genes:\n")

print(
  as.data.frame(c9_genes),
  row.names=FALSE
)

# ==========================================================================
# 9B. C12 NEURAL / SCHWANN RETURN AUDIT
# ==========================================================================

cat("\n============================================================\n")
cat("CLUSTER 12 — NEURAL/SCHWANN RETURN AUDIT\n")
cat("============================================================\n")

c12_cells <- colnames(x)[
  as.character(x$immune_cluster_06)=="12"
]

c12_md <- x@meta.data[c12_cells,,drop=FALSE]

cat("Cells:",length(c12_cells),"\n")

cat("\nFrozen source L1:\n")
print(
  sort(
    table(c12_md$L1),
    decreasing=TRUE
  )
)

cat("\nPatient support:\n")
print(
  sort(
    table(c12_md$patient_id),
    decreasing=TRUE
  )
)

c12_genes <- gene_ev %>%
  filter(
    cluster=="12",
    gene %in% c(
      "PTPRC","LST1",
      "CD3D","CD3E","TRAC","LYZ",
      "PMP2","MPZ","PLP1",
      "SOX10","S100B","CDH19"
    )
  ) %>%
  arrange(desc(mean_expression))

cat("\nKey C12 lineage genes:\n")

print(
  as.data.frame(c12_genes),
  row.names=FALSE
)

# ==========================================================================
# 9C. C9/C12 QC / DOUBLET COMPARISON
# ==========================================================================

doublet_summary <- md %>%
  mutate(
    audit_group=case_when(
      cluster=="9" ~ "C9_epithelial_candidate",
      cluster=="12" ~ "C12_Schwann_candidate",
      TRUE ~ "other_immune_universe"
    )
  ) %>%
  group_by(audit_group) %>%
  summarise(
    n_cells=n(),

    median_doublet_score=
      median(doublet_score,na.rm=TRUE),

    q25_doublet_score=
      as.numeric(
        quantile(
          doublet_score,
          0.25,
          na.rm=TRUE
        )
      ),

    q75_doublet_score=
      as.numeric(
        quantile(
          doublet_score,
          0.75,
          na.rm=TRUE
        )
      ),

    median_nCount_RNA=
      median(nCount_RNA,na.rm=TRUE),

    median_nFeature_RNA=
      median(nFeature_RNA,na.rm=TRUE),

    median_percent_mt=
      median(percent.mt,na.rm=TRUE),

    .groups="drop"
  )

cat("\n=== C9/C12 QC / DOUBLET COMPARISON ===\n")

print(
  as.data.frame(doublet_summary),
  row.names=FALSE
)

write_csv(
  doublet_summary,
  file.path(
    TDIR,
    "08E_C9_C12_doublet_QC_summary.csv"
  )
)

# ==========================================================================
# 9D. COMPLETE CLUSTER ADJUDICATION REPORT
# ==========================================================================

cat("\n============================================================\n")
cat("COMPLETE IMMUNE ADJUDICATION MATRIX\n")
cat("============================================================\n")

for(cl in cluster_levels){

  q <- qc %>% filter(cluster==cl)

  cat("\n------------------------------------------------------------\n")
  cat(
    "CLUSTER",cl,
    "| n =",q$n_cells,
    "| patients =",q$n_patients,
    "| max patient =",q$dominant_patient,
    sprintf("(%.1f%%)",100*q$dominant_patient_fraction),
    "\n"
  )

  ident <- ev %>%
    filter(
      cluster==cl,
      class=="identity"
    ) %>%
    arrange(desc(programme_z)) %>%
    select(
      programme,
      programme_z,
      mean_programme,
      pct_positive
    )

  nonimmune <- ev %>%
    filter(
      cluster==cl,
      class=="nonimmune"
    ) %>%
    arrange(desc(programme_z)) %>%
    select(
      programme,
      programme_z,
      mean_programme,
      pct_positive
    )

  states <- ev %>%
    filter(
      cluster==cl,
      class=="state"
    ) %>%
    arrange(desc(programme_z)) %>%
    select(
      programme,
      programme_z,
      mean_programme,
      pct_positive
    )

  cat("\nTOP IMMUNE PROGRAMMES:\n")
  print(
    as.data.frame(head(ident,8)),
    row.names=FALSE
  )

  cat("\nNON-IMMUNE PROGRAMMES:\n")
  print(
    as.data.frame(nonimmune),
    row.names=FALSE
  )

  cat("\nSTATE PROGRAMMES:\n")
  print(
    as.data.frame(states),
    row.names=FALSE
  )

  keygenes <- gene_ev %>%
    filter(cluster==cl) %>%
    arrange(desc(mean_expression)) %>%
    slice_head(n=20)

  cat("\nTOP PANEL GENES:\n")
  print(
    as.data.frame(keygenes),
    row.names=FALSE
  )

  dp <- disposition %>%
    filter(cluster==cl)

  cat(
    "\nPROVISIONAL DISPOSITION:",
    dp$provisional_disposition,
    "\n"
  )

  cat(
    "Leukocyte core:",
    sprintf(
      "z=%.2f | mean=%.3f | positive=%.1f%%",
      dp$leukocyte_core_z,
      dp$leukocyte_core_mean,
      dp$leukocyte_core_pct
    ),
    "\n"
  )

  cat(
    "Best immune:",
    dp$top_immune_programme,
    sprintf("(z=%.2f)",dp$top_immune_z),
    "\n"
  )

  cat(
    "Best non-immune:",
    dp$top_nonimmune_programme,
    sprintf("(z=%.2f)",dp$top_nonimmune_z),
    "\n"
  )

  cat(
    "Immune minus non-immune z:",
    sprintf("%.2f",dp$immune_minus_nonimmune_z),
    "\n"
  )
}

# ==========================================================================
# 10. FINAL SUMMARY
# ==========================================================================

cat("\nSTEP 10/10 — final evidence report\n")

disposition_print <- disposition %>%
  mutate(
    cluster_numeric=as.numeric(cluster)
  ) %>%
  arrange(cluster_numeric) %>%
  select(
    cluster,
    n_cells,
    n_patients,
    n_HC,
    n_AA,
    n_FFA,
    n_batches,
    dominant_patient,
    dominant_patient_fraction,
    leukocyte_core_z,
    leukocyte_core_mean,
    leukocyte_core_pct,
    top_immune_programme,
    top_immune_z,
    top_nonimmune_programme,
    top_nonimmune_z,
    immune_minus_nonimmune_z,
    provisional_disposition
  )

cat("\n============================================================\n")
cat("PROVISIONAL DISPOSITION SUMMARY\n")
cat("============================================================\n")

print(
  as.data.frame(disposition_print),
  row.names=FALSE
)

write_csv(
  disposition_print,
  file.path(
    TDIR,
    "08E_final_evidence_summary.csv"
  )
)

# ==========================================================================
# FORMAL INVARIANTS
# ==========================================================================

cat("\n=== FORMAL INVARIANTS ===\n")

checks <- data.frame(
  check=c(
    "cells_4762",
    "genes_23412",
    "clusters_17",
    "cluster_counts_sum_4762",
    "all_clusters_present",
    "disposition_rows_17",
    "C9_present",
    "C12_present",
    "no_cells_modified"
  ),
  pass=c(
    ncol(x)==EXPECTED_CELLS,
    nrow(x)==EXPECTED_GENES,
    nrow(qc)==17,
    sum(qc$n_cells)==EXPECTED_CELLS,
    all(cluster_levels %in% qc$cluster),
    nrow(disposition)==17,
    "9" %in% qc$cluster,
    "12" %in% qc$cluster,
    ncol(x)==EXPECTED_CELLS
  ),
  stringsAsFactors=FALSE
)

print(
  checks,
  row.names=FALSE
)

write_csv(
  checks,
  file.path(
    TDIR,
    "08E_formal_checks.csv"
  )
)

stopifnot(all(checks$pass))

# ==========================================================================
# AUDIT TRAIL
# ==========================================================================

audit <- c(
  "PIPELINE_V2 — 08E immune identity/state/contamination audit",
  paste0("Generated: ",Sys.time()),
  "Input: immune_v2_08B_reclustered.rds",
  "Working resolution: 0.6",
  paste0("Cells: ",ncol(x)),
  paste0("Genes: ",nrow(x)),
  "Clusters: 17",
  "Identity, state and non-immune programmes evaluated separately.",
  "Leukocyte backbone explicitly audited.",
  "Epithelial evidence explicitly audited.",
  "Neural/Schwann evidence explicitly audited.",
  "Fibroblast and endothelial evidence explicitly audited.",
  "C9 and C12 evaluated at cluster and cell level.",
  "C9/C12 patient-level lineage support exported.",
  "QC and scDblFinder score retained as orthogonal evidence.",
  "Disease was not used to assign identity.",
  "Provisional disposition categories are evidence flags only.",
  "No cells reassigned.",
  "No cells removed.",
  "No immune identities frozen.",
  "No immune states frozen.",
  "No immune resolution frozen.",
  "No global L1 object modified."
)

writeLines(
  audit,
  file.path(
    TDIR,
    "08E_audit.txt"
  )
)

cat("\n============================================================\n")
cat("08E COMPLETE\n")
cat("============================================================\n")
cat("Cells:",ncol(x),"\n")
cat("Clusters audited:",nrow(qc),"\n")
cat("No cells reassigned.\n")
cat("No cells removed.\n")
cat("No immune identities frozen.\n")
cat("No immune states frozen.\n")
cat("No L1 object modified.\n")
cat("PASS — identity/state/contamination evidence generated.\n")
