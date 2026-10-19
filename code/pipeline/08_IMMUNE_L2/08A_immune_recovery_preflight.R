suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(tidyr)
  library(readr)
})

cat("\n============================================================\n")
cat("PIPELINE_V2 — 08A IMMUNE RECOVERY PREFLIGHT\n")
cat("============================================================\n\n")

IN <- "PIPELINE_V2/objects/annotated/global_atlas_v2_L1_frozen.rds"

OUT <- "PIPELINE_V2/results/annotation/08_IMMUNE_L2/08A_recovery_preflight"
TDIR <- file.path(OUT,"tables")
dir.create(TDIR, recursive=TRUE, showWarnings=FALSE)

x <- readRDS(IN)
DefaultAssay(x) <- "RNA"

stopifnot(
  ncol(x)==82893,
  nrow(x)==23412,
  "L1" %in% colnames(x@meta.data),
  length(unique(x$patient_id))==16
)

cat("Global cells:",ncol(x),"\n")
cat("Genes:",nrow(x),"\n")
cat("Patients:",length(unique(x$patient_id)),"\n")

# ============================================================
# 1. FROZEN L1 INVENTORY
# ============================================================

cat("\nSTEP 1/8 — frozen L1 inventory\n")

l1_counts <- x@meta.data %>%
  count(L1,name="n_cells") %>%
  mutate(pct=100*n_cells/sum(n_cells)) %>%
  arrange(desc(n_cells))

cat("\n=== FROZEN L1 COUNTS ===\n")
print(l1_counts,row.names=FALSE)

write_csv(
  l1_counts,
  file.path(TDIR,"08A_frozen_L1_counts.csv")
)

expected_immune <- c("T_NK","Myeloid_DC")

stopifnot(
  all(expected_immune %in% unique(as.character(x$L1)))
)

canonical_immune_cells <- colnames(x)[
  as.character(x$L1) %in% expected_immune
]

cat("\nCanonical broad immune cells:",
    length(canonical_immune_cells),"\n")

stopifnot(length(canonical_immune_cells)==4762)

# ============================================================
# 2. IMMUNE MARKER PROGRAMMES
#
# These are recovery/evidence programmes, not annotations.
# ============================================================

cat("\nSTEP 2/8 — immune marker coverage\n")

programmes <- list(

  Pan_leukocyte=c(
    "PTPRC","LST1","CD37","CORO1A","LCP1"
  ),

  Pan_T=c(
    "CD3D","CD3E","CD3G",
    "TRAC","TRBC1","TRBC2"
  ),

  CD4_T=c(
    "CD4","IL7R","LTB",
    "MAL","TCF7","LEF1"
  ),

  CD8_T=c(
    "CD8A","CD8B","CCL5",
    "GZMK","CTSW"
  ),

  NK=c(
    "NKG7","GNLY","KLRD1",
    "KLRF1","FCGR3A","PRF1"
  ),

  Cytotoxic=c(
    "NKG7","GNLY","PRF1",
    "GZMB","GZMH","CTSW"
  ),

  Treg=c(
    "FOXP3","IL2RA","CTLA4",
    "TIGIT","IKZF2"
  ),

  B_cell=c(
    "CD79A","CD79B","MS4A1",
    "CD37","CD74","HLA-DRA"
  ),

  Plasma_cell=c(
    "MZB1","JCHAIN","SDC1",
    "DERL3","TNFRSF17","IGHG1"
  ),

  Myeloid_core=c(
    "LYZ","LST1","FCER1G",
    "TYROBP","CTSS","AIF1"
  ),

  Monocyte=c(
    "FCN1","S100A8","S100A9",
    "CTSS","LYZ","VCAN"
  ),

  Macrophage=c(
    "C1QA","C1QB","C1QC",
    "APOE","CD163","MRC1"
  ),

  cDC2=c(
    "CD1C","CLEC10A","FCER1A",
    "CD1E","HLA-DQA1"
  ),

  cDC1=c(
    "CLEC9A","XCR1","CADM1",
    "BATF3","IRF8"
  ),

  DC_mature=c(
    "LAMP3","CCR7","FSCN1",
    "CD83","CCL19"
  ),

  pDC=c(
    "GZMB","GZMB","JCHAIN",
    "TCF4","IL3RA","CLEC4C"
  ),

  Mast=c(
    "TPSAB1","TPSB2","CPA3",
    "KIT","MS4A2","HDC"
  ),

  IFN_response=c(
    "ISG15","IFIT1","IFIT2",
    "IFIT3","MX1","OAS1"
  ),

  Cycling=c(
    "MKI67","TOP2A","UBE2C",
    "RRM2","CENPF"
  )
)

# remove duplicated genes within programmes
programmes <- lapply(programmes,unique)

coverage <- bind_rows(
  lapply(names(programmes),function(p){

    requested <- programmes[[p]]
    present <- intersect(requested,rownames(x))

    data.frame(
      programme=p,
      n_present=length(present),
      n_requested=length(requested),
      coverage=length(present)/length(requested),
      missing=paste(
        setdiff(requested,present),
        collapse=";"
      )
    )
  })
)

cat("\n=== IMMUNE MARKER COVERAGE ===\n")
print(coverage,row.names=FALSE)

write_csv(
  coverage,
  file.path(TDIR,"08A_immune_marker_coverage.csv")
)

# ============================================================
# 3. PROGRAMME SCORES ACROSS ALL 82,893 CELLS
#
# Critical: search globally, not only inside frozen immune L1.
# ============================================================

cat("\nSTEP 3/8 — global immune programme evidence\n")

genes <- intersect(
  unique(unlist(programmes)),
  rownames(x)
)

mat <- LayerData(
  x,
  assay="RNA",
  layer="data"
)[genes,,drop=FALSE]

score_list <- lapply(
  programmes,
  function(g){

    gg <- intersect(g,rownames(mat))

    if(length(gg)==0){
      return(rep(NA_real_,ncol(x)))
    }

    Matrix::colMeans(
      mat[gg,,drop=FALSE]
    )
  }
)

score_mat <- do.call(cbind,score_list)
rownames(score_mat) <- colnames(x)

# ============================================================
# 4. PROGRAMME EVIDENCE BY FROZEN L1
# ============================================================

cat("\nSTEP 4/8 — programme evidence by frozen L1\n")

md <- x@meta.data
md$cell <- rownames(md)

l1_programme_rows <- list()

for(l1 in unique(as.character(md$L1))){

  cells <- md$cell[as.character(md$L1)==l1]

  for(p in names(programmes)){

    v <- score_mat[cells,p]

    l1_programme_rows[[length(l1_programme_rows)+1]] <-
      data.frame(
        L1=l1,
        programme=p,
        n_cells=length(cells),
        mean_programme=mean(v,na.rm=TRUE),
        median_programme=median(v,na.rm=TRUE),
        pct_positive=100*mean(v>0,na.rm=TRUE)
      )
  }
}

l1_programme <- bind_rows(l1_programme_rows)

l1_programme <- l1_programme %>%
  group_by(programme) %>%
  mutate(
    z={
      s <- sd(mean_programme,na.rm=TRUE)
      if(is.finite(s) && s>0)
        as.numeric(scale(mean_programme))
      else
        rep(0,n())
    }
  ) %>%
  ungroup()

write_csv(
  l1_programme,
  file.path(TDIR,"08A_programme_evidence_by_L1.csv")
)

cat("\n=== TOP L1 FOR EACH IMMUNE PROGRAMME ===\n")

top_l1 <- l1_programme %>%
  group_by(programme) %>%
  arrange(desc(z),.by_group=TRUE) %>%
  slice_head(n=3) %>%
  ungroup()

print(
  as.data.frame(
    top_l1 %>%
      select(
        programme,L1,z,
        mean_programme,pct_positive
      )
  ),
  row.names=FALSE
)

write_csv(
  top_l1,
  file.path(TDIR,"08A_top_L1_by_programme.csv")
)

# ============================================================
# 5. GLOBAL CLUSTER IMMUNE EVIDENCE
# ============================================================

cat("\nSTEP 5/8 — global cluster immune evidence\n")

GC <- "RNA_snn_res.0.6"
stopifnot(GC %in% colnames(md))

md$global_cluster <- as.character(md[[GC]])

cluster_rows <- list()

for(cl in sort(unique(md$global_cluster))){

  cells <- md$cell[md$global_cluster==cl]

  for(p in names(programmes)){

    v <- score_mat[cells,p]

    cluster_rows[[length(cluster_rows)+1]] <-
      data.frame(
        cluster=cl,
        programme=p,
        n_cells=length(cells),
        frozen_L1=paste(
          sort(unique(as.character(md$L1[md$global_cluster==cl]))),
          collapse=";"
        ),
        mean_programme=mean(v,na.rm=TRUE),
        pct_positive=100*mean(v>0,na.rm=TRUE)
      )
  }
}

cluster_ev <- bind_rows(cluster_rows) %>%
  group_by(programme) %>%
  mutate(
    z={
      s <- sd(mean_programme,na.rm=TRUE)
      if(is.finite(s) && s>0)
        as.numeric(scale(mean_programme))
      else
        rep(0,n())
    }
  ) %>%
  ungroup()

write_csv(
  cluster_ev,
  file.path(TDIR,"08A_global_cluster_immune_evidence.csv")
)

# ============================================================
# 6. TARGETED SEARCH FOR B / PLASMA / MAST OUTSIDE
#    CURRENT IMMUNE L1
#
# Evidence only. Thresholds are deliberately descriptive.
# ============================================================

cat("\nSTEP 6/8 — searching for overlooked immune populations\n")

target_programmes <- c(
  "B_cell",
  "Plasma_cell",
  "Mast",
  "cDC1",
  "cDC2",
  "pDC"
)

candidate_rows <- list()

for(p in target_programmes){

  z <- cluster_ev %>%
    filter(programme==p) %>%
    arrange(desc(z))

  candidate_rows[[length(candidate_rows)+1]] <-
    z %>%
    mutate(
      currently_broad_immune =
        grepl(
          "T_NK|Myeloid_DC",
          frozen_L1
        )
    )
}

candidates <- bind_rows(candidate_rows)

write_csv(
  candidates,
  file.path(TDIR,"08A_target_population_cluster_scan.csv")
)

cat("\n=== TARGET POPULATION CLUSTER SCAN — TOP 5 EACH ===\n")

candidate_top <- candidates %>%
  group_by(programme) %>%
  arrange(desc(z),.by_group=TRUE) %>%
  slice_head(n=5) %>%
  ungroup()

print(
  as.data.frame(
    candidate_top %>%
      select(
        programme,
        cluster,
        frozen_L1,
        z,
        mean_programme,
        pct_positive,
        currently_broad_immune
      )
  ),
  row.names=FALSE
)

# ============================================================
# 7. PATIENT SUPPORT OF CURRENT BROAD IMMUNE COMPARTMENT
# ============================================================

cat("\nSTEP 7/8 — current broad immune patient support\n")

immune_md <- md[
  md$cell %in% canonical_immune_cells,
  ,
  drop=FALSE
]

immune_by_patient <- immune_md %>%
  count(
    L1,
    patient_id,
    disease,
    batch,
    name="n_cells"
  ) %>%
  arrange(L1,disease,desc(n_cells))

write_csv(
  immune_by_patient,
  file.path(TDIR,"08A_current_immune_by_patient.csv")
)

immune_summary <- immune_md %>%
  count(
    patient_id,
    disease,
    batch,
    name="n_cells"
  ) %>%
  summarise(
    total_cells=sum(n_cells),
    n_patients=n_distinct(patient_id),
    n_HC=n_distinct(patient_id[disease=="HC"]),
    n_AA=n_distinct(patient_id[disease=="AA"]),
    n_FFA=n_distinct(patient_id[disease=="FFA"]),
    n_batches=n_distinct(batch),
    min_cells_patient=min(n_cells),
    median_cells_patient=median(n_cells),
    max_cells_patient=max(n_cells)
  )

cat("\n=== CURRENT BROAD IMMUNE SUPPORT ===\n")
print(immune_summary,row.names=FALSE)

cat("\n=== CURRENT BROAD IMMUNE CELLS BY PATIENT ===\n")

immune_total_patient <- immune_md %>%
  count(
    patient_id,
    disease,
    batch,
    name="n_cells"
  ) %>%
  arrange(disease,desc(n_cells))

print(
  immune_total_patient,
  row.names=FALSE
)

write_csv(
  immune_total_patient,
  file.path(TDIR,"08A_current_immune_total_by_patient.csv")
)

# ============================================================
# 8. AUDIT
# ============================================================

cat("\nSTEP 8/8 — audit\n")

stopifnot(
  length(canonical_immune_cells)==4762,
  sum(immune_total_patient$n_cells)==4762
)

audit <- c(
  "PIPELINE_V2 — 08A immune recovery preflight",
  paste0("Generated: ",Sys.time()),
  "Input: global_atlas_v2_L1_frozen.rds",
  "Global cells: 82893",
  "Genes: 23412",
  "Current broad immune definition: frozen L1 T_NK + Myeloid_DC",
  paste0(
    "Current broad immune cells: ",
    length(canonical_immune_cells)
  ),
  "Global immune programme scan performed across all 82893 cells",
  "Targeted recovery audit includes T, NK, B, plasma, myeloid, monocyte, macrophage, DC, pDC and mast programmes",
  "Programme scores are evidence only and are not cell labels",
  "No L1 identities modified",
  "No cells removed",
  "No reclustering performed",
  "No object written"
)

writeLines(
  audit,
  file.path(TDIR,"08A_audit.txt")
)

cat("\n============================================================\n")
cat("08A COMPLETE\n")
cat("============================================================\n")
cat("Current broad immune cells:",length(canonical_immune_cells),"\n")
cat("No cells reassigned.\n")
cat("No immune subatlas created.\n")
cat("PASS — recovery evidence generated; immune universe remains unfrozen.\n")
