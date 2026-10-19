suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(readr)
})

cat("\n============================================================\n")
cat("PIPELINE_V2 — 07D EPITHELIAL MARKER COMPARISON\n")
cat("============================================================\n\n")

IN <- "PIPELINE_V2/objects/annotated/epithelial_v2_07B_reclustered.rds"
OUT <- "PIPELINE_V2/results/annotation/07_EPITHELIAL_L2/07D_markers"
dir.create(OUT,recursive=TRUE,showWarnings=FALSE)

x <- readRDS(IN)
DefaultAssay(x) <- "RNA"

stopifnot(
  ncol(x)==63893,
  "epi_snn_res.0.4" %in% colnames(x@meta.data),
  "epi_snn_res.0.5" %in% colnames(x@meta.data)
)

# Seurat v5 layer sanity
x <- JoinLayers(x,assay="RNA")

cat("Cells:",ncol(x),"\n")
cat("Genes:",nrow(x),"\n\n")

# ============================================================
# 1. RESOLUTION 0.4 MARKERS
# ============================================================

cat("STEP 1/5 — markers resolution 0.4\n")

Idents(x) <- factor(x$epi_snn_res.0.4)

m04 <- FindAllMarkers(
  x,
  assay="RNA",
  only.pos=TRUE,
  min.pct=0.10,
  logfc.threshold=0.25,
  test.use="wilcox",
  verbose=FALSE
)

write_csv(
  m04,
  file.path(OUT,"07D_markers_res0.4_all.csv")
)

top04 <- m04 %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC),.by_group=TRUE) %>%
  slice_head(n=15) %>%
  ungroup()

write_csv(
  top04,
  file.path(OUT,"07D_markers_res0.4_top15.csv")
)

# ============================================================
# 2. RESOLUTION 0.5 MARKERS
# ============================================================

cat("STEP 2/5 — markers resolution 0.5\n")

Idents(x) <- factor(x$epi_snn_res.0.5)

m05 <- FindAllMarkers(
  x,
  assay="RNA",
  only.pos=TRUE,
  min.pct=0.10,
  logfc.threshold=0.25,
  test.use="wilcox",
  verbose=FALSE
)

write_csv(
  m05,
  file.path(OUT,"07D_markers_res0.5_all.csv")
)

top05 <- m05 %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC),.by_group=TRUE) %>%
  slice_head(n=15) %>%
  ungroup()

write_csv(
  top05,
  file.path(OUT,"07D_markers_res0.5_top15.csv")
)

# ============================================================
# 3. DIRECT TEST OF CLEAN SPLIT: C15 vs C18 at res 0.5
# ============================================================

cat("STEP 3/5 — direct C15 vs C18 split markers\n")

Idents(x) <- factor(x$epi_snn_res.0.5)

split_15_vs_18 <- FindMarkers(
  x,
  ident.1="15",
  ident.2="18",
  assay="RNA",
  min.pct=0.10,
  logfc.threshold=0.10,
  test.use="wilcox",
  verbose=FALSE
)

split_15_vs_18$gene <- rownames(split_15_vs_18)

split_15_vs_18 <- split_15_vs_18 %>%
  arrange(desc(abs(avg_log2FC)))

write_csv(
  split_15_vs_18,
  file.path(OUT,"07D_res0.5_C15_vs_C18.csv")
)

# ============================================================
# 4. SUSPICIOUS CLUSTER QC
# ============================================================

cat("STEP 4/5 — suspicious cluster audit\n")

qcvars <- intersect(
  c(
    "nCount_RNA",
    "nFeature_RNA",
    "percent.mt",
    "doublet_score"
  ),
  colnames(x@meta.data)
)

audit_clusters <- bind_rows(
  lapply(
    list(
      list(res="0.4",cluster="14"),
      list(res="0.4",cluster="16"),
      list(res="0.5",cluster="11"),
      list(res="0.5",cluster="14")
    ),
    function(z){

      cc <- paste0("epi_snn_res.",z$res)

      md <- x@meta.data %>%
        filter(
          as.character(.data[[cc]])==z$cluster
        )

      base <- data.frame(
        resolution=z$res,
        cluster=z$cluster,
        n_cells=nrow(md),
        n_patients=n_distinct(md$patient_id),
        dominant_patient=
          names(sort(table(md$patient_id),decreasing=TRUE))[1],
        dominant_patient_fraction=
          max(table(md$patient_id))/nrow(md),
        dominant_disease=
          names(sort(table(md$disease),decreasing=TRUE))[1],
        dominant_disease_fraction=
          max(table(md$disease))/nrow(md)
      )

      for(v in qcvars){
        base[[paste0(v,"_median")]] <-
          median(md[[v]],na.rm=TRUE)
      }

      base
    }
  )
)

write_csv(
  audit_clusters,
  file.path(OUT,"07D_suspicious_cluster_QC.csv")
)

cat("\n=== SUSPICIOUS CLUSTERS ===\n")
print(audit_clusters,row.names=FALSE)

# ============================================================
# 5. COMPACT MARKER REPORT
# ============================================================

cat("\nSTEP 5/5 — compact marker report\n")

cat("\n================ RESOLUTION 0.4 ================\n")

for(cl in sort(unique(as.integer(as.character(top04$cluster))))){

  genes <- top04 %>%
    filter(as.character(cluster)==as.character(cl)) %>%
    pull(gene)

  cat(
    sprintf(
      "C%-2s | %s\n",
      cl,
      paste(genes,collapse=", ")
    )
  )
}

cat("\n================ RESOLUTION 0.5 ================\n")

for(cl in sort(unique(as.integer(as.character(top05$cluster))))){

  genes <- top05 %>%
    filter(as.character(cluster)==as.character(cl)) %>%
    pull(gene)

  cat(
    sprintf(
      "C%-2s | %s\n",
      cl,
      paste(genes,collapse=", ")
    )
  )
}

cat("\n================ C15 vs C18 ================\n")

print(
  as.data.frame(
    split_15_vs_18 %>%
      select(
        gene,
        avg_log2FC,
        pct.1,
        pct.2,
        p_val_adj
      ) %>%
      slice_head(n=40)
  ),
  row.names=FALSE
)

audit <- c(
  "PIPELINE_V2 — 07D epithelial marker comparison",
  paste0("Generated: ",Sys.time()),
  "Cells: 63893",
  "Expression: SoupX-corrected RNA",
  "Compared candidate resolutions: 0.4 and 0.5",
  "Marker testing is exploratory annotation evidence",
  "Direct split test: res0.5 C15 vs C18",
  "Suspicious clusters audited: 0.4 C14/C16; 0.5 C11/C14",
  "No resolution selected",
  "No L2 labels assigned"
)

writeLines(
  audit,
  file.path(OUT,"07D_audit.txt")
)

cat("\n============================================================\n")
cat("07D COMPLETE\n")
cat("============================================================\n")
cat("No resolution selected.\n")
cat("No L2 labels assigned.\n")
cat("PASS — marker evidence generated.\n")
