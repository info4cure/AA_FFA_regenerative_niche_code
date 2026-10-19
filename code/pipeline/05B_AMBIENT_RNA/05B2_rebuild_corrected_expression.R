# =============================================================================
# PIPELINE_V2 — 05B2_rebuild_corrected_expression.R
# Rebuild a SoupX-corrected EXPRESSION object while preserving the exact
# canonical atlas cells, metadata, cluster assignments and reductions.
#
# No PCA/Harmony/neighbors/clustering/UMAP is recalculated.
# =============================================================================
suppressPackageStartupMessages({library(Seurat);library(Matrix)})
options(stringsAsFactors=FALSE);set.seed(20260924)
source("PIPELINE_V2/config/paths.R")

cat("\n============================================================\n")
cat("PIPELINE_V2 — 05B2 REBUILD CORRECTED EXPRESSION\n")
cat("============================================================\n\n")

ATLAS_FILE<-file.path(V2_FROZEN,"global_atlas_v2_frozen.rds")
CORR_DIR<-file.path(V2_ROOT,"objects","ambient_corrected","per_sample")
OUT_DIR<-file.path(V2_ROOT,"objects","ambient_corrected")
QC_DIR<-file.path(V2_ROOT,"results","QC","ambient_RNA")
OUT_FILE<-file.path(OUT_DIR,"global_atlas_v2_SoupX_expression.rds")
dir.create(OUT_DIR,recursive=TRUE,showWarnings=FALSE);dir.create(QC_DIR,recursive=TRUE,showWarnings=FALSE)

atlas<-readRDS(ATLAS_FILE);DefaultAssay(atlas)<-"RNA"
stopifnot(ncol(atlas)==82893,nrow(atlas)==23412)
canonical_genes<-rownames(atlas);canonical_cells<-colnames(atlas)
patients<-unique(as.character(atlas$patient_id))
cat("Canonical cells:",length(canonical_cells),"\nCanonical genes:",length(canonical_genes),"\nPatients:",length(patients),"\n")

mats<-vector("list",length(patients));names(mats)<-patients;audit<-list()
for(i in seq_along(patients)){
  pt<-patients[i];f<-file.path(CORR_DIR,paste0(pt,"_SoupX_corrected_singlets.rds"))
  if(!file.exists(f))stop("Missing corrected matrix: ",f)
  z<-readRDS(f);expected_cells<-canonical_cells[atlas$patient_id==pt]
  if(!identical(colnames(z),expected_cells))stop(pt,": corrected cell IDs/order do not match canonical atlas.")
  missing_genes<-setdiff(canonical_genes,rownames(z))
  if(length(missing_genes))stop(pt,": missing ",length(missing_genes)," canonical genes.")
  z<-z[canonical_genes,,drop=FALSE]
  mats[[pt]]<-z
  audit[[pt]]<-data.frame(patient_id=pt,cells=ncol(z),genes=nrow(z),cell_order_exact=identical(colnames(z),expected_cells),canonical_genes_complete=identical(rownames(z),canonical_genes))
  cat(sprintf("[%02d/%02d] %-7s cells=%5d genes=%d\n",i,length(patients),pt,ncol(z),nrow(z)))
}
corr<-do.call(cbind,mats)
# Restore canonical global order explicitly.
corr<-corr[,canonical_cells,drop=FALSE]
stopifnot(identical(colnames(corr),canonical_cells),identical(rownames(corr),canonical_genes),ncol(corr)==82893,nrow(corr)==23412)

# Build corrected RNA assay object, then restore canonical metadata/reductions.
obj<-CreateSeuratObject(counts=corr,meta.data=atlas@meta.data,project="AA_FFA_HC_V2_SoupX",min.cells=0,min.features=0)
obj<-NormalizeData(obj,normalization.method="LogNormalize",scale.factor=10000,verbose=FALSE)

# Preserve reductions exactly; they are coordinates from canonical atlas,
# not recomputed from corrected expression.
for(red in Reductions(atlas)) obj[[red]]<-atlas[[red]]

# Preserve graph-independent cluster metadata and active identity.
Idents(obj)<-factor(as.character(atlas$RNA_snn_res.0.6),levels=as.character(0:26))
obj$expression_source<-"SoupX_corrected"
obj$geometry_source<-"global_atlas_v2_frozen"
obj$ambient_correction<-"SoupX_1.6.2_per_sample"

# Hard invariants.
counts2<-LayerData(obj,assay="RNA",layer="counts");data2<-LayerData(obj,assay="RNA",layer="data")
stopifnot(
 nrow(counts2)==23412,nrow(data2)==23412,ncol(obj)==82893,
 identical(rownames(counts2),canonical_genes),
 identical(rownames(data2),canonical_genes),
 identical(colnames(obj),canonical_cells),
 identical(obj@meta.data$patient_id,atlas@meta.data$patient_id),
 identical(as.character(obj@meta.data$RNA_snn_res.0.6),as.character(atlas@meta.data$RNA_snn_res.0.6))
)
for(red in Reductions(atlas)){
 stopifnot(identical(rownames(Embeddings(obj,red)),rownames(Embeddings(atlas,red))))
 stopifnot(isTRUE(all.equal(Embeddings(obj,red),Embeddings(atlas,red),tolerance=0)))
}

saveRDS(obj,OUT_FILE,compress=FALSE)
a<-do.call(rbind,audit);rownames(a)<-NULL;write.csv(a,file.path(QC_DIR,"05B2_rebuild_per_sample_audit.csv"),row.names=FALSE)
writeLines(c(
 "PIPELINE_V2 — 05B2 corrected expression audit",
 paste0("Generated: ",format(Sys.time(),"%Y-%m-%d %H:%M:%S")),
 paste0("Cells: ",ncol(obj)),
 paste0("Genes counts: ",nrow(counts2)),
 paste0("Genes LogNormalized data: ",nrow(data2)),
 "Cell IDs/order: EXACT canonical match",
 "Patient metadata: EXACT canonical match",
 "RNA_snn_res.0.6: EXACT canonical match",
 paste0("Reductions preserved: ",paste(Reductions(obj),collapse=", ")),
 "Reductions were copied, NOT recalculated.",
 "No HVG/PCA/Harmony/neighbors/clustering/UMAP was recalculated.",
 paste0("Output: ",OUT_FILE)
),file.path(QC_DIR,"05B2_rebuild_corrected_expression_audit.txt"))

cat("\n============================================================\n05B2 COMPLETE\n============================================================\n")
cat("Cells:",ncol(obj),"\nCounts genes:",nrow(counts2),"\nLogNormalized genes:",nrow(data2),"\n")
cat("Cell order exact:",identical(colnames(obj),canonical_cells),"\n")
cat("Cluster metadata exact:",identical(as.character(obj$RNA_snn_res.0.6),as.character(atlas$RNA_snn_res.0.6)),"\n")
cat("Output:\n ",OUT_FILE,"\n\nPASS\n")
