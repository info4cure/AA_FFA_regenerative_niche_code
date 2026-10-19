# =============================================================================
# PIPELINE_V2 — 05B1_SoupX_per_sample.R
# Per-sample ambient RNA correction with SoupX.
#
# READ-ONLY inputs:
#   raw_feature_bc_matrix + filtered_feature_bc_matrix for each of 16 samples
#   global_atlas_v2_frozen.rds (only for exact singlet barcode selection)
#
# OUTPUT:
#   corrected counts for EXACT V2 singlets, one RDS per sample
#   contamination/QC summaries
#
# This script DOES NOT modify, recluster, or overwrite the frozen atlas.
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(SoupX)
})

options(stringsAsFactors=FALSE)
set.seed(20260924)
source("PIPELINE_V2/config/paths.R")

cat("\n============================================================\n")
cat("PIPELINE_V2 — 05B1 SoupX PER SAMPLE\n")
cat("============================================================\n\n")

ATLAS <- file.path(V2_FROZEN,"global_atlas_v2_frozen.rds")
MANIFEST <- file.path(V2_FROZEN,"scrnaseq_provenance_manifest.csv")
OUT_ROOT <- file.path(V2_ROOT,"objects","ambient_corrected","per_sample")
QC_ROOT <- file.path(V2_ROOT,"results","QC","ambient_RNA")
dir.create(OUT_ROOT,recursive=TRUE,showWarnings=FALSE)
dir.create(QC_ROOT,recursive=TRUE,showWarnings=FALSE)

if(!file.exists(ATLAS)) stop("Missing frozen atlas: ",ATLAS)
if(!file.exists(MANIFEST)) stop("Missing provenance manifest: ",MANIFEST)

atlas <- readRDS(ATLAS)
m <- read.csv(MANIFEST,stringsAsFactors=FALSE)
stopifnot(ncol(atlas)==82893,nrow(m)==16)

cat("SoupX version:",as.character(packageVersion("SoupX")),"\n")
cat("Atlas cells:",ncol(atlas),"\nSamples:",nrow(m),"\n\n")

summary_list <- list()

for(i in seq_len(nrow(m))){
  src <- m$source_sample_id[i]
  patient <- m$patient_id[i]
  batch <- m$batch[i]

  cat(sprintf("[%02d/16] %-7s %-12s\n",i,patient,batch))

  base <- if(batch=="historical"){
    file.path("01_RAW_DATA/scrnaseq/01_iScience_GSE314972/cellranger_links",src)
  } else {
    file.path("01_RAW_DATA/scrnaseq/02_CABIMER_2026/cellranger",src)
  }

  raw_dir <- file.path(base,"raw_feature_bc_matrix")
  filt_dir <- file.path(base,"filtered_feature_bc_matrix")
  if(!dir.exists(raw_dir)||!dir.exists(filt_dir)) stop("Missing raw/filtered matrix for ",patient)

  raw <- Read10X(raw_dir,gene.column=2,unique.features=TRUE)
  filt <- Read10X(filt_dir,gene.column=2,unique.features=TRUE)
  if(is.list(raw)) raw <- raw[["Gene Expression"]]
  if(is.list(filt)) filt <- filt[["Gene Expression"]]

  common <- intersect(rownames(raw),rownames(filt))
  raw <- raw[common,,drop=FALSE]
  filt <- filt[common,,drop=FALSE]

  # Exact V2 singlet barcodes for this patient.
  atlas_full <- colnames(atlas)[atlas$patient_id==patient]
  atlas_bc <- sub(paste0("^",patient,"_"),"",atlas_full)

  if(anyDuplicated(atlas_bc)) stop("Duplicated stripped atlas barcodes: ",patient)
  missing <- setdiff(atlas_bc,colnames(filt))
  if(length(missing)>0) stop(patient,": ",length(missing)," atlas singlets missing from filtered matrix")

  # SoupX channel: raw droplets + Cell Ranger filtered cells.
  sc <- SoupChannel(tod=raw,toc=filt,calcSoupProfile=TRUE)

  # SoupX autoEstCont needs clusters. Generate sample-internal clusters solely
  # for contamination estimation; these NEVER replace atlas clusters.
  tmp <- CreateSeuratObject(counts=filt,min.cells=0,min.features=0)
  tmp <- NormalizeData(tmp,verbose=FALSE)
  tmp <- FindVariableFeatures(tmp,nfeatures=min(2000,nrow(tmp)),verbose=FALSE)
  tmp <- ScaleData(tmp,features=VariableFeatures(tmp),verbose=FALSE)
  tmp <- RunPCA(tmp,features=VariableFeatures(tmp),npcs=30,seed.use=20260924,verbose=FALSE)
  tmp <- FindNeighbors(tmp,dims=1:20,verbose=FALSE)
  tmp <- FindClusters(tmp,resolution=0.5,random.seed=20260924,verbose=FALSE)

  cl <- setNames(as.character(Idents(tmp)),colnames(tmp))
  sc <- setClusters(sc,cl)

  # Estimate contamination and correct counts.
  sc <- autoEstCont(sc,doPlot=FALSE)
  rho <- sc$metaData$rho
  rho_sample <- unique(round(rho,12))
  if(length(rho_sample)!=1) rho_sample <- median(rho,na.rm=TRUE)

  corrected <- adjustCounts(sc,roundToInt=TRUE)

  # Keep EXACT V2 singlets, preserving atlas barcode order.
  corrected <- corrected[,atlas_bc,drop=FALSE]
  colnames(corrected) <- atlas_full

  if(ncol(corrected)!=length(atlas_full)) stop("Corrected singlet count mismatch: ",patient)

  # Save full corrected feature universe for this sample.
  out_file <- file.path(OUT_ROOT,paste0(patient,"_SoupX_corrected_singlets.rds"))
  saveRDS(corrected,out_file,compress=FALSE)

  pre <- filt[,atlas_bc,drop=FALSE]
  pre_umi <- sum(pre)
  post_umi <- sum(corrected)

  summary_list[[patient]] <- data.frame(
    source_sample_id=src,
    patient_id=patient,
    disease=unique(atlas$disease[atlas$patient_id==patient]),
    batch=batch,
    raw_droplets=ncol(raw),
    filtered_cells=ncol(filt),
    atlas_singlets=length(atlas_bc),
    genes=nrow(corrected),
    rho=as.numeric(rho_sample[1]),
    singlet_UMI_pre=pre_umi,
    singlet_UMI_post=post_umi,
    UMI_removed=pre_umi-post_umi,
    pct_UMI_removed=100*(pre_umi-post_umi)/pre_umi,
    corrected_file=out_file,
    stringsAsFactors=FALSE
  )

  cat(sprintf(
    "   raw=%d filtered=%d singlets=%d rho=%.4f UMI removed=%.2f%%\n",
    ncol(raw),ncol(filt),length(atlas_bc),as.numeric(rho_sample[1]),
    100*(pre_umi-post_umi)/pre_umi
  ))

  rm(raw,filt,sc,corrected,tmp,pre,cl)
  gc(verbose=FALSE)
}

tab <- do.call(rbind,summary_list)
rownames(tab) <- NULL
write.csv(tab,file.path(QC_ROOT,"05B1_SoupX_per_sample_summary.csv"),row.names=FALSE)

stopifnot(sum(tab$atlas_singlets)==82893)
if(any(!is.finite(tab$rho))) stop("Non-finite SoupX contamination estimate detected.")

writeLines(c(
  "PIPELINE_V2 — 05B1 SoupX per-sample audit",
  paste0("Generated: ",format(Sys.time(),"%Y-%m-%d %H:%M:%S")),
  paste0("SoupX: ",packageVersion("SoupX")),
  "Correction performed independently for 16 samples.",
  paste0("Total exact V2 singlets retained: ",sum(tab$atlas_singlets)),
  paste0("Median rho: ",format(median(tab$rho),digits=4)),
  paste0("Range rho: ",format(min(tab$rho),digits=4)," - ",format(max(tab$rho),digits=4)),
  "Sample-internal clustering was used only by autoEstCont.",
  "Frozen atlas clustering/UMAP/identities were not modified.",
  "Corrected matrices retain exact V2 singlet cell IDs."
),file.path(QC_ROOT,"05B1_SoupX_audit.txt"))

cat("\n============================================================\n")
cat("05B1 COMPLETE\n")
cat("============================================================\n")
print(tab[,c("patient_id","disease","batch","atlas_singlets","rho","pct_UMI_removed")],row.names=FALSE)
cat("\nTotal singlets:",sum(tab$atlas_singlets),"\n")
cat("Median rho:",median(tab$rho),"\n")
cat("Output:",OUT_ROOT,"\n")
cat("\nPASS — frozen atlas unchanged.\n")
