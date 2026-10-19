# PIPELINE_V2 — 05_integration.R
suppressPackageStartupMessages({library(Seurat);library(harmony);library(Matrix)})
options(stringsAsFactors=FALSE); source("PIPELINE_V2/config/paths.R")
cat("\n============================================================\nPIPELINE_V2 — 05 INTEGRATION\n============================================================\n\n")
IN_RDS<-file.path(V2_CHECKPOINTS,"04_singlet_objects.rds"); OUT_RDS<-file.path(V2_CHECKPOINTS,"05_integrated_atlas.rds"); PCA_RDS<-file.path(V2_CHECKPOINTS,"05_PCA_checkpoint.rds")
COMMON_GENES_FILE<-file.path(V2_TABLES,"05_common_genes.txt"); CLUSTER_TABLE<-file.path(V2_TABLES,"05_cluster_sizes.csv"); QC_DIR<-file.path(V2_QC,"integration")
dir.create(V2_CHECKPOINTS,recursive=TRUE,showWarnings=FALSE);dir.create(V2_TABLES,recursive=TRUE,showWarnings=FALSE);dir.create(QC_DIR,recursive=TRUE,showWarnings=FALSE)
OUT_AUDIT<-file.path(QC_DIR,"05_integration_audit.txt"); CLUSTER_COMPARISON_FILE<-file.path(QC_DIR,"05_historical_cluster_comparison.csv"); CLUSTER_CONTINGENCY_FILE<-file.path(QC_DIR,"05_cluster_contingency.csv")
PARAM_FILE<-file.path(V2_ROOT,"config","integration_parameters.csv"); RES_FILE<-file.path(V2_ROOT,"config","clustering_resolutions.csv")
HISTORICAL_ATLAS<-file.path(LEGACY_SCRNA,"03_integration","A03_AA_FFA_HC_de_novo_atlas.rds")
EXPECTED_CELLS<-82893L;EXPECTED_GENES<-23412L;EXPECTED_PATIENTS<-16L;EXPECTED_CLUSTERS<-27L
EXPECTED_DISEASE<-c(AA=36131L,FFA=24085L,HC=22677L);EXPECTED_BATCH<-c(CABIMER2026=48016L,historical=34877L)
EXPECTED_CLUSTER_SIZES<-c("0"=12372L,"1"=6459L,"2"=5794L,"3"=5311L,"4"=4808L,"5"=4771L,"6"=4304L,"7"=4194L,"8"=3831L,"9"=3308L,"10"=3247L,"11"=2858L,"12"=2607L,"13"=2356L,"14"=2121L,"15"=1989L,"16"=1719L,"17"=1588L,"18"=1556L,"19"=1454L,"20"=1388L,"21"=1273L,"22"=1230L,"23"=1192L,"24"=651L,"25"=330L,"26"=182L)
if(!file.exists(PARAM_FILE))stop("Missing parameter file: ",PARAM_FILE);if(!file.exists(RES_FILE))stop("Missing resolution file: ",RES_FILE)
params<-read.csv(PARAM_FILE,stringsAsFactors=FALSE);get_param<-function(name){x<-params$value[params$parameter==name];if(length(x)!=1)stop("Parameter must occur exactly once: ",name);x[[1]]}
NORMALIZATION_METHOD<-get_param("normalization_method");SCALE_FACTOR<-as.numeric(get_param("scale_factor"));VARIABLE_METHOD<-get_param("variable_feature_method");N_HVG<-as.integer(get_param("n_variable_features"));N_PCS<-as.integer(get_param("n_pcs"));HARMONY_DIMS<-as.integer(get_param("harmony_dims"));NEIGHBOR_DIMS<-as.integer(get_param("neighbors_dims"));UMAP_DIMS<-as.integer(get_param("umap_dims"));PRIMARY_RES<-as.numeric(get_param("primary_resolution"));SEED<-as.integer(get_param("seed"));HARMONY_COVARIATE<-get_param("harmony_covariate")
resolutions<-as.numeric(read.csv(RES_FILE,stringsAsFactors=FALSE)$resolution);set.seed(SEED)
packages<-c("Seurat","SeuratObject","harmony","uwot");versions<-vapply(packages,function(pkg)as.character(packageVersion(pkg)),character(1));cat("Software:\n");print(c(R=R.version.string,versions));cat("\nParameters:\n");print(params,row.names=FALSE);cat("\nClustering resolutions:\n");print(resolutions)
adjusted_rand_index<-function(x,y){if(length(x)!=length(y))stop("ARI vectors differ");tab<-table(x,y);n<-sum(tab);c2<-function(z)z*(z-1)/2;s<-sum(c2(tab));r<-sum(c2(rowSums(tab)));c<-sum(c2(colSums(tab)));t<-c2(n);e<-r*c/t;m<-.5*(r+c);d<-m-e;if(d==0)return(ifelse(s==e,1,0));(s-e)/d}
if(!file.exists(IN_RDS))stop("Missing V2 singlet checkpoint: ",IN_RDS);objs<-readRDS(IN_RDS);if(length(objs)!=EXPECTED_PATIENTS)stop("Expected 16 sample objects");total_cells<-sum(vapply(objs,ncol,numeric(1)));if(total_cells!=EXPECTED_CELLS)stop("Expected 82,893 singlets");cat("\nSamples:",length(objs),"\nCells:",total_cells,"\n")
cat("\nSTEP 1/8 — harmonizing count layers\n");for(nm in names(objs)){DefaultAssay(objs[[nm]])<-"RNA";objs[[nm]][["RNA"]]<-JoinLayers(objs[[nm]][["RNA"]])}
gene_sets<-lapply(objs,function(z)rownames(LayerData(z,assay="RNA",layer="counts")));common_genes<-Reduce(intersect,gene_sets);cat("Common genes across all 16 samples:",length(common_genes),"\n");if(length(common_genes)!=EXPECTED_GENES)stop("Expected 23,412 common genes");writeLines(common_genes,COMMON_GENES_FILE)
cat("STEP 2/8 — rebuilding common-gene objects\n");clean<-vector("list",length(objs));names(clean)<-names(objs)
for(nm in names(objs)){z<-objs[[nm]];counts<-LayerData(z,assay="RNA",layer="counts")[common_genes,,drop=FALSE];clean[[nm]]<-CreateSeuratObject(counts=counts,meta.data=z@meta.data,project=nm,min.cells=0,min.features=0);rm(z,counts);gc(verbose=FALSE)}
rm(objs,gene_sets);gc(verbose=FALSE)
cat("STEP 3/8 — merging 16 samples\n");joint<-merge(clean[[1]],y=clean[-1],project="AA_FFA_HC_V2");rm(clean);gc(verbose=FALSE);joint[["RNA"]]<-JoinLayers(joint[["RNA"]])
stopifnot(ncol(joint)==EXPECTED_CELLS,nrow(joint)==EXPECTED_GENES,length(unique(joint$patient_id))==EXPECTED_PATIENTS)
cat("Merged cells:",ncol(joint),"\nMerged genes:",nrow(joint),"\nPatients:",length(unique(joint$patient_id)),"\n")
dc<-table(joint$disease);bc<-table(joint$batch);cat("\nDisease:\n");print(dc);cat("\nBatch:\n");print(bc)
for(nm in names(EXPECTED_DISEASE))if(is.na(dc[nm])||as.integer(dc[nm])!=EXPECTED_DISEASE[[nm]])stop("Disease mismatch: ",nm)
for(nm in names(EXPECTED_BATCH))if(is.na(bc[nm])||as.integer(bc[nm])!=EXPECTED_BATCH[[nm]])stop("Batch mismatch: ",nm)
cat("\nSTEP 4/8 — normalization and HVGs\n");joint<-NormalizeData(joint,normalization.method=NORMALIZATION_METHOD,scale.factor=SCALE_FACTOR,verbose=FALSE);joint<-FindVariableFeatures(joint,selection.method=VARIABLE_METHOD,nfeatures=N_HVG,verbose=FALSE);if(length(VariableFeatures(joint))!=N_HVG)stop("HVG mismatch");joint<-ScaleData(joint,features=VariableFeatures(joint),verbose=FALSE)
cat("STEP 5/8 — PCA\n");joint<-RunPCA(joint,features=VariableFeatures(joint),npcs=N_PCS,seed.use=SEED,verbose=FALSE);if(ncol(Embeddings(joint,"pca"))!=N_PCS)stop("PCA dimension mismatch")
joint<-RunUMAP(joint,reduction="pca",dims=seq_len(UMAP_DIMS),reduction.name="umap_uncorrected",reduction.key="UNC_",seed.use=SEED,verbose=FALSE);saveRDS(joint,PCA_RDS,compress=FALSE)
cat("STEP 6/8 — Harmony batch correction\n");if(HARMONY_COVARIATE!="batch")stop("Historical reproduction requires batch only")
joint<-RunHarmony(joint,group.by.vars=HARMONY_COVARIATE,reduction.use="pca",dims.use=seq_len(HARMONY_DIMS),reduction.save="harmony",verbose=FALSE);if(ncol(Embeddings(joint,"harmony"))<HARMONY_DIMS)stop("Harmony dimension mismatch")
cat("STEP 7/8 — graph and clustering\n");joint<-FindNeighbors(joint,reduction="harmony",dims=seq_len(NEIGHBOR_DIMS),verbose=FALSE);for(res in resolutions)joint<-FindClusters(joint,resolution=res,random.seed=SEED,verbose=FALSE)
primary_col<-paste0("RNA_snn_res.",format(PRIMARY_RES,scientific=FALSE,trim=TRUE));if(!primary_col%in%colnames(joint@meta.data))stop("Missing primary cluster column: ",primary_col);Idents(joint)<-primary_col;n_clusters<-length(unique(as.character(Idents(joint))));cat("Primary clusters:",n_clusters,"\n");if(n_clusters!=EXPECTED_CLUSTERS)stop("Expected 27 clusters")
cat("STEP 8/8 — Harmony UMAP\n");joint<-RunUMAP(joint,reduction="harmony",dims=seq_len(UMAP_DIMS),reduction.name="umap_harmony",reduction.key="HAR_",seed.use=SEED,verbose=FALSE)
cluster_sizes<-as.data.frame(table(cluster=as.character(Idents(joint))),stringsAsFactors=FALSE);names(cluster_sizes)[2]<-"n_cells";cluster_sizes$cluster<-as.character(cluster_sizes$cluster);cluster_sizes<-cluster_sizes[order(as.integer(cluster_sizes$cluster)),,drop=FALSE];write.csv(cluster_sizes,CLUSTER_TABLE,row.names=FALSE);cat("\nV2 cluster sizes:\n");print(cluster_sizes,row.names=FALSE)
observed_sizes<-setNames(as.integer(cluster_sizes$n_cells),cluster_sizes$cluster);raw_size_vector_exact<-setequal(names(observed_sizes),names(EXPECTED_CLUSTER_SIZES))&&identical(as.integer(observed_sizes[names(EXPECTED_CLUSTER_SIZES)]),as.integer(EXPECTED_CLUSTER_SIZES));cat("\nHistorical raw cluster-size vector exact:",raw_size_vector_exact,"\n")
historical_available<-file.exists(HISTORICAL_ATLAS);cell_identity_exact<-NA;cluster_label_exact<-NA;cluster_partition_ari<-NA_real_;cluster_partition_exact<-NA
if(historical_available){cat("\nLoading historical A03 atlas for comparison...\n");hist<-readRDS(HISTORICAL_ATLAS);cell_identity_exact<-identical(sort(colnames(joint)),sort(colnames(hist)));cat("Cell identities exact:",cell_identity_exact,"\n");if(!cell_identity_exact)stop("V2 cells differ from historical A03")
shared_cells<-colnames(joint);hist_col<-"RNA_snn_res.0.6";if(!hist_col%in%colnames(hist@meta.data))stop("Historical primary clustering missing")
new_cluster<-as.character(joint@meta.data[shared_cells,primary_col]);old_cluster<-as.character(hist@meta.data[shared_cells,hist_col]);cluster_label_exact<-identical(new_cluster,old_cluster);cluster_partition_ari<-adjusted_rand_index(old_cluster,new_cluster);cluster_partition_exact<-isTRUE(all.equal(cluster_partition_ari,1,tolerance=1e-12))
comparison<-data.frame(cell_id=shared_cells,historical_cluster=old_cluster,v2_cluster=new_cluster,exact_label_match=old_cluster==new_cluster);write.csv(comparison,CLUSTER_COMPARISON_FILE,row.names=FALSE)
contingency<-table(historical=old_cluster,V2=new_cluster);write.csv(as.data.frame.matrix(contingency),CLUSTER_CONTINGENCY_FILE,row.names=TRUE)
cat("Raw numerical cluster labels exact:",cluster_label_exact,"\nHistorical vs V2 ARI:",sprintf("%.12f",cluster_partition_ari),"\nCluster partition exact (ARI=1):",cluster_partition_exact,"\n");rm(hist);gc(verbose=FALSE)
}else warning("Historical A03 not found; historical comparison skipped.")
if(historical_available&&!isTRUE(cluster_partition_exact))stop("Historical and V2 cluster partitions differ (ARI != 1). Do not continue.")

# =============================================================================
# FULL RNA EXPRESSION UNIVERSE INVARIANT
#
# HVGs define atlas geometry; they DO NOT define the biological gene universe.
# Full counts and LogNormalized expression are retained for pseudobulk,
# tape-strip RNA-seq mapping and proteomic projection.
# =============================================================================

rna_counts <- LayerData(
  joint,
  assay = "RNA",
  layer = "counts"
)

rna_data <- LayerData(
  joint,
  assay = "RNA",
  layer = "data"
)

rna_scaled <- LayerData(
  joint,
  assay = "RNA",
  layer = "scale.data"
)

hvg <- VariableFeatures(joint)

stopifnot(
  nrow(rna_counts) == EXPECTED_GENES,
  nrow(rna_data) == EXPECTED_GENES,
  nrow(rna_scaled) == N_HVG,
  length(hvg) == N_HVG,
  identical(
    rownames(rna_counts),
    rownames(rna_data)
  ),
  all(hvg %in% rownames(rna_counts)),
  setequal(
    rownames(rna_scaled),
    hvg
  )
)

cat("\nFULL RNA EXPRESSION UNIVERSE\n")
cat("Counts:        ", nrow(rna_counts), " genes\n")
cat("LogNormalized: ", nrow(rna_data), " genes\n")
cat("Scaled/HVG:    ", nrow(rna_scaled), " genes\n")
cat("Status: PASS\n")

rm(
  rna_counts,
  rna_data,
  rna_scaled,
  hvg
)

gc(verbose = FALSE)

saveRDS(joint,OUT_RDS,compress=FALSE)
audit<-c("PIPELINE_V2 — 05 integration audit",paste0("Generated: ",format(Sys.time(),"%Y-%m-%d %H:%M:%S")),"",paste0("R: ",R.version.string),paste0(names(versions),": ",versions),"",paste0("Seed: ",SEED),paste0("Cells: ",ncol(joint)),paste0("Genes: ",nrow(joint)),paste0("HVGs: ",length(VariableFeatures(joint))),paste0("PCA dimensions: ",ncol(Embeddings(joint,"pca"))),paste0("Harmony covariate: ",HARMONY_COVARIATE),paste0("Harmony dims used: 1:",HARMONY_DIMS),paste0("Primary resolution: ",PRIMARY_RES),paste0("Primary clusters: ",n_clusters),"",paste0("Historical atlas available: ",historical_available),paste0("Cell identities exact: ",cell_identity_exact),paste0("Raw cluster labels exact: ",cluster_label_exact),paste0("Raw cluster-size vector exact: ",raw_size_vector_exact),paste0("Historical vs V2 ARI: ",cluster_partition_ari),paste0("Partition exact (ARI=1): ",cluster_partition_exact),"",paste0("PCA checkpoint: ",PCA_RDS),paste0("Validated atlas: ",OUT_RDS),paste0("Common genes: ",COMMON_GENES_FILE),paste0("Cluster sizes: ",CLUSTER_TABLE))
writeLines(audit,OUT_AUDIT)
cat("\n============================================================\nV2 INTEGRATION COMPLETE\n============================================================\n")
cat("Cells:",ncol(joint),"\nGenes:",nrow(joint),"\nHVGs:",length(VariableFeatures(joint)),"\nClusters:",n_clusters,"\n")
if(historical_available){cat("Historical cell identities exact:",cell_identity_exact,"\n");cat("Historical raw labels exact:",cluster_label_exact,"\n");cat("Historical partition ARI:",sprintf("%.12f",cluster_partition_ari),"\n")}
cat("\nValidated atlas:\n  ",OUT_RDS,"\nAudit:\n  ",OUT_AUDIT,"\n")
cat("\nPASS\n")
