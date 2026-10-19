# =============================================================================
# PIPELINE_V2 — 06A_cluster_markers.R
# Unbiased evidence for L1 annotation of the 27 global clusters.
# No L1 assignment; no FibConsensus/BJD signatures; descriptive markers only.
# =============================================================================
suppressPackageStartupMessages({library(Seurat);library(ggplot2)})
options(stringsAsFactors=FALSE); set.seed(20260924)
source("PIPELINE_V2/config/paths.R")

cat("\n============================================================\nPIPELINE_V2 — 06A GLOBAL CLUSTER MARKERS\n============================================================\n\n")

IN_RDS<-file.path(V2_FROZEN,"global_atlas_v2_frozen.rds")
OUT_ROOT<-file.path(V2_ROOT,"results","annotation","06_GLOBAL_L1")
MARKER_DIR<-file.path(OUT_ROOT,"markers");EVIDENCE_DIR<-file.path(OUT_ROOT,"evidence");FIG_DIR<-file.path(OUT_ROOT,"figures")
dir.create(MARKER_DIR,recursive=TRUE,showWarnings=FALSE);dir.create(EVIDENCE_DIR,recursive=TRUE,showWarnings=FALSE);dir.create(FIG_DIR,recursive=TRUE,showWarnings=FALSE)

ALL_MARKERS_FILE<-file.path(MARKER_DIR,"06A_all_positive_markers.csv")
TOP20_FILE<-file.path(MARKER_DIR,"06A_top20_markers.csv")
TOP50_FILE<-file.path(MARKER_DIR,"06A_top50_markers.csv")
PANEL_FILE<-file.path(EVIDENCE_DIR,"06A_L1_canonical_marker_panel.csv")
PANEL_COVERAGE_FILE<-file.path(EVIDENCE_DIR,"06A_L1_marker_coverage.csv")
PATIENT_COMP_FILE<-file.path(EVIDENCE_DIR,"06A_cluster_patient_composition.csv")
DISEASE_COMP_FILE<-file.path(EVIDENCE_DIR,"06A_cluster_disease_composition.csv")
BATCH_COMP_FILE<-file.path(EVIDENCE_DIR,"06A_cluster_batch_composition.csv")
CLUSTER_SIZE_FILE<-file.path(EVIDENCE_DIR,"06A_cluster_sizes.csv")
AUDIT_FILE<-file.path(EVIDENCE_DIR,"06A_audit.txt")

if(!file.exists(IN_RDS))stop("Frozen global atlas not found: ",IN_RDS)
x<-readRDS(IN_RDS);DefaultAssay(x)<-"RNA";CLUSTER_COL<-"RNA_snn_res.0.6"
if(!CLUSTER_COL%in%colnames(x@meta.data))stop("Missing clustering column: ",CLUSTER_COL)
stopifnot(ncol(x)==82893,nrow(x)==23412,length(unique(x$patient_id))==16)
x$global_cluster<-factor(as.character(x@meta.data[[CLUSTER_COL]]),levels=as.character(0:26));Idents(x)<-"global_cluster";stopifnot(length(levels(Idents(x)))==27)
cat("Cells:",ncol(x),"\nGenes:",nrow(x),"\nPatients:",length(unique(x$patient_id)),"\nClusters:",length(levels(Idents(x))),"\n\n")

rna_counts<-LayerData(x,assay="RNA",layer="counts");rna_data<-LayerData(x,assay="RNA",layer="data")
stopifnot(nrow(rna_counts)==23412,nrow(rna_data)==23412,identical(rownames(rna_counts),rownames(rna_data)));rm(rna_counts,rna_data);gc(verbose=FALSE)

cluster_sizes<-as.data.frame(table(global_cluster=x$global_cluster),stringsAsFactors=FALSE);names(cluster_sizes)[2]<-"n_cells";cluster_sizes$global_cluster<-as.integer(as.character(cluster_sizes$global_cluster));cluster_sizes<-cluster_sizes[order(cluster_sizes$global_cluster),];write.csv(cluster_sizes,CLUSTER_SIZE_FILE,row.names=FALSE)

marker_panel<-list(
Keratinocyte=c("KRT5","KRT14","KRT15","KRT17","KRT1","KRT10"),
Fibroblast=c("COL1A1","COL1A2","COL3A1","DCN","LUM","PDGFRA"),
Vascular_endothelial=c("PECAM1","VWF","KDR","CDH5","EMCN","ENG"),
Lymphatic_endothelial=c("PROX1","PDPN","FLT4","CCL21","LYVE1"),
Pericyte=c("RGS5","CSPG4","MCAM","NOTCH3","PDGFRB","ABCC9"),
Smooth_muscle=c("ACTA2","TAGLN","MYH11","CNN1","MYL9"),
Lymphocyte=c("PTPRC","CD3D","CD3E","TRBC1","TRBC2","NKG7","GNLY"),
Myeloid_DC=c("PTPRC","LST1","TYROBP","FCER1G","HLA-DRA","CD74","CTSS"),
Melanocyte=c("PMEL","MLANA","TYR","DCT","MITF"),
Sebaceous_gland=c("PPARG","PLIN2","KRT79","MGST1","FASN"),
Sweat_gland=c("AQP5","DCD","PIP","KRT8","KRT18","KRT19"),
Neural_Schwann=c("S100B","SOX10","PLP1","MPZ","S100A1")
)
panel_long<-do.call(rbind,lapply(names(marker_panel),function(lineage)data.frame(L1_candidate=lineage,gene=marker_panel[[lineage]],stringsAsFactors=FALSE)));rownames(panel_long)<-NULL;panel_long$present<-panel_long$gene%in%rownames(x);write.csv(panel_long,PANEL_FILE,row.names=FALSE)
coverage<-do.call(rbind,lapply(split(panel_long,panel_long$L1_candidate),function(z)data.frame(L1_candidate=z$L1_candidate[1],n_present=sum(z$present),n_total=nrow(z),coverage=mean(z$present))));rownames(coverage)<-NULL;write.csv(coverage,PANEL_COVERAGE_FILE,row.names=FALSE);cat("Canonical marker coverage:\n");print(coverage,row.names=FALSE)

cat("\nFinding positive markers across 27 clusters...\n")
markers<-FindAllMarkers(x,assay="RNA",only.pos=TRUE,min.pct=0.10,logfc.threshold=0.25,test.use="wilcox",verbose=TRUE)
markers$cluster<-as.integer(as.character(markers$cluster));markers<-markers[order(markers$cluster,-markers$avg_log2FC),];write.csv(markers,ALL_MARKERS_FILE,row.names=FALSE)
get_top_n<-function(df,n){out<-do.call(rbind,lapply(split(df,df$cluster),function(z){z<-z[order(-z$avg_log2FC,z$p_val_adj),,drop=FALSE];head(z,n)}));rownames(out)<-NULL;out}
top20<-get_top_n(markers,20);top50<-get_top_n(markers,50);write.csv(top20,TOP20_FILE,row.names=FALSE);write.csv(top50,TOP50_FILE,row.names=FALSE)

md<-x@meta.data;md$global_cluster<-factor(as.character(md[[CLUSTER_COL]]),levels=as.character(0:26))
patient_comp<-as.data.frame(table(global_cluster=md$global_cluster,patient_id=md$patient_id),stringsAsFactors=FALSE);names(patient_comp)[3]<-"n_cells";pt<-aggregate(n_cells~patient_id,patient_comp,sum);names(pt)[2]<-"patient_total_cells";patient_comp<-merge(patient_comp,pt,by="patient_id",all.x=TRUE);patient_comp$fraction_of_patient<-patient_comp$n_cells/patient_comp$patient_total_cells;patient_comp$global_cluster<-as.integer(as.character(patient_comp$global_cluster));patient_comp<-patient_comp[order(patient_comp$global_cluster,patient_comp$patient_id),];write.csv(patient_comp,PATIENT_COMP_FILE,row.names=FALSE)
disease_comp<-as.data.frame(table(global_cluster=md$global_cluster,disease=md$disease),stringsAsFactors=FALSE);names(disease_comp)[3]<-"n_cells";dt<-aggregate(n_cells~disease,disease_comp,sum);names(dt)[2]<-"disease_total_cells";disease_comp<-merge(disease_comp,dt,by="disease",all.x=TRUE);disease_comp$fraction_of_disease<-disease_comp$n_cells/disease_comp$disease_total_cells;disease_comp$global_cluster<-as.integer(as.character(disease_comp$global_cluster));disease_comp<-disease_comp[order(disease_comp$global_cluster,disease_comp$disease),];write.csv(disease_comp,DISEASE_COMP_FILE,row.names=FALSE)
batch_comp<-as.data.frame(table(global_cluster=md$global_cluster,batch=md$batch),stringsAsFactors=FALSE);names(batch_comp)[3]<-"n_cells";bt<-aggregate(n_cells~batch,batch_comp,sum);names(bt)[2]<-"batch_total_cells";batch_comp<-merge(batch_comp,bt,by="batch",all.x=TRUE);batch_comp$fraction_of_batch<-batch_comp$n_cells/batch_comp$batch_total_cells;batch_comp$global_cluster<-as.integer(as.character(batch_comp$global_cluster));batch_comp<-batch_comp[order(batch_comp$global_cluster,batch_comp$batch),];write.csv(batch_comp,BATCH_COMP_FILE,row.names=FALSE)

p_umap<-DimPlot(x,reduction="umap_harmony",group.by="global_cluster",label=TRUE,repel=TRUE,label.size=4,raster=TRUE,pt.size=.25)+ggtitle("Global atlas - computational clusters")+theme_void(base_size=12)+theme(plot.title=element_text(face="bold",hjust=.5),legend.position="none")
ggsave(file.path(FIG_DIR,"06A_global_27clusters_UMAP.pdf"),p_umap,width=8,height=7,units="in",bg="white");ggsave(file.path(FIG_DIR,"06A_global_27clusters_UMAP.png"),p_umap,width=8,height=7,units="in",dpi=600,bg="white")
panel_present<-panel_long[panel_long$present,,drop=FALSE];dot_genes<-unique(panel_present$gene)
p_dot<-DotPlot(x,features=dot_genes,assay="RNA",group.by="global_cluster",dot.scale=5)+RotatedAxis()+labs(x=NULL,y="Global cluster",color="Average\nexpression",size="% expressed")+theme_classic(base_size=9)+theme(axis.text.x=element_text(angle=55,hjust=1,vjust=1,color="black",size=7),axis.text.y=element_text(color="black",size=8))
ggsave(file.path(FIG_DIR,"06A_L1_canonical_marker_DotPlot.pdf"),p_dot,width=17,height=8,units="in",bg="white");ggsave(file.path(FIG_DIR,"06A_L1_canonical_marker_DotPlot.png"),p_dot,width=17,height=8,units="in",dpi=600,bg="white")

audit<-c("PIPELINE_V2 — 06A global cluster marker audit",paste0("Generated: ",format(Sys.time(),"%Y-%m-%d %H:%M:%S")),paste0("Input: ",IN_RDS),paste0("Cells: ",ncol(x)),paste0("Genes: ",nrow(x)),paste0("Patients: ",length(unique(x$patient_id))),paste0("Clusters: ",length(levels(Idents(x)))),"","FindAllMarkers: wilcox; only.pos=TRUE; min.pct=0.10; logfc.threshold=0.25",paste0("Positive marker rows: ",nrow(markers)),paste0("Canonical markers requested: ",nrow(panel_long)),paste0("Canonical markers present: ",sum(panel_long$present)),"","Marker p-values are descriptive cell-level evidence only; not disease inference.","No L1 identity assigned.","No FibConsensus signatures used.","No follicular/BJD signatures used.")
writeLines(audit,AUDIT_FILE)
cat("\n============================================================\n06A COMPLETE\n============================================================\n")
cat("Positive marker rows:",nrow(markers),"\nTop20:",TOP20_FILE,"\nTop50:",TOP50_FILE,"\nEvidence:",EVIDENCE_DIR,"\nFigures:",FIG_DIR,"\n\nPASS\n")
