# =============================================================================
# PIPELINE_V2 — 06A4_KRT_ambient_contamination_audit.R
# Quantify keratin signal outside epithelial clusters without modifying atlas.
# =============================================================================
suppressPackageStartupMessages({library(Seurat);library(ggplot2)})
options(stringsAsFactors=FALSE);set.seed(20260924);source("PIPELINE_V2/config/paths.R")
cat("\n=== 06A4 KRT AMBIENT CONTAMINATION AUDIT ===\n")
IN<-file.path(V2_FROZEN,"global_atlas_v2_frozen.rds")
OUT<-file.path(V2_ROOT,"results","annotation","06_GLOBAL_L1","KRT_ambient_audit")
TD<-file.path(OUT,"tables");FD<-file.path(OUT,"figures");dir.create(TD,recursive=TRUE,showWarnings=FALSE);dir.create(FD,recursive=TRUE,showWarnings=FALSE)
x<-readRDS(IN);DefaultAssay(x)<-"RNA";CC<-"RNA_snn_res.0.6";x$global_cluster<-factor(as.character(x@meta.data[[CC]]),levels=as.character(0:26))
stopifnot(ncol(x)==82893,nrow(x)==23412)
KRT<-c("KRT5","KRT14","KRT15","KRT17","KRT1","KRT10")
KRT<-KRT[KRT%in%rownames(x)];stopifnot(length(KRT)>=5)
counts<-LayerData(x,assay="RNA",layer="counts");dat<-LayerData(x,assay="RNA",layer="data")
det<-t(as.matrix(counts[KRT,,drop=FALSE]>0));expr<-t(as.matrix(dat[KRT,,drop=FALSE]))
n_krt<-rowSums(det);krt_sum<-rowSums(expr);krt_mean<-rowMeans(expr)
cell<-data.frame(cell_id=colnames(x),cluster=as.character(x$global_cluster),patient_id=x$patient_id,disease=x$disease,batch=x$batch,n_KRT_detected=n_krt,KRT_sum=krt_sum,KRT_mean=krt_mean,nCount_RNA=x$nCount_RNA,nFeature_RNA=x$nFeature_RNA,doublet_score=x$doublet_score,stringsAsFactors=FALSE)
for(g in KRT){cell[[paste0(g,"_det")]]<-det[,g];cell[[paste0(g,"_expr")]]<-expr[,g]}
con<-gzfile(file.path(TD,"06A4_cell_KRT_metrics.csv.gz"),"wt");write.csv(cell,con,row.names=FALSE);close(con)

summ<-do.call(rbind,lapply(0:26,function(cl){z<-cell[cell$cluster==as.character(cl),];r<-data.frame(cluster=cl,n_cells=nrow(z),pct_KRT_ge1=100*mean(z$n_KRT_detected>=1),pct_KRT_ge2=100*mean(z$n_KRT_detected>=2),pct_KRT_ge3=100*mean(z$n_KRT_detected>=3),pct_KRT_ge4=100*mean(z$n_KRT_detected>=4),mean_KRT_program=mean(z$KRT_mean),median_KRT_program=median(z$KRT_mean),median_nFeature=median(z$nFeature_RNA),median_doublet_score=median(z$doublet_score,na.rm=TRUE));for(g in KRT)r[[paste0("pct_",g,"_positive")]]<-100*mean(z[[paste0(g,"_det")]]);r}))
write.csv(summ,file.path(TD,"06A4_cluster_KRT_summary.csv"),row.names=FALSE)

# Per patient × cluster: crucial to distinguish generalized ambient signal from sample-specific effects.
pc<-do.call(rbind,lapply(split(cell,list(cell$patient_id,cell$cluster),drop=TRUE),function(z)data.frame(patient_id=z$patient_id[1],disease=z$disease[1],batch=z$batch[1],cluster=as.integer(z$cluster[1]),n_cells=nrow(z),pct_KRT_ge1=100*mean(z$n_KRT_detected>=1),pct_KRT_ge3=100*mean(z$n_KRT_detected>=3),mean_KRT_program=mean(z$KRT_mean),median_doublet_score=median(z$doublet_score,na.rm=TRUE))))
write.csv(pc,file.path(TD,"06A4_patient_cluster_KRT_summary.csv"),row.names=FALSE)

# Reference non-epithelial clusters based on strong independent marker evidence.
non_epi<-c(4,9,11,18,19,21,25,26)
write.csv(summ[summ$cluster%in%non_epi,],file.path(TD,"06A4_non_epithelial_KRT_summary.csv"),row.names=FALSE)

# Relationship of multi-KRT signal to doublet score / complexity.
cell$KRT_category<-cut(cell$n_KRT_detected,breaks=c(-Inf,0,1,2,Inf),labels=c("0","1","2","3+"))
qc<-aggregate(cbind(nFeature_RNA,doublet_score)~KRT_category,data=cell,FUN=function(v)median(v,na.rm=TRUE))
write.csv(qc,file.path(TD,"06A4_KRT_burden_QC_relationship.csv"),row.names=FALSE)

# Heatmap: % cells positive for each KRT by cluster.
long<-do.call(rbind,lapply(KRT,function(g)data.frame(cluster=summ$cluster,gene=g,pct_positive=summ[[paste0("pct_",g,"_positive")]])))
long$cluster<-factor(long$cluster,levels=rev(0:26))
p1<-ggplot(long,aes(gene,cluster,fill=pct_positive))+geom_tile()+scale_fill_gradient(low="white",high="black")+labs(x=NULL,y="Global cluster",fill="% positive",title="Keratin detection across global clusters")+theme_classic(base_size=10)
ggsave(file.path(FD,"06A4_KRT_detection_heatmap.pdf"),p1,width=7,height=8,bg="white");ggsave(file.path(FD,"06A4_KRT_detection_heatmap.png"),p1,width=7,height=8,dpi=600,bg="white")

# Burden plot: increasingly stringent multi-KRT detection.
bur<-rbind(data.frame(cluster=summ$cluster,criterion=">=1 KRT",pct=summ$pct_KRT_ge1),data.frame(cluster=summ$cluster,criterion=">=2 KRT",pct=summ$pct_KRT_ge2),data.frame(cluster=summ$cluster,criterion=">=3 KRT",pct=summ$pct_KRT_ge3),data.frame(cluster=summ$cluster,criterion=">=4 KRT",pct=summ$pct_KRT_ge4))
bur$cluster<-factor(bur$cluster,levels=0:26);bur$criterion<-factor(bur$criterion,levels=c(">=1 KRT",">=2 KRT",">=3 KRT",">=4 KRT"))
p2<-ggplot(bur,aes(cluster,pct,group=criterion,linetype=criterion))+geom_line(linewidth=.45)+geom_point(size=1)+labs(x="Global cluster",y="% cells",linetype=NULL,title="Increasingly stringent multi-keratin detection")+theme_classic(base_size=9)
ggsave(file.path(FD,"06A4_multiKRT_burden.pdf"),p2,width=10,height=4.5,bg="white");ggsave(file.path(FD,"06A4_multiKRT_burden.png"),p2,width=10,height=4.5,dpi=600,bg="white")

# Non-epithelial clusters only: KRT programme distribution.
dne<-cell[cell$cluster%in%as.character(non_epi),];dne$cluster<-factor(dne$cluster,levels=as.character(non_epi))
p3<-ggplot(dne,aes(cluster,KRT_mean))+geom_violin(scale="width",trim=TRUE,fill="grey80",color="grey20",linewidth=.25)+labs(x="Non-epithelial cluster",y="Mean LogNormalized KRT expression",title="Keratin signal in independently non-epithelial clusters")+theme_classic(base_size=10)
ggsave(file.path(FD,"06A4_non_epithelial_KRT_violin.pdf"),p3,width=8,height=4.5,bg="white");ggsave(file.path(FD,"06A4_non_epithelial_KRT_violin.png"),p3,width=8,height=4.5,dpi=600,bg="white")

# Batch comparison within independently non-epithelial cells.
batch<-aggregate(cbind(n_KRT_detected,KRT_mean)~batch,data=dne,FUN=mean)
write.csv(batch,file.path(TD,"06A4_non_epithelial_batch_KRT.csv"),row.names=FALSE)

writeLines(c("PIPELINE_V2 — 06A4 KRT ambient contamination audit",paste0("Input: ",IN),paste0("Cells: ",ncol(x)),paste0("KRT panel: ",paste(KRT,collapse=", ")),"Detection defined from RNA/counts > 0; programme magnitude from RNA/data LogNormalized.","Non-epithelial reference clusters: 4,9,11,18,19,21,25,26.","No cell filtering or annotation performed. This is a sensitivity/audit analysis."),file.path(TD,"06A4_audit.txt"))
cat("\n=== NON-EPITHELIAL REFERENCE CLUSTERS ===\n");print(summ[summ$cluster%in%non_epi,c("cluster","n_cells","pct_KRT_ge1","pct_KRT_ge2","pct_KRT_ge3","pct_KRT_ge4","mean_KRT_program","median_doublet_score")],row.names=FALSE)
cat("\n=== KRT BURDEN vs QC ===\n");print(qc,row.names=FALSE)
cat("\nPASS — audit only; frozen atlas unchanged.\n")
