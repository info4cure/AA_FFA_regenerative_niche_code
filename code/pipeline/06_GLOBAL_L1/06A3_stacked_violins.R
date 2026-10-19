# =============================================================================
# PIPELINE_V2 — 06A3_stacked_violins.R
# Publication-oriented stacked violins for global L1 identity/purity review.
# READ-ONLY: frozen atlas is never modified.
# =============================================================================
suppressPackageStartupMessages({library(Seurat);library(ggplot2);library(patchwork)})
options(stringsAsFactors=FALSE);set.seed(20260924)
source("PIPELINE_V2/config/paths.R")

cat("\n=== 06A3 STACKED VIOLINS ===\n")
IN<-file.path(V2_FROZEN,"global_atlas_v2_frozen.rds")
OUT<-file.path(V2_ROOT,"results","annotation","06_GLOBAL_L1","stacked_violins")
dir.create(OUT,recursive=TRUE,showWarnings=FALSE)

x<-readRDS(IN);DefaultAssay(x)<-"RNA";CC<-"RNA_snn_res.0.6"
x$global_cluster<-factor(as.character(x@meta.data[[CC]]),levels=as.character(0:26))
Idents(x)<-"global_cluster"
stopifnot(ncol(x)==82893,nrow(x)==23412)

# Broad L1 + selected cluster-discriminating markers.
panel<-list(
  Epithelial=c("KRT5","KRT14","KRT15","KRT17","KRT1","KRT10"),
  Follicular=c("KRT75","KRT83","KRT85","TCHH","LGR5","FGF18"),
  Sebaceous=c("PPARG","AWAT2","MOGAT1"),
  Sweat_gland=c("AQP5","DCD","PIP","MUC7","FOXI1"),
  Fibroblast=c("COL1A1","COL1A2","DCN","LUM","PDGFRA","SFRP2"),
  Mural=c("RGS5","CSPG4","PDGFRB","ACTG2","CNN1","MYH11"),
  Vascular_EC=c("PECAM1","VWF","ACKR1","EMCN"),
  Lymphatic_EC=c("CCL21","PROX1","LYVE1","FLT4"),
  T_NK=c("CD3D","CD3E","NKG7","GNLY"),
  Myeloid_DC=c("LST1","LYZ","CD1C","FCER1G"),
  Melanocyte=c("PMEL","MLANA","TYR","DCT"),
  Neural_Schwann=c("S100B","PLP1","SOX10")
)
genes<-unique(unlist(panel,use.names=FALSE));present<-genes[genes%in%rownames(x)];missing<-setdiff(genes,present)
write.csv(data.frame(gene=genes,present=genes%in%present),file.path(OUT,"06A3_marker_coverage.csv"),row.names=FALSE)
if(length(missing))cat("Missing:",paste(missing,collapse=", "),"\n")
cat("Markers present:",length(present),"/",length(genes),"\n")

# Pull only requested LogNormalized genes. FetchData keeps memory manageable.
df<-FetchData(x,vars=c("global_cluster",present),layer="data")
df$global_cluster<-factor(df$global_cluster,levels=as.character(0:26))

# Quantitative cluster × gene means for downstream review.
mean_rows<-lapply(present,function(g){
  a<-aggregate(df[[g]],list(global_cluster=df$global_cluster),mean)
  names(a)<-c("global_cluster","mean_log_normalized")
  a$gene<-g;a
})
means<-do.call(rbind,mean_rows)
write.csv(means,file.path(OUT,"06A3_cluster_gene_mean_expression.csv"),row.names=FALSE)

# Chatomics-like stacked violin:
# Each gene gets its own panel and y-scale. Clusters remain columns.
make_violin<-function(g,show_x=FALSE){
  d<-data.frame(global_cluster=df$global_cluster,expression=df[[g]])
  ggplot(d,aes(x=global_cluster,y=expression,group=global_cluster))+
    geom_violin(scale="width",trim=TRUE,linewidth=.18,fill="grey75",color="grey20")+
    labs(x=NULL,y=g)+
    theme_classic(base_size=7)+
    theme(
      axis.title.y=element_text(size=7,face="bold",angle=0,vjust=.5,hjust=1,margin=margin(r=5)),
      axis.text.y=element_blank(),axis.ticks.y=element_blank(),axis.line.y=element_blank(),
      axis.text.x=if(show_x) element_text(size=6,color="black") else element_blank(),
      axis.ticks.x=if(show_x) element_line(linewidth=.2) else element_blank(),
      axis.line.x=element_line(linewidth=.25),
      plot.margin=margin(0,3,0,32)
    )
}

plots<-lapply(seq_along(present),function(i)make_violin(present[i],show_x=i==length(present)))
stack<-wrap_plots(plots,ncol=1)+plot_annotation(
  title="Global atlas — stacked marker expression by computational cluster",
  theme=theme(plot.title=element_text(face="bold",size=11,hjust=.5))
)

# Large audit version.
nG<-length(present)
ggsave(file.path(OUT,"06A3_stacked_violins_full.pdf"),stack,width=11,height=max(14,nG*.34),units="in",bg="white",limitsize=FALSE)
ggsave(file.path(OUT,"06A3_stacked_violins_full.png"),stack,width=11,height=max(14,nG*.34),units="in",dpi=300,bg="white",limitsize=FALSE)

# Compact high-specificity panel for rapid L1 review.
compact<-c("KRT14","KRT10","KRT83","TCHH","AWAT2","DCD","PIP",
           "COL1A1","LUM","RGS5","MYH11","PECAM1","ACKR1","CCL21","LYVE1",
           "CD3D","NKG7","LST1","CD1C","MLANA","TYR","PLP1","SOX10")
compact<-compact[compact%in%present]
plots2<-lapply(seq_along(compact),function(i)make_violin(compact[i],show_x=i==length(compact)))
stack2<-wrap_plots(plots2,ncol=1)+plot_annotation(
  title="Global atlas — compact L1 marker panel",
  theme=theme(plot.title=element_text(face="bold",size=11,hjust=.5))
)
ggsave(file.path(OUT,"06A3_stacked_violins_compact.pdf"),stack2,width=11,height=max(9,length(compact)*.36),units="in",bg="white",limitsize=FALSE)
ggsave(file.path(OUT,"06A3_stacked_violins_compact.png"),stack2,width=11,height=max(9,length(compact)*.36),units="in",dpi=600,bg="white",limitsize=FALSE)

# Additional focused incompatibility panels: visualization only.
focus_sets<-list(
  Epithelial_vs_immune=c("KRT5","KRT14","KRT17","CD3D","NKG7","LST1","FCER1G"),
  Epithelial_vs_stroma=c("KRT5","KRT14","KRT17","COL1A1","DCN","LUM","PECAM1","VWF"),
  Mural=c("RGS5","CSPG4","PDGFRB","ACTG2","CNN1","MYH11")
)
for(nm in names(focus_sets)){
  gg<-focus_sets[[nm]][focus_sets[[nm]]%in%present]
  pp<-lapply(seq_along(gg),function(i)make_violin(gg[i],show_x=i==length(gg)))
  obj<-wrap_plots(pp,ncol=1)+plot_annotation(title=gsub("_"," ",nm))
  ggsave(file.path(OUT,paste0("06A3_focus_",nm,".pdf")),obj,width=11,height=max(4,length(gg)*.48),bg="white",limitsize=FALSE)
}

writeLines(c(
  "PIPELINE_V2 — 06A3 stacked violin audit",
  paste0("Input: ",IN),
  paste0("Cells: ",ncol(x)),
  paste0("Genes in atlas: ",nrow(x)),
  paste0("Requested markers: ",length(genes)),
  paste0("Present markers: ",length(present)),
  "Expression layer: RNA/data (LogNormalized full-gene expression)",
  "Clusters: RNA_snn_res.0.6, ordered 0:26",
  "No cells removed; no annotation changed.",
  "FibConsensus and BJD follicular signatures intentionally deferred to L2."
),file.path(OUT,"06A3_audit.txt"))

cat("\nOutputs:",OUT,"\n")
cat("PASS — frozen atlas unchanged.\n")
