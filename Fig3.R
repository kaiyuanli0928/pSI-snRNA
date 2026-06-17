# Fig3.R
# Load Packages
library(Seurat)  
library(multtest)
library(dplyr)
library(ggplot2)
library(patchwork)
library(SeuratData)
library(tidyverse)
library(Rcpp)
library(harmony)
library(ggstatsplot)
library(EnhancedVolcano)
library(dplyr)


# load GABAergic neuron cluster rData
load("Result/multiFiltGABA_pbmc.rData")
#Analysis select DGE gene
marker<-c("Hcn1","Hcn2","Hcn3","Hcn4")

# featureplot
p <- FeaturePlot(pbmc_filt,features=marker,ncol=2,pt.size=0.01,label=T,cols = c("lightgrey","red"),split.by ="orig.ident2")
ggsave(p,file="GABAResult/HCN Marker gene.pdf",width=8,height=14)

# VlnPlot
p <- VlnPlot(pbmc_filt,features=marker,pt.size=0.05,stack=T,flip=T,fill.by = "ident",adjust=0.8,cols=plan,split.by ="orig.ident2")+NoLegend()
ggsave(p,file="GABAResult/HCNvlnplot_pro.pdf",width=16,height=8)


# select gene positive cell
expr <- LayerData(pbmc_filt, layer = "counts")
gene_expression <- expr %>% .["Hcn1",] %>% as.data.frame()
colnames(gene_expression) <- "Hcn1"
gene_expression$cell <- rownames(gene_expression)
# positive cell
gene_expression_sel_pos <- gene_expression[which(gene_expression$Hcn1>0),]
pbmc_filt_pos <- pbmc_filt[,rownames(gene_expression_sel_pos)]
# select negative cell
gene_expression_sel_neg <- gene_expression[which(gene_expression$Hcn1==0),]  
pbmc_filt_neg <- pbmc_filt[,rownames(gene_expression_sel_neg)]

# add idents OE
sample_info<-c("OE","OE") 
Header<-colnames(pbmc_filt_pos)
project<-lapply(1:length(sample_info),function(x){
  sample_name<-sample_info[x]
  pattern<-paste0("_",x,"$")
  Index<-grep(pattern,Header,value=F)
  info<-c(rep(sample_name,length(Index)))
  return(info)
})
project=unlist(project)
pbmc_filt_pos[['orig.ident3']]<-project
head(pbmc_filt_pos)


# add idents KO
sample_info<-c("KO","KO") 
Header<-colnames(pbmc_filt_neg)
project<-lapply(1:length(sample_info),function(x){
  sample_name<-sample_info[x]
  pattern<-paste0("_",x,"$")
  Index<-grep(pattern,Header,value=F)
  info<-c(rep(sample_name,length(Index)))
  return(info)
})
project=unlist(project)
pbmc_filt_neg[['orig.ident3']]<-project
head(pbmc_filt_neg)

# Merge
pbmc_filt3<- merge(pbmc_filt_pos, y=c(pbmc_filt_neg))
head(pbmc_filt3)
tail(pbmc_filt3)
table(pbmc_filt3@meta.data$orig.ident3)
DefaultAssay(pbmc_filt3) <- "RNA"   #注意要加入

#DGE
different <- FindMarkers(
  pbmc_filt3,
  ident.1 = "OE",
  ident.2 = "KO",
  logfc.threshold = 0.1,  
  min.pct = 0.2,           
  test.use = "MAST"        
)
write.csv(different,file="GABAResult/DGE_HCN_OEKO_MAST.csv")

# GO
library(org.Mm.eg.db)
gene=rownames(different[different$avg_log2FC>0,])
gene=as.character(na.omit(AnnotationDbi::select(org.Mm.eg.db,
                                                   keys = gene_up,
                                                   columns = "ENTREZID",
                                                   keytype = "SYMBOL")[,2]))

go<-enrichGO(gene = gene, OrgDb = org.Mm.eg.db, ont="BP", pAdjustMethod = "BH",
                pvalueCutoff = 0.9,qvalueCutoff =0.9,readable = TRUE)

dotplot(go)
ggsave(file="GABAResult/go_HCN_dotplot.pdf",width=10,height=8)













