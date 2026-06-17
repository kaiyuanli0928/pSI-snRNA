# FigS2.R
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

# load rData
load("Result/multiFilt_pbmc.rData")

#umap
umap1<-DimPlot(pbmc_filt,reduction = "umap",group.by ="orig.ident2",pt.size=1.2)
umap2<-DimPlot(pbmc_filt,reduction = "umap",pt.size=1.2,label = TRUE)
p<-umap1 + umap2
ggsave(file="Result/umap6.pdf",dpi=300,width=14,height=6)
#umap
umap3<-DimPlot(pbmc_filt,reduction = "umap",split.by ="orig.ident2",pt.size=1.2,label = TRUE)
ggsave(file="Result/umapsplit.pdf",dpi=300,width=14,height=6)


# Correlation analysis
table(pbmc_filt$seurat_clusters)
av <- AverageExpression(pbmc_filt, group.by = "seurat_clusters", assays = "RNA")
av <- av[[1]]
av <- as.matrix(av) 
head(av)

cg <- names(tail(sort(apply(av, 1, sd)), 1000))
View(av[cg, ])
View(cor(av[cg, ], method = "spearman"))

p<-pheatmap::pheatmap(cor(av[cg,],method="spearman")) 
ggsave(p,file="Result/correlation_pheatmap.pdf",width=8,height=6)


# Cluster Rename
new.cluster.ids <- c("/","/") #ident cluster name
names(new.cluster.ids)<- levels(pbmc_filt)
#  names(new.cluster.ids)
pbmc_filt <- RenameIdents(pbmc_filt, new.cluster.ids) 
DimPlot(pbmc_filt, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

# change cluster name and cluster sequence 
my_levels <- c( "/","/")
#change cluster sequence 
my_levels <- c( "/","/",)
# Relevel object@timepoint
levels(pbmc_filt) <- my_levels

# DotPlot
#Cluster Marker gene
marker<-c("/")
p <- DotPlot(pbmc_filt,features=marker,dot.scale = 10,cols = c("blue","red"),col.min = 0)
ggsave(p,file="Result/Cluster Marker_dotplot_pro.pdf",width=13,height=3.5)


# DGE
pbmc_filt <- JoinLayers(pbmc_filt, assay = "RNA")
Idents(pbmc_filt)<-"orig.ident2" 
different <- FindMarkers(pbmc_filt,
                         ident.1 = "F",
                         ident.2 = "M",
                         logfc.threshold = 0.1,  
                         min.pct = 0.2,           
                         test.use = "MAST"        
)
write.csv(different,file="Result/DGE.csv")


# Volcano plot
key_genes <- c("/")

res <- different
EnhancedVolcano(
  toptable = res,
  lab = ifelse(res$gene %in% key_genes, res$gene, " "),  # 只标指定基因
  x = "log2FoldChange",
  y = "padj",
  col = c("grey", "blue", "blue", "red"),
  pCutoff = 0.05,
  FCcutoff = 0.5,
  labSize = 4.0,
  max.overlaps = Inf,
  xlim = c(-4,4), #限制x轴的范围
  ylim = c(-20,350) #限制y轴的范围
)
ggsave(file="Result/volcano",width=8,height=12)



#select Glut cluster
cells <- WhichCells(object = pbmc_filt, idents = c("/"))
pbmc_filt1 <- subset(pbmc_filt,cells=cells)
pbmc_filt1
DefaultAssay(pbmc_filt1) <- "RNA" 

#add Glu idents
sample_info<-c("Glu","Glu") 
Header<-colnames(pbmc_filt1)
project<-lapply(1:length(sample_info),function(x){
  sample_name<-sample_info[x]
  pattern<-paste0("_",x,"$")
  Index<-grep(pattern,Header,value=F)
  info<-c(rep(sample_name,length(Index)))
  return(info)
})
project=unlist(project)
pbmc_filt1[['orig.ident4']]<-project
head(pbmc_filt1)


# select GABA cluster
cells <- WhichCells(object = pbmc_filt, idents = c("/"))
cells  #cells可以查看细胞的barcode
pbmc_filt2 <- subset(pbmc_filt,cells=cells)
pbmc_filt2
#数据intergrate
DefaultAssay(pbmc_filt2) <- "RNA"  

# add GABA idents
sample_info<-c("GABA","GABA")
Header<-colnames(pbmc_filt2)
project<-lapply(1:length(sample_info),function(x){
  sample_name<-sample_info[x]
  pattern<-paste0("_",x,"$")
  Index<-grep(pattern,Header,value=F)
  info<-c(rep(sample_name,length(Index)))
  return(info)
})
project=unlist(project)
pbmc_filt2[['orig.ident3']]<-project
head(pbmc_filt2)

# merge
pbmc_filt3<- merge(pbmc_filt1, y=c(pbmc_filt2))
head(pbmc_filt3)
tail(pbmc_filt3)
table(pbmc_filt3@meta.data$orig.ident3)
DefaultAssay(pbmc_filt3) <- "RNA"   #注意要加入


pbmc_filt3 <- NormalizeData(pbmc_filt3,normalization.method = "LogNormalize", scale.factor = 10000)
pbmc_filt3 <- FindVariableFeatures(object = pbmc_filt3,selection.method = "vst",nfeatures = 2000)
top10 <- head(VariableFeatures(pbmc_filt3), 10)
plot1 <- VariableFeaturePlot(pbmc_filt3)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
CombinePlots(plots = list(plot1, plot2))
ggsave(file="GABAResult/VariableFeature_GABA_Glu.pdf",dpi=300,width=12,height=6)

all.genes <- rownames(pbmc_filt3)
pbmc_filt3 <- ScaleData(object = pbmc_filt3,features =all.genes) 
pbmc_filt3 <- RunPCA(object = pbmc_filt3,features = VariableFeatures(object = pbmc_filt3)) 
print(pbmc_filt3[["pca"]], dims = 1:5, nfeatures = 5)
#可视化
DimPlot(pbmc_filt3, reduction = "pca")
DimHeatmap(pbmc_filt3, dims = 1:6, cells = 500, balanced = TRUE)
pbmc_filt3 <- JackStraw(pbmc_filt3, num.replicate = 100)
pbmc_filt3 <- ScoreJackStraw(pbmc_filt3, dims = 1:20)
plot1<-JackStrawPlot(pbmc_filt3, dims = 1:15)
plot2<-ElbowPlot(pbmc_filt3)
CombinePlots(plots = list(plot1, plot2))
ggsave(file="GABAResult/JackStrawPlot_ElbowPlot_GABA_Glu.pdf",dpi=300,width=12,height=6)

#harmony
pbmc_filt3 = pbmc_filt3 %>% RunHarmony("orig.ident", plot_convergence = TRUE)
pbmc_filt3 <- pbmc_filt3 %>%
  RunUMAP(reduction = "harmony", dims = 1:16) %>%   
  FindNeighbors(reduction = "harmony", dims = 1:16) %>%
  FindClusters(resolution = 0.4) %>%  
  identity()

save(pbmc_filt3, file = "GABAResult/sub_pbmc.rData")


# DGE
pbmc_filt3 <- JoinLayers(pbmc_filt3, assay = "RNA")
Idents(pbmc_filt3)<-"orig.ident2" 
different <- FindMarkers(pbmc_filt3,
                         ident.1 = "F",
                         ident.2 = "M",
                         logfc.threshold = 0.1,  
                         min.pct = 0.2,           
                         test.use = "MAST"        
)
write.csv(different,file="GABAResult/DGE_Glut_GABA.csv")


diff_order <- different[order(-different$avg_log2FC),] 
write.csv(diff_order,file="GABAResult/diff_order_Glu_GABA.csv")
select.gene <- subset(diff_order, p_val < 3.88E-50) 
write.csv(select.gene,file="GABAResult/diff_order2_Glu_GABA.csv")
select.gene <- row.names(select.gene)[1:25] 

p <- DoHeatmap(pbmc_filt3, size = 3, angle = -50, hjust=0.8,features = select.gene)+
  scale_fill_gradientn(colors = c("blue4","blue3","blue1","white","red","darkred"))
p
ggsave(file="GABAResult/Heatmap1_GABA_Glu.pdf",dpi=100,height=4,width=12)
ggsave(p, filename="GABAResult/Heatmap_GABA_Glu.emf", width=10,
       height=4, units=c("cm")) 


diff_order <- different[order(different$avg_log2FC),] 
select.gene <- subset(diff_order, p_val < 3.88E-50)
select.gene <- row.names(select.gene)[1:25]

p <- DoHeatmap(pbmc_filt3, size = 3, angle = -50, hjust=0.8,features = select.gene)+
  scale_fill_gradientn(colors = c("blue4","blue3","blue1","white","red","darkred"))
p
ggsave(file="GABAResult/Heatmap2_GABA_Glu.pdf",dpi=100,height=4,width=12)
ggsave(p, filename="GABAResult/Heatmap2_GABA_Glu.emf", width=10,
       height=4, units=c("cm")) 


#Analysis select DGE gene
marker<-c("/")
# VlnPlot
p <- VlnPlot(pbmc_filt,features=marker,pt.size=0.05,stack=T,flip=T,fill.by = "ident",adjust=0.8,cols=plan,split.by ="orig.ident2")+NoLegend()
ggsave(p,file="GABAResult/Glut_GABA_vlnplot_pro.pdf",width=16,height=8)


#Analysis select DGE gene
marker<-c("/")
# Ridge plot
p <-RidgePlot(pbmc_filt, features = marker, ncol = 1,group.by ="orig.ident2")
ggsave(p,file="GABAResult/Agg_gene_RidgePlot.pdf",width=8,height=15)


#Analysis select DGE gene
marker<-c("/")
# FeaturePlot
p <- FeaturePlot(pbmc_filt,features=marker,ncol=2,pt.size=0.01,label=T,cols = c("lightgrey","red"),split.by ="orig.ident2")
ggsave(p,file="GABAResult/F_M Marker_gene.pdf",width=8,height=14)






