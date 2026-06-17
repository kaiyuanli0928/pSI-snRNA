# Fig1.R
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

#set way
setwd("")
dir()
 # LoadFiles
M001 <- Read10X(data.dir = "")
# Filter
M001 <- CreateSeuratObject(counts = M001, project = "M1", min.cells = 3, min.features = 200)

#set way
setwd("")
dir()
# LoadFiles
F001 <- Read10X(data.dir = "")
# Filter
F001 <- CreateSeuratObject(counts = F001, project = "F1", min.cells = 3, min.features = 200)

#set way
setwd("")
dir()
# LoadFiles
M002 <- Read10X(data.dir = "")
# Filter
M002 <- CreateSeuratObject(counts = M002, project = "M2", min.cells = 3, min.features = 200)

#set way
setwd("")
dir()
# LoadFiles
F002 <- Read10X(data.dir = "")
# Filter
F002 <- CreateSeuratObject(counts = F002, project = "F2", min.cells = 3, min.features = 200)

# Merge
pbmc<- merge(M001,y=c(F001, M002, F002))


sample_info<-c("M","F","M","F")
Header<-colnames(pbmc)
project<-lapply(1:length(sample_info),function(x){
  sample_name<-sample_info[x]
  pattern<-paste0("_",x,"$")
  Index<-grep(pattern,Header,value=F)
  info<-c(rep(sample_name,length(Index)))
  return(info)
})
project=unlist(project)
pbmc[['orig.ident2']]<-project


pbmc[["percent.mito"]] <-PercentageFeatureSet(pbmc, pattern = "^MT-")


# Set save path
setwd("")
###创建输出目录
if(!dir.exists("Result")){dir.create("Result")}
dir()

p<-VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3,pt.size=0,group.by = "orig.ident")
ggsave(file="Result/ngene_numi_pmito_vlnplot3.pdf",dpi=300,height=6,width=8)
p<-VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3,pt.size=0,group.by = "orig.ident2")
ggsave(file="Result/ngene_numi_pmito_vlnplot4.pdf",dpi=300,height=6,width=8)

pbmc_filt <- subset(pbmc,subset = nFeature_RNA>950 & nFeature_RNA<11500 & nCount_RNA>-Inf & nCount_RNA<Inf & percent.mito<15)
pbmc_filt

p<-VlnPlot(pbmc_filt, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3,group.by = "orig.ident")
ggsave(file="Result/filter_ngene_numi_pmito_vlnplot.pdf",dpi=300,height=10,width=8)

p<-VlnPlot(pbmc_filt, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3,group.by = "orig.ident2")
ggsave(file="Result/filter_ngene_numi_pmito_vlnplot2.pdf",dpi=300,height=10,width=8)

pbmc_filt <- NormalizeData(pbmc_filt,normalization.method = "LogNormalize", scale.factor = 10000)
pbmc_filt <- FindVariableFeatures(object = pbmc_filt,selection.method = "vst",nfeatures = 2000)
top10 <- head(VariableFeatures(pbmc_filt), 10)
plot1 <- VariableFeaturePlot(pbmc_filt)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
p <- plot1 + plot2
ggsave(file="Result/VariableFeature.pdf",dpi=300,width=12,height=6)

all.genes <- rownames(pbmc_filt)
pbmc_filt <- ScaleData(object = pbmc_filt,features =all.genes) 

pbmc_filt <- RunPCA(object = pbmc_filt)###PCA降维
# Examine and visualize PCA results a few different ways
print(pbmc_filt[["pca"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(pbmc_filt, dims = 1:2, reduction = "pca")
DimPlot(pbmc_filt, reduction = "pca")
#DimHeatmap
DimHeatmap(pbmc_filt, dims = 1:10, cells = 500, balanced = TRUE)

#JackStraw
pbmc_filt <- JackStraw(pbmc_filt, num.replicate = 100) 
pbmc_filt <- ScoreJackStraw(pbmc_filt, dims = 1:20)
JackStrawPlot(pbmc_filt,dims = 1:20)
ggsave(file="Result/JackStraw.pdf",dpi=300,height=6,width=8)
#ElbowPlot
ElbowPlot(pbmc_filt)
ggsave(file="Result/ElbowPlot.pdf",dpi=300,height=6,width=8)

#harmony
pbmc_filt = pbmc_filt %>% RunHarmony("orig.ident", plot_convergence = TRUE)

pbmc_filt <- pbmc_filt %>%
  RunUMAP(reduction = "harmony", dims = 1:12) %>%   
  FindNeighbors(reduction = "harmony", dims = 1:12) %>%
  FindClusters(resolution = 0.6) %>% 
  identity()

#umap
umap1<-DimPlot(pbmc_filt,reduction = "umap",group.by ="orig.ident2",pt.size=1.2)
umap2<-DimPlot(pbmc_filt,reduction = "umap",pt.size=1.2,label = TRUE)
p<-umap1 + umap2
ggsave(file="Result/umap6.pdf",dpi=300,width=14,height=6)
#umap
umap3<-DimPlot(pbmc_filt,reduction = "umap",split.by ="orig.ident2",pt.size=1.2,label = TRUE)
ggsave(file="Result/umapsplit.pdf",dpi=300,width=14,height=6)

# Save R Data 
save(pbmc_filt, file = "Result/multiFilt_pbmc.rData")


# Cluster statistic
count_table <- table(pbmc_filt@meta.data$seurat_clusters, pbmc_filt@meta.data$orig.ident2)
count_table
write.csv(count_table,file="Result/ncounts_cluster.csv")
#### 可视化
pn<-count_table %>%  as.data.frame() %>% ggbarstats(x = Var2, y = Var1, counts = Freq)
ggsave(pn,file="Result/ncounts_cluster.pdf",width=16,height=8)


# Feature plot
marker<-c("/")
p <- FeaturePlot(pbmc_filt,features=marker,ncol=2,pt.size=0.05,cols = c("lightgrey","#E41A1C"))
ggsave(p,file="Result/Marker gene.pdf",width=8,height=10)


#Print Veen
library(ggvenn)
#select gene
Thy1_pos <- WhichCells(pbmc_filt, expression = Thy1 > 0)
Camk2a_pos <- WhichCells(pbmc_filt, expression = Camk2a > 0)
Gad2_pos <- WhichCells(pbmc_filt, expression = Gad2 > 0)
Slc32a1_pos <- WhichCells(pbmc_filt, expression = Slc32a1 > 0)

venn_data <- list(Thy1 = Thy1_pos,Camk2a = Camk2a_pos,Gad2 = Gad2_pos,Slc32a1 = Slc32a1_pos)

p <- ggvenn(venn_data,
             fill_color = c('#31C53F','#F68282','#D4D915','#ff9a36','#E63863'),  # Nature期刊常用配色
             stroke_size = 0.8,
             set_name_size = 5,
             text_size = 4.5,
             show_percentage = TRUE, 
             digits = 1) +
  labs(title = "neuron Expression Overlap",
       subtitle = paste("Total cells:", ncol(pbmc_filt))) +
  theme_void() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 12)
  )
ggsave(p,file="Result/Gene_Venn.pdf",width=16,height=20)


# Set cluster colors
plan <- FetchData(pbmc_filt, vars = c("umap_1", "umap_2", "seurat_clusters"))
# set plan
plan <- c("#FF0000", "#00FFFF", "#008000", "#A52A2A","#FFA500","#5F9EA0",  "#5F9EA0","#FFA500", "#FF1493") 
#Define cluster colors based on their characteristics

umap_coords <- Embeddings(pbmc_filt, reduction = "umap")
plot_data <- data.frame(
  umap_1 = umap_coords[, 1],         # X axis
  umap_2 = umap_coords[, 2],         # Y axis
  seurat_clusters = pbmc_filt$seurat_clusters 
)
p <- ggplot(plot_data, aes(x = umap_1, y = umap_2, color = seurat_clusters)) +
  geom_point(size = 0.5, alpha = 0.8) +
  scale_color_manual(values = plan) +  
  theme_classic() +
  labs(
    title = "Custom Cluster Colors with ggplot2",
    x = "UMAP_1",
    y = "UMAP_2",
    color = "Cluster" ) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "right")
ggsave(p,file="Result/cell ident.pdf",width=7,height=6)





