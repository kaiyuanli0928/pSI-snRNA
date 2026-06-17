# Fig2.R
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
# select cluster __ neuron cluster
cells <- WhichCells(object = pbmc_filt, idents = c("/"))  #select neuron cluster

pbmc_filt <- subset(pbmc_filt,cells=cells)
pbmc_filt
#数据intergrate
DefaultAssay(pbmc_filt) <- "RNA"   #注意要加入

#JackStraw
pbmc_filt <- JackStraw(pbmc_filt, num.replicate = 100) 
pbmc_filt <- ScoreJackStraw(pbmc_filt, dims = 1:20)
JackStrawPlot(pbmc_filt,dims = 1:20)
ggsave(file="NeuronResult/JackStraw.pdf",dpi=300,height=6,width=8)
#ElbowPlot
ElbowPlot(pbmc_filt)
ggsave(file="NeuronResult/ElbowPlot.pdf",dpi=300,height=6,width=8)

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
ggsave(file="NeuronResult/umap6.pdf",dpi=300,width=14,height=6)
#umap
umap3<-DimPlot(pbmc_filt,reduction = "umap",split.by ="orig.ident2",pt.size=1.2,label = TRUE)
ggsave(file="NeuronResult/umapsplit.pdf",dpi=300,width=14,height=6)

# Save neuron R Data 
save(pbmc_filt, file = "Result/multiFiltneuron_pbmc.rData")


# Set cluster colors
plan <- FetchData(pbmc_filt, vars = c("umap_1", "umap_2", "seurat_clusters"))
# set plan
plan <- c("#FF0000",  "#00FFFF")   #Define cluster colors based on their characteristics

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
ggsave(p,file="NeuronResult/cell ident.pdf",width=7,height=6)



# load neuron cluster rData
load("Result/multiFiltneuron_pbmc.rData")
# select cluster __ GABAergic neuron cluster
cells <- WhichCells(object = pbmc_filt, idents = c("/")) #select GABAergic neuron cluster

pbmc_filt <- subset(pbmc_filt,cells=cells)
pbmc_filt
#数据intergrate
DefaultAssay(pbmc_filt) <- "RNA"   #注意要加入

#JackStraw
pbmc_filt <- JackStraw(pbmc_filt, num.replicate = 100) 
pbmc_filt <- ScoreJackStraw(pbmc_filt, dims = 1:20)
JackStrawPlot(pbmc_filt,dims = 1:20)
ggsave(file="NeuronResult/JackStraw.pdf",dpi=300,height=6,width=8)
#ElbowPlot
ElbowPlot(pbmc_filt)
ggsave(file="NeuronResult/ElbowPlot.pdf",dpi=300,height=6,width=8)

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
ggsave(file="NeuronResult/umap6.pdf",dpi=300,width=14,height=6)
#umap
umap3<-DimPlot(pbmc_filt,reduction = "umap",split.by ="orig.ident2",pt.size=1.2,label = TRUE)
ggsave(file="NeuronResult/umapsplit.pdf",dpi=300,width=14,height=6)

# Save GABAergic neuron R Data 
save(pbmc_filt, file = "Result/multiFiltGABA_pbmc.rData")


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
write.csv(different,file="GABAResult/DGE.csv")


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
ggsave(file="GABAResult/volcano",width=8,height=12)


#Heatmap_ DGE in male and female GABAergic neuron
diff_order <- different[order(-different$avg_log2FC),] 
write.csv(diff_order,file="GABAResult/diff_order1_GABA_FM.csv")
select.gene <- subset(diff_order, p_val < 3.88E-50) 
write.csv(select.gene,file="GABAResult/diff_order2_GABA_FM.csv")
select.gene <- row.names(select.gene)[1:25] 

# pbmc_filt cell ID
cell_ids <- Cells(pbmc_filt) 
shuffled_cells <- sample(cell_ids)

DoHeatmap(pbmc_filt, features = select.gene)
p <- DoHeatmap(pbmc_filt, size = 3, angle = -50, hjust=0.8,
               features = select.gene,                     
               cells = shuffled_cells,                     
               group.colors = c(F = "red", M = "cyan")
)+  
  scale_fill_gradientn(colors = c("blue4","blue3","blue1","white","red","darkred"))
p
ggsave(file="GABAResult/Heatmap1.pdf",dpi=100,height=4,width=12)
ggsave(p, filename="GABAResult/Heatmap1.emf", width=10,
       height=4, units=c("cm")) 


####
diff_order <- different[order(different$avg_log2FC),] 
select.gene <- subset(diff_order, p_val < 3.88E-50)
select.gene <- row.names(select.gene)[1:25]

# pbmc_filt cell ID
cell_ids <- Cells(pbmc_filt) 
shuffled_cells <- sample(cell_ids)  

DoHeatmap(pbmc_filt, features = select.gene)
p <- DoHeatmap(pbmc_filt, size = 3, angle = -50, hjust=0.8,
               features = select.gene,                    
               cells = shuffled_cells,                      
               group.colors = c(F = "red", M = "cyan")
)+
  scale_fill_gradientn(colors = c("blue4","blue3","blue1","white","red","darkred"))
p
ggsave(file="GABAResult/Heatmap2.pdf",dpi=100,height=4,width=12)
ggsave(p, filename="GABAResult/Heatmap2.emf", width=10,
       height=4, units=c("cm")) 


# GO term
library(org.Mm.eg.db)
gene=rownames(different[different$avg_log2FC<0,])
gene=as.character(na.omit(AnnotationDbi::select(org.Mm.eg.db,
                                                     keys = gene,
                                                     columns = "ENTREZID",
                                                     keytype = "SYMBOL")[,2]))

go<-enrichGO(gene = gene, OrgDb = org.Mm.eg.db, ont="BP", pAdjustMethod = "BH",
                   pvalueCutoff = 0.9,qvalueCutoff =0.9,readable = TRUE)

dotplot(go)
ggsave(file="GABAResult/GO_dotplot.pdf",width=10,height=8)




