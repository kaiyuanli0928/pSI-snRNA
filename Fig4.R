# Fig6+FigS11.R
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
marker<-c("Ar")

# featureplot
p <- FeaturePlot(pbmc_filt,features=marker,ncol=2,pt.size=0.01,label=T,cols = c("lightgrey","red"),split.by ="orig.ident2")
ggsave(p,file="GABAResult/Ar Marker gene.pdf",width=8,height=14)