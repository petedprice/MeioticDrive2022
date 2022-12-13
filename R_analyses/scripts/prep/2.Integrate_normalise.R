#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

###LIBRARIES ----
library(Seurat)
library(tidyverse)
library(Matrix)
library(scales)
library(cowplot)
library(RCurl)
library(stringr)
library(ggpubr)

#Load data
#load("outdata/RData/filtered_seurat.RData")
load(args[1]) # path to filtered seurat RData
                          
# split object into a list by sample
split_seurat <- SplitObject(filtered_seurat, split.by = "sample")

#normalise data
split_seurat <- lapply(split_seurat, NormalizeData, verbose = TRUE)

# Score cells for cell cycle
#seurat_phase <- CellCycleScoring(split_seurat, 
#                                 g2m.features = g2m_genes, 
#                                 s.features = s_genes)


#View cell cycle scores and phases assigned to cells                                 
#View(seurat_phase@meta.data)    
#split_seurat[[i]] <- lapply(split_seurat, CellCycleScoring, g2m.features=g2m_genes, s.features=s_genes) # Need to score for cell-cycle markers


#SCT or normalise data (Compare)
split_seurat <- lapply(split_seurat, SCTransform, vars.to.regress = 'mitoRatio') #may potentially have to regress out cell cycle 


#prep data for integration 
features <- SelectIntegrationFeatures(object.list = split_seurat)
split_seurat <- lapply(split_seurat, FindVariableFeatures, 
                       selection.method = "vst", nfeatures = 2000)
split_seurat <- PrepSCTIntegration(object.list = split_seurat, 
                                   anchor.features = features)
anchors <- FindIntegrationAnchors(object.list = split_seurat, 
                                  anchor.features = features, 
                                  normalization.method = "SCT")


#Integrate data 
seurat_integrated <- IntegrateData(anchorset = anchors, 
                                   normalization.method = "SCT")
save(seurat_integrated, file = "outdata/RData/integrated_seurat.RData")
load("outdata/RData/integrated_seurat.RData")
seurat_integrated <- FindVariableFeatures(seurat_integrated, 
                                          selection.method = "vst",
                                          nfeatures = 2000, 
                                          verbose = FALSE)

seurat_integrated <- ScaleData(seurat_integrated)

seurat_integrated <- RunPCA(object = seurat_integrated)
seurat_integrated <- RunUMAP(seurat_integrated, 
                             dims = 1:40,
                             reduction = "pca")
seurat_integrated <- FindNeighbors(object = seurat_integrated, 
                                   dims = 1:40)
seurat_integrated <- FindClusters(object = seurat_integrated,
                                  resolution = c(0.4, 0.6, 0.8, 1.0, 1.4))
# Visualization
#p1 <- DimPlot(seurat_integrated, reduction = "umap", group.by = "treatment")
#p2 <- DimPlot(seurat_integrated, reduction = "umap", label = TRUE, repel = TRUE)
#ggsave("plots/UMAP1.pdf", p2 + p1, width = 15, height = 10)

#DefaultAssay(seurat_integrated) <- "RNA"
#nk.markers <- FindConservedMarkers(seurat_integrated, ident.1 = 6, grouping.var = "treatment", verbose = FALSE)
#head(nk.markers)
########################### ORTHOLOGS ---------
load("outdata/RData/orthologs.RData")
markers <- filter(orthologs_testis, is.na(Cluster) == FALSE) %>% 
  filter(sub("^gene-", "", TDel_GID) %in% rownames(seurat_integrated))





plot_func <- function(cluster, mk_df = markers){
  print(cluster)
  mks <- filter(mk_df, Cluster == cluster)
  mks2 <- str_split(mks$TDel_GID, "gene-", simplify = TRUE)[,2]
  size = length(mks2) * 1.5
  f <- FeaturePlot(seurat_integrated, features = mks2, min.cutoff = "q10")
  ggsave(paste("plots/", cluster, "_feature_plot.pdf", sep = ""), f, height = size, width = 1.5* size)
}

lapply(unique(markers$Cluster), plot_func, mk_df = markers)


ggsave("plots/featureplot_markers.pdf", plots, height = 30, width = 30)


f <- FeaturePlot(seurat_integrated, features = marker_genes[51:80], min.cutoff = "q9")
ggsave("plots/featureplot_markers.pdf", f, height = 30, width = 30)
d <- DimPlot(seurat_integrated, reduction = "umap", split.by = "treatment", height = 25, width = 25)
ggsave("plots/del.pdf", d)

