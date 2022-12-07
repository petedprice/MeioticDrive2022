###LIBRARIES ----
library(SingleCellExperiment)
library(Seurat)
library(tidyverse)
library(Matrix)
library(scales)
library(cowplot)
library(RCurl)
library(stringr)
library(ggpubr)

#Load data
load("outdata/RData/filtered_seurat.RData")

                            

# split object into a list by sample
split_seurat <- SplitObject(filtered_seurat, split.by = "sample")

#normalise data
split_seurat <- lapply(split_seurat, NormalizeData, verbose = TRUE)

# Score cells for cell cycle
seurat_phase <- CellCycleScoring(split_seurat, 
                                 g2m.features = g2m_genes, 
                                 s.features = s_genes)

# View cell cycle scores and phases assigned to cells                                 
View(seurat_phase@meta.data)    
#split_seurat[[i]] <- lapply(split_seurat, CellCycleScoring, g2m.features=g2m_genes, s.features=s_genes) # Need to score for cell-cycle markers

#SCT transform data 
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
