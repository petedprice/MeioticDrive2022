#!/usr/bin/env Rscript

###### SETTING UP INPUT COMMANDS ----
library("optparse")
option_list = list(
  make_option(c("-d", "--path_to_seurat_object"), type="character", default=".", 
              help="path to where you have the stored your seurat object", metavar="character"),
  make_option(c("-o", "--output_path"), type="character", default=".", 
              help="where you want to save your output plots and RData files", metavar="character"),
  make_option(c("-t", "--threads"), type="numeric", default=1, 
              help="number of threads for parallelising", metavar="numeric"),
  make_option(c("-s", "--samples"), type="character", default="all", 
              help="path to dataframe containing samples (see format on github)", metavar="character"), 
  make_option(c("-c", "--cellcycle"), type="character", default="data/cell_cycle_markers_complete.csv", 
              help="path to dataframe containing cell cycle markers", metavar="character")
  
)


opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
opt$cellcycle <- "indata/markers/cell_cycle_markers_complete.csv"
cellcycle <- read.table(opt$cellcycle)
cellcycle$TDel_GID <- gsub("gene-", "", cellcycle$TDel_GID)

if (is.null(opt$path_to_seurat_object)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file)", call.=FALSE)
}

###LIBRARIES ----
library(Seurat)
library(tidyverse)
library(Matrix)
library(scales)
library(cowplot)
library(RCurl)
library(stringr)
library(ggpubr)
library(future)
library(future.apply)

#Load data and parsing commands
output_path <- opt$output_path
load(opt$path_to_seurat_object) # path to filtered seurat RData

#make folders
outdatapath = paste(output_path, "/outdata", sep = "")
dir.create(outdatapath, showWarnings = F, recursive = T)
plotpath = paste(output_path, "/plots/", sep = "")
dir.create(plotpath, showWarnings = F, recursive = T)

#parallelise
plan("multicore", workers = opt$threads)
options(future.globals.maxSize = 8000 * 1024^5)

# split object into a list by sample
filtered_seurat <- CellCycleScoring(filtered_seurat, 
                             g2m.features = cellcycle$TDel_GID[cellcycle$phase == "G2/M"], 
                             s.features = cellcycle$TDel_GID[cellcycle$phase == "S"])
split_seurat <- SplitObject(filtered_seurat, split.by = "sample")

if (opt$samples != 'all'){
  print("sample removing")
  keep_samples <- read.table(opt$samples)[,1]
  split_seurat <- split_seurat[keep_samples]
  print(paste("keeping samples ", keep_samples, sep = ""))
}



#SCT normalize the data (SCTransform also accounts for sequencing depth)
print("SCTransform")
split_seurat <- future_lapply(split_seurat, SCTransform, vars.to.regress = 
                                 c("mitoRatio","nUMI","S.Score","G2M.Score")) #may potentially have to regress out cell cycle 


#prep data for integration 
# Identify variable features for integrating
print("SelectIntegrationFeatures")
features <- SelectIntegrationFeatures(object.list = split_seurat, nfeatures = 3000)

#Prepossessing step necessary if SCT transformed
print("PrepSCTIntegration")

split_seurat <- PrepSCTIntegration(object.list = split_seurat, 
                                   anchor.features = features)

#Find anchors that link datasets
print("FindIntegrationAnchors")
anchors <- FindIntegrationAnchors(object.list = split_seurat, 
                                  anchor.features = features, 
                                  normalization.method = "SCT")

#Initial plot making for comparisons to before integration 

# Perform PCA
print("first unintegrated UMAPS")
remerged <- Reduce(merge, split_seurat)
remerged <- RunPCA(object = remerged, features = features)
remerged <- RunUMAP(remerged, 
                             dims = 1:30,
                             reduction = "pca")


pdf("controled_plot.pdf")
DimPlot(remerged, group.by = 'Phase', split.by = "Phase")
dev.off()

remerged <- FindNeighbors(remerged, dims = 1:30, verbose = FALSE)
remerged <- FindClusters(remerged, verbose = FALSE)

save(remerged, file = paste(outdatapath, "/remerged.RData", sep = ""))

fs_PCA1 <- DimPlot(remerged,
                   split.by = "sample")
ggsave(filename = paste(plotpath, "fs_sample_PCA.pdf", sep = ""))
fs_PCA2 <- DimPlot(remerged,
                   split.by = "treatment")
ggsave(filename = paste(plotpath, "fs_treatment_PCA.pdf", sep = ""))


#Integrate data 
print("integrating")
seurat_integrated <- IntegrateData(anchorset = anchors, 
                                   normalization.method = "SCT", 
                                   features.to.integrate = unique(unlist(lapply(split_seurat, rownames))))

#Save data
print("saving data")
save(split_seurat, seurat_integrated, anchors, features, file = paste(outdatapath, "/integrated_seurat.RData", sep = ""))

