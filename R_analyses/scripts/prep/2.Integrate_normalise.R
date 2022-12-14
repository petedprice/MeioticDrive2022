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
              help="path to dataframe containing samples (see format on github)", metavar="character")
  
)


opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

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
split_seurat <- SplitObject(filtered_seurat, split.by = "sample")

if (opt$samples != 'all'){
  print("sample removing")
  keep_samples <- read.table(opt$samples)[,1]
  split_seurat <- split_seurat[keep_samples]
}

split_seurat <- lapply(split_seurat, SCTransform, vars.to.regress = 'mitoRatio') #may potentially have to regress out cell cycle 
split_seurat <- future.apply::future_lapply(split_seurat, SCTransform, vars.to.regress = 'mitoRatio') #may potentially have to regress out cell cycle 

#SCT normalise the data

#prep data for integration 
# Identify variable features for integrating
features <- SelectIntegrationFeatures(object.list = split_seurat, nfeatures = 3000)

#Preprosssesing step neccesary if SCT transformed
split_seurat <- PrepSCTIntegration(object.list = split_seurat, 
                                   anchor.features = features)


#Find anchors that link datasets
anchors <- FindIntegrationAnchors(object.list = split_seurat, 
                                  anchor.features = features, 
                                  normalization.method = "SCT")

#Initial plot making for comparisons to before integration 


# Perform PCA
remerged <- Reduce(merge, split_seurat)
remerged <- RunPCA(object = remerged, features = features)
save(remerged, file = paste(outdatapath, "/remerged.RData", sep = ""))
fs_PCA1 <- PCAPlot(remerged,
                   split.by = "sample")
ggsave(filename = paste(plotpath, "fs_sample_PCA.pdf", sep = ""))
fs_PCA2 <- PCAPlot(remerged,
                   split.by = "treatment")
ggsave(filename = paste(plotpath, "fs_treatment_PCA.pdf", sep = ""))


#Integrate data 
seurat_integrated <- IntegrateData(anchorset = anchors, 
                                   normalization.method = "SCT")

#Save data
save(split_seurat, seurat_integrated, file = paste(outdatapath, "/integrated_seurat.RData", sep = ""))

