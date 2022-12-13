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

##### Data and custom scripts ----
path_to_MD2022 <- args[1]
#datapath = paste("indata/cellranger/filtered", sep = "")
datapath = args[2]
output_path = args[3]

files <- list.files(datapath, pattern = "stalkie_dros")
#source("scripts/Bits_and_bobs/Usefull_functions.R")
source(paste(path_to_MD2022, "/R_analyses/scripts/Bits_and_bobs/Usefull_functions.R", sep = ""))

#filtering thresholds 
filt = "filtered"
ftr = 200

#READING IN COUNT MATRICES and creating sample variables. 
obj_list <- list()
for (file in files){
  print(file)
  
  sample <- str_split(file, "_", simplify = TRUE)[,1] %>% 
    str_split("-", simplify = TRUE)   
  sample <- sample[,2]
  
  seurat_data <- Read10X(data.dir = paste(datapath, "/", file, sep = ""))
  
  seurat_obj <- CreateSeuratObject(counts = seurat_data, 
                                   min.features = ftr, ### FILTERING CELLS WITH LESS THAN 200 GENES
                                   project = sample, 
                                   min.cells = 3) #Might be able to include min.cells = 3 for keeping a gene rather than doing this later
  seurat_obj <- RenameCells(seurat_obj, add.cell.id = sample)
  seurat_obj@meta.data$sample <- sample
  seurat_obj@meta.data$treatment <- substr(sample, 1, 2)
  obj_list <- c(obj_list, seurat_obj)
}

merged_seurat <- reduce(obj_list, merge)
merged_seurat$log10GenesPerUMI <- log10(merged_seurat$nFeature_RNA) / log10(merged_seurat$nCount_RNA) # Calculating the number of features per UMI 


#CALCULATING PROPORTION OF READS MAPPING TO MITOCHONDRIAL GENOME (x100 as originally a percentage)
merged_seurat$mitoRatio <- PercentageFeatureSet(object = merged_seurat, pattern = "mt")/ 100 
metadata <- merged_seurat@meta.data # Seperately creating metadata dataframe to save metrics etc without risking affecting seurat object
metadata$cells <- rownames(metadata) #add cell IDs to metadat
# Rename columns
metadata <- metadata %>%
  dplyr::rename(seq_folder = orig.ident,
                nUMI = nCount_RNA,
                nGene = nFeature_RNA)
merged_seurat@meta.data <- metadata #Save the more complete metadat to the seurat object

#FILTERING DATA
filtered_seurat <- subset(x = merged_seurat, 
                            #(log10GenesPerUMI > 0) & # Can be dying cells or simple cell types such as blood cells
                            (mitoRatio < 0.05))
metadata_clean <- filtered_seurat@meta.data

outdatapath = paste(output_path, "/outdata", sep = "")
dir.create(outdatapath, showWarnings = F, recursive = T)
save(filtered_seurat, metadata_clean, paste(outdatapath, "/filtered_seurat.RData", sep = ""))

plotpath = paste(output_path, "/plots", sep = "")
dir.create(plotpath, showWarnings = T, recursive = T)
make_plots_function(metadata_clean, plotpath = plotpath)




