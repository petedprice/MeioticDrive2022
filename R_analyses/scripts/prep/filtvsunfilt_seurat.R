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

##### DATA -----
sum_data_fin <- data.frame(
  sample = NA,
  ncells = NA,
  mean_nUMI = NA, 
  mean_nGene = NA, 
  mean_complexity = NA,
  treatment = NA, 
  filt = NA, 
  cell_filt = NA
)

for (filt in c("filtered", "unfiltered")){
  datapath = paste("indata/cellranger/", filt, sep = "")
  files <- list.files(datapath, pattern = "stalkie_dros")
  
  #READING IN COUNT MATRICES and creating sample variables. 
  obj_list <- list()
  for (file in files){
    print(file)
    
    sample <- str_split(file, "_", simplify = TRUE)[,1] %>% 
      str_split("-", simplify = TRUE)   
    sample <- sample[,2]
    
    seurat_data <- Read10X(data.dir = paste(datapath, "/", file, sep = ""))
    seurat_obj <- CreateSeuratObject(counts = seurat_data, 
                                     min.features = 200, ### FILTERING CELLS WITH LESS THAN 200 GENES
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
  for (f in c(100,200)){
    filtered_seurat <- subset(x = merged_seurat, 
                              subset= (nUMI >= 500) & 
                                (nGene >= f) & 
                                #(log10GenesPerUMI > 0) & # Can be dying cells or simple cell types such as blood cells
                                (mitoRatio < 0.05))
    metadata_clean <- filtered_seurat@meta.data
    sum_data <- metadata_clean %>% 
      group_by(sample) %>% 
      summarise(
        ncells = n(), 
        mean_nUMI = mean(nUMI), 
        mean_nGene = mean(nGene), 
        mean_complexity = mean(log10GenesPerUMI), 
        treatment = max(treatment)
        )
    sum_data$filt <- filt
    sum_data$cell_filt <- f
    
    sum_data_fin <- rbind(sum_data_fin, sum_data)
    
    assign(paste(filt, "metadata", sep = "_"), metadata_clean)
    assign(paste(filt, "seurat", sep = "_"), filtered_seurat)
    
    p1 <- ggplot(metadata_clean, aes(x = sample, fill = treatment)) + geom_bar() + 
      labs(y = "Number of cells")
    
    p2 <- ggplot(sum_data, aes(x = mean_nUMI, y = ncells, colour = treatment, shape = treatment)) + 
      geom_point() +
      labs(y = "", x = "mean(number of UMIs)")
    
    save_plot <- ggarrange(p1, p2, common.legend = TRUE, legend = 'bottom')
    save_plot <- annotate_figure(save_plot, top = text_grob(paste("QC for ", filt, " with cell filter", f, " data", sep = ""),
                                          color = "Black", face = "bold", size = 14))
    ggsave(paste("plots/", filt, "_cell_filter", f, "_QC.pdf", sep = ""), save_plot, width = 10, height = 6)
  }
}
sum_data_fin <- sum_data_fin[-which(is.na(sum_data_fin$sample)),]
write.csv(sum_data_fin, "outdata/QC/paramater_checks.csv", quote = FALSE, row.names = FALSE)

#save(filtered_metadata, filtered_seurat, 
#     file="outdata/RData/filtered_seurat.RData")

