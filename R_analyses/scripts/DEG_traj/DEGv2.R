## based off http://bioconductor.org/books/3.15/OSCA.multisample/multi-sample-comparisons.html#putting-it-all-together
# load libraries and functions
library(dplyr)
library(Seurat)
library(HGNChelper)
library(openxlsx)
library(SCINA)
library(ggpubr)
#library(scater)
library(scuttle)
library(edgeR)
library(statmod)
library(scran)
library(ggpubr)
library(gridExtra)
library(tidyverse)
library(data.table)


load("data/RData/testis.Cluster_marker_seurat.RData")
ortholog_table <- read.table("data/ortholog_table.txt")

new_names <- ortholog_table$Dros_GID %>% 
  setNames(gsub("gene-", "", ortholog_table$TDel_GID))
new_names[is.na(new_names) == TRUE] <- ortholog_table$product[is.na(new_names)==T]
new_names[new_names == "missing"] <- ortholog_table$TDel_GID[new_names == "missing"]


##### BIOCONDUCTOR EDGER2 -----
## PLOT FUNCTION ---- 
make_pca_function <- function(x, ortholog_table){
  md <- metadata(x)
  mds_data <- md$y %>% 
    cpm(log = TRUE) %>% 
    plotMDS()
  
  PCA <- data.frame(x = mds_data$x, 
                      y = mds_data$y, 
                      names = md$y$samples$seq_folder,
                      treatment = md$fit$samples$treatment) %>% 
    ggplot(aes(x = x, y = y, fill = names, label = names, colour = treatment))  +
    geom_text(hjust="inward", vjust="inward") +
    labs(x = "PC2", y = "PC1", title = md$y$samples$customclassif[1]) 
  pdf(paste("plots/DEG/PCA_", md$y$samples$customclassif[1], ".pdf", sep = ""))
  print(PCA) 
  dev.off()
  x <- x[is.na(x$logFC) == FALSE,] %>%
    as.data.frame()
  
  
  x$geneID <- rownames(x)
  x$gene <- do.call(recode, c(list(rownames(x)), new_names))
  
  top_p <- x$FDR[order(x$FDR)][11]
  top_p2 <- x$FDR[order(x$FDR)][21]
  
  Volcano <- x %>% 
    ggplot(aes(x = logFC, y = -log10(FDR), label = gene,
               colour = (abs(logFC) > 1 &  FDR < 0.05))) + 
    geom_point() + geom_text(data=subset(x, FDR < top_p & abs(logFC) > 1),
                             aes(logFC,-log10(FDR),label=gene), 
                             hjust="inward", vjust="inward") +
    labs(y = "-log10 Pvalue (FDR adjusted)", x = "log2 fold-change", 
         title = paste(md$y$samples$customclassif[1], "(positive is ST bias)"))+ 
    theme(legend.position = "none")  + 
    geom_hline(yintercept = -log10(0.05), color = 'grey', linetype = 'dashed') + 
    geom_vline(xintercept = c(-1,1), color = 'grey', linetype = 'dashed')
  pdf(paste("plots/DEG/DEG_", md$y$samples$customclassif[1], ".pdf", sep = ""))
  print(Volcano) 
  dev.off()
  ST_bias_genes <- filter(x, logFC > 1, FDR < top_p2)[,c('gene', 'geneID')] 
  SR_bias_genes <- filter(x, logFC < -1, FDR < top_p2)[,c('gene', 'geneID')]
  output <- lst(PCA,Volcano, ST_bias_genes, SR_bias_genes)

  return(output)

}


sce <- as.SingleCellExperiment(seurat_marker, assay = "RNA")


#### BULK RNASEQ ACROSS WHOLE TISSUE ----
whole_tissue <- aggregateAcrossCells(sce, id=colData(sce)[,c("sample")])
de.results_wt <- pseudoBulkDGE(whole_tissue, 
                            label='bulk',
                            condition = whole_tissue$treatment,
                            design=~treatment,
                            coef = "treatmentst",
)
metadata(de.results_wt$bulk)$y$samples$customclassif <- "bulk"

bulk_DEG_output <- make_pca_function(de.results_wt$bulk)


### PSEUDOBULK INDIVIDUAL CELL TYPES ----
agg_cell_types <- aggregateAcrossCells(sce, id=colData(sce)[,c("customclassif", "sample")])
agg_cell_types <- agg_cell_types[,agg_cell_types$ncells >= 10]


y <- DGEList(counts(agg_cell_types), samples=colData(agg_cell_types))
keep <- filterByExpr(y, group=agg_cell_types$treatment)
y <- y[keep,]
mds_data <-cpm(y, log = TRUE) %>% 
  plotMDS()
plot <- data.frame(x = mds_data$x, 
                   y = mds_data$y, 
                   names = y$samples$seq_folder,
                   treatment = y$samples$treatment, 
                   cell_type = y$samples$customclassif) %>% 
  ggplot(aes(x = x, y = y, label = cell_type, colour = treatment))  +
  geom_text(hjust="inward", vjust="inward") +
  labs(x = "PC2", y = "PC1", title = "ind cell types")

pdf("plots/DEG/PCA_pseudobulk_ind_cell_types.pdf")
plot 
dev.off()

de.results_act <- pseudoBulkDGE(agg_cell_types, 
                            label=agg_cell_types$customclassif,
                            condition = agg_cell_types$treatment,
                            design=~treatment,
                            coef = "treatmentst"
)



cell_type_DEG_output <- lapply(de.results_act, make_pca_function)

PCAs <- lapply(cell_type_DEG_output, function(x)(return(x$PCA)))
Volcanos <- lapply(cell_type_DEG_output, function(x)(return(x$Volcano)))

pdf("plots/DEG/PCA_compiled.pdf", width = 14, height = 10)
ggarrange(plotlist =  PCAs, common.legend = T)
dev.off()


pdf("plots/DEG/Volcanos_compiled.pdf", width = 30, height = 22)
ggarrange(plotlist =Volcanos, common.legend = T)
dev.off()

save(cell_type_DEG_output, de.results_act, file = "data/DEG.RData")

read.table("data")