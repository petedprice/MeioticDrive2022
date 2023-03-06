# load libraries and functions
library(dplyr)
library(Seurat)
library(HGNChelper)
library(openxlsx)
library(SCINA)
library(ggpubr)
library(scater)
library(scuttle)
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R"); source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")
# load gene set preparation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
# load cell type annotation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")
load("data/RData/integrated_seurat.RData")
load("outdata/RData/ortholog_table.RData")


#### FUNCTIONS ----
swap_names <- function(x, tab, srt){
  names <- unlist(lapply(x, function(g)(return(tab$TDel_GID[tab$Dros_GID == g])))) %>% 
    gsub(pattern = "gene-", replacement = "")
  return(intersect(names, rownames(srt)))
}

check_subset <- function(cell, ml){
  matches <- lapply(ml, function(x)(return(prod(cell %in% x)))) %>% 
    unlist()
  if (sum(matches) > 1) {
    return(FALSE)
  } else {
    return(TRUE)
  }
}

DEG_func <- function(c, seurat_obj, ortholog_table, s = 'customclassif'){
  if (s == "customclassif"){
    ss_obj <- subset(seurat_obj, customclassif == c)
  } else if (s == "scina_labels") {
    ss_obj <- subset(seurat_obj, scina_labels == c)
  }
  if (length(unique(ss_obj$treatment)) < 2){
    return(NULL)
  } else {
    t1 <- unique(seurat_obj$treatment)[1]
    t2 <- unique(seurat_obj$treatment)[2]
    DEG_test <- FindMarkers(ss_obj, ident.1 = t1, 
                            ident.2 = t2, 
                            logfc.threshold = log(2), 
                            min.pct = 0.5)
    DEG_test$DEG_comp <- c
    DEG_test$gene <- rownames(DEG_test)
    rownames(DEG_test) <- NULL
    return(DEG_test)
  }
}


markers <- readxl::read_excel("indata/markers/elife2019/elife-47138-supp1-v1.xlsx", col_names = TRUE) %>% 
  dplyr::select("Gene", "Cluster")


markerslist <- lapply(unique(markers$Cluster), function(x)(return(markers$Gene[markers$Cluster == x])))
names(markerslist) <- unique(markers$Cluster)
nmarkerslist <- sapply(markerslist, swap_names, tab = ortholog_table, srt = seurat_integrated)
gs_list <- list()
gs_list$gs_positive <- nmarkerslist


#### SC_TYPE MARKERS ----
# get cell-type by cell matrix
es.max = sctype_score(scRNAseqData = seurat_integrated[["integrated"]]@scale.data, scaled = TRUE, 
                      gs = gs_list$gs_positive, gs2 = gs_list$gs_negative) 

# merge by cluster
cL_resutls = do.call("rbind", lapply(unique(seurat_integrated@meta.data$seurat_clusters), function(cl){
  es.max.cl = sort(rowSums(es.max[ ,rownames(seurat_integrated@meta.data[seurat_integrated@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
  head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(seurat_integrated@meta.data$seurat_clusters==cl)), 10)
}))
sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  

# set low-confident (low ScType score) clusters to "unknown"
sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"
print(sctype_scores[,1:3])


seurat_integrated@meta.data$customclassif = ""
for(j in unique(sctype_scores$cluster)){
  cl_type = sctype_scores[sctype_scores$cluster==j,]; 
  seurat_integrated@meta.data$customclassif[seurat_integrated@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
}

#### SCINA ----
keep_or_not <- lapply(nmarkerslist, check_subset, ml = nmarkerslist) %>% 
  unlist()
SCINA_markerslist <- nmarkerslist[keep_or_not]


scina.data <- as.data.frame(seurat_integrated@assays$integrated[,]) 

results = SCINA(scina.data, SCINA_markerslist, 
                max_iter = 1, convergence_n = 10, 
                convergence_rate = 0.999, sensitivity_cutoff = 0.9, 
                rm_overlap=FALSE, allow_unknown=TRUE, log_file='SCINA.log')
seurat_integrated$scina_labels <- results$cell_labels


##### BIOCONDUCTOR DGESEQ2 -----
library(edgeR)
library(statmod)
library(scran)
library(ggpubr)
library(gridExtra)
sce <- as.SingleCellExperiment(seurat_integrated, assay = "RNA")

#### BULK RNASEQ ACROSS WHOLE TISSUE ----
whole_tissue <- aggregateAcrossCells(sce, id=colData(sce)[,c("sample")])
de.results_wt <- pseudoBulkDGE(whole_tissue, 
                            label='bulk',
                            condition = whole_tissue$treatment,
                            design=~treatment,
                            coef = "treatmentst",
)

mds_data <- metadata(de.results_wt$bulk)$y %>% 
  cpm(log = T) %>% 
  plotMDS()

plot <- data.frame(x = mds_data$x, 
                        y = mds_data$y, 
                        names = colnames(mds_data$distance.matrix.squared),
                        treatment = metadata(de.results_wt$bulk)$fit$samples$treatment) %>% 
  ggplot(aes(x = x, y = y, fill = names, label = names, colour = treatment))  +
  geom_text(hjust="inward", vjust="inward") +
  labs(x = "PC2", y = "PC1", title = "whole tissue pseudobulk")

pdf("plots/DEG/PCA_pseudobulk_WT.pdf")
plot 
dev.off()


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
                            coef = "treatmentst",
)

make_pca_function <- function(x){
  md <- metadata(x)
  mds_data <- md$y %>% 
    cpm(log = TRUE) %>% 
    plotMDS()
  
  plot1 <- data.frame(x = mds_data$x, 
                     y = mds_data$y, 
                     names = md$y$samples$seq_folder,
                     treatment = md$fit$samples$treatment) %>% 
    ggplot(aes(x = x, y = y, fill = names, label = names, colour = treatment))  +
    geom_text(hjust="inward", vjust="inward") +
    labs(x = "PC2", y = "PC1", title = md$y$samples$customclassif[1]) 
  pdf(paste("plots/DEG/PCA_", md$y$samples$customclassif[1], ".pdf", sep = ""))
  print(plot1) 
  dev.off()
  x <- as.data.frame(x)
  x$gene <- rownames(x)
  plot2 <- x %>% 
    ggplot(aes(x = logFC, y = -log10(PValue), label = gene,
                 colour = (abs(logFC) > 2 &  PValue < 0.05))) + 
    geom_point() + geom_text(data=subset(x, -log10(PValue) > 7.5 & abs(logFC) > 2),
                             aes(logFC,-log10(PValue),label=gene), 
                             hjust="inward", vjust="inward") +
    labs(y = "-log10 Pvalue", x = "log2 fold-change", 
         title = md$y$samples$customclassif[1]) + 
    theme(legend.position = "none")  
  pdf(paste("plots/DEG/DEG_", md$y$samples$customclassif[1], ".pdf", sep = ""))
  print(plot2) 
  dev.off()
  return(list(plot1,plot2))

}

plots_list <- lapply(de.results_act, make_pca_function)

pdf("plots/DEG/PCA_compiled.pdf", width = 14, height = 10)
ggarrange(plotlist = plots_list, common.legend = T)
dev.off()

agt_cov <- cor(counts(agg_cell_types))
colnames(agt_cov) <- paste(agg_cell_types$sample, agg_cell_types$scina_labels)
rownames(agt_cov) <- paste(agg_cell_types$sample, agg_cell_types$scina_labels)

dev.off()
pdf("plots/heatmap_scina.pdf")
heatmap(agt_cov)
dev.off()
heatmap.2(agt_cov)


summed
label = "Cyst"
#current <- summed[,label==summed$




y <- DGEList(counts(summed.filt), samples=colData(summed.filt))
discarded <- summed$ncells < 10
y <- y[,!discarded]
summary(discarded)

keep <- filterByExpr(y, group=summed$treatment)
y <- y[keep,]
summary(keep)

y <- calcNormFactors(y)
y$samples

par(mfrow=c(2,4))
for (i in seq_len(ncol(y))) {
  plotMD(y, column=i)
}
dev.off()
plotMDS(cpm(y, log=TRUE), 
        col= as.numeric(factor(y$samples$treatment)))


design <- model.matrix(~factor(treatment), y$samples)
y <- estimateDisp(y, design)
summary(y$trended.dispersion)
plotBCV(y)

fit <- glmQLFit(y, design, robust=TRUE)
summary(fit$var.prior)

