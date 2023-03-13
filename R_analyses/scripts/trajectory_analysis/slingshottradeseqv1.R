devtools::install_github("statOmics/tradeSeq")

library(tradeSeq)
library(RColorBrewer)
library(SingleCellExperiment)
library(slingshot)
library(Seurat)
library(tidyverse)
library(devtools)
install_github('kstreet13/bioc2020trajectories')
library(pheatmap)
load("data/RData/integrated_seurat.RData")
# run cell-type IDing
siss <- subset(seurat_integrated, customclassif != 'Unknown' &
                 integrated_snn_res.0.4 %in% c('0',1,2,3,4,5,6,7,8,9,10,11))

dim(siss)
dim(seurat_integrated)

siss$new_clusters <- siss$customclassif
#siss$new_clusters[siss$customclassif %in% c("Mature spermatids", "Early spermatids", 
#                                            "Early spermatocytes")] <- "germ"
DimPlot(siss, group.by = 'new_clusters', reduction = 'umap', split.by = 'treatment')

sce_trad <- as.SingleCellExperiment(siss, assay = "RNA")

shuffle <- sample(ncol(sce_trad))
layout(matrix(1:2, nrow = 1))
par(mar = c(4.5,4,1,1))

plot(reducedDims(sce_trad)$UMAP[shuffle, ],
     asp = 1, pch = 16, xlab = "UMAP-1", ylab = "UMAP-2",
     col = alpha(c(1:2)[factor(colData(sce_trad)$treatment)][shuffle], alpha = .5))
legend("topright", pch = 16, col = 1:2, bty = "n", 
       legend = levels(factor(colData(sce_trad)$treatment)))
cols = 1:length(unique(colData(sce_trad)$new_clusters))
plot(reducedDims(sce_trad)$UMAP[shuffle, ], asp = 1, pch = 16, xlab = "UMAP-1", ylab = "UMAP-2", 
     col = alpha(cols[factor(colData(sce_trad)$new_clusters)][shuffle], alpha = .5))
legend("topright", pch = 16, col = cols, bty = "n", legend = levels(factor(colData(sce_trad)$new_clusters)))

scores <- bioc2020trajectories::imbalance_score(
  rd = reducedDims(sce_trad)$UMAP, 
  cl = colData(sce_trad)$treatment,
  k = 20, smooth = 40)


grad <- viridis::plasma(10, begin = 0, end = 1)
names(grad) <- levels(cut(scores$scaled_scores, breaks = 10))
plot(reducedDims(sce_trad)$UMAP, col = grad[cut(scores$scaled_scores, breaks = 10)],
     asp = 1, pch = 16, xlab = "UMAP-1", ylab = "UMAP-2", cex = .8)
legend("topleft", legend = names(grad), col = grad, pch = 16, bty = "n", cex = 2 / 3)


sce_trad <- slingshot(sce_trad, reducedDim = 'UMAP', clusterLabels = sce_trad$new_clusters,
                 start.clus = "Epithelial cells", approx_points = 500)

colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
plotcol <- colors[cut(sce_trad$slingPseudotime_1, breaks=100)]
plot(reducedDims(sce_trad)$UMAP, col = plotcol)
lines(SlingshotDataSet(sce_trad))

data.frame(cluster = sce_trad$new_clusters, pt = sce_trad$slingPseudotime_1, 
           treatment = sce_trad$treatment) %>% 
  ggplot(aes(x = pt, fill = cluster)) + geom_density(alpha = 0.2) + 
  facet_grid(.~treatment)

data.frame(cluster = sce_trad$new_clusters, pt = sce_trad$slingPseudotime_1, 
           treatment = sce_trad$treatment) %>% 
  ggplot(aes(x = pt, fill = treatment)) + geom_density(alpha = 0.2)


pseudotime <- slingPseudotime(sce_trad, na = FALSE)
cellWeights <- slingCurveWeights(sce_trad)
gois <- lapply(cell_type_DEG_output, function(x)(return(c(x$ST_bias_genes, x$SR_bias_genes)))) %>% 
  unlist()
gois <- c(gois, sample(rownames(sce), 81), "LOC119669221") %>% unique()
#gois <- cell_type_DEG_output$`Mature spermatids`$ST_bias_genes


icMat <- evaluateK(counts = as.matrix(assays(sce_trad)$counts),
                   pseudotime = pseudotime,
                   cellWeights = cellWeights,
                   conditions = factor(colData(sce_trad)$treatment),
                   nGenes = 300,
                   k = 3:7, parallel=F)


sce_trad_ss <- fitGAM(counts = counts(sce_trad)[gois,], 
              pseudotime = pseudotime, 
              cellWeights = cellWeights,
              conditions = factor(colData(sce_trad)$treatment),
              nknots = 5, parallel=F)

rowData(sce_trad_ss)$assocRes <- associationTest(sce_trad_ss, lineages = TRUE, l2fc = log2(0.5))
assocRes <- rowData(sce_trad_ss)$assocRes
assocRes <- assocRes[is.na(assocRes$waldStat) == F,]
sr_genes <-  rownames(assocRes)[
  which(p.adjust(assocRes$pvalue_lineage1_conditionsr, "fdr") <= 0.05)
]
st_genes <-  rownames(assocRes)[
  which(p.adjust(assocRes$pvalue_lineage1_conditionst, "fdr") <= 0.05)
]

yhatSmooth <- predictSmooth(sce_trad_ss, gene = st_genes, nPoints = 50, tidy = FALSE)
heatSmooth <- pheatmap(t(scale(t(yhatSmooth[, 1:50]))),
                       cluster_cols = FALSE,
                       show_rownames = T, 
                       show_colnames = T)
 

View(filter(ortholog_table, TDel_GID %in% gois))
plotSmoothers(sce_trad_ss, assays(sce_trad_ss)$counts, gene = "LOC119669221", alpha = 1, border = TRUE) + ggtitle("RpL40")


condRes <- conditionTest(sce_trad_ss, l2fc = log2(2))
condRes$padj <- p.adjust(condRes$pvalue, "fdr")
mean(condRes$padj <= 0.05, na.rm = TRUE)

conditionGenes <- rownames(condRes)[condRes$padj <= 0.05]
conditionGenes <- conditionGenes[!is.na(conditionGenes)]

oo <- order(condRes$waldStat, decreasing = TRUE)

# most significant gene
plotSmoothers(sce_trad_ss, assays(sce_trad_ss)$counts,
              gene = rownames(assays(sce_trad_ss)$counts)[oo[1]],
              alpha = 1, border = TRUE)

# least significant gene
plotSmoothers(sce_trad_ss, assays(sce_trad_ss)$counts,
              gene = rownames(assays(sce_trad_ss)$counts)[oo[nrow(sce_trad_ss)]],
              alpha = 1, border = TRUE)






### based on mean smoother
yhatSmooth <- predictSmooth(sce_trad_ss, gene = conditionGenes, nPoints = 50, tidy = FALSE)
yhatSmoothScaled <- t(scale(t(yhatSmooth)))
heatSmooth_TGF <- pheatmap(yhatSmoothScaled[, 51:100],
                           cluster_cols = FALSE,
                           show_rownames = F, show_colnames = T, main = "Standard", legend = FALSE,
                           silent = TRUE
)

matchingHeatmap_mock <- pheatmap(yhatSmoothScaled[heatSmooth_TGF$tree_row$order, 1:50],
                                 cluster_cols = FALSE, cluster_rows = FALSE,
                                 show_rownames = T, show_colnames = T, main = "Drive",
                                 legend = FALSE, silent = TRUE
)

grid.arrange(heatSmooth_TGF[[4]], matchingHeatmap_mock[[4]], ncol = 2)
