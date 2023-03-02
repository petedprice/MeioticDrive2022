library(tradeSeq)
library(RColorBrewer)
library(SingleCellExperiment)
library(slingshot)

sce <- as.SingleCellExperiment(seurat_integrated, assay = "RNA")

shuffle <- sample(ncol(sce))
layout(matrix(1:2, nrow = 1))
par(mar = c(4.5,4,1,1))

plot(reducedDims(sce)$UMAP[shuffle, ],
     asp = 1, pch = 16, xlab = "UMAP-1", ylab = "UMAP-2",
     col = alpha(c(1:10)[factor(colData(sce)$customclassif)][shuffle], alpha = .5))
legend("topright", pch = 16, col = 1:10, bty = "n", 
       legend = levels(factor(colData(sce)$customclassif)))
plot(reducedDims(sce)$UMAP[shuffle, ], asp = 1, pch = 16, xlab = "UMAP-1", ylab = "UMAP-2", 
     col = alpha(c(3, 4)[factor(colData(sce)$treatment)][shuffle], alpha = .5))
legend("topright", pch = 16, col = 3:4, bty = "n", legend = levels(factor(colData(sce)$treatment)))


sce <- slingshot(sce, reducedDim = 'UMAP', clusterLabels = colData(sce)$customclassif,
                 start.clus = "GSC, Early spermatogonia")

colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
plotcol <- colors[cut(sce$slingPseudotime_1, breaks=100)]

plot(reducedDims(sce)$UMAP, col = plotcol, pch=16, asp = 1)
lines(SlingshotDataSet(sce), lwd=2, col='black')


sce <- fitGAM(sce)
ATres <- associationTest(sce)
ti

counts <- seurat_integrated@assays$RNA@counts
sds <- SlingshotDataSet(sce)

icMat <- evaluateK(counts = counts, sds = sds, k = 3:10, 
                   nGenes = 200, verbose = T)

set.seed(7)
pseudotime <- slingPseudotime(sds, na = FALSE)
cellWeights <- slingCurveWeights(sds)

sce <- fitGAM(counts = counts, pseudotime = pseudotime, cellWeights = cellWeights,
              nknots = 6, verbose = FALSE)

table(rowData(sce)$tradeSeq$converged)

assoRes <- associationTest(sce)
head(assoRes)

startRes <- startVsEndTest(sce)

oStart <- order(startRes$waldStat, decreasing = TRUE)
sigGeneStart <- names(sce)[oStart[3]]
plotSmoothers(sce, counts, gene = sigGeneStart)








pal <- c(RColorBrewer::brewer.pal(9, "Set1"), RColorBrewer::brewer.pal(8, "Set2"))


dimred <- seurat_integrated@reductions$pca@cell.embeddings
clustering <- seurat_integrated$seurat_clusters
map <- setNames(1:length(unique(clustering)), unique(clustering))
clustering <- map[clustering] %>% as.numeric()

counts <- as.matrix(seurat_integrated@assays$RNA@counts[seurat_integrated@assays$RNA@var.features, ])
lineages <- SlingshotDataSet(getLineages(data = dimred, clusterLabels = clustering, 
                                         start.clus = 6))
plot(dimred[, 1:2], col = pal[as.numeric(clustering)], cex = 0.5, pch = 16)
lines(lineages, lwd = 3, col = "black")
for (i in unique(clustering)) {
  text(mean(dimred[clustering == i, 1]), mean(dimred[clustering == i, 2]), labels = i, font = 2, col = 'white')
}

curves <- getCurves(lineages, approx_points = 300, thresh = 0.01, stretch = 0.8, allow.breaks = FALSE, shrink = 0.99)
plot(dimred, col = pal[clustering], asp = 1, pch = 16)
lines(SlingshotDataSet(curves), lwd = 3, col = "black")

pdf("del2.pdf")
DoHeatmap(seurat_integrated, features = seurat_integrated@assays$integrated@var.features[1:200], 
          group.by = 'scina_labels')
dev.off()

