BiocManager::install("slingshot")
library(slingshot)

sce <- as.SingleCellExperiment(seurat_integrated, assay = "RNA")
sce <- slingshot(sce, reducedDim = 'UMAP', clusterLabels = colData(sce)$seurat_clusters,
                 start.clus = 'inner', approx_points = 150)

ks.test(slingPseudotime(sce)[colData(sce)$treatment == "sr", 1],
        slingPseudotime(sce)[colData(sce)$treatment == "st", 1])


plot(reducedDims(sce)$UMAP, col = sce$slingPseudotime_5,
     asp = 1, pch = 16, xlab = "UMAP-1", ylab = "UMAP-2", cex = .8)
legend("topleft", legend = names(grad), col = grad, pch = 16, bty = "n", cex = 2 / 3)


colData(sce)
assocRes <- rowData(sce)$
mockGenes <-  rownames(assocRes)[
  which(p.adjust(assocRes$pvalue_lineage1_conditionMock, "fdr") <= 0.05)
]
tgfbGenes <-  rownames(assocRes)[
  which(p.adjust(assocRes$, "fdr") <= 0.05)
]

length(mockGenes)
