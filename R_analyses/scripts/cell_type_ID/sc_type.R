# load libraries and functions
library(dplyr)
library(Seurat)
library(HGNChelper)
library(openxlsx)
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R"); source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")


# load gene set preparation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
# load cell type annotation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")


load("indata/RData/integrated_seurat.RData")
load("outdata/RData/ortholog_table.RData")
seurat_integrated <- RunPCA(object = seurat_integrated)
seurat_integrated <- RunUMAP(seurat_integrated, 
                             dims = 1:40,
                             reduction = "pca")

seurat_integrated <- FindNeighbors(object = seurat_integrated, 
                                   dims = 1:40)
seurat_integrated <- FindClusters(object = seurat_integrated,
                                  resolution = 0.4)


markers <- readxl::read_excel("indata/markers/elife2019/elife-47138-supp1-v1.xlsx", col_names = TRUE) %>% 
  dplyr::select("Gene", "Cluster")
markerslist <- lapply(unique(markers$Cluster), function(x)(return(markers$Gene[markers$Cluster == x])))
names(markerslist) <- unique(markers$Cluster)
nmarkerslist <- sapply(markerslist, swap_names, tab = ortholog_table, srt = seurat_integrated)
gs_list <- list()
gs_list$gs_positive <- nmarkerslist





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

d <- DimPlot(seurat_integrated, reduction = "umap", label = TRUE, repel = TRUE, 
             group.by = 'customclassif', split.by = 'treatment')        
pdf("plots/Cell_types/sc_type/sc_type_cell_types_UMAP.pdf", width = 22, height = 8)
d
dev.off()


