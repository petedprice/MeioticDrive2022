###LIBRARIES ----
library(Seurat)
library(tidyverse)
library(Matrix)
library(scales)
library(cowplot)
library(RCurl)
library(stringr)
library(ggpubr)
library(devtools)
devtools::install_github('satijalab/seurat-wrappers')
library(SeuratWrappers)
if (!requireNamespace("BiocManager"))
  install.packages("BiocManager")
install.packages('SCINA')
library('SCINA')
BiocManager::install()
BiocManager::install(c("monocle"))
# Next install a few more dependencies
BiocManager::install(c('DelayedArray', 'DelayedMatrixStats', 'org.Hs.eg.db', 'org.Mm.eg.db'))

load("indata/RData/integrated_seurat.RData")
load("outdata/RData/ortholog_table.RData")

markers <- readxl::read_excel("indata/markers/flyatlas_dros_markers.xlsx", col_names = TRUE) %>% 
  filter(Tissue == "testis") %>% 
  dplyr::select("Annotation", "Broad annotation", "Top Computed Markers (no CR, no CG)")

markers2 <- readxl::read_excel("indata/markers/elife2019/elife-47138-supp1-v1.xlsx", col_names = TRUE) %>% 
  dplyr::select("Gene", "Cluster")
markerslist3 <- lapply(unique(markers2$Cluster), function(x)(return(markers2$Gene[markers2$Cluster == x])))
names(markerslist3) <- unique(markers2$Cluster)

markerslist <- apply(markers, 1, function(x)(
  str_split(x[3], ", ")))
markerslist <- lapply(markerslist, unlist)
names(markerslist) <- markers$`Broad annotation`
markerslist2 <- list()
markerslist2$MRS <- markerslist[names(markerslist) == "male reproductive system"] %>% 
  unlist() %>% 
  unique()
markerslist2$MGC <- markerslist[names(markerslist) == "male germline cell"] %>% 
  unlist()%>% 
  unique()
markerslist2$SPC <- markerslist[names(markerslist) == "somatic precursor cell"] %>% 
  unlist()%>% 
  unique()

swap_names <- function(x, tab, srt){
  names <- unlist(lapply(x, function(g)(return(tab$TDel_GID[tab$Dros_GID == g])))) %>% 
    gsub(pattern = "gene-", replacement = "")
  return(intersect(names, rownames(srt)))
}



nmarkerslist <- sapply(markerslist3, swap_names, tab = ortholog_table, srt = seurat_integrated)
nmarkerslist <- nmarkerslist[-which(sapply(nmarkerslist, length) == 0)]
nmarkerslist <- nmarkerslist[-which(duplicated(nmarkerslist))]

check_subset <- function(cell, ml){
  matches <- lapply(ml, function(x)(return(prod(cell %in% x)))) %>% 
    unlist()
  if (sum(matches) > 1) {
    return(FALSE)
  } else {
    return(TRUE)
  }
}

keep_or_not <- lapply(nmarkerslist, check_subset, ml = nmarkerslist) %>% 
  unlist()
nmarkerslist <- nmarkerslist[keep_or_not]


#DimPlot(seurat_integrated, group.by = "ident")

scina.data <- as.data.frame(seurat_integrated@assays$integrated[,]) 

results = SCINA(scina.data, nmarkerslist, 
                max_iter = 1, convergence_n = 10, 
                convergence_rate = 0.999, sensitivity_cutoff = 0.9, 
                rm_overlap=FALSE, allow_unknown=TRUE, log_file='SCINA.log')

seurat_integrated$scina_labels <- results$cell_labels
d <- DimPlot(seurat_integrated, split.by= "scina_labels", group.by = "scina_labels")
pdf("plots/cell_types_UMAP.pdf", width = 28, height = 8)
d
dev.off()

DimPlot(seurat_integrated, group.by = "scina_labels", shape.by = "scina_labels", 
        raster = FALSE)

tibble(
  cluster = seurat_integrated$seurat_clusters,
  cell_type = seurat_integrated$scina_labels
) %>%
  group_by(cluster,cell_type) %>%
  count() %>%
  group_by(cluster) %>%
  mutate(
    percent=(100*n)/sum(n)
  ) %>%
  ungroup() %>%
  mutate(
    cluster=paste("Cluster",cluster)
  ) %>%
  ggplot(aes(x="",y=percent, fill=cell_type)) +
  geom_col(width=1) +
  coord_polar("y", start=0) +
  facet_wrap(vars(cluster)) +  
  theme(axis.text.x=element_blank()) +
  xlab(NULL) +
  ylab(NULL)


