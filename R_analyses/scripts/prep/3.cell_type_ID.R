#!/usr/bin/env Rscript
##### LIBRARIES -------
library("optparse")

#SC_TYPE FUNCTIONS
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")

###### SETTING UP INPUT COMMANDS ----
option_list = list(
  make_option(c("-d", "--path_to_seurat_object"), type="character", default=".", 
              help="path to where you have the stored your seurat object", metavar="character"),
  make_option(c("-o", "--output_path"), type="character", default=".", 
              help="where you want to save your output plots and RData files", metavar="character"),
  make_option(c("-t", "--threads"), type="numeric", default=1, 
              help="number of threads for parallelising", metavar="numeric"),
  make_option(c("-l", "--ortholog_table"), type="character", default="/home/bop20pp/software/MeioticDrive2022/R_analyses/data/ortholog_table.txt", 
              help="path to dataframe containing ortholog information", metavar="character"),
  make_option(c("-s", "--marker_source"), type="character", default="testis.Cluster", 
              help="which markers to use for classifying cell types. Column name from ortholog_table.", metavar="character")
  
)


opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if (is.null(opt$path_to_seurat_object)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file)", call.=FALSE)
}

##### FUNCTIONS ----
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


#### DETLETE -----
opt$path_to_seurat_object <- "data/RData/integrated_seurat.RData"
opt$ortholog_table <- "data/ortholog_table.txt"
###################

##### LOADING DATA ----
load(opt$path_to_seurat_object)
ortholog_table <- read.table(opt$ortholog_table)
marker_source <- opt$marker_source

clusters <- unique(ortholog_table[,marker_source])
clusters <- clusters[is.na(clusters) == F]

markerslist <- lapply(clusters, function(x)(return(ortholog_table$TDel_GID[ortholog_table[,marker_source] == x])))
markerslist <- lapply(markerslist, function(x)(return(x[which(is.na(x) == F)])))
names(markerslist) <- clusters

gs_list <- list()
gs_list$gs_positive <- markerslist

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

DimPlot(seurat_integrated, group.by = 'customclassif', 
        label = T)


#### SCINA ----
keep_or_not <- lapply(markerslist, check_subset, ml = markerslist) %>% 
  unlist()
SCINA_markerslist <- markerslist[keep_or_not]


scina.data <- as.data.frame(seurat_integrated@assays$integrated[,]) 

results = SCINA(scina.data, SCINA_markerslist, 
                max_iter = 1, convergence_n = 10, 
                convergence_rate = 0.999, sensitivity_cutoff = 0.9, 
                rm_overlap=FALSE, allow_unknown=TRUE, log_file='SCINA.log')
seurat_integrated$scina_labels <- results$cell_labels
DimPlot(seurat_integrated, group.by = 'scina_labels', label = T)








