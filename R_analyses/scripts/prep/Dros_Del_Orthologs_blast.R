library(GenomicFeatures)
library("readxl")
library(stringr)
library(dplyr)
library(tidyr)

dros_txdb_coords <- makeTxDbFromGFF("indata/ref_files/dros.gtf", format = 'gff')
k <- keys(dros_txdb_coords, keytype = "GENEID")
dros_txdf <- AnnotationDbi::select(dros_txdb_coords, keys = k,  columns = c("TXNAME", "TXCHROM"), keytype = "GENEID")
colnames(dros_txdf) <- c("Dros_GID", "Dros_TXN", "Dros_CHR")



TDel_txdb_coords <- makeTxDbFromGFF("indata/ref_files/stalkie.gtf", format = 'gtf')
k <- keys(TDel_txdb_coords, keytype = "GENEID")
TDel_txdf <- AnnotationDbi::select(TDel_txdb_coords, keys = k,  columns = c("TXNAME", "TXCHROM"), keytype = "GENEID")
colnames(TDel_txdf) <- c("TDel_GID", "TDel_TXN", "TDel_CHR")

orths_t0 <- read.table("indata/orthologs/top_hits/top_hits.txt", 
                        sep = ',', header = FALSE)

colnames(orths_t0) <- c("TDel_TXN", "Dros_TXN")
orths_t0[,1] <- str_split(orths_t0[,1], "stalkie_", simplify = TRUE)[,2]
orths_t0[,2] <- str_split(orths_t0[,2], "dros_rna-", simplify = TRUE)[,2]
orths_t0[,2] <- str_split(orths_t0[,2], "[)]", simplify = TRUE)[,1]
orths_t1 <- merge(orths_t0, TDel_txdf, by.x = 'TDel_TXN', by.y = 'TDel_TXN')
orths_t2 <- merge(orths_t1, dros_txdf, by.x = 'Dros_TXN', 
               by.y = 'Dros_TXN')
##### NEXT THING TO DO IS LINK IN THE MARKER GENES 
markers <- read_excel("indata/markers/elife2019/elife-47138-supp1-v1.xlsx")
markers2 <- read_excel("indata/markers/flyatlas_dros_markers.xlsx") %>% 
  filter(Tissue == "testis")

markers2list <- apply(markers2, 1, function(x)(return(
  data.frame(
    Cluster = rep(x[2]), 
    gene = t(str_split(x[15], ",", simplify = TRUE)))
  )))

markers2df <- bind_rows(markers2list)
orthologs_testis <- merge(orths_t2, markers, by.x = 'Dros_GID', by.y = 'Gene', all.x = TRUE)[,c(3,4,5,6,1,2,7:12)]
#or 
intersect(unique(markers2df$gene), orths_t2$Dros_GID)
orthologs_dros_atlas <- merge(orths_t2, markers2df, by.x = 'Dros_GID', by.y = 'gene', all.x = TRUE)[,c(3,4,5,6,1,2,7)]
save(orthologs_testis, orthologs_dros_atlas, file = "outdata/RData/orthologs.RData")
