for this to work you'll need the following packages installed for R
library(Seurat)
library(tidyverse)
library(Matrix)
library(scales)
library(cowplot)
library(RCurl)
library(stringr)
library(ggpubr)


Functions what they do and what not to use 

### NORMALISATION ----

SCTransform: 
- Normalisation method that in a single command replaces NormalizeData(), ScaleData(), and FindVariableFeatures().
- Transformed data will be available in the SCT assay, which is set as the default after running sctransform





#### Identifying cell markers for clusters ----
FindAllMarkers
- Finds markers that distinguish cell type clusters 




#### Functions not to use ----
ScaleData, NormaliseData, FindVariableFeatures: 
- do not use if you have used SCTransform

