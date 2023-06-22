
# library(SeuratDisk) didn't work
# https://github.com/satijalab/seurat/issues/3414
if(F){
  Sys.getenv("GITHUB_PAT")
  Sys.unsetenv("GITHUB_PAT")
  Sys.getenv("GITHUB_PAT")
  remotes::install_github("mojaveazure/seurat-disk")
  library(SeuratDisk)
}


BiocManager::install("HDF5Array")
BiocManager::install("rhdf5")
# library(zellkonverter) works
#https://rdrr.io/github/theislab/zellkonverter/man/readH5AD.html
#BiocManager::install("zellkonverter")
library(zellkonverter)
library(SummarizedExperiment)
sce<-readH5AD('./GSE127832_wing/GSM3639558_DGE_Drop-Seq_wing.disc_rep1_Lib1.h5ad',use_hdf5 = TRUE)
class(assay(sce))
#[1] "DelayedMatrix"
#attr(,"package")
#[1] "DelayedArray"
sce

dim(assay(sce))
df.expr=assay(sce)
df.expr[1:3,1:3]
dim(df.expr) #10941   150

#https://bioconductor.org/packages/release/bioc/html/SingleCellExperiment.html

# Interoperability between single-cell object formats
#https://satijalab.org/seurat/articles/conversion_vignette.html
library(scater)
library(Seurat)
library(SeuratDisk)
library(SeuratData)
library(patchwork)

#########################################################################################################
library(zellkonverter)
library(SingleCellExperiment)
library(scater);library(scran)
library(ggplot2);library(gridExtra);
library(tidyverse);
library(zellkonverter)
library(SummarizedExperiment)
library(viridis)
library(ggpubr)

## read in data (https://hongjielilab.shinyapps.io/AFCA/)
sce=readH5AD('adata_headBody_S_v1.0.h5ad') # 22966 31001  
sce #SingleCellExperiment,  15992 566254 
cell.meta=colData(sce)
colnames(cell.meta)
table(cell.meta$sex,cell.meta$age)
head(cell.meta)
#table(cell.meta$afca_annotation_broad,cell.meta$afca_annotation)
length(table(cell.meta$afca_annotation)) # 163 cell types



