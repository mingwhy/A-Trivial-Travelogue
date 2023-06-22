
###############################################################
## read in data
library(reticulate)
#Sys.setenv(RETICULATE_PYTHON = "/Users/mingyang/anaconda3/envs/RCFGL/bin/python")
#py_config()
sc=import('scanpy')  #py_install('scanpy')
#ad=sc$read_h5ad('/Users/mingyang/Documents/Data_worm_aing/ad_worm_aging.h5ad') #https://scanpy.readthedocs.io/en/stable/generated/scanpy.read_h5ad.html
ad=sc$read_h5ad('/Users/mingyang/Documents/Data_worm_aing/ad_worm_aging_umi.h5ad') #https://scanpy.readthedocs.io/en/stable/generated/scanpy.read_h5ad.html
cell.meta=ad$obs
gene.meta=ad$var

ad$X$shape # 47423 20305
class(ad$X)
mat=ad$X$toarray()
class(mat) #"matrix" "array" 
dim(mat) #  47423 20305 
dim(cell.meta) #47423

# create a SingleCellExperiment
library(zellkonverter)
library(SingleCellExperiment)
library(scater);library(scran)
library(ggplot2);library(gridExtra)
library(tidyverse);
library(SummarizedExperiment)
#https://www.bioconductor.org/packages/devel/bioc/vignettes/SingleCellExperiment/inst/doc/intro.html
sce <- SingleCellExperiment( list(counts=as(t(mat), "sparseMatrix")), #list(counts=mat)),
                             colData=cell.meta,
                             rowData=gene.meta);
sce
writeH5AD(sce, 'sce_worm.h5ad')
