
library(zellkonverter)
library(SingleCellExperiment)
library(scater);library(scran)
library(ggplot2);library(gridExtra);
library(tidyverse);
library(zellkonverter)
library(SummarizedExperiment)
library(viridis)
library(ggpubr)


###############################################################
## read in data
sce=readH5AD('~/Documents/Data_AgingFlyCellAtlas/adata_headBody_S_v1.0.h5ad') # 22966 31001 
#sce=readH5AD('~/Documents/Data_AgingFlyCellAtlas/adata_head_S_v1.0.h5ad') # 22966 31001 
#sce=readH5AD('~/Documents/Data_AgingFlyCellAtlas/adata_Body_S_v1.0.h5ad') # 15992 276273
sce #SingleCellExperiment,  
assayNames(sce) #'X'
cell.meta=colData(sce)
colnames(cell.meta)
table(cell.meta$tissue)
#body   head 
#276273 289981 
table(cell.meta$sex,cell.meta$age)
head(cell.meta)
#table(cell.meta$afca_annotation_broad,cell.meta$afca_annotation)
length(table(cell.meta$afca_annotation)) # 163 cell types, body=73
length(table(cell.meta$afca_annotation_broad)) #17 broad cell classes
#grep('head',cell.meta$afca_annotation,value = T)

##################################################################
## reverse back to umi.count matrix
mat=assay(sce,'X')
mat[1:3,1:3]

# test on 10 cells
test.mat=mat[,1:10]
cell.meta[1:10,]$total_counts
emat<-(exp(test.mat)-1)
colSums(emat) #all the same size, 1e4
tmp=emat[,1]/10^4*cell.meta[1:10,]$total_counts[1]
table(tmp)

# apply to all cells (may take time and exhaust memory)
mat=exp(assay(sce,'X'))-1
lib.size=cell.meta$total_counts
umi.mat=((exp(mat)-1)/10^4 ) * rep(lib.size, rep(nrow(mat),length(lib.size))) 
#https://stackoverflow.com/questions/17080099/fastest-way-to-multiply-matrix-columns-with-vector-elements-in-r



