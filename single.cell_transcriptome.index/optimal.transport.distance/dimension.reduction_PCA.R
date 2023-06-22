

library(SingleCellExperiment)
library(SummarizedExperiment)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(Seurat)

# for fast PCA on large matrix: https://github.com/bwlewis/irlba
library(irlba)
#setClass("scaled_matrix", contains="matrix", slots=c(scale="numeric"))
#setMethod("%*%", signature(x="scaled_matrix", y="numeric"), function(x ,y) x@.Data %*% (y / x@scale))
#setMethod("%*%", signature(x="numeric", y="scaled_matrix"), function(x ,y) (x %*% y@.Data) / y@scale)
########################################################
## read in dataset
#BiocManager::install("zellkonverter")
library(zellkonverter)
library(SummarizedExperiment)

## read in data
sce<-readH5AD('~/Documents/Data_mouse_aging_atlas//TMS.gene.data_final/tabula-muris-senis-facs-official-raw-obj.h5ad');       
class(assay(sce))
sce

dim(assay(sce))
df.expr=assay(sce)
df.expr[1:10,1:10] #umi count
dim(df.expr) #22966 110824

## sample meta information
cell.meta=colData(sce)
head(cell.meta)
dim(cell.meta) #110824     13
colnames(cell.meta)


table(cell.meta$age) #6 age group: 1, 3, 18, 21, 24,30month. fac: 3,18,21,24 month
table(cell.meta$tissue) #23 tissue
length(table(cell.meta$cell_ontology_class)) #120 cell types
unique(paste(cell.meta$tissue,cell.meta$cell_ontology_class)) #207 unique #consistent with `elife-62293-supp1-v2.xlsx`
x=as.data.frame(cell.meta) %>% group_by(tissue) %>% summarise(n.cell.type=length(unique(cell_ontology_class)))
x$n.cell.type

tissue_cell.type=paste(cell.meta$tissue,cell.meta$cell_ontology_class,sep=':')
cell.meta$tissue_cell.type=tissue_cell.type

table(cell.meta[cell.meta$age=='21m',]$tissue_cell.type) #remove these cells
i=cell.meta$age!='21m'
cell.meta=cell.meta[i,]
df.expr=df.expr[,i]
cell.meta$age=as.character(cell.meta$age)
cell.meta$raw.age=cell.meta$age;

cell.meta[cell.meta$age=='3m',]$age<-'young';
cell.meta[cell.meta$age!='young',]$age<-'old';
dim(cell.meta) # 110096     14
dim(df.expr) # 22966 110096

#######################################################################
##filter tc to cell.types which contain>=100 cells in both young and old
x=as.data.frame(cell.meta) %>% group_by(tissue_cell.type,age) %>% summarise(n=n())
tcs=names(which(table(x[x$n>=100,]$tissue_cell.type)==2)) #76 tc
#x=as.data.frame(cell.meta) %>% group_by(tissue_cell.type,raw.age) %>% summarise(n=n())
#tcs=names(which(table(x[x$n>=20,]$tissue_cell.type)==3)) #115 tc

df.expr=df.expr[,cell.meta$tissue_cell.type %in% tcs]
cell.meta=cell.meta[cell.meta$tissue_cell.type %in% tcs,]
dim(df.expr);dim(cell.meta) # 22966 100817 for 115tc; 22966 94573 for 76 tc

###############################################################################################
## calculate optimal.transport distance between cell populations
all.tcs=sort(unique(cell.meta$tissue_cell.type))

#for(tc in all.tcs[1:4]){
for(tc in all.tcs){
  expression_matrix=df.expr[,cell.meta$tissue_cell.type==tc]
  
  # gene filter
  geneFilter <- apply(expression_matrix,1,function(x){
    sum(x >= 3) >= 10
  })
  expression_matrix <- expression_matrix[geneFilter, ]
  cell_metadata<-cell.meta[cell.meta$tissue_cell.type==tc,]
  
  
  ## create single cell obj to perform PCA
  obj=CreateSeuratObject(expression_matrix)
  obj=NormalizeData(obj)
  obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000)
  var.genes=VariableFeatures(obj)
  
  #obj <- ScaleData(obj, features = rownames(obj)) #for PCA
  obj <- ScaleData(obj, features = var.genes)
  ncell=ncol(obj);
  obj <- RunPCA(obj, features = var.genes,npcs = 200,verbose = F)
  dim(obj@reductions$pca@cell.embeddings) #ncell x #PC
  
  # PC variance explained (https://github.com/satijalab/seurat/issues/982)
  mat <- Seurat::GetAssayData(obj, assay = "RNA", slot = "scale.data")
  mat=mat[var.genes,] 
  pca <- obj[["pca"]]
  # Get the total variance:
  total_variance <- sum(matrixStats::rowVars(mat))
  total_variance; #gene.var has already been scaled, show equal to #gene
  eigValues = (pca@stdev)^2  ## EigenValues
  varExplained = eigValues / total_variance
  cat(tc,',varExplained',sum(varExplained),'\n')
  
  if(F){
    A=obj@assays$RNA@scale.data;
    dim(A) #ngene x ncell
    n1=sum(cell_metadata$age=='young')
    n2=sum(cell_metadata$age=='old')
    S=irlba(A, nv=min(n1,n2,200))
    #length(S$d) #eigenvalue
    dim(S$v) #ncell X pc1..100
    expr.pc=S$v
  }
 
 
  cat('tc',tc,'is done\n')
}



