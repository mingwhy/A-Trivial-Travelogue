
library(zellkonverter)
library(SingleCellExperiment)
library(scater);library(scran)
library(ggplot2);library(gridExtra)
library(tidyverse);
library(SummarizedExperiment)
library(ggpointdensity)
library(viridis)
library(patchwork); #for plot_annotation
library(scMerge)

#########################################################
## subsample cell number to be equal in both age gropus
## read in data (a "SingleCellExperiment" object)
if(!file.exists('select.tc_subsampled.h5ad')){
  sce=readH5AD('select.tc.h5ad') # 22966 31001 
  sce=sce[,sce$age %in% c('3m','24m')]
  sce #22966 24443 
  
  (pick.cell.types=as.character(unique(sce$tissue_cell.type))) #39 tc
  
  ## sub-sample cell: min(3,18,24), if 18m contain<50cell, discard this age group
  cell.meta=colData(sce)
  as.data.frame(cell.meta) %>% group_by(tissue_cell.type,age,sex) %>% summarise(n=n())
  
  min.ncell=50
  (tcs=unique(sce$tissue_cell.type))
  sce.subs<-lapply(tcs,function(tc){
    sce.tmp=sce[,sce$tissue_cell.type==tc]
    meta.tmp=colData(sce.tmp)
    x=names(which(table(meta.tmp$age)>=min.ncell))
    sce.tmp=sce.tmp[,sce.tmp$age %in% x]
    sce.tmp$age=factor(sce.tmp$age,levels=x)
    meta.tmp=colData(sce.tmp)
    n=min(table(meta.tmp$age))
    cat(as.character(tc),x,n,'\n')
    tmp=lapply(levels(sce.tmp$age),function(age.i){
      x=sce.tmp[,sce.tmp$age==age.i]
      set.seed(1224)
      i=sample(1:ncol(x),n,replace = F)
      x[,i]
    })
    sce.sub<-Reduce(`cbind`,tmp)
    sce.sub
  })
  lapply(sce.subs,dim)
  sce<-Reduce(cbind,sce.subs)
  assayNames(sce)[1] <- "counts" 
  sce #22966 17442
  writeH5AD(sce, 'select.tc_subsampled.h5ad')
}

sce=readH5AD('select.tc_subsampled.h5ad')
(pick.cell.types=as.character(unique(sce$tissue_cell.type))) #39 tc
table(sce$age,sce$tissue_cell.type) #same #cell for each tc across age groups

############################################################
## filter gene: expr in 10% cell per cell type (in all ages)
if(T){
  sce.filter<-lapply(pick.cell.types,function(tc){
    tmp=sce[,sce$tissue_cell.type==tc]
    mat=assay(tmp,'counts')
    #i=which(Matrix::rowSums(mat>0)> ncol(mat)*0.1)
    i=which(Matrix::rowSums(mat>0)!=0)
    tmp[i,]
  })
  names(sce.filter)<-pick.cell.types;
}

############################################################
## subsample cell read count
library(countland) #https://github.com/shchurch/countland/blob/main/tutorials_and_vignettes/R_tutorials_and_vignettes/vignette-tutorial.Rmd

set.seed(84095) # set random seed for reproducibility
minCellCounts<-lapply(sce.filter,function(tmp){
  m=assay(tmp,'counts')
  
  C <- countland(m)
  C <- Subsample(C,cell_counts='min')
  # preview the resulting subsampled count matrix
  #C@subsample[1:10,1:10]
  return(C)
})

Matrix::colSums(minCellCounts[[1]]@subsample)
minCellCounts[[1]]@subsample[1:3,1:3]
names(minCellCounts)<-names(sce.filter)

# create single cell object
sce.minCellCounts<-lapply(names(minCellCounts),function(x){
  tmp=minCellCounts[[x]]@subsample
  n.expr.gene=Matrix::colSums(tmp>0)
  #cat(x);print(summary(n.expr.gene))
  mat=tmp[,n.expr.gene>=100] #expr at least 100 genes
  cell.info=colData(sce.filter[[x]])
  cell.info=cell.info[n.expr.gene>=100,]
  cat(x);print(table(cell.info$age)) #all contain >=50cells per age per tc
  #sce <- SingleCellExperiment(list(counts=minCellCounts[[x]]@subsample),colData=colData(sce.filter[[x]]) )
  sce <- SingleCellExperiment(list(counts=mat),colData=cell.info) 
  sce 
})
names(sce.minCellCounts)<-names(minCellCounts)
saveRDS(sce.minCellCounts, file='sce_minCellCounts.rds')

sce.minCellCounts=readRDS('sce_minCellCounts.rds')
lapply(sce.minCellCounts,function(i) table(i$age))
# all contain >=55 cells, some tc have unequal #cell at 3m vs 24m


##############################################################
# save a h5ad format for python use 
if(T){
sce.minCellCounts=readRDS('sce_minCellCounts.rds')
min(unlist(lapply(sce.minCellCounts,function(i) table(i$age)))) #55 cells

genes.union=sort(unique(unlist(lapply(sce.minCellCounts,rownames))))
length(genes.union)
tc.names=sort(names(sce.minCellCounts))

sce.minCellCounts=lapply(tc.names,function(tc){
  tmp=sce.minCellCounts[[tc]]
  na.genes=genes.union[!genes.union %in% rownames(tmp)]
  na.mat=Matrix::Matrix(0,nrow=length(na.genes),ncol=ncol(tmp))
  rownames(na.mat)=na.genes
  new.mat=rbind(na.mat,assay(tmp,'counts'))
  new.mat=new.mat[genes.union,] #all mat have the same rownames order
  sce <- SingleCellExperiment(list(counts=new.mat),
                              colData=colData(tmp) )
  sce
})
sce=Reduce(`cbind`,sce.minCellCounts)
#sce=SingleCellExperiment::cbind(sce.minCellCounts[[1]],sce.minCellCounts[[2]])
sce #22837 17244 

writeH5AD(sce, 'sce_minCellCounts.h5ad')

# run PCA
sce=readH5AD('sce_minCellCounts.h5ad') # 22966 31001 
sce #22837 17244 
assayNames(sce) #"counts" 
logcounts(sce)<-log2(counts(sce)+1)
assayNames(sce) #"counts"    "logcounts"

pca_out<-scater::runPCA(sce,ncomponents=50,ntop=2000, 
         scale=FALSE,exprs_values='logcounts') #https://rdrr.io/github/davismcc/scater/man/runPCA.html
pca_out # a single cell object
reducedDimNames(pca_out) #https://bioconductor.org/packages/devel/bioc/vignettes/SingleCellExperiment/inst/doc/intro.html#3_Adding_low-dimensional_representations
#head(reducedDim(pca_out, "PCA")[,1:2])
colnames(colData(pca_out))[14]='tissue_cell_type'
writeH5AD(pca_out, 'sce_minCellCounts_PCA50.h5ad')
}


############################################################
## find commonly expressed genes (>=10% cells per age group)
if(F){
sce.minCellCounts=readRDS('sce_minCellCounts.rds')
tc.names=names(sce.minCellCounts)

## find shared expr genes
sce.shared<-lapply(tc.names,function(tc){
  sce=sce.minCellCounts[[tc]]
  count=assay(sce,'counts')
  
  # gene expr in >=10% cell in each age
  cell.meta=colData(sce)
  ages=unique(cell.meta$age)
  include.genes<-sapply(ages,function(age){
    tmp=count[,cell.meta$age==age]
    i=Matrix::rowSums(tmp>0)
    i>=ncol(tmp)*0.1
  })
  i=Matrix::rowSums(include.genes)==length(ages)
  sce=sce[i,]
  cat(tc,sum(i),'genes kept\n')
  return(sce)
})
names(sce.shared)=tc.names
sapply(sce.shared,dim)
saveRDS(sce.shared,'sce_minCellCounts_shared.rds')
}

