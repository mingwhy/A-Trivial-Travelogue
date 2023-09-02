library(Matrix)
#https://stackoverflow.com/questions/51467276/how-to-find-the-column-means-for-a-sparse-matrix-excluding-0-values
# Note: rows or columsn with all 0 are dropped
colMeans_drop0 <- function (dgCMat) {
  nnz_per_col <- diff(dgCMat@p)
  ColInd <- rep.int(1:ncol(dgCMat), nnz_per_col)
  sapply(split(dgCMat@x, ColInd), mean)
}
rowMeans_drop0 <- function (dgCMat) {
  RowInd <- dgCMat@i + 1
  sapply(split(dgCMat@x, RowInd), mean)
}
# number of non-zero: https://stackoverflow.com/questions/51560456/r-package-matrix-get-number-of-non-zero-entries-per-rows-columns-of-a-sparse
##########
library(SingleCellExperiment)
library(zellkonverter)
library(scater);library(scran)
library(ggplot2);library(gridExtra);
library(tidyverse);
library(SummarizedExperiment)

###############################################################
## read in kegg module info
gene.meta=readRDS('~/Documents/Data_AgingFlyCellAtlas/AFCA_gene.id.rds')
#gene.meta=readRDS('~/Documents/Data_AgingFlyCellAtlas/AFCA_gene.meta.rds')
gene.meta=gene.meta[gene.meta$FLYBASE!='',]
head(gene.meta)
dup.names=names(which(table(gene.meta$original.id)>1))

###############################################################
## read in data
if(F){ #only run once
  sce0=readH5AD('~/Documents/Data_AgingFlyCellAtlas/adata_headBody_S_v1.0.h5ad') # 15992  566254  
  #sce0=readH5AD('~/Documents/Data_AgingFlyCellAtlas/adata_head_S_v1.0.h5ad') # 15992 289981 
  #sce0=readH5AD('~/Documents/Data_AgingFlyCellAtlas/adata_body_S_v1.0.h5ad') # 15992 276273 
  sce0 #SingleCellExperiment 
  summary(sce0$n_genes_by_counts) 
  summary(sce0$total_counts)
  summary(sce0$pct_counts_mt) #all below 5%
  table(sce0$tissue)
  sce0=sce0[,-grep('unannotated',ignore.case = T,sce0$afca_annotation)] #remove the 'unannotation' 
  
  sce=sce0[,sce0$sex=='female' & sce0$n_genes_by_counts>=200 & sce0$total_counts>=1000 & sce0$pct_counts_mt<5]
  sce$afca_annotation=droplevels(sce$afca_annotation)
  sce$cell.type=sce$afca_annotation;
  sce$tc=paste(sce$tissue,sce$cell.type,sep=';')
  cell.meta=colData(sce)
  #saveRDS(cell.meta,'female_filtered.cell.meta.rds')
  #saveRDS(cell.meta,'female_head_filtered.cell.meta.rds')
  #saveRDS(cell.meta,'female_body_filtered.cell.meta.rds')
  
  #cell.meta$tc=paste(cell.meta$tissue,cell.meta$afca_annotation,sep=';')
  length(table(cell.meta$cell.type))  #152
  
  # select tc with >=100 cell in all 4 ages
  x=as.data.frame(cell.meta) %>% group_by(tissue,age,cell.type) %>% summarise(ncell=n()) %>% filter(ncell>=100)
  x1=x %>% spread(age,ncell)
  dim(x1) #63
  nage=apply(x1,1,function(i) sum(!is.na(i[-c(1,2)])) )
  pick.cell.type=x1[nage==4,]
  dim(pick.cell.type) #27 (16 in head or 11 in body)
  
  pick.cell.type$tc=paste(pick.cell.type$tissue,pick.cell.type$cell.type,sep=';')
  
  keep.sce=sce[,sce$tc %in% pick.cell.type$tc]
  writeH5AD(keep.sce,'sce_filteredBy_ncell100.h5ad')
}

######################################################################
## calculate mean gene expr score per cell type per age

sce=readH5AD('sce_filteredBy_ncell100.h5ad')
#sce=keep.sce

assayNames(sce)='counts'
sce=logNormCounts(sce)
assayNames(sce) # "counts"    "logcounts"

out.file='log1p_female_gene.mean.expr.rds';
unique(sce$tc)

if(!file.exists(out.file)){
  out<-lapply(as.character(unique(sce$tc)),function(tc){
    sce.tc=sce[,sce$tc==tc]
    expr.m=assay(sce.tc,'logcounts')
    #expr.m[expr.m==0]=NA #remove zero expr gene, transform NA here too slow
    #class(expr.m) # "dgCMatrix"
    
    x=lapply(unique(sce.tc$age),function(age){
      tmp=expr.m[,sce.tc$age==age]
      tmp1=tmp[rowSums(tmp)!=0,]
      #gene.na=apply(tmp,1,function(i) sum(!is.na(i)))
      #gene.na2=apply(tmp,1,function(i) sum(i!=0))
      gene.na=tabulate(tmp1@i + 1L, nrow(tmp1)) ### nnz per row,number of non-zeros
      exprs=rowMeans_drop0(tmp)
      #length(exprs);length(gene.na);nrow(tmp1)
      #exprs=Matrix::rowMeans(tmp[,sce.tc$age==age],na.rm=T)
      #exprs[gene.na<50]=NA #gene which expr in less than 50 cells, do not estimate mean.expr  
      #exprs[gene.na<50]=0 #gene which expr in less than 50 cells, do not estimate mean.expr  
      names(exprs)=rownames(tmp1)
      exprs
    })
    cat('tc',tc,'is done\n')
    #dfx=as.data.frame(Reduce(`cbind`,x))
    names(x)=paste(tc,unique(sce.tc$age),sep=';') #tissue;cell.type;age
    return(x)
  })
  #out0=as.data.frame(Reduce(`cbind`,out))
  #out1=out0[!Matrix::rowSums(out0)==0,]
  #genes=rownames(out1)
  #out2=cbind(genes,out1)
  #data.table::fwrite(out2,file=out.file)
  
  names(out)=as.character(unique(sce$tc))
  saveRDS(out,out.file)
}

##############################################
##############################################
##############################################
## load FruitMat GEM
library(R.matlab) #in matlab `load(Fruitfly-GEM.mat)`
gem=readMat('./Fruitfly-GEM-main/model/Fruitfly-GEM.mat')
ifly=gem[[1]][,,1]
names(ifly)
dim(ifly$S) # 8135mz X 11898 fluxes (how many overlapped with metabolome assay)
genes=unlist(ifly$genes)
length(genes) #1753 genes
sum(genes %in% gene.meta$SYMBOL) #1653

## after iMat is done
res=R.matlab::readMat("results_rep1/iMat_fly.body.mat") 
#res=R.matlab::readMat("results_rep1/iMat_fly.head.mat") 
experiment=res$experiment[,,1]
names(experiment)
experiment$name
cell.type_per.age=unlist(experiment$conditions)
f=experiment$fluxes.exp.all
datf=lapply(f,'[[',1)
sapply(datf,dim)
datf=as.data.frame(Reduce(`cbind`,datf))
summary(datf[,2])
apply(datf,2,summary)
colnames(datf)=cell.type_per.age

# remove `unannotated` label
datf=datf[,-grep('unannotated', cell.type_per.age)]
cell.type_per.age=cell.type_per.age[-grep('unannotated', cell.type_per.age)]

#
dim(datf) #11898 reactions

#http://www.sthda.com/english/articles/31-principal-component-methods-in-r-practical-guide/112-pca-principal-component-analysis-essentials/
library("FactoMineR")
library("factoextra")
res.pca <- PCA(t(datf), graph = FALSE)
eig.val <- get_eigenvalue(res.pca)
eig.val
fviz_eig(res.pca, addlabels = TRUE, ylim = c(0, 50))
#fviz_pca_ind(res.pca)
grp.cell.type=unlist(lapply(strsplit(cell.type_per.age,';'),'[[',2))
grp.age=unlist(lapply(strsplit(cell.type_per.age,';'),'[[',3))
fviz_pca_ind(res.pca, geom='point',pointshape = 21,
             pointsize = 2.5,
             #col.ind = grp.cell.type,
             #fill.ind  = grp.age) # color by groups
             fill.ind  = grp.cell.type) # color by groups
