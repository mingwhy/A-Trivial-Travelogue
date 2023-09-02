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
##############################################
## load FruitMat GEM
library(R.matlab) #in matlab `load(Fruitfly-GEM.mat)`
gem=readMat('~/Documents/bioinfo_software/GEM_metabolic.models/Fruitfly-GEM-main/model/Fruitfly-GEM.mat')
# run GEM_convertOldStyleModel.m to get `rules` in GEM fields, used in exprs2fluxes.R
rules=readMat('rules.mat')
length(rules)

ifly=gem[[1]][,,1]
length(ifly)
lapply(ifly,dim)
names(ifly)
dim(ifly$S) # 8135mz X 11898 fluxes (how many overlapped with metabolome assay)
#genes=unlist(ifly$genes)
#length(genes) #1753 genes
ifly$rules=rules$rules
length(ifly$subSystems)
ifly$subSystems=unlist(ifly$subSystems)

ifly$genes=unlist(ifly$genes)
length(ifly$genes)
length(ifly$rules)
ifly$ub=ifly$ub[,1]
ifly$lb=ifly$lb[,1]
ifly$rowub=rep(0,nrow(ifly$S)) #used in iMat, check with human data `data("recon1")`, all 0
ifly$rowlb=rep(0,nrow(ifly$S))

##
#not do this: ifly$rxnNames=unlist(ifly$rxnNames) #unlist would drop NULL automatically
ifly$rxnNames=as.character(ifly$rxnNames)
length(unlist(ifly$rxns)) #11898
ifly$rxns=unlist(ifly$rxns)
ifly$mets=unlist(ifly$mets) #all members non-empty, checked length
ifly$metNames=unlist(ifly$metNames)


##############################################
library(SingleCellExperiment)
library(zellkonverter)
library(scater);library(scran)
library(ggplot2);library(gridExtra);
library(tidyverse);
library(SummarizedExperiment)

###############
## read in kegg module info
gene.meta=readRDS('~/Documents/Data_AgingFlyCellAtlas/AFCA_gene.id.rds')
#gene.meta=readRDS('~/Documents/Data_AgingFlyCellAtlas/AFCA_gene.meta.rds')
gene.meta=gene.meta[gene.meta$FLYBASE!='',]
head(gene.meta)
dup.names=names(which(table(gene.meta$original.id)>1))

###############################################################
## read in data
sce=readH5AD('sce_filteredBy_ncell100.h5ad')

assayNames(sce)='counts'
sce=logNormCounts(sce)
assayNames(sce) # "counts"    "logcounts"

nc <- 4L # number of CPUs/cores to use for parallelization
library(gembox);
source('src_exprs2fluxes.R')
model=ifly

out.dir='imat_out_shuffle';
dir.create(out.dir)
nrep=10;
source('src_gembox_fluxVec.R')

for(tc in unique(sce$tc)){
  
  #test on 'adult brain perineurial glial cell'
  #tc='head;adult brain perineurial glial cell'
  #tc='head;gamma Kenyon cell'
  id=gsub('\\/','\\-',tc) #save filename can not contain '/'
  out.file=paste0(out.dir,'/shuffle_',id,'.rds');
  out.file; #tissue; cell.type; .rds
  
  if(!file.exists(out.file)){
    sce.tc=sce[,sce$tc==tc]
    dim(sce.tc)
    expr.m=assay(sce.tc,'logcounts')
    true.age=sce.tc$age
    
    if(F){
      tmp=expr.m[,true.age==5]
      tmp1=tmp[rowSums(tmp)!=0,]
      gene.na=tabulate(tmp1@i + 1L, nrow(tmp1)) ### nnz per row,number of non-zeros
      exprs=rowMeans_drop0(tmp)  
      #exprs[gene.na<50]=0 #gene which expr in less than 50 cells, do not estimate mean.expr  
      names(exprs)=rownames(tmp1)
      obs.exprs=exprs
    }
    
    res=c()
    for(irep in 1:nrep){
      sce.tc$age=sample(true.age,ncol(sce.tc),replace = F)
      x=sapply(unique(sce.tc$age),function(age){
        tmp=expr.m[,sce.tc$age==age]
        tmp1=tmp[rowSums(tmp)!=0,]
        gene.na=tabulate(tmp1@i + 1L, nrow(tmp1)) ### nnz per row,number of non-zeros
        exprs=rowMeans_drop0(tmp)  
        #exprs[gene.na<50]=0 #gene which expr in less than 50 cells, do not estimate mean.expr  
        names(exprs)=rownames(tmp1)
        #tmp.names=intersect(names(exprs),names(obs.exprs))
        #cor(obs.exprs[tmp.names],exprs[tmp.names],use = 'complete')
        #exprs=exprs[exprs!=0] #handled in src_geombox_fluxVec.R
        
        tmp=gembox_fluxVec(exprs=exprs,model=model) #GEM model
        #summary(tmp)
        return(tmp)
      })
      colnames(x)=paste0('rep',irep,'_','age',unique(sce.tc$age))
      res=cbind(res,x)
    }
    #colnames(res)
    colnames(res)=paste0(colnames(res),'_',tc) #rep; age; tc
    saveRDS(res,out.file)
  }
}
  


