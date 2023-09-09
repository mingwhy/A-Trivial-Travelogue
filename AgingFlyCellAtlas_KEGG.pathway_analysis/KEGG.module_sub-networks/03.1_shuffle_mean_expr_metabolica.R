
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
############################################

library(ggplot2);library(gridExtra)
library(tidyverse);library(Matrix)
#https://www.preprints.org/manuscript/202101.0280/v1
source('propagationFuncs.R')

load("./files/all_producers.RData")
length(all_producers) #1559, each one <=> one metabolite
names(all_producers)
load("./files/rxn2gene.RData")
dim(rxn2gene) #df,6420 obs. of  4 variables, reaction and its associated enzymes
head(rxn2gene) #3rd column: Entrez 
length(unique(rxn2gene$Reaction)) #1901

# end product of each (sub)module
mz.produced=lapply(1:length(all_producers),function(i){
  all_producers[[i]]$isProduced
})
table(unlist(lapply(mz.produced,length))) #1559,each sub-graph only produce one met
mz.produced=unlist(mz.produced)
tail(sort(table(mz.produced))) #there are metabolites which are end product of multiple module
length(unique(mz.produced)) #1270
which(mz.produced=='cpd:C00388')
all_producers[[718]] #sub-network which produce Histamine

# read in fly-human orthologs (00_process_h5ad_AFCA.R)
orthologs=readRDS('~/Documents/Data_AgingFlyCellAtlas/human_fly_ortholog_entrezid.rds')
colnames(orthologs)
sum(orthologs$human_entrezid=='')
sum(is.na(orthologs$human_entrezid)) #1558
orthologs=orthologs[!is.na(orthologs$human_entrezid),]
dim(orthologs) # 24085    10
gene.order=unique(orthologs$human_entrezid)
length(gene.order)#11310

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

out.dir='shuffle_out';
dir.create(out.dir)

for(irep in 21:100){

  out.file=paste0(out.dir,'/rep_',irep,'.rds');
  out.file; #tissue; cell.type; .rds
  if(file.exists(out.file)){next}
  
  out<-lapply(as.character(unique(sce$tc)),function(tc){
    sce.tc=sce[,sce$tc==tc]
    expr.m=assay(sce.tc,'logcounts')
    #expr.m[expr.m==0]=NA #remove zero expr gene, transform NA here too slow
    #class(expr.m) # "dgCMatrix"
    true.age=sce.tc$age
    pseudo.age=sample(true.age,length(true.age),replace = F)
    sce.tc$age=pseudo.age;
    
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
  
  names(out)=as.character(unique(sce$tc))
  saveRDS(out,out.file)
  
}
  


