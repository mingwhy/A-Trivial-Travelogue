
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
library(ggplot2);library(gridExtra)
library(tidyverse)
# from mz to gene, or from gene to mz from `00_Fruitfly-GEM_metabolite_neighbors.R`
mz_gene_list=readRDS('mz_gene_list.rds')
names(mz_gene_list)
########################################################
## read in gene symbol mapping
gene.meta=readRDS('~/Documents/Data_AgingFlyCellAtlas/AFCA_gene.id.rds')
#gene.meta=readRDS('~/Documents/Data_AgingFlyCellAtlas/AFCA_gene.meta.rds')
gene.meta=gene.meta[gene.meta$FLYBASE!='',]
head(gene.meta)
dup.names=names(which(table(gene.meta$original.id)>1))

##########################################################
## read in mz measurement data for trajectory comparison
load('~/Documents/aging_metabolism/batch_block_corrected_mzdata')
ls() #'dat' 'mzs'
head(dat)
dim(dat) #181  93
dat[1:3,1:8]
sort(unique(dat$AgeNum)) #4 10 24 45 69 80 days
colnames(dat)[1:8]
metabolome.dat=dat;
sample.info=metabolome.dat[,1:6]

load('~/Documents/aging_metabolism/age-associated.metabolites.for.Ming')
ls() #"cmat",age-independent covariance between metabolites 
# "for.Ming",  mz ~ age betas 
mz.age.betas=for.Ming
dim(mz.age.betas) #86

######################################################################
## calculate mean gene expr score per cell type per age
## only run once!!
library(SingleCellExperiment)
library(zellkonverter)
library(scater);library(scran)
library(ggplot2);library(gridExtra);
library(tidyverse);
library(SummarizedExperiment)

sce=readH5AD('sce_filteredBy_ncell100.h5ad')

assayNames(sce)='counts'
if(!file.exists('sizeFactors_sce.rds')){
  #Pooling normalization
  #https://bioconductor.org/packages/devel/bioc/vignettes/scuttle/inst/doc/norm.html#3_Pooling_normalization
  clusters <- quickCluster(sce)
  sce <- computePooledFactors(sce, clusters=clusters)
  summary(sizeFactors(sce))
  saveRDS(sizeFactors(sce),'sizeFactors_sce.rds')
}
#size.factors=readRDS('sizeFactors_sce.rds')
#length(size.factors);dim(sce)
#sce=logNormCounts(sce,size.factors=size.factors,pseudo.count =1)#https://rdrr.io/github/LTLA/scuttle/man/normalizeCounts.html

sce=logNormCounts(sce);
assayNames(sce) # "counts"    "logcounts"

out.dir='shuffle_out';
dir.create(out.dir)

bootstrap=20; #each age label permutation, only took one sample
prop.cell=0.8; #sample 80%

for(irep in 1:50){
  
  out.file=paste0(out.dir,'/rep_',irep,'.rds');
  out.file; #tissue; cell.type; .rds
  if(file.exists(out.file)){next}

  # shuffle age label, then assign genes mz according to mz_gene_list
  all.mz.out=list();
  
  for(tc in as.character(unique(sce$tc))){
    sce.tc=sce[,sce$tc==tc]
    expr.m=assay(sce.tc,'logcounts')
    
    true.age=sce.tc$age
    pseudo.age=sample(true.age,length(true.age),replace = F)
    sce.tc$age=pseudo.age;
    
    bootstrap.samples=list()
    for(i.bootstrap in 1:bootstrap){
      
      x=lapply(unique(sce.tc$age),function(age){
        sce.tc.age=expr.m[,sce.tc$age==age,drop=F]
        tmp=sce.tc.age[,sample(1:ncol(sce.tc.age),ceiling(ncol(sce.tc.age)*prop.cell),replace = F),drop=F]
        tmp1=tmp[rowSums(tmp)!=0,] # 1xncell matrix, rowsum==0
        #gene.na=tabulate(tmp1@i + 1L, nrow(tmp1)) ### nnz per row,number of non-zeros
        exprs=rowMeans_drop0(tmp)
        names(exprs)=rownames(tmp1) #only non-zero genes expr values would be returned
        exprs
      })
      #sapply(x,length) 
      names(x)=unique(sce.tc$age)
      bootstrap.samples[[i.bootstrap]]<-x;
    }
    
    # finish cell type, start mz one by one
    for(mz in names(mz_gene_list)){
      cat(mz,'\n')
      # obs metabolome data
      mz.name=mz.age.betas[mz.age.betas$KEGGid==mz,]$mz
      
      # extract this mz associated gene expr from scRNA-seq
      df.g=mz_gene_list[[mz]]
      genes=c(df.g$from,df.g$to)
      genes=gsub("^C[0-9]{5}", "", genes)
      genes=genes[genes!='']  
      genes=unlist(lapply(genes,function(i){unlist(strsplit(i,'or'))}))
      genes=gsub('^\\s+|\\s+$','',genes)
      genes=unique(genes)
      
      one.mz.out=c()
      for(i.bootstrap in 1:length(bootstrap.samples)){
        x=bootstrap.samples[[i.bootstrap]]
        for(age in names(x)){
          tmp=x[[age]][genes];
          names(tmp)=genes
          tmp[is.na(tmp)]=0;
          one.mz.out=rbind(one.mz.out,c(tc,age,tmp))
        }
      }       
      colnames(one.mz.out)[1:2]=c('cell.type','age')
      tmp=reshape::melt(as.data.frame(one.mz.out),id=c("cell.type","age"))
      df.one.mz.one.tc=data.frame(age=tmp$age,expr.values=tmp$value,cell.type=tmp$cell.type,gene=tmp$variable)
    
      all.mz.out[[mz]]=rbind(all.mz.out[[mz]],df.one.mz.one.tc)
    }
  }
  
  length(all.mz.out) #number of metabolites
  tmp=all.mz.out[[1]]
  unique(table(tmp$cell.type) )
  length(unique(tmp$age)) * bootstrap * length(unique(tmp$gene))
  
  saveRDS(all.mz.out, out.file);
  cat('irep',irep,'done\n')
}

