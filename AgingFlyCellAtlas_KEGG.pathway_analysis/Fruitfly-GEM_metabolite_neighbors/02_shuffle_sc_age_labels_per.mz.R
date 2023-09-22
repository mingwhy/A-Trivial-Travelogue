
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

out.dir='shuffle_out1';
dir.create(out.dir)

bootstrap=1; #each age label permutation, only took one sample
prop.cell=0.8; #sample 80%

for(irep in 1:100){
  
  out.file=paste0(out.dir,'/rep_',irep,'.rds');
  out.file; #tissue; cell.type; .rds
  if(file.exists(out.file)){next}
  
  all.mz.out=list();
  
  for(mz in names(mz_gene_list)){
    
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
    
    #sce.sub=sce[intersect(genes,rownames(sce)),,drop=F]
    sce.sub=sce[intersect(genes,rownames(sce)),]
    
    one.mz.out=list()
    
    for(gene.name in rownames(sce.sub)){
      out<-lapply(as.character(unique(sce.sub$tc)),function(tc){
        sce.tc=sce.sub[gene.name,sce.sub$tc==tc]
        expr.m=assay(sce.tc,'logcounts')
        
        true.age=sce.tc$age
        pseudo.age=sample(true.age,length(true.age),replace = F)
        sce.tc$age=pseudo.age;
        
        x=lapply(unique(sce.tc$age),function(age){
          sce.tc.age=expr.m[,sce.tc$age==age,drop=F]
          expr.values=c() 
          while(length(expr.values)<bootstrap){
            tmp=sce.tc.age[,sample(1:ncol(sce.tc.age),ceiling(ncol(sce.tc.age)*prop.cell),replace = F),drop=F]
            tmp1=tmp[rowSums(tmp)!=0,] # 1xncell matrix, rowsum==0
            if(nrow(tmp1)==0){
              expr.values=c(expr.values,0)
              next
            } #sampled values are all 0
            #gene.na=tabulate(tmp1@i + 1L, nrow(tmp1)) ### nnz per row,number of non-zeros
            exprs=rowMeans_drop0(tmp)
            names(exprs)=rownames(tmp1) #only non-zero genes expr values would be returned
            expr.values=c(expr.values,exprs)
          }
          data.frame(age=age,expr.values=expr.values)
        })
        
        df=as.data.frame(Reduce(`rbind`,x))
        df$cell.type=tc
        cat('cell.type',tc,'is done\n')
        return(df)
      })
      
      df.out=as.data.frame(Reduce(`rbind`,out))
      df.out$gene=gene.name;
      one.mz.out[[gene.name]]=df.out
    }
    df.one.mz=as.data.frame(Reduce(`rbind`,one.mz.out))
    all.mz.out[[mz]]=df.one.mz
  }
  
  saveRDS(all.mz.out, out.file);
  
}

