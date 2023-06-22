
library(ggplot2);library(gridExtra)
library(dplyr)
library(ggpubr)
#install.packages('Seurat')
library(Seurat)
library(irlba) # for fast PCA on large matrix: https://github.com/bwlewis/irlba
library(RcppML) #for fast NMF
library(SingleCellExperiment)
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
## use liver cell and perform DR (dimension reduction)
library(grDevices);library(RColorBrewer)
library(uwot); #for UMAP
#library(scran) #for clustering cells
#library(slingshot, quietly = FALSE)
#https://bioconductor.org/packages/devel/bioc/vignettes/slingshot/inst/doc/vignette.html
#full quantile normalization, a well-established method which forces each cell to have the same distribution of expression values.
FQnorm <- function(counts){
  rk <- apply(counts,2,rank,ties.method='min')
  counts.sort <- apply(counts,2,sort)
  refdist <- apply(counts.sort,1,median)
  norm <- apply(rk,2,function(r){ refdist[r] })
  rownames(norm) <- rownames(counts)
  return(norm)
}

table(cell.meta$tissue_cell.type,cell.meta$age)
all.tcs=sort(unique(cell.meta$tissue_cell.type))
#my.colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
plotcol <- brewer.pal(8,'Dark2')

for(tc in all.tcs){
  expression_matrix=df.expr[,cell.meta$tissue_cell.type==tc]
  gene.names <- rownames(expression_matrix)
  cell.names <- colnames(expression_matrix)
  
  cell_metadata<-cell.meta[cell.meta$tissue_cell.type==tc,]
  rownames(cell_metadata)=cell.names
  
  sce <- SingleCellExperiment(assays = List(counts = expression_matrix))
  # gene filter
  geneFilter <- apply(assays(sce)$counts,1,function(x){
    sum(x >= 3) >= 10
  })
  sce <- sce[geneFilter, ]
  
  # quantile normalize
  #assays(sce)$FQnorm <- FQnorm(assays(sce)$counts)
  # log-transformed normalize (https://bioconductor.org/packages/devel/bioc/vignettes/SingleCellExperiment/inst/doc/intro.html#3_Adding_low-dimensional_representations)
  counts <- assay(sce, "counts")
  libsizes <- colSums(counts)
  size.factors <- libsizes/mean(libsizes)
  logcounts(sce) <- log2(t(t(counts)/size.factors) + 1)
  
  assayNames(sce)
  # "counts"    "FQnorm"    "logcounts"
  #A=log1p(assays(sce)$norm); #ngene x ncell
  A=assays(sce)$logcounts; #ngene x ncell
  
  # pca
  #pca <- prcomp(t(A), scale. = FALSE) 
  pca <- prcomp(t(A), scale. = TRUE) #may be slow for big data
  rd1 <- pca$x[,1:2]
  if(F){
    n1=sum(cell_metadata$age=='young')
    n2=sum(cell_metadata$age=='old')
    S=irlba(A, nv=min(n1+n2,200))
    dim(S$v) #ncell X pc1..100
    rd1 <- S$v[,1:2]
  }
  #plot(rd1, col = rgb(0,0,0,.5), pch=16, asp = 1)
  
  # umap
  #rd2 <- uwot::umap(t(log1p(assays(sce)$norm)))
  rd2 <- uwot::umap(t(assays(sce)$logcounts))
  colnames(rd2) <- c('UMAP1', 'UMAP2')
  #plot(rd2, col = rgb(0,0,0,.5), pch=16, asp = 1)
  
  # NMF
  nmf_model<-RcppML::nmf(A,k=200,tol=1e-5)
  dim(A) # ngene x ncell
  dim(nmf_model@w) # ngene x k
  dim(nmf_model@h) # k x ncell
  length(nmf_model@d) # k
  rd3=t(nmf_model@h[c(1,2),])
  plot(rd3, col = rgb(0,0,0,.5), pch=16, asp = 1)
  #add both dimensionality reductions to the SingleCellExperiment object
  #reducedDims(sce) <- SimpleList(PCA = rd1, UMAP = rd2, NMF=rd3)
  R_mse <- mean((A - nmf_model$w %*% Diagonal(x = nmf_model$d) %*% nmf_model$h)^2)
  R_mse
  
  tmp=as.data.frame(cbind(rd1,rd2,rd3))
  tmp$age=cell_metadata$age
  p1<-ggplot(tmp,aes(x=PC1,y=PC2,col=age))+geom_point()+theme_classic()+scale_color_manual(values=plotcol)+ggtitle(tc)
  p2<-ggplot(tmp,aes(x=UMAP1,y=UMAP2,col=age))+geom_point()+theme_classic()+scale_color_manual(values=plotcol)+ggtitle(tc)
  grid.arrange(p1,p2,ncol=2)
  #par(mfrow=c(1,2),mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
  #plot(rd1, col = plotcol[as.numeric(factor(cell_metadata$age))], pch=16, asp = 1)
  #legend(0,0,legend=unique(cell_metadata$age),col=plotcol[1:2],pch=16, horiz=TRUE, cex=0.5,bty='n')
  #plot(rd2, col =  plotcol[as.numeric(factor(cell_metadata$age))], pch=16, asp = 1)
  
  
}


