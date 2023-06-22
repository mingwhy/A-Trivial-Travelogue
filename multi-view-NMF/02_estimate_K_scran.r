
## for each cell type, select effectively expressed genes
## estimate the number of latent factors
library(Seurat)
library(ggplot2);library(gridExtra);
options(stringsAsFactors = F)

## save data in small size and compatible format between R and Matlab
#https://github.com/HenrikBengtsson/R.matlab
#BiocManager::install('R.matlab')
library(R.matlab)
#source('src_EstNumModule.R')
source('src_estimate_K_RMT.R')

## read in processed wholebrain data
out.file=paste0('sc_scran_100_k_est.out.txt');
file="../single.cell_datasets/fly.brain.atlas/wholebrain_filtered_valid.rds";

#out.file=paste0('sn_scran_100_k_est.out.txt'); #the min #cell to be included
#file="../single.cell_datasets/FCA_head/whole_head_filtered_valid.rds";

dat=readRDS(file);
df.expr=dat@assays$RNA@counts #12094 features across 56192 samples within 1 assay 
colnames(dat@meta.data)
unique(dat@meta.data$sex)

## separate female and male, select cell clusters that contain>=500cells in both sexes
min.cell=500;max.cell=Inf;
i=apply(table(dat$sex,dat$annotation),2,function(i)sum(i>=min.cell & i<max.cell))
pick.cell.clusters=names(which(i==2)) #79
pick.cell.clusters

# cell clusters with known cell type annotations
#pick.cell.clusters=pick.cell.clusters[grep('[A-Z]',pick.cell.clusters)]
pick.cell.clusters #41


## if you want to use scran to normalize data
library(scran)
scran_norm<-function(expr.mat){
  sce <- SingleCellExperiment(list(counts=expr.mat),
                              #colData=DataFrame(cell.type=cell.type),
                              rowData=DataFrame(gene=rownames(expr.mat)) )
  #sce
  clusters <- quickCluster(sce)
  sce <- computeSumFactors(sce, clusters=clusters)
  #summary(sizeFactors(sce))
  sce <- logNormCounts(sce)
  log.sce=sce@assays@data$logcounts
  #dim(expr.mat);dim(log.sce);
  return(log.sce)
}


## effective.expressed.gene >=1umi in max(5,ncell*10%) cells

out=c();
for(i.sex in c('male','female')){
  for(i.cluster in pick.cell.clusters){
    mat=df.expr[,dat$sex==i.sex & dat$annotation==i.cluster]
    #gene.filter=Matrix::rowSums(mat>0) > max(5,ncol(mat)*0.1)
    gene.filter=Matrix::rowSums(mat>0) > max(10,ncol(mat)*0.1)
    mat=as.matrix(mat[gene.filter,])
    mat=scran_norm(mat);
    
    ncell=ncol(mat)
    #system.time(EstNumModule(data.m=mat)) #coop is faster than Rfast
    #i.cluster2=gsub('\\/','\\.',i.cluster)#there is _ and -, use . to replace /
    k=estimate_K(data.m=mat)
    out=rbind(out,c(i.sex,i.cluster,ncell,k))
    cat(i.sex,i.cluster,'is done\n')
  }
}
out1=as.data.frame(out)
colnames(out1)=c('sex','cluster','ncell','k')
write.table(out1,out.file,sep='\t',quote=F,row.names = F)

max(as.numeric(out1$k)) 
max(out1[out1$ncell>=500,]$k)
#ncell=100, sn=25, sc=17

