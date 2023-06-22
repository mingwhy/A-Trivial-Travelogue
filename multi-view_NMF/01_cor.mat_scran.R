
library(Seurat)
library(ggplot2);library(gridExtra);
options(stringsAsFactors = F)

#out.dir1='brain_scRNA-seq/';
#out.dir1='brain_snRNA-seq_n5c0.05/';
out.dir1='brain_snRNA-seq_n10c0.1/';
if(!dir.exists(out.dir1)) dir.create(out.dir1)

out.dir=paste0(out.dir1,'/pearson_coexpr/');
if(!dir.exists(out.dir)) dir.create(out.dir)

## save data in small size and compatible format between R and Matlab
#https://github.com/HenrikBengtsson/R.matlab
#BiocManager::install('R.matlab')
library(R.matlab)

## read in processed wholebrain data
if(length(grep('sn',out.dir1))==1){
  file="../single.cell_datasets/FCA_head/whole_head_filtered_valid.rds";
}else{
  file="../single.cell_datasets/fly.brain.atlas/wholebrain_filtered_valid.rds";
}

dat=readRDS(file);
#dat=NormalizeData(dat);
#df.umi=dat@assays$RNA@counts #12094 features across 56192 samples within 1 assay 
df.expr=dat@assays$RNA@data #logNormal

unique(dat@meta.data$sex)
table(dat@meta.data$annotation)

## separate female and male, select cell clusters that contain>=500cells in both sexes
min.cell=100;max.cell=Inf;
i=apply(table(dat$sex,dat$annotation),2,function(i)sum(i>=min.cell & i<max.cell))
pick.cell.clusters=names(which(i==2)) #79
pick.cell.clusters #sn,54; sc,60.

dat=subset(dat,annotation %in% pick.cell.clusters)
df.expr=dat@assays$RNA@counts;
dim(df.expr); #sn:12602 50588

sparsity=sum(df.expr==0)/nrow(df.expr)/ncol(df.expr)
cat(out.dir,sparsity,'\n') #sn=0.9537,sc=0.9055

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
#system.time({df.expr.norm=scran_norm(df.expr)}) #15min
#saveRDS(df.expr.norm,file=paste0(out.dir1,'/expr_scran.norm.rds'))

## effective.expressed.gene >=1umi in max(10,ncell*10%) cells
pdf('check_pearson.cor.mat_hist.pdf',useDingbats = T,width = 14)
par(mfrow=c(4,4),mai=c(1,1,1,0.5),mar=c(2.5,1,1.5,0.5))

genes=list();
for(i.sex in c('male','female')){
  for(i.cluster in pick.cell.clusters){
    mat=df.expr[,dat$sex==i.sex & dat$annotation==i.cluster]
    
    #gene.filter=Matrix::rowSums(mat>0) > max(5,ncol(mat)*0.05)
    gene.filter=Matrix::rowSums(mat>0) > max(10,ncol(mat)*0.1)
    
    mat=as.matrix(mat[gene.filter,])
    mat=scran_norm(mat);
    
    ncell=ncol(mat)
    #set.seed(2049)
    #sample.mat=mat[,sample(1:ncol(mat),min.cell,replace = F)]
    #cor.mat=Rfast::cora(t(sample.mat)) #cor.mat
    cor.mat=coop::tpcor(mat) #cor.mat
    
    i.cluster2=gsub('\\/','\\.',i.cluster)#there is _ and -, use . to replace /
    genes[[paste0(i.sex,i.cluster2,sep='__')]]=rownames(cor.mat)
    
    out.file=paste0(out.dir,'/ncell_',ncell,'_',i.sex,'_',i.cluster2,'_pearson.rds')
    saveRDS(cor.mat,file=out.file)
    #filename=paste0(out.dir,'/',i.sex,'_',i.cluster2,'_perb.mat')
    #writeMat(filename, cell.type=cor.mat)
    
    x=cor.mat[upper.tri(cor.mat)]
    hist(x,main=paste0(i.sex,i.cluster,sep=','),xlab='')
  }
}
dev.off()


c.genes=unlist(genes)
common.genes=names(which(table(c.genes)==length(genes)))
all.genes=unique(c.genes)
length(all.genes) #5801
length(common.genes) #914

saveRDS(list(genes=genes,common.genes=common.genes,all.genes=all.genes),
        file=paste0(out.dir1,'all_common_genes.rds'))

