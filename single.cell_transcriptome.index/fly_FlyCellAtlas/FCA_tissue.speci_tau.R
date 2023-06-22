library(tidyverse)
library(Seurat)

##########################################################
## read in gene meta info for FCA
gene.meta=data.table::fread('~/Documents/Data_fly_FCA/gene.meta.txt')
dim(gene.meta) #16276     11
gene.meta=as.data.frame(gene.meta)
rownames(gene.meta)=gene.meta$current_symbol
dim(gene.meta) #16276    11

#####################################################################
## read in FCA data, average gene expr per tissue
folders=folders=Sys.glob('~/Documents/Data_fly_FCA/*FCA*')
#folders=folders[c(1,7)] #antenna and body
#folders=folders[c(2,9)]  #gut and fat.body

tissue.gene.mean.expr=list()
for(folder in folders){
  tissue=strsplit(basename(folder),'\\_')[[1]][2]
  input.file=Sys.glob(paste0(folder,'/*valid.rds'))
  cat('begin',tissue,input.file,'\n')
  
  ## read in data
  #tissue='gut'
  #input.file = "../../single.cell_datasets/FCA_gut/whole_gut_filtered_valid.rds"
  sc = readRDS(input.file)
  dim(sc)
  
  # antenna and body has female, male, mix, three class labels
  # only keep cell types with both female and male cells profiled
  sc=sc[,sc$sex!='mix' & sc$annotation!='unannotated']
  meta.info=sc@meta.data
  
  sc=NormalizeData(sc)
  expr.mat=sc@assays$RNA@data
  genes=Matrix::rowMeans(expr.mat)  
  tissue.gene.mean.expr[[tissue]]=genes;
}

sapply(tissue.gene.mean.expr,length); #different tissues have different number of genes
names(tissue.gene.mean.expr);
tissue.gene.mean.expr[[1]]
saveRDS(tissue.gene.mean.expr,'tissue.gene.mean.expr.rds')

## different tissue, total expr abundance
sapply(tissue.gene.mean.expr,sum) #between 1500~2700 (testis)

genes=unique(unlist(sapply(tissue.gene.mean.expr,names)))
length(genes) #16276
df.expr=matrix(0,length(genes),length(tissue.gene.mean.expr))
rownames(df.expr)=genes
colnames(df.expr)=names(tissue.gene.mean.expr)
for(i in 1:length(tissue.gene.mean.expr)){
  x=tissue.gene.mean.expr[[i]]
  df.expr[names(x),i]=x
}
dim(df.expr) # 16276    17
df.expr[1:3,1:3]

## compute spatial (tissue) specificity
tc.names=colnames(df.expr)
Xmax=tc.names[apply(df.expr,1,which.max)]#which tc has the max
length(Xmax) #16276 genes

gene.max.tissue=apply(df.expr,1,which.max)
gene.max.tissue=tc.names[gene.max.tissue]
gene.max=apply(df.expr,1,max)
summary(gene.max) #no 0
n=ncol(df.expr)
df.expr2=1-df.expr/gene.max
tau=Matrix::rowSums(df.expr2)/(n-1)
hist(tau) #0,HK. 1,high specificity

head(gene.meta)
df.tau=as.data.frame(tau)
df.tau$SYMBOL=names(tau)
df.tau$tissue=gene.max.tissue
tmp=merge(df.tau, gene.meta)
dim(tmp) #16276 genes

tmp$group=cut(tmp$tau,breaks=seq(0,1,0.05),include.lowest = T)
table(tmp$group)

pdf('tissue.gene.mean.expr_tau.pdf',useDingbats = T,width = 16)
print( ggplot(subset(tmp,tau>=0.15),aes(x=group,fill=tissue))+
  geom_bar()+theme_classic()
)
dev.off()

data.table::fwrite(tmp,file='tissue.gene.mean.expr_tau.txt')
