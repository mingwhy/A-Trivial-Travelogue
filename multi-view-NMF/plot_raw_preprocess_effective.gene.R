library(SingleCellExperiment)
library(SCopeLoomR)
library(Seurat)
library(Seurat)
library(ggplot2);library(gridExtra);
options(stringsAsFactors = F)

dataset='sc';
#dataset='sn';

## read in raw data
if(dataset=='sn'){
  # data download from https://flycellatlas.org/
  loom_path <- '../single.cell_datasets/FCA_head/s_fca_biohub_head_10x.loom';
  loom <- open_loom(loom_path, mode="r+")
  
  cell.annotation.all=get_cell_annotation(loom)
  dim(cell.annotation.all) #11788   435
  
  labels=colnames(cell.annotation.all)
  tmp=cell.annotation.all[,-grep('TrackRegulonsAUC|MotifRegulonsAUC',labels)]
  colnames(tmp)
  table(tmp$batch_id)
  sort(table(tmp$batch))
  sort(table(tmp$id))  #only two samples
  sort(table(tmp$batch))==sort(table(tmp$id)) #batch <=> id
  
  cell.annotation.df=tmp
  colnames(cell.annotation.df)
  head(cell.annotation.df$annotation)
  table(cell.annotation.df$annotation)
  
  #genes=get_genes(loom)
  #length(genes) #13056 genes
  #raw <- get_dgem(loom)
  #raw[1:5,1:5]
  #dim(raw) #13056 gene by 100527 cell
  close_loom(loom)  
  x=table(cell.annotation.df$sex,cell.annotation.df$annotation)
  df.raw=reshape2::melt(x)
  colnames(df.raw)=c('sex','cell.cluster','ncell')
}else{
  ## cell type info
  df.cell=read.table("../single.cell_datasets/fly.brain.atlas/GSE107451_DGRP-551_w1118_WholeBrain_57k_Metadata.tsv.gz",as.is=T,header=T);
  table(df.cell$annotation);table(df.cell$sex)
  #sum(df.cell$annotation=='Hsp') #668 cells annotated to 'Hsp'
  # remove them due to personal email communication, they are 'stressed' cells.
  x=table(df.cell$sex,df.cell$annotation)
  df.raw=reshape2::melt(x)
  colnames(df.raw)=c('sex','cell.cluster','ncell')
}
dim(df.raw)
df.sub=df.raw[df.raw$sex=='female',]
x=df.sub[order(df.sub$ncell),]$cell.cluster
df.raw$cell.cluster=factor(df.raw$cell.cluster,levels=x)
p1=ggplot(df.raw,aes(x=cell.cluster,y=ncell))+
  geom_bar(stat='identity',fill='white',col='blue')+
  facet_wrap(.~sex)+scale_y_log10()+theme_classic()+
  coord_flip()
if(dataset=='sn'){p1=p1+ggtitle('Single-nuclei dataset')
}else{p1=p1+ggtitle('Single-cell dataset') }


## read in processed wholebrain data
if(dataset=='sn'){
  file="../single.cell_datasets/FCA_head/whole_head_filtered_valid.rds";
}else{
  file="../single.cell_datasets/fly.brain.atlas/wholebrain_filtered_valid.rds";
}

dat=readRDS(file);
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

df=as.numeric();
genes=list();
for(i.sex in c('male','female')){
  for(i.cluster in pick.cell.clusters){
    mat=df.expr[,dat$sex==i.sex & dat$annotation==i.cluster]
    #gene.filter=Matrix::rowSums(mat>0) > max(5,ncol(mat)*0.05)
    gene.filter=Matrix::rowSums(mat>0) > max(10,ncol(mat)*0.1)
    mat=as.matrix(mat[gene.filter,])
    ncell=ncol(mat);
    ngene=sum(gene.filter)
    df=rbind(df,c(i.sex,i.cluster,ncell,ngene))
    genes[[paste0(i.sex,i.cluster,sep='__')]]=rownames(mat)
  }
}
df=as.data.frame(df)
colnames(df)=c('sex','cell.cluster','ncell','ngene')
df$ncell=as.numeric(df$ncell)
df$ngene=as.numeric(df$ngene)

df.sub=df[df$sex=='female',]
x=df.sub[order(df.sub$ncell),]$cell.cluster
df$cell.cluster=factor(df$cell.cluster,levels=x)
p2=ggplot(df,aes(x=cell.cluster,y=ncell))+geom_bar(stat='identity')+
  facet_wrap(.~sex)+scale_y_log10()+theme_classic()+
  coord_flip()
if(dataset=='sn'){
  p2=p2+ggtitle('Single-nuclei dataset, after pre-processing')
}else{
  p2=p2+ggtitle('Single-cell dataset, after pre-processing')
}


p3=ggplot(df,aes(x=cell.cluster,y=ngene))+geom_bar(stat='identity')+
  facet_wrap(.~sex)+scale_y_log10()+theme_classic()+
  coord_flip()
if(dataset=='sn'){p3=p3+ggtitle('Single-nuclei dataset, after pre-processing')
}else{p3=p3+ggtitle('Single-cell dataset, after pre-processing') }

## genes express xx cell clusters
all.genes=unlist(genes)
x=table(table(all.genes))
dfx=as.data.frame(x)

colnames(dfx)=c('#cell.cluster one gene expressed',ylab='#gene')
p4=ggplot(dfx,aes(x=dfx[,1],y=dfx[,2]))+geom_bar(stat='identity')+
  theme_classic()+
  xlab('#cell.cluster one gene is effectively expressed')+
  ylab('#gene')+coord_flip()
  
if(dataset=='sn'){p4=p4+ggtitle('Single-nuclei dataset, gene exprssion breadth')
}else{p4=p4+ggtitle('Single-cell dataset, gene exprssion breadth')}

pdf(paste0('dataset_',dataset,'.pdf'),useDingbats = T,height = 9,width = 18)
grid.arrange(p1+theme(axis.text.y = element_text(size=5))+
               geom_bar(stat='identity',fill='white',col='blue'),
             p2+geom_bar(stat='identity',fill='white',col='blue'),
             p3+ylab('#expressd gene')+geom_bar(stat='identity',fill='white',col='blue'),
             p4+theme(axis.text.y = element_text(size=5))+geom_bar(stat='identity',fill='white',col='blue'),
             ncol=4)
dev.off()

