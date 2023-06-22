
library(Seurat)
library(tidyverse)
library(Matrix);
library(ggplot2);library(gridExtra)
library(RColorBrewer)
source("./gcl_r_version.R")

## read in pathway data
## download dme-kegg.ID.df.txt and kegg-flygenes.rds from 
# https://github.com/mingwhy/bioinfo_homemade_tools/tree/main/KEGG.decompose
kegg.info=read.table('dme-kegg.ID.df.txt',as.is=T,sep="\t",header=T)
head(kegg.info)
dim(kegg.info) #137 pathway

kegg.genes=readRDS('kegg-flygenes.rds')
names(kegg.genes)
length(kegg.genes) #137 pathways and all their fly genes
sapply(kegg.genes,nrow)
length(unique(unlist(sapply(kegg.genes,'[',,2)))) 
#3249 unique fly genes
sum(sapply(kegg.genes,nrow)>=10); #113 pathway contain >=10 genes
sum(sapply(kegg.genes,nrow)>=20); #97 pathway contain >=20 genes

dme.id=names(kegg.genes)
n.gene=sapply(kegg.genes,nrow)
sum(kegg.info$kegg.id==dme.id) #137
df.kegg = data.frame(dme.name=kegg.info$kegg.name,dme.id=dme.id,n.gene=n.gene)
df.kegg=df.kegg[order(df.kegg$n.gene),]
head(df.kegg) #plot kegg size(ngene below)
df.kegg$dme.name=factor(df.kegg$dme.name,levels=df.kegg$dme.name)
kegg.ngene.plot=ggplot(df.kegg,aes(x=dme.name,y=n.gene))+geom_bar(stat='identity')+
  coord_flip()+theme_bw()+geom_text(label=df.kegg$n.gene,hjust=-0.5)+
  ylim(0,1300)+
  expand_limits(x=0)
#kegg.ngene.plot


## read in processed wholebrain data (repo: fly.brain.core_coexpr.net)
# https://github.com/mingwhy/fly.brain.core_coexpr.net
file="../processed_data/wholebrain_filtered.rds";
dat=readRDS(file);
colnames(dat@meta.data)
unique(dat@meta.data$sex)
unique(dat@meta.data$Age)
sum(dat@meta.data$annotation=='Hsp')

df.expr0=dat@assays$RNA@counts 

#dat <- NormalizeData(dat, normalization.method = "LogNormalize", scale.factor = 10000)
#df.expr0=dat@assays$RNA@data

## keep female cells
#sex='female'
#table(dat@meta.data$sex)
#i=dat@meta.data$sex==sex;
#sum(i) #28850 female cells
#df.expr=df.expr0[,i]
#df.meta=dat@meta.data[i,]
df.expr=df.expr0;
df.meta=dat@meta.data;

dim(df.expr)
dim(df.meta) #56192 
sum(colnames(df.expr)==rownames(df.meta))

#########################################
# a matrix of GCL: pathway X fly.cell.age
(ages=sort(unique(df.meta$Age)))

sub.sample=20;
num_divisions=20; #num_divisions in subsetting genes into 2 sets when calcualte GCL 
#outfile=paste0(sex,'_GCL_',num_divisions,'_sub',sub.sample,'.rds');
(outfile=paste0('pathway_GCL_div',num_divisions,'_min',sub.sample,'.rds'));

gcl.per.pathway=list()
if(!file.exists(outfile)){
  all.gene.names=rownames(df.expr)
  for(age.i in ages){
    pick=df.meta[df.meta$Age==age.i,]
    if(nrow(pick)<sub.sample){ # min cell 
      next
    }
    df=df.expr[,rownames(pick)];
    dim(df); #gene by cell
    #if(ncol(df)>sub.sample){
    #  set.seed(123)
    #  df=df[,sample(1:ncol(df),sub.sample,replace=F)]
    #}
    result=list()
    for(pathway.id in names(kegg.genes)){
      pathway.name=kegg.info[kegg.info$kegg.id==pathway.id,]$kegg.name
      genes=kegg.genes[[pathway.id]]$SYMBOL
      if(sum(genes %in% all.gene.names)<20){next} #minimal gene per pathway
      df1=df[all.gene.names %in% genes,]
      out=gcl(df1,num_divisions)
      result[[pathway.name]]=out
      cat(pathway.id,'at age',age.i,'is done\n')
    }
    gcl.per.pathway[[as.character(age.i)]]<-result;
  }
  saveRDS(gcl.per.pathway,outfile)
}  
    
   
gcl.per.pathway=readRDS(outfile)
sapply(gcl.per.pathway,length)

cell.type=rep(names(gcl.per.pathway),sapply(gcl.per.pathway,length))
df=data.frame(gcl=unlist(gcl.per.pathway))
df$cell.type=cell.type

x=df %>% group_by(cell.type) %>% summarise(median=median(gcl))
x=x[order(x$median),]
df$cell.type=factor(df$cell.type,levels=x$cell.type)

pdf(paste0(sex,'GCL_ranked_cell.type_',num_divisions,'_sub',sub.sample,'.pdf'),width=20)
ggplot(df,aes(x=cell.type,y=gcl))+geom_boxplot(outlier.shape = NA)+
  geom_jitter(size=0.001)+theme_bw()+
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1)
  )
dev.off()
