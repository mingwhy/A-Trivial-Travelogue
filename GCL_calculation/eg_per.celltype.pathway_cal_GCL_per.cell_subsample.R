
library(Seurat)
library(tidyverse)
library(Matrix);
library(ggplot2);library(gridExtra)
library(RColorBrewer)
source("./gcl_r_version.R")

## read in pathway data
## download dme-kegg.ID.df.txt and kegg-flygenes.rds from 
# https://github.com/mingwhy/bioinfo_homemade_tools/tree/main/KEGG.decompose
kegg.info=read.table('../gene.age_per.pathway/dme-kegg.ID.df.txt',as.is=T,sep="\t",header=T)
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


# plot #cell over cell types, mix sex, pick cell.type>=#200 cells
cell.type=df.meta$annotation
cell.type.n=table(cell.type)
cell.type.status=data.frame(cell.type=names(cell.type.n),number.of.cells=as.numeric(cell.type.n),
                            status=ifelse(cell.type.n>=200,1,0))
head(cell.type.status)

# unknown cluster VS unkonw cell type
cell.type.status$cell.type.known='no';
cell.type.status[grep('[a-zA-Z]',cell.type.status$cell.type),]$cell.type.known='yes'
table(cell.type.status$status,cell.type.status$cell.type.known)

df=cell.type.status[grep('[a-zA-Z]',cell.type.status$cell.type),] #74 known cell type 
pick.cell.type=cell.type.status[cell.type.status$cell.type.known=='yes' & 
                                  cell.type.status$number.of.cells>=200,]
dim(pick.cell.type) #33 cell types
sum(pick.cell.type$number.of.cells) # 17499 cells in total

cell.types=as.character(pick.cell.type$cell.type)
tmp=df.meta[df.meta$annotation %in% cell.types,]
tmp1=table(tmp$annotation,tmp$Age)
length(tmp1)
sum(tmp1>=20) #203 items. the original study, min.cell=22

###########################################################
# cal GCL for each pathway in each cell type individually
(ages=sort(unique(df.meta$Age)))

sub.sample=20;
num_divisions=50; #num_divisions in subsetting genes into 2 sets when calcualte GCL 
#outfile=paste0(sex,'_GCL_',num_divisions,'_sub',sub.sample,'.rds');

dir='./cell.type_pathway_div50/';
if(!dir.exists(dir)){dir.create(dir)}

gcl_cell.type_pathway=list();
all.gene.names=rownames(df.expr)
for(cell.type in cell.types){
  cell.type2=stringr::str_replace_all(cell.type, '\\/', '\\-')
  (filename=paste0(dir,'celltype_',cell.type2,'.rds'))
  if(!file.exists(filename)){  
    df0=df.expr[,df.meta$annotation==cell.type];
    dim(df0); #gene by cell
    #if(ncol(df)>sub.sample){
    #  set.seed(123)
    #  df=df[,sample(1:ncol(df),sub.sample,replace=F)]
    #}
    cell.type.gcl=list();
    for(age.i in ages){
      pick=df.meta[df.meta$Age==age.i & df.meta$annotation==cell.type,]
      if(nrow(pick)<sub.sample){ # min cell 
        next
      }
      df=df0[,rownames(pick)];
      dim(df); #gene by cell
      #if(ncol(df)>sub.sample){
      #  set.seed(123)
      #  df=df[,sample(1:ncol(df),sub.sample,replace=F)]
      #}
      result=list();
      for(pathway.id in names(kegg.genes)){
        pathway.name=kegg.info[kegg.info$kegg.id==pathway.id,]$kegg.name
        genes=kegg.genes[[pathway.id]]$SYMBOL
        if(sum(genes %in% all.gene.names)<20){next} #minimal gene per pathway
        df1=df[all.gene.names %in% genes,]
        out=gcl(df1,num_divisions)
        result[[pathway.name]]=out
        cat(cell.type,pathway.id,'at age',age.i,'is done\n')
      }
      cell.type.gcl[[as.character(age.i)]]<-result;
    }
    saveRDS(cell.type.gcl,filename)
  }
  cell.type.gcl=readRDS(filename)
  gcl_cell.type_pathway[[cell.type]]=cell.type.gcl
}  

for(cell.type.name in names(gcl_cell.type_pathway)){
  cell.type2=stringr::str_replace_all(cell.type.name, '\\/', '\\-')
  x=gcl_cell.type_pathway[[cell.type.name]]
  df.list=list()
  for(age in names(x)){
    df=data.frame(pathway=rep(names(x[[age]]),sapply(x[[age]],length)),gcl=unlist(x[[age]]))
    df$age=rep(age,nrow(df))
    df.list[[age]]=df
  }
  df=Reduce(`rbind`,df.list)
  df[is.infinite(df$gcl),]
  
  df$age=as.numeric(df$age)
  df$age=factor(df$age,levels=sort(unique(df$age)))
  pdf(paste0(dir,'celltype_',cell.type2,'.pdf'),width=20,height = 20)
  print(ggplot(df,aes(x=age,y=gcl))+geom_boxplot(outlier.shape = NA)+
    geom_jitter(size=0.001)+theme_bw()+
      facet_wrap(~pathway,ncol=5)+
      theme(
        panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1)
      )
  )
  dev.off()
}
