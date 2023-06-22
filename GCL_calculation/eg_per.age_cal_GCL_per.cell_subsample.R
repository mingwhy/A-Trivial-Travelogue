
library(Seurat)
library(tidyverse)
library(Matrix);
library(ggplot2);library(gridExtra)
library(RColorBrewer)
source("./gcl_r_version.R")

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
sex='female'
table(dat@meta.data$sex)
i=dat@meta.data$sex==sex;
sum(i) #28850 female cells
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

sub.sample=20;
num_divisions=50; #num_divisions in subsetting genes into 2 sets when calcualte GCL 
#outfile=paste0(sex,'_GCL_',num_divisions,'_sub',sub.sample,'.rds');
(outfile=paste0('GCL_div',num_divisions,'_min',sub.sample,'.rds'));

if(!file.exists(outfile)){
  ## compute cell type specific E-R anticorrelation and plot against #detected.genes
  gcl.per.cell.type=list()
  for(cell.type in cell.types){ 
    cell.type2=stringr::str_replace_all(cell.type, '\\/', '\\-')
    
    df=df.expr[,df.meta$annotation==cell.type];
    dim(df); #gene by cell
    #if(ncol(df)>sub.sample){
    #  set.seed(123)
    #  df=df[,sample(1:ncol(df),sub.sample,replace=F)]
    #}
    result=list()
    for(age.i in unique(df.meta[df.meta$annotation==cell.type,]$Age)){
      pick=df.meta[df.meta$annotation==cell.type &df.meta$Age==age.i,]
      if(nrow(pick)<sub.sample){
        #result[[age.i]]=NA;
        next
      }
      df1=df[,rownames(pick)]
      out=gcl(df1,num_divisions)
      result[[as.character(age.i)]]=out
      cat(cell.type,'at age',age.i,'is done\n')
    }
    gcl.per.cell.type[[cell.type]]<-result;
  }
  saveRDS(gcl.per.cell.type,outfile)
}  
    
   
gcl.per.cell.type=readRDS(outfile)

sapply(gcl.per.cell.type,length)
gcl_cell.type=list()
for(cell.type.name in names(gcl.per.cell.type)){
  x=gcl.per.cell.type[[cell.type.name]]
  ages=rep(as.numeric(names(x)),sapply(x,length))
  gcl=unlist(x)
  x1=data.frame(cell.type=rep(cell.type.name,length(gcl)),
                ages=ages,gcl=gcl)
  gcl_cell.type[[cell.type.name]]=x1
}
df.gcl_cell.type=Reduce(`rbind`,gcl_cell.type)
dim(df.gcl_cell.type)

df.gcl_cell.type$ages=factor(df.gcl_cell.type$ages,levels=sort(unique(df.gcl_cell.type$ages)))

pdf(paste0('per.age_GCL_div',num_divisions,'_min',sub.sample,'.pdf'),width=12)
ggplot(df.gcl_cell.type,aes(x=ages,y=gcl))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(size=0.001)+theme_bw()+
  facet_wrap(~cell.type,ncol=5)+
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1)
  )
dev.off()
