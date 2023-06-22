#stem cell, long-lived cell, short-lived cell.
#x: longevity, y: magnitude of aging 

# calculate magnitude of aging for each tc (young-old / mean(young-young, old-old))

library(ggplot2);library(gridExtra)
library(dplyr);library(ggpubr)
library(Seurat);
library(grDevices);library(RColorBrewer)
library(transport)
plotcol <- brewer.pal(8,'Dark2')

all.OT.distances=readRDS('nmf_73tc_OT.dist.rds')
length(all.OT.distances)
names(all.OT.distances)

dist.norm<-lapply(names(all.OT.distances),function(tc){
  x=all.OT.distances[[tc]];
  ds=sapply(1:300,function(i){
    x1=x[c(i,i+300,i+600),2]
    d=x1[3]/max(x1[2],x1[1])
    d
  })
  #median(ds)
  ds
})
names(dist.norm)=names(all.OT.distances)
df=data.frame(OT.dist.norm=unlist(dist.norm),tc=rep(names(dist.norm),sapply(dist.norm,length)))
head(df)
x=df %>% group_by(tc) %>% summarise(median.dist=median(OT.dist.norm))
x=x[order(x$median.dist),]
df$tc=factor(df$tc,levels=x$tc)

ggplot(df,aes(x=tc,y=OT.dist.norm))+geom_violin()+geom_jitter(size=0.2)+
  theme_classic()+coord_flip()


## read in cell function annotation to color labels
## more cell type annotation data (downstream analysis: https://github.com/czbiohub/tabula-muris-senis/blob/master/2_aging_signature/README.md)
meta2=data.table::fread('~/Documents/Data_mouse_aging_atlas/TMS.gene.data_final/annotation_data/cell_ontology_class_functional_annotation.073020.tsv')
dim(meta2) #237
head(meta2)
colnames(meta2)

table(meta2$tissue); #except 'Heart_and_Aorta', 23 tissues in total, consistent with 2021 elife
table(meta2$`cell category`) # 6 cell functional classes

meta2[meta2$`cell category`=='stem cell/progenitor;muscle cell',] #check elife paper fig3A
meta2[meta2$`cell category`=='parenchymal;epithelial',]

meta2$cell.category=sapply(strsplit(meta2$`cell category`,';'),'[',1)
table(meta2$cell.category)  # 6 cell functional classes

table(meta2[grep('Brain',meta2$tissue),]$cell.category)

meta2$tc=paste(meta2$tissue,meta2$cell_ontology_class,sep=':')
unique(meta2$tc) #237


dfc=merge(df,meta2)
head(dfc)
pdf('nmf_73tc_OT.dist.norm.pdf',useDingbats = T,height = 8,width = 7)
print( ggplot(dfc,aes(x=tc,y=OT.dist.norm,col=cell.category))+geom_violin()+geom_jitter(size=0.2)+
  theme_classic()+coord_flip()+scale_color_manual(values=sample(plotcol)) )
dev.off()

tmp=dfc[!duplicated(dfc$tc),]
table(tmp$cell.category)

print( ggplot(subset(dfc,cell.category=='immune'),aes(x=tc,y=OT.dist.norm,col=cell.category))+geom_violin()+geom_jitter(size=0.2)+
         theme_classic()+coord_flip()+scale_color_manual(values=plotcol) )

table(tmp$tissue)
print( ggplot(subset(dfc,tissue %in% c('Heart','Large_Intestine','Pancreas','Liver')),
              aes(x=tc,y=OT.dist.norm,col=tissue))+geom_violin()+geom_jitter(size=0.2)+
         theme_classic()+coord_flip()+scale_color_manual(values=plotcol) )
