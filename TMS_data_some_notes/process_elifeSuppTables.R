
library(ggplot2)
library(tidyverse)
library(Seurat)
## read in elife-62293-supp1-v2.xlsx
library(readxl)
df=read_xlsx('~/Documents/Data_mouse_aging_atlas/SuppTables1-3/elife-62293-supp1-v2.xlsx',
             sheet='all FACS tissue-cell')
colnames(df)
dim(df) #207 tc
head(df)
unique(df$tissue) #23
unique(df$cell_ontology_class) #120
tc=paste(df$tissue,df$cell_ontology_class,sep=':')
length(unique(tc)) #207
df$tc=tc;

min.cell=20
x=df[df$`3m`>=min.cell & df$`18m`>=min.cell & df$`24m`>=min.cell,]
dim(x) #115 tc
table(x$tissue)

## lifespan data
df.lifespan=read_xlsx('~/Documents/Data_mouse_aging_atlas/SuppTables1-3/elife-62293-supp1-v2.xlsx',
             sheet='76 TMS FACS tissue-cell')
colnames(df.lifespan)
dim(df.lifespan) #76 tc
head(df.lifespan)
sum(is.na(df.lifespan$binary_lifespan)) #48

df.lifespan$tc=paste(df.lifespan$tissue,df.lifespan$cell_ontology_class,sep=':')
sum(df.lifespan$tc %in% df$tc) #76 all exist

## cell functional category data
## more cell type annotation data (downstream analysis: https://github.com/czbiohub/tabula-muris-senis/blob/master/2_aging_signature/README.md)
meta2=data.table::fread('~/Documents/Data_mouse_aging_atlas/TMS.gene.data_final/annotation_data/cell_ontology_class_functional_annotation.073020.tsv')
dim(meta2) #237
head(meta2)
colnames(meta2)
meta2$tc=paste(meta2$tissue,meta2$cell_ontology_class,sep=':')

sum(df$tc %in% meta2$tc) #all 207 exist
table(meta2$dataset)

meta2=meta2[meta2$tc %in% df$tc,]
table(meta2$tissue)
table(meta2$`cell category`)

meta2[meta2$`cell category`=='stem cell/progenitor;muscle cell',] #check elife paper fig3A
meta2[meta2$`cell category`=='parenchymal;epithelial',]

meta2$cell.category=sapply(strsplit(meta2$`cell category`,';'),'[',1)
table(meta2$cell.category)  # 6 cell functional classes

table(meta2[grep('Brain',meta2$tissue),]$cell.category)

x=df[df$`3m`>=min.cell & df$`18m`>=min.cell & df$`24m`>=min.cell,]
meta2=meta2[meta2$tc %in% x$tc,]
table(meta2$tissue)
table(meta2$cell.category)
tmp=meta2 %>% group_by(tissue,cell.category) %>% summarise(n=n())

(mycol=RColorBrewer::brewer.pal(8,"Dark2"))
pdf('tc115_distribution.pdf',useDingbats = T,width = 16)
ggplot(tmp,aes(x=tissue,y=n,fill=cell.category))+
  geom_bar(stat='identity')+theme_classic()+ylab('Number of cell.types')+xlab('Tissue')+
  scale_fill_manual(values=mycol)+
  theme(axis.text.x = element_text(angle=45,hjust=0.5,vjust = 0.6))
dev.off()


