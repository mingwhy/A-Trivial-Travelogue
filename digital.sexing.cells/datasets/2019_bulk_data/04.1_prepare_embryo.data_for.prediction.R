
library(ggplot2)
library(gridExtra);
library(Seurat);
########################################################
## process fly embryo data, read in gene chro info
df.gene.table=data.table::fread('../../../single.cell_datasets/embryo/gene.meta_embryo.txt',header=T,sep='\t')
head(df.gene.table)
table(df.gene.table$LOCATION_ARM) #15 Y-chromosome genes

df.gene.table=df.gene.table[df.gene.table$LOCATION_ARM %in% c('2R','3R','3L','4','2L','X','Y'),]
table(df.gene.table$LOCATION_ARM)

df.gene.table=as.data.frame(df.gene.table)
rownames(df.gene.table)=df.gene.table$current_symbol

########################################################
## read in dataset
#dat=readRDS('../../single.cell_datasets/FCA_wing/whole_wing_filtered_valid.rds')
#dat=readRDS('../../single.cell_datasets/FCA_heart/whole_heart_filtered_valid.rds')
tissue='embryo'
dat=readRDS('../../../single.cell_datasets/embryo/whole_embryo_filtered_valid.rds')
#dat=subset(dat,annotation=="spermatocyte 3")

mat <- dat@assays$RNA@counts
gene.names <- rownames(mat)
cell.names <- colnames(mat)
#cell.types<-dat$annotation
#sex.labels<-dat$sex
sex='unknown'
cell.type='unknown';

##
dim(mat) #9710 gene x 1297 cells
length(gene.names)
length(cell.names)
#table(cell.types)
#table(sex.labels)
rownames(mat)=gene.names
colnames(mat)=cell.names
mat[1:3,1:3]

umi.mat=as.matrix(mat)
# filter out some genes
if(F){
  filter <- rowSums(umi.mat>2)>=5
  table(filter) 
  umi.mat<- umi.mat[filter,]
  dim(umi.mat) #5779 1297
}
if(F){
  seurat.obj=CreateSeuratObject(count=umi.mat)
  seurat.obj=NormalizeData(seurat.obj)
  umi.mat=seurat.obj@assays$RNA@data #normalized data
}
cell.size=Matrix::colSums(umi.mat)
tmp=sapply(1:ncol(umi.mat), function(i) umi.mat[,i]/cell.size[i] )
umi.mat=log(tmp*10^6+1,base=2) #log2TPM
colnames(umi.mat)=colnames(mat)

########################################################
# use gene chr info to filter genes
# (some genes annotated to non A or X chr, discard these genes)
i=rownames(umi.mat) %in% df.gene.table$current_symbol
sum(i);#9705
umi.mat=umi.mat[i,]
dim(umi.mat) #9705 1297
########################################################
pick.genes1=df.gene.table[grep('msl|mof|mle|roX',df.gene.table$SYMBOL),]$SYMBOL
pick.genes2=df.gene.table[df.gene.table$LOCATION_ARM=='Y',]$SYMBOL
pick.genes3=c('Sxl','run','sisA')

male.cells.by.chrY= names(which(Matrix::colSums(umi.mat[c(pick.genes2),])>0)) #222 cells
length(male.cells.by.chrY) #222
male.cells.by.roX2= names(which(Matrix::colSums(umi.mat['lncRNA:roX2',,drop=F])>0)) #222 cells
length(male.cells.by.roX2) #79
male.cells=unique(c(male.cells.by.roX2,male.cells.by.chrY))
length(male.cells) #281


# due to gene name with '-' are not allowed in `multiclassPairs`
original.gene.names=rownames(umi.mat)
new.gene.names=gsub('-','_',original.gene.names)
rownames(umi.mat)=new.gene.names

sample.meta=data.frame(barcode=colnames(umi.mat),sex='unknown')
rownames(sample.meta)=sample.meta$barcode
#sample.meta[male.cells,]$sex='male'
sample.meta[male.cells.by.chrY,]$sex='male.Y'
sample.meta[male.cells.by.roX2,]$sex='male.roX2'
saveRDS(list(umi.mat=umi.mat,sample.meta=sample.meta),'prepared_embryo.data.rds')





