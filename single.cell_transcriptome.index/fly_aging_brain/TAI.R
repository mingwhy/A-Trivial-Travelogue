
library(Seurat)
library(ggplot2)
##
age=as.data.frame(data.table::fread('./gene.age/Drosophila_melanogaster.csv'))
rownames(age)=age$ensembl_id
age[duplicated(age$ensembl_id),]
dim(age) #13931 genes

age[grep('>',age$gene_age),]
max(as.numeric(age$gene_age),na.rm=T) #1934
length(grep('>',age$gene_age)) #790, remove these genes
age=age[-grep('>',age$gene_age),]
dim(age) #9806
age$gene_age=as.numeric(age$gene_age)

##
dat=readRDS('~/Documents/single.cell_datasets/fly.brain.atlas/wholebrain_filtered_valid.rds')
dim(dat) # 12094 56192
## 
gene.meta=data.table::fread('~/Documents/single.cell_datasets/fly.brain.atlas/gene.meta_brain.txt')
dim(gene.meta) #12094    11
gene.meta=as.data.frame(gene.meta)
rownames(gene.meta)=gene.meta$current_symbol
gene.meta=gene.meta[rownames(dat),]
dim(gene.meta) #12094    11

## overlap 
sum(gene.meta$current_symbol %in% rownames(dat))#12094
sum(gene.meta$FBID_KEY %in% age$ensembl_id) #10596

##
dat=NormalizeData(dat)
#expr.mat=dat@assays$RNA@counts
expr.mat=dat@assays$RNA@data; #check out 2013 plos genetics paper
dim(expr.mat);sum(rownames(expr.mat)==rownames(gene.meta)) #9806
rownames(expr.mat)=gene.meta$FBID_KEY

expr.mat=expr.mat[rownames(expr.mat) %in% age$ensembl_id,]
dim(expr.mat) #9806 56192

overlap.genes=intersect(rownames(expr.mat),age$ensembl_id)
age=age[overlap.genes,]
dim(age); #9806   10

## apply to all cells
cell.meta=dat@meta.data
unique(cell.meta$annotation)
unique(cell.meta$sex)
unique(cell.meta$Age)

dim(expr.mat);dim(age);sum(age$ensembl_id==rownames(expr.mat))

x1=as.numeric(age$gene_age %*% expr.mat)
#sum(age$gene_age * expr.mat[,1]);x1[1]
x2=as.numeric(Matrix::colSums(expr.mat))
index.per.cell=x1/x2;

cell.meta$TI=index.per.cell;
cell.meta$Age=factor(cell.meta$Age,levels=sort(as.numeric(unique(cell.meta$Age))))

cell.type='Ensheathing_glia'
test1=subset(cell.meta,annotation==cell.type);
ggplot(test1,aes(x=Age,y=TI))+geom_violin()+geom_jitter(size=0.2)+
  ggtitle(cell.type)+
  stat_summary(fun.y=median, geom="point", size=2, color="red")+
  theme_classic()

summary(cell.meta$nGene) #number of expr genes, minimal 275, max 4405

plots=lapply(unique(cell.meta$annotation),function(cell.type){
  test1=subset(cell.meta,annotation==cell.type);
  ggplot(test1,aes(x=Age,y=TI))+geom_violin()+geom_jitter(size=0.2)+
    ggtitle(cell.type)+
    stat_summary(fun.y=median, geom="point", size=2, color="red")+
    theme_classic()
})
length(plots)
plots[[10]]

pdf('test_TAI.pdf',useDingbats = T)
for(i in seq(1,length(plots),9)){
  gridExtra::grid.arrange(grobs=plots[i:min((i+8),length(plots))],ncol=3)
}
dev.off();


