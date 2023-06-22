
library(Seurat)
library(ggplot2)
## PPI
ppi=data.table::fread('./PPI/ensembl_fly_ppi-cutoff400_link.txt')
dim(ppi);#854428
head(ppi)
# is there any duplicate rows
x=apply(ppi[,c(4,5)],1,function(i) {paste(sort(i),collapse = '-')})
head(x);length(x)
head(x[duplicated(x)])
x[x=="FBgn0031082-FBgn0031085"]
ppi[which(x=="FBgn0031082-FBgn0031085"),]

duplicate.rows=duplicated(x)
unique.ppi=ppi[!duplicate.rows,]
dim(unique.ppi) #427214

gene.degree=as.data.frame(table(c(unique.ppi$Ensembl_gene1, unique.ppi$Ensembl_gene2)))
head(gene.degree)
colnames(gene.degree)=c('id','degree')
dim(gene.degree) #11989
rownames(gene.degree)=gene.degree$id
summary(gene.degree$degree)
#gene.degree$degree=log10(gene.degree$degree) #may be the same issue as 'some high express genes'.
#summary(gene.degree$degree) #make min and max less than 5 fold 

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
sum(gene.meta$FBID_KEY %in% gene.degree$id) #9497

##
dat=NormalizeData(dat)
#expr.mat=dat@assays$RNA@counts
expr.mat=dat@assays$RNA@data; #check out 2013 plos genetics paper
dim(expr.mat);sum(rownames(expr.mat)==rownames(gene.meta)) #9806
rownames(expr.mat)=gene.meta$FBID_KEY

expr.mat=expr.mat[rownames(expr.mat) %in% gene.degree$id,]
dim(expr.mat) #9497 56192

overlap.genes=intersect(rownames(expr.mat),gene.degree$id)
gene.degree=gene.degree[overlap.genes,]
dim(gene.degree); #9497

## apply to all cells
cell.meta=dat@meta.data
unique(cell.meta$annotation)
unique(cell.meta$sex)
unique(cell.meta$Age)

dim(expr.mat);dim(gene.degree);sum(gene.degree$id==rownames(expr.mat))

x1=as.numeric(gene.degree$degree %*% expr.mat)
#sum(gene.degree$degree * expr.mat[,1]);x1[1]
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

pdf('test_TPI.pdf',useDingbats = T)
for(i in seq(1,length(plots),9)){
  gridExtra::grid.arrange(grobs=plots[i:min((i+8),length(plots))],ncol=3)
}
dev.off();


