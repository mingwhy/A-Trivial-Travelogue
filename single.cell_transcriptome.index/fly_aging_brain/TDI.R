
library(Seurat)
library(ggplot2)
##
dnds=as.data.frame(data.table::fread('../gene.properties/dNdS/Dmel_Dsim_analysis_results_flydivas_v1.2.txt'))
head(dnds)
dnds[duplicated(dnds$id),]

dim(dnds);#10765
str(dnds)
summary(dnds$omega) 
sum(dnds$omega>2,na.rm=T) #16
dnds=dnds[dnds$omega<=2,]
dim(dnds)#10749    10
rownames(dnds)=dnds$id;

##
dat=readRDS('~/Documents/Data_fly_FCA/fly.brain.atlas/wholebrain_filtered_valid.rds')
dim(dat) # 12094 56192
## 
gene.meta=data.table::fread('~/Documents/Data_fly_FCA/fly.brain.atlas/gene.meta_brain.txt')
dim(gene.meta) #12094    11
gene.meta=as.data.frame(gene.meta)
rownames(gene.meta)=gene.meta$current_symbol
gene.meta=gene.meta[rownames(dat),]
dim(gene.meta) #12094    11

## overlap 
sum(gene.meta$current_symbol %in% rownames(dat))#12094
sum(gene.meta$FBID_KEY %in% dnds$id) #8409

##
dat=NormalizeData(dat)
#expr.mat=dat@assays$RNA@counts
expr.mat=dat@assays$RNA@data; #check out 2013 plos genetics paper
dim(expr.mat);sum(rownames(expr.mat)==rownames(gene.meta)) #12094
rownames(expr.mat)=gene.meta$FBID_KEY
# filter low-expr genes out, expr in >=5% cells
#i= Matrix::rowSums(expr.mat>0) > (0.05*ncol(expr.mat))
#sum(i) #4706 genes
#expr.mat=expr.mat[i,]

expr.mat=expr.mat[rownames(expr.mat) %in% dnds$id,]
dim(expr.mat) # 8409 56192

overlap.genes=intersect(rownames(expr.mat),dnds$id)
dnds=dnds[overlap.genes,]
dim(dnds); #8409   10
min(dnds$omega) #no 0

## apply to all cells
cell.meta=dat@meta.data
unique(cell.meta$annotation)
unique(cell.meta$sex)
unique(cell.meta$Age)

dim(expr.mat);dim(dnds);sum(dnds$id==rownames(expr.mat))

x1=as.numeric(dnds$omega %*% expr.mat)
#sum(dnds$omega * expr.mat[,1]);x1[1]
x2=as.numeric(Matrix::colSums(expr.mat))
index.per.cell=x1/x2;

cell.meta$TDI=index.per.cell;
cell.meta$Age=factor(cell.meta$Age,levels=sort(as.numeric(unique(cell.meta$Age))))

cell.type='Ensheathing_glia'
test1=subset(cell.meta,annotation==cell.type);
ggplot(test1,aes(x=Age,y=TDI))+geom_violin()+geom_jitter(size=0.2)+
  ggtitle(cell.type)+
  stat_summary(fun.y=median, geom="point", size=2, color="red")+
  theme_classic()

summary(cell.meta$nGene) #number of expr genes, minimal 275, max 4405

plots=lapply(unique(cell.meta$annotation),function(cell.type){
  test1=subset(cell.meta,annotation==cell.type);
  ggplot(test1,aes(x=Age,y=TDI))+geom_violin()+geom_jitter(size=0.2)+
    ggtitle(cell.type)+
    stat_summary(fun.y=median, geom="point", size=2, color="red")+
    theme_classic()
})
length(plots)
plots[[10]]

pdf('test_TDI.pdf',useDingbats = T)
for(i in seq(1,length(plots),9)){
  gridExtra::grid.arrange(grobs=plots[i:min((i+8),length(plots))],ncol=3)
}
dev.off();


##plot gene number
ngene=Matrix::colSums(expr.mat>0)
summary(ngene) #209 ~ 3658
cell.meta$ngene=ngene

plots=lapply(unique(cell.meta$annotation),function(cell.type){
  test1=subset(cell.meta,annotation==cell.type);
  ggplot(test1,aes(x=Age,y=ngene))+geom_violin()+geom_jitter(size=0.2)+
    ggtitle(cell.type)+
    stat_summary(fun.y=median, geom="point", size=2, color="red")+
    theme_classic()
})
length(plots)
plots[[10]]

pdf('test_nGene.pdf',useDingbats = T)
for(i in seq(1,length(plots),9)){
  gridExtra::grid.arrange(grobs=plots[i:min((i+8),length(plots))],ncol=3)
}
dev.off();

