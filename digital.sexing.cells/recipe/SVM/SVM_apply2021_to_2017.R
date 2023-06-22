
library(ggplot2)
library(gridExtra)
library(Seurat)
library(e1071)
library(caTools)
###############################
classifier=readRDS('svm_2021sc_train80samples.rds')
classifier.gene.names=colnames(classifier$SV)
length(classifier.gene.names) #6073 genes

###############################
tissue='embryo'
dat=readRDS('~/Documents/Data_fly_FCA/embryo/whole_embryo_filtered_valid.rds')
dat #9710 features across 1297 samples within 1 assay 
mat <- dat@assays$RNA@counts
gene.names <- rownames(mat)
cell.names <- colnames(mat)

gene.meta=data.table::fread('~/Documents/Data_fly_FCA/embryo/gene.meta_embryo.txt')
gene.meta=as.data.frame(gene.meta)
head(gene.meta)
dim(gene.meta) #9710

sum(gene.names %in% gene.meta$current_symbol)
rownames(gene.meta)<-gene.meta$current_symbol

rownames(mat)<-gene.meta[gene.names,]$FBID_KEY
gene.names=rownames(mat)

################################################
overlap.genes=gene.names[gene.names %in% classifier.gene.names]
length(overlap.genes) #5729

inp=(matrix(0,nrow=ncol(mat),ncol=length(classifier.gene.names)))
rownames(inp)=colnames(mat)
colnames(inp)=classifier.gene.names

inp[,overlap.genes]<-t(as.matrix(mat[overlap.genes,]))

test=scale(inp); 
dim(test) #sample by feature
sum(is.nan(test))
test[is.nan(test)]=0

pred_2017= predict(classifier, test,probability = TRUE)
table(pred_2017) #all correct
#embryoFemale   embryoMale 
#879          418 

########################################
dat #9710 features across 1297 samples within 1 assay 
mat <- dat@assays$RNA@counts
gene.names <- rownames(mat)
cell.names <- colnames(mat)

y.mat=mat[gene.meta[gene.meta$LOCATION_ARM=='Y',]$current_symbol,]
dim(y.mat) #15 y genes
y.expr=Matrix::colSums(y.mat)
table(y.expr)
#y.expr
#0    1    2    3    4    6 
#1075  171   36    7    6    2 

y.cells=names(which(Matrix::colSums(y.mat)>0))

table(pred_2017[y.cells])
#embryoFemale   embryoMale 
#135           87 

pred_2017[which(y.expr>2)]

##
preds_probs=attr(pred_2017,'probabilities')
head(preds_probs)
table(ifelse(preds_probs[,1]>0.5,'female','male'))

preds_probs=as.data.frame(as.matrix(preds_probs))
preds_probs$larger_prob=apply(preds_probs[,c(1,2)],1,max)
head(preds_probs)
preds_probs$pred_sex=ifelse(preds_probs[,1]>0.5,'embryoFemale','embryoMale')

preds_probs[y.cells,]


