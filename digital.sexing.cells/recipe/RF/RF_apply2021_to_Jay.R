
library(ggplot2)
library(Seurat)
library(gridExtra)
library(tidyverse)
library(rsample)      # data splitting 
library(randomForest) # basic implementation
library(ranger)       # a faster implementation of randomForest
library(caret)        # an aggregator package for performing many machine learning models
library(h2o)          # an extremely fast java-based platform
###############################
classifier=readRDS('rf_ranger_normData_500tree_80train.rds');
#classifier=readRDS('rf_ranger_normData_1000tree_80train.rds');
classifier.gene.names=names(classifier$variable.importance)
length(classifier.gene.names) #6073 genes


###############################
# read in gene.meta info which contains chr info
gene.meta=data.table::fread('../validated_15392genes.txt')
y.genes=gene.meta[gene.meta$chromosome_name=='Y',]$SYMBOL
length(y.genes) #34 y genes
gene.meta=as.data.frame(gene.meta)
rownames(gene.meta)=gene.meta$SYMBOL

# read in embryo data by pred time window
(files=Sys.glob('../RNA_seurat_object/pred_windows/*.rds'))

start.time=Sys.time()
out=c()
for(file in files){
  #x=readRDS('./RNA_seurat_object/pred_windows/18-20hrs_finished_processing.rds')
  x=readRDS(file)
  x=NormalizeData(x)
  expr.mat=x@assays$RNA@data
  
  tmp.gene.names=rownames(expr.mat)
  rownames(expr.mat)=gene.meta[tmp.gene.names,]$FLYBASE
  tmp.gene.names=rownames(expr.mat)
  
  overlap.genes=tmp.gene.names[tmp.gene.names %in% classifier.gene.names]
  length(overlap.genes) #6071
  
  inp=(matrix(0,nrow=ncol(expr.mat),ncol=length(classifier.gene.names)))
  rownames(inp)=colnames(expr.mat)
  colnames(inp)=classifier.gene.names
  
  inp[,overlap.genes]<-t(as.matrix(expr.mat[overlap.genes,]))
  
  #test=scale(inp); #do not scale in RF
  test=inp;
  dim(test) #sample by feature
  sum(is.nan(test))
  test[is.nan(test)]=0
  
  preds = predict(classifier, test,probability = TRUE)
  head(preds$predictions)##probabilities
  
  time=gsub('_finished_processing.rds','',basename(file))
  out[[time]]=preds;
  cat('time',time,'finished\n')
}
#saveRDS(out,'RF_normData_1000tree_apply2021_to_Jay.out.rds')
saveRDS(out,'RF_normData_500tree_apply2021_to_Jay.out.rds')

end.time=Sys.time()
print(end.time-start.time) #4min

names(out)
out=out[order(as.numeric(lapply(strsplit(names(out),'\\-'),'[[',1)))]

#pdf('result_RF_normData_1000tree_apply2021_to_Jay.out.pdf',useDingbats = T)
pdf('result_RF_normData_500tree_apply2021_to_Jay.out.pdf',useDingbats = T)
par(mfrow=c(3,4))
for(i in 1:length(out)){
  preds=out[[i]]
  names(preds)
  head(preds$predictions)
  larger.probs=apply(preds$predictions,1,max)
  hist(larger.probs,main=names(out)[i],xlim=c(0.5,1))
}
dev.off();

predict.labels=ifelse(preds$predictions[,1]>0.5,'embryoFemale','embryoMale')
table(predict.labels)

lapply(out,function(i){
  predict.labels=ifelse(i$predictions[,1]>0.5,'embryoFemale','embryoMale')
  table(predict.labels)
})
      