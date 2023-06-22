
library(ggplot2)
library(Seurat)
library(gridExtra)
library(tidyverse)
library(rsample)      # data splitting 
library(randomForest) # basic implementation
library(ranger)       # a faster implementation of randomForest
library(caret)        # an aggregator package for performing many machine learning models
library(h2o)          # an extremely fast java-based platform
library(multiclassPairs)
###############################
classifier=readRDS('sc_classifier_train80.rds')

filtered_genes=readRDS('sc_classifier_filtered_genes.rds')
length(filtered_genes$OnevsrestScheme$filtered_genes$embryoFemale)
length(filtered_genes$OnevsrestScheme$filtered_genes$embryoMale)
genes=unique(c(filtered_genes$OnevsrestScheme$filtered_genes$embryoFemale,filtered_genes$OnevsrestScheme$filtered_genes$embryoMale))
length(genes) #6011

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
  expr.mat=x@assays$RNA@data
  
  test=as.matrix(expr.mat) #gene by cell matrix
  class(test)
  rownames(test)=gene.meta[rownames(test),]$FLYBASE #change to FBgn 
  
  # apply on the testing data
  results_test <- predict_one_vs_rest_TSP(classifier = classifier,
                                          Data = test,
                                          tolerate_missed_genes = TRUE,
                                          weighted_votes = TRUE,
                                          classes = c("embryoFemale",'embryoMale'),
                                          verbose = TRUE)
  table(results_test$max_score)
  head(results_test)
  
  
  time=gsub('_finished_processing.rds','',basename(file))
  out[[time]]=results_test;
  cat('time',time,'finished\n')
}
saveRDS(out,'TSP_apply2021sc_to_Jay.out.rds')

end.time=Sys.time()
print(end.time-start.time) #7min

out=readRDS('TSP_apply2021sc_to_Jay.out.rds')
names(out)
out=out[order(as.numeric(lapply(strsplit(names(out),'\\-'),'[[',1)))]

pdf('test.pdf',useDingbats = T)
par(mfrow=c(3,4))
for(i in 1:length(out)){
  preds=out[[i]]
  #head(preds)
  larger.probs=apply(preds[,c(1,2)],1,max)
  hist(larger.probs,main=names(out)[i],xlim=c(0.5,1))
  #p=hist(larger.probs,main=names(out)[i],xlim=c(0.5,1),type='N')
  #p$freq=p$counts/sum(p$counts)
  #p$counts=p$freq;
  #plot(p,main=names(out)[i],xlim=c(0.5,1))  
}
dev.off();


lapply(out,function(i){
  preds=i[,c(1,2)]
  #preds=attr(i,'probabilities')
  predict.labels=ifelse(preds[,1]>0.5,'embryoFemale','embryoMale')
  table(predict.labels)
})
