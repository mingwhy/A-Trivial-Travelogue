
library(ggplot2)
library(gridExtra)

########################################################################################
## read in gene expression data and sample meta information 
expr.mat=data.table::fread('../external_data/2019_paper_reproduce.result/gene_by_sample_log2TPM.txt')
expr.mat=as.matrix(expr.mat,rownames=1)

sample.meta=data.table::fread('../external_data/2019_paper_reproduce.result/sample.meta_sex.label.txt')

dim(expr.mat) #8934, 54
dim(sample.meta) #54, 4
sum(colnames(expr.mat)==sample.meta$GSM.id) #54


# pick female and male samples
expr.mat=as.data.frame(expr.mat)
sample.meta=as.data.frame(sample.meta)
expr.mat=expr.mat[,sample.meta$cluster!='PB']
sample.meta=sample.meta[sample.meta$cluster!='PB',]

dim(expr.mat) #gene by sample
dim(sample.meta) #sample by attributes

gene.names.2019=rownames(expr.mat)

## read in 2011 data to subset genes
dat.2011=data.table::fread('../external_data/2011_paper_data/2011_log2.RPKM.txt')
dim(dat.2011) #12353    24

gene.names=dat.2011$gene.name
overlap.genes=intersect(gene.names.2019,gene.names)
length(overlap.genes) #8215

expr.mat=expr.mat[overlap.genes,]
dim(expr.mat) #gene by sample
dim(sample.meta) #sample by attributes

data = data.frame(t(expr.mat), y =sample.meta$cluster)
head(data)
tail(colnames(data)) #last column is sex label

# create training and valida
library(e1071)
library(caTools)
gene.names=colnames(data)
length(unique(gene.names))# 8935
dim(data) # 49 8935

dataset=data
dim(dataset) #49 sample x 5 feature gene
#dataset=data[,-3] #remove roX2
colnames(dataset)

# use all 2019 data to construct classifier Predicting the Test set results
if(F){
  dim(dataset)
  scaled.dataset= scale(dataset[-label.column])
  scaled.dataset=as.data.frame(scaled.dataset)
  scaled.dataset=cbind(scaled.dataset,dataset$y)
  colnames(scaled.dataset)[label.column]='y'
  scaled.dataset$y=as.factor(scaled.dataset$y)
  classifier = svm(formula = y ~ .,
                   data = scaled.dataset,
                   type = 'C-classification',
                   kernel = 'linear')
  y_pred = predict(classifier, newdata = scaled.dataset)
  table(scaled.dataset[,label.column],y_pred)
  classifier
  #saveRDS(classifier,'svm_classifier_2019_all49samples.rds')
}

if(F){
  
  set.seed(321) 
  split = sample.split(dataset$y, SplitRatio = 0.8)
  training_set = subset(dataset, split == TRUE)
  test_set = subset(dataset, split == FALSE)
  (label.column=ncol(dataset)) #5nd column is the cluster label
  
  # Feature Scaling for train and test separately
  training_set[-label.column] = scale(training_set[-label.column])
  test_set[-label.column] = scale(test_set[-label.column])
  training_set$y=factor(training_set$y)
  classifier = svm(formula = y ~ .,
                   data = training_set,
                   type = 'C-classification',
                   probability = TRUE,
                   kernel = 'linear')
          #gamma=0.05,kernel = 'radial', cost=10)
  y_pred = predict(classifier, newdata = test_set[-label.column])
  #y_pred = predict(classifier, newdata = test_set[-label.column],prob=T)
  table(test_set[,label.column],y_pred)
  #classifier
  #saveRDS(classifier,'svm_classifier_2019_train34samples.rds')
  saveRDS(classifier,'svm_classifier_2019_train80samples.rds')
}
classifier=readRDS('svm_classifier_2019_train80samples.rds')
#https://stackoverflow.com/questions/34781495/how-to-find-important-factors-in-support-vector-machine
w<-t(classifier$coefs) %*% classifier$SV # weight vectors
w <- apply(w, 2, function(v){sqrt(sum(v^2))})  # weight
w <- sort(w, decreasing = T)
length(w) #8934
head(w)

library(org.Dm.eg.db)
dim(dataset) #49 x 5
top.genes=names(w)[1:10]
x=AnnotationDbi::select(org.Dm.eg.db,keys=top.genes,
                        keytype='FLYBASE',columns=c('SYMBOL'))
x$SYMBOL  

gene.names.2019=names(w)
length(gene.names.2019) #8215
#####################################################
## apply RF model to other bulk dataset
# gene name match (FBgn), expr.mat format
dat.2015=data.table::fread('../external_data/2015_paper_data/2015_FPKM_normalization_6mel.txt')

gene.names=dat.2015$gene.name
dat.2015=as.data.frame(dat.2015[,-1])
rownames(dat.2015)=gene.names
dim(dat.2015) #6003    6
dat.2015[1:3,1:3]

sum(gene.names %in% gene.names.2019) #5011

inp=(matrix(0,nrow=ncol(dat.2015),ncol=length(gene.names.2019)))
rownames(inp)=colnames(dat.2015)
colnames(inp)=gene.names.2019
overlap.genes=gene.names[gene.names %in% gene.names.2019]
inp[,overlap.genes]<-t(dat.2015[overlap.genes,])

test=scale(inp); 
dim(test) #sample by feature
test[is.nan(test)]=0
pred_2015 = predict(classifier, test)
pred_2015 #all correct

##########################################
## read in 2011 data
dat.2011=data.table::fread('../external_data/2011_paper_data/2011_log2.RPKM.txt')
dim(dat.2011) #12353    24

gene.names=dat.2011$gene.name
dat.2011=as.data.frame(dat.2011[,-1])
rownames(dat.2011)=gene.names
dim(dat.2011) #6003    6
dat.2011[1:3,1:3]
#dat.2011=dat.2011[,grep('13|14',colnames(dat.2011))] #limit to stage14
apply(dat.2011,2,function(i) sum(is.na(i)))
# remove M14D M14D_r2 
ncol(dat.2011)
dat.2011=dat.2011[,1:22]
sum(gene.names %in% gene.names.2019) #8215

inp=(matrix(0,nrow=ncol(dat.2011),ncol=length(gene.names.2019)))
rownames(inp)=colnames(dat.2011)
colnames(inp)=gene.names.2019
overlap.genes=gene.names[gene.names %in% gene.names.2019]
inp[,overlap.genes]<-t(dat.2011[overlap.genes,])

test=scale(inp); 
dim(test) #sample by feature
test[is.nan(test)]=0
pred_2011 = predict(classifier, test)
pred_2011 #all correct
table(substr(names(pred_2011),0,1),pred_2011)
# late stage are more likely to be correctly labelled

##########################################
## read in 2020 data
dat_2020=data.table::fread('../external_data/2020_paper_data/expression-noise-across-fly-embryogenesis-master/150samples_retained.txt')
dim(dat_2020) #8004  151

gene.names=dat_2020$Ensembl.Gene.ID
dat_2020=as.data.frame(dat_2020[,-1])
rownames(dat_2020)=gene.names
dim(dat_2020) #8004  150
dat_2020[1:3,1:3]

sum(gene.names %in% gene.names.2019) #6530

inp=(matrix(0,nrow=ncol(dat_2020),ncol=length(gene.names.2019)))
rownames(inp)=colnames(dat_2020)
colnames(inp)=gene.names.2019
overlap.genes=gene.names[gene.names %in% gene.names.2019]
inp[,overlap.genes]<-t(dat_2020[overlap.genes,])

test=scale(inp); 
dim(test) #sample by feature
test[is.nan(test)]=0
pred_2020 = predict(classifier,test,prob=TRUE)
table(pred_2020)
# female   male 
#  84     66
#probplot(pred_2020)
pred_2020_probs=attr(pred_2020,'probabilities')
head(pred_2020_probs)
table(ifelse(pred_2020_probs[,1]>0.5,'female','male'))
# female   male 
#  84     66

pred_2020_probs=as.data.frame(as.matrix(pred_2020_probs))
table(substr(rownames(pred_2020_probs),0,2))
pred_2020_probs$stage=substr(rownames(pred_2020_probs),0,2)
pred_2020_probs$larger_prob=apply(pred_2020_probs[,c(1,2)],1,max)
head(pred_2020_probs)
pred_2020_probs$pred_sex=ifelse(pred_2020_probs[,1]>0.5,'female','male')
pred_2020_probs$pred_sex=factor(pred_2020_probs$pred_sex)
ggplot(pred_2020_probs,aes(x=stage,y=larger_prob,colour=pred_sex))+
  geom_jitter(aes(colour=pred_sex),size=0.5,
         position=position_jitterdodge(dodge.width = 1))+
  geom_violin(fill=NA)+
  theme_classic()
table(pred_2020_probs$pred_sex,pred_2020_probs$stage)

