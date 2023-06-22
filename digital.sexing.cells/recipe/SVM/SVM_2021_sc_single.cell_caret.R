
library(ggplot2)
library(gridExtra)
library(Seurat)

dat.both=readRDS('~/Documents/to_be_removed/sc_sex.differences/embryo_sex.cells/integrated.sexed.samples_seurat.obj.rds')
grep("Sxl|msl-2",rownames(dat.both))
gene.names=rownames(dat.both)

########################################################
## read in gene chro info
df.gene.table=data.table::fread('~/Documents/Data_fly_FCA/embryo_germline/gene.meta_embryo.txt',header=T,sep='\t')
head(df.gene.table)
table(df.gene.table$LOCATION_ARM) #15 Y-chromosome genes

df.gene.table=df.gene.table[df.gene.table$LOCATION_ARM %in% c('2R','3R','3L','4','2L','X','Y'),]
table(df.gene.table$LOCATION_ARM)

df.gene.table=as.data.frame(df.gene.table)
rownames(df.gene.table)=df.gene.table$current_symbol
sum(gene.names %in% df.gene.table$submitted_item) #11740
length(gene.names) #12185

rownames(df.gene.table)=df.gene.table$submitted_item
overlap.genes=gene.names[gene.names %in% df.gene.table$submitted_item]
dat.both=dat.both[overlap.genes,]
dat.both


###########################################################
expr.mat=dat.both@assays$RNA@counts

# gene filter
## calculate two set of gene: HVG, 
## (binary or continuous) mutual information about sex label
## filter gene: expr in 5% cell
if(F){
dat.both<-NormalizeData(dat.both)
dat.both <- FindVariableFeatures(dat.both, selection.method = "vst", 
                          nfeatures = 2000)
features=VariableFeatures(dat.both)
expr.mat=expr.mat[features,]
}
if(T){
  i=Matrix::rowSums(expr.mat>0)
  expr.mat=expr.mat[i>0.05*ncol(expr.mat),]
}

#change symbol to FBgn
dim(expr.mat) #6369 17142, 2000 17142 if HVG
#rownames(expr.mat)<-df.gene.table[overlap.genes,]$FBID_KEY
rownames(expr.mat)<-df.gene.table[rownames(expr.mat),]$FBID_KEY
expr.mat[1:3,1:3]
sex=dat.both$sex

## create data frame for SVM
data = data.frame(t(as.matrix(expr.mat)), y =sex)
head(data)
tail(colnames(data)) #last column is sex label

# create training and valida
library(e1071)
library(caTools)
gene.names=rownames(expr.mat)
length(unique(gene.names))# 6073, or 2000 if HVG
dim(data) # 17142  6074 (gene+sex.label.column)

dataset=data
dim(dataset) #17142  6074
colnames(dataset)


start.time=Sys.time();
set.seed(321) 
split = sample.split(dataset$y, SplitRatio = 0.8)
training_set = subset(dataset, split == TRUE)
test_set = subset(dataset, split == FALSE)
(label.column=ncol(dataset)) #5nd column is the cluster label
dim(training_set) #13714  6074
dim(test_set) #3428 6074
table(training_set$y)
table(test_set$y)

# Feature Scaling for train and test separately
training_set[-label.column] = scale(training_set[-label.column])
test_set[-label.column] = scale(test_set[-label.column])
training_set$y=factor(training_set$y)

# use caret to perform SVM https://rpubs.com/uky994/593668
library(caret)
library(doMC)
registerDoMC()

#By default caret builds the SVM linear classifier using C = 1.
# Set up Repeated k-fold Cross Validation
train_control <- trainControl(method="repeatedcv", number=5, repeats=1)
classifier = train(#formula = y ~ ., #not so sure which one, form or formula
                   form = y ~ ., #
                   data = dataset,
                   method = "svmLinear",
                   trControl = train_control,
                   preProcess = c("center","scale"),
                   probability = TRUE);
                   #tuneGrid = expand.grid(C = seq(0, 2, length = 20))); #compute SVM for a grid values of C and choose automatically the final model for predictions:

                   
  y_pred = predict(classifier, newdata = test_set[-label.column])
  #y_pred = predict(classifier, newdata = test_set[-label.column],prob=T)
  table(test_set[,label.column],y_pred)
  #classifier
  #saveRDS(classifier,'svm_classifier_2019_train34samples.rds')
  #saveRDS(classifier,'svm_2021sc_HVG2000_train80samples.rds')
  saveRDS(classifier,'svm_2021sc_train80samples.rds')
  end.time=Sys.time();
  print(end.time-start.time)
  #Time difference of 1.4 hours, 80%
  #Time difference of 2.1 hours, 70%
}

classifier=readRDS('svm_2021sc_train80samples.rds')
#classifier=readRDS('svm_2021sc_train70samples.rds')
#classifier=readRDS('svm_2021sc_HVG2000_train80samples.rds')
# 1hr for 12k cell in 6k genes
#https://stackoverflow.com/questions/34781495/how-to-find-important-factors-in-support-vector-machine
#https://stats.stackexchange.com/questions/39243/how-does-one-interpret-svm-feature-weights
#classifier
#Number of Support Vectors:  2349 sample vetors
dim(classifier$coefs) #2349 sample vectors x 1 (n.class-1) 
dim(classifier$SV) #2349 x 6073 original gene features
w<-t(classifier$coefs) %*% classifier$SV # weight vectors
dim(w) #1 x 6073 original gene features
w <- apply(w, 2, function(v){sqrt(sum(v^2))})  # weight
w <- sort(w, decreasing = T)
length(w) #8934
head(w,30) #top 30 contain, Sxl, msl-2, roX1, and roX2


library(org.Dm.eg.db)
dim(dataset) #49 x 5
#top.genes=names(w)[1:30]
top.genes=names(w)[1:100]
x=AnnotationDbi::select(org.Dm.eg.db,keys=top.genes,
                        keytype='FLYBASE',columns=c('SYMBOL'))
x$SYMBOL  
sort(x$SYMBOL)
table(df.gene.table[x$SYMBOL,]$LOCATION_ARM)
x$SYMBOL[grep('Sxl|roX|msl',ignore.case = T,x$SYMBOL)] #in top30

## prediction also takes time
dim(test_set) # 3428 6074, cell by gene
start.time=Sys.time();
preds = predict(classifier, test_set[,-label.column],probability = T)
end.time=Sys.time();
print(end.time-start.time) #Time difference of 1.116371 mins


table(preds) 
table(test_set[,label.column],preds)
#              preds
#              embryoFemale embryoMale
#embryoFemale         2007         70
#embryoMale             83       1268

preds_probs=attr(preds,'probabilities')
head(preds_probs)
table(ifelse(preds_probs[,1]>0.5,'female','male'))

preds_probs=as.data.frame(as.matrix(preds_probs))
preds_probs$larger_prob=apply(preds_probs[,c(1,2)],1,max)
head(preds_probs)
preds_probs$pred_sex=ifelse(preds_probs[,1]>0.5,'embryoFemale','embryoMale')
preds_probs$original_sex= test_set[,label.column]

min(preds_probs$larger_prob)
tmp=subset(preds_probs,larger_prob>0.6)
table(tmp$pred_sex,tmp$original_sex)

tmp=preds_probs[preds_probs$pred_sex!=preds_probs$original_sex,]

par(mfrow=c(1,2))
hist(tmp$larger_prob)
hist(preds_probs$larger_prob) #not very satisfying

par(mfrow=c(2,1))
hist(tmp$larger_prob,main='larger prob from wrong labelled cases',xlab='')
hist(preds_probs$larger_prob,main='larger prob from all cases',xlab='')

##############################################################
library(caret)

y <- as.factor(test_set[,label.column]) # factor of positive / negative cases
predictions <- as.factor(preds_probs$pred_sex) # factor of predictions
# comfusion matrix
cm <- confusionMatrix(predictions, reference =y)
cm$byClass

(precision <- posPredValue(predictions, y))
(recall <- sensitivity(predictions, y))
(F1 <- (2 * precision * recall) / (precision + recall))

#https://stackoverflow.com/questions/8499361/easy-way-of-counting-precision-recall-and-f1-score-in-r
library (ROCR);
dim(test_set); #3428 cell by 6074 gene
y=test_set[,label.column]=='embryoFemale'# logical array of positive / negative cases
predictions <- preds_probs$embryoFemale # numeric,array of predictions， 1nd <=> embryoFemale

pred <- prediction(predictions, y);

# Recall-Precision curve             
RP.perf <- performance(pred, "prec", "rec");
plot (RP.perf);

# ROC curve
ROC.perf <- performance(pred, "tpr", "fpr");
plot (ROC.perf);

# ROC area under the curve
auc.tmp <- performance(pred,"auc");
(auc <- as.numeric(auc.tmp@y.values))

# F1-score performance(pred,"f") gives a vector of F1-scores 
performance(pred,"f")


