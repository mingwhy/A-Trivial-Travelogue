
library(ggplot2)
library(Seurat)
library(gridExtra)
library(tidyverse)
library(rsample)      # data splitting 
library(randomForest) # basic implementation
library(ranger)       # a faster implementation of randomForest
library(caret)        # an aggregator package for performing many machine learning models
library(h2o)          # an extremely fast java-based platform

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
input.dat=as.data.frame(t(expr.mat))
input.dat[1:3,1:3]
input.dat$sex=sample.meta$cluster
input.dat$sex=factor(input.dat$sex)

gene.names.2019=rownames(expr.mat)
head(gene.names.2019)
# create training and validation data 
set.seed(123)
valid_split <- initial_split(input.dat, .8)
ames_train <- training(valid_split)
ames_test  <- testing(valid_split)
dim(ames_train) # 39 8935
dim(ames_test) # 10 8935
features=rownames(expr.mat)
tail(features)

system.time(
  rf_ranger<- ranger(
    formula=sex ~ ., 
    data = ames_train, 
    #num.trees = 500,
    importance = 'impurity',
    mtry = floor(length(features) / 3),
    probability = TRUE
  )
)

rf_ranger
rf_ranger$variable.importance %>% 
  tidy() %>%
  dplyr::arrange(desc(x)) %>%
  dplyr::top_n(25) %>%
  ggplot(aes(reorder(names, x), x)) +
  geom_col() +
  coord_flip() +
  ggtitle("Top 25 important variables")


tmp=rf_ranger$variable.importance %>% 
  tidy() %>%
  dplyr::arrange(desc(x)) %>%
  dplyr::top_n(30)
library(AnnotationDbi);library(org.Dm.eg.db)
tmp1=AnnotationDbi::select(org.Dm.eg.db,keys=tmp$names,keytype = 'FLYBASE',columns=c('SYMBOL') )
head(tmp1)

#pred_ranger = predict(rf_ranger, ames_test)
pred_ranger = predict(rf_ranger, ames_test[,-ncol(ames_test)])
head(pred_ranger$predictions) #probabilities
predict.labels=ifelse(pred_ranger$predictions[,1]>0.5,'female','male')
table(predict.labels,ames_test$sex)
#predict.labels female male
#embryoFemale      5    0
#embryoMale        0    5

library(caret)

y <- as.factor(ames_test$sex) # factor of positive / negative cases
predictions <- as.factor(predict.labels) # factor of predictions
# comfusion matrix
cm <- confusionMatrix(predictions, reference =y)
cm$byClass

(precision <- posPredValue(predictions, y))
(recall <- sensitivity(predictions, y))
(F1 <- (2 * precision * recall) / (precision + recall))

#https://stackoverflow.com/questions/8499361/easy-way-of-counting-precision-recall-and-f1-score-in-r
library (ROCR);
y=ames_test$sex=='female'# logical array of positive / negative cases
predictions <- pred_ranger$predictions[,1] # array of predictions， 1nd <=> embryoFemale

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

#####################################################
## apply RF model to other bulk dataset
# gene name match (FBgn), expr.mat format
dat.2015=data.table::fread('../external_data/2015_paper_data/2015_FPKM_normalization_6mel.txt')
                           
gene.names=dat.2015$gene.name
dat.2015=as.data.frame(dat.2015[,-1])
rownames(dat.2015)=gene.names
dim(dat.2015) #6003    6
dat.2015[1:3,1:3]

length(names(rf_ranger$variable.importance)) #8934
sum(gene.names %in% names(rf_ranger$variable.importance)) #5011

inp=(matrix(0,nrow=ncol(dat.2015),ncol=length(names(rf_ranger$variable.importance))))
rownames(inp)=colnames(dat.2015)
colnames(inp)=names(rf_ranger$variable.importance)
overlap.genes=gene.names[gene.names %in% names(rf_ranger$variable.importance)]
inp[,overlap.genes]<-t(dat.2015[overlap.genes,])

test=inp;
dim(test)
pred_2015 = predict(rf_ranger, test)
pred_2015$predictions #all predict to be males
