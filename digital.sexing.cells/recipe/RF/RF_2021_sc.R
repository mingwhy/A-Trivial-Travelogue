
library(ggplot2)
library(Seurat)
library(gridExtra)
library(tidyverse)
library(rsample)      # data splitting 
library(randomForest) # basic implementation
library(ranger)       # a faster implementation of randomForest
library(caret)        # an aggregator package for performing many machine learning models
library(h2o)          # an extremely fast java-based platform

########################################################
## read in gene chro info
df.gene.table=data.table::fread('../../single.cell_datasets/embryo_germline/gene.meta_embryo.txt',header=T,sep='\t')
head(df.gene.table)
df.gene.table=df.gene.table[df.gene.table$LOCATION_ARM %in% c('2R','3R','3L','4','2L','X'),]
table(df.gene.table$LOCATION_ARM)

df.gene.table=as.data.frame(df.gene.table)
rownames(df.gene.table)=df.gene.table$current_symbol
#rownames(df.gene.table)=df.gene.table$FBID_KEY
########################################################################################
## read in gene expression data and sample meta information 

dat.both=readRDS('../../single.cell_sex.differences/embryo_sex.cells/integrated.sexed.samples_seurat.obj.rds')
grep("Sxl|msl-2",rownames(dat.both))

dim(dat.both) #12185 x 17142
sum(rownames(dat.both) %in% df.gene.table$submitted_item)#11740
dat.both=dat.both[rownames(dat.both) %in% df.gene.table$submitted_item,]

i=match(rownames(dat.both),df.gene.table$submitted_item)
sum(df.gene.table[i,]$submitted_item==rownames(dat.both)) #11740


###########################################################
dat.both=NormalizeData(dat.both)
#expr.mat=as.matrix(dat.both@assays$RNA@counts)
expr.mat=as.matrix(dat.both@assays$RNA@data) #use normalized data

# gene filter
## calculate two set of gene: HVG, 
## (binary or continuous) mutual information about sex label
## filter gene: expr in 5% cell
if(F){
  dat.both <- FindVariableFeatures(dat.both, selection.method = "vst", 
                                   nfeatures = 2000)
  features=VariableFeatures(dat.both)
  expr.mat=expr.mat[features,]
  rownames(expr.mat)=df.gene.table[features,]$validated_id 
}
if(T){
  i=Matrix::rowSums(expr.mat>0)
  expr.mat=expr.mat[i>0.05*ncol(expr.mat),]
  rownames(expr.mat)=df.gene.table[rownames(expr.mat),]$validated_id #for ranger to tolerate feature names
}

head(rownames(expr.mat)); #use FBgn 
sample.meta=dat.both@meta.data
  
dim(expr.mat) # 11740 gene X 17142cell
dim(sample.meta) #17142cell  4
head(sample.meta)

## filter out low expressed genes: express in 5% cells
input.dat=expr.mat;
features=rownames(input.dat)
dim(input.dat) #6073 17142

input.dat=as.data.frame(t(input.dat))
input.dat[1:3,1:3]
input.dat$sex=sample.meta$sex
input.dat$sex=factor(input.dat$sex)

##############################################
## use ranger on all cells
if(F){
system.time(
  rf_ranger<- ranger(
                  formula=sex ~ ., 
                  data = input.dat, 
                  num.trees = 500,
                  importance = 'impurity',
                  mtry = floor(length(features) / 3)
            )
)
#   user   system  elapsed 
# 7578.153   16.446  402.549 
rf_ranger
model=rf_ranger
rf_ranger$prediction.error

# feature importance in random forest
#https://www.statistik.uni-dortmund.de/useR-2008/slides/Strobl+Zeileis.pdf
# Gini importance
model$importance
importance(model,type=2)
top.pick=names( importance(model,type=2)[1:20])
rownames(df.gene.table)=df.gene.table$FBID_KEY
df.gene.table[top.pick,]

# permutation importance
#mean decrease in classification accuracy after permuting Xj over all trees
importance(model,type=1)
top.pick=names( importance(model,type=1)[1:20])
rownames(df.gene.table)=df.gene.table$FBID_KEY
df.gene.table[top.pick,]
}
##############################################
# training and testing, prediction confidence

# create training and validation data 
set.seed(321)
#valid_split <- initial_split(input.dat, .8)
#ames_train <- training(valid_split)
#ames_test  <- testing(valid_split)
dataset=input.dat
split = caTools::sample.split(dataset$sex, SplitRatio = 0.8)
ames_train = subset(dataset, split == TRUE)
ames_test = subset(dataset, split == FALSE)
dim(ames_train) #13714  6074 (gene + sex label)
dim(ames_test) #3428 6074 (gene + sex label)
table(ames_train$sex) #8309         5405 
table(ames_test$sex) # 2077         1351 

if(T){
  system.time(
    rf_ranger<- ranger(
      formula=sex ~ ., 
      data = ames_train, 
      num.trees = 500, #1000 tree, 1.5hr, 500 tree, 41min
      importance = 'impurity',
      mtry = floor(length(features) / 3),
      probability = TRUE
    )
  )
  #saveRDS(rf_ranger,file='rf_ranger_500tree_80train.rds')
  #saveRDS(rf_ranger,file='rf_ranger_1000tree_80train.rds')
  saveRDS(rf_ranger,file='rf_ranger_normData_500tree_80train.rds')
  #saveRDS(rf_ranger,file='rf_ranger_normData_1000tree_80train.rds')
}

rf_ranger=readRDS('rf_ranger_normData_500tree_80train.rds')
#rf_ranger=readRDS('rf_ranger_normData_1000tree_80train.rds')
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
  dplyr::top_n(100)
df.gene.table1=df.gene.table;
rownames(df.gene.table1)=df.gene.table1$FBID_KEY
tmp1=df.gene.table1[tmp$names,]
head(tmp1)
grep('Sxl|roX|msl',ignore.case = T,tmp1$SYMBOL) #34 79 91,"lncRNA:roX1" "lncRNA:roX2" "Sxl" 
tmp1$SYMBOL[grep('Sxl|roX|msl',ignore.case = T,tmp1$SYMBOL)] #in top100


dim(ames_test) #3428 6074, cell by gene
#pred_ranger = predict(rf_ranger, ames_test)
pred_ranger = predict(rf_ranger, ames_test[,-ncol(ames_test)])
head(pred_ranger$predictions) #probabilities
predict.labels=ifelse(pred_ranger$predictions[,1]>0.5,'embryoFemale','embryoMale')
table(ames_test$sex,predict.labels)
#              predict.labels
#           embryoFemale embryoMale
#embryoFemale         2035         42
#embryoMale            142       1209

tmp=as.data.frame(pred_ranger$predictions)
tmp$pred=predict.labels;
tmp$real=ames_test$sex
tmp$larger.prob=apply(tmp[,c(1,2)],1,max)
head(tmp)
sub.tmp=tmp[tmp$pred!=tmp$real,]

par(mfrow=c(2,1))
hist(sub.tmp$larger.prob,main='larger prob from wrong labelled cases',xlab='')
hist(tmp$larger.prob,main='larger prob from all cases',xlab='')

########################################################
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
y=ames_test$sex=='embryoFemale'# logical array of positive / negative cases
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
