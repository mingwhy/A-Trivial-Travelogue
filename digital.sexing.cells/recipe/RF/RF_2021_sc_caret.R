
library(ggplot2)
library(Seurat)
library(gridExtra)
library(tidyverse)
library(randomForest) 
library(caret)        # an aggregator package for performing many machine learning models


########################################################
## read in gene chro info
df.gene.table=data.table::fread('../../single.cell_datasets/embryo_germline/gene.meta_embryo.txt',header=T,sep='\t')
head(df.gene.table)
df.gene.table=df.gene.table[df.gene.table$LOCATION_ARM %in% c('2R','3R','3L','4','2L','X'),]
table(df.gene.table$LOCATION_ARM)

df.gene.table=as.data.frame(df.gene.table)
#rownames(df.gene.table)=df.gene.table$current_symbol
rownames(df.gene.table)=df.gene.table$FBID_KEY

########################################################################################
## read in gene expression data and sample meta information 

dat.both=readRDS('../../single.cell_sex.differences/embryo_sex.cells/integrated.sexed.samples_seurat.obj.rds')
grep("Sxl|msl-2",rownames(dat.both))

dim(dat.both) #12185 x 17142
sum(rownames(dat.both) %in% df.gene.table$submitted_item)#11740
dat.both=dat.both[rownames(dat.both) %in% df.gene.table$submitted_item,]

i=match(rownames(dat.both),df.gene.table$submitted_item)
sum(df.gene.table[i,]$submitted_item==rownames(dat.both)) #11740

expr.mat=as.matrix(dat.both@assays$RNA@counts)
rownames(expr.mat)=df.gene.table[i,]$validated_id #for ranger to tolerate feature names

sample.meta=dat.both@meta.data

dim(expr.mat) # 11740 gene X 17142cell
dim(sample.meta) #17142cell  4
head(sample.meta)

## filter out low expressed genes: express in 5% cells
i=Matrix::rowSums(expr.mat>0)
input.dat=expr.mat[i>ncol(expr.mat)*0.05,]
features=rownames(input.dat)
dim(input.dat) #6073 17142

input.dat=as.data.frame(t(input.dat))
input.dat[1:3,1:3]
input.dat$sex=sample.meta$sex
input.dat$sex=factor(input.dat$sex)

floor(length(features) / 3)

##############################################
# create training and validation data 
set.seed(321)
dataset=input.dat
split = caTools::sample.split(dataset$sex, SplitRatio = 0.8)
ames_train = subset(dataset, split == TRUE)
ames_test = subset(dataset, split == FALSE)
dim(ames_train) #13714  6074
dim(ames_test) #3428 6074
table(ames_train$sex)
table(ames_test$sex)

########################################################
# tune RF model (https://rpubs.com/phamdinhkhanh/389752)

#5 folds repeat 1 times
control <- trainControl(method='repeatedcv', 
                        number=5, 
                        repeats=1,search='grid')

#Metric compare model is Accuracy
metric <- "Accuracy"
set.seed(123)
tail(colnames(ames_train)) #xxx, 'sex'

# tune by grid search
if(T){
  #create tunegrid with 15 values from 1:15 for mtry to tunning model. Our train function will change number of entry variable at each split according to tunegrid. 
  ncol(ames_train)/3 #mtry~nfeature/3
  tunegrid <- expand.grid(.mtry=c(200,500,1000,1500,2000,2500))
  system.time(
    rf_gridsearch <- train(sex~., data=ames_train, method="rf", metric=metric, tuneGrid=tunegrid, trControl=control)
  )
  print(rf_gridsearch)
  plot(rf_gridsearch)
}

if(F){
  # tune by tools 
  #https://stackoverflow.com/questions/61151082/tunerf-vs-caret-tunning-for-random-forest
  set.seed(123)
  bestMtry <- tuneRF(ames_train[,-ncol(ames_train)],y=ames_train$sex, stepFactor = 1.5, improve = 1e-5, ntree = 500)
}


#Number randomely variable selected is mtry
mtry <- sqrt(ncol(x))
tunegrid <- expand.grid(.mtry=mtry)
rf_default <- train(Class~., 
                    data=dataset, 
                    method='rf', 
                    metric='Accuracy', 
                    tuneGrid=tunegrid, 
                    trControl=control)
print(rf_default)


