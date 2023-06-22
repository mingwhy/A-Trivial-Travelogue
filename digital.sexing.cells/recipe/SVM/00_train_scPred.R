
library("scPred") #https://github.com/powellgenomicslab/scPred
library("Seurat")
library("magrittr")


library(ggplot2)
library(gridExtra)
library(Seurat)

dat.both=readRDS('~/Documents/Data_fly_FCA/embryo_germline/embryo_sex.cells/integrated.sexed.samples_seurat.obj.rds')
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

dim(dat.both) # 11740 17142


########################################################
## 80/20 split to train a model
if(!file.exists('scPred_train80.rds')){
  #reference <- dat.both
  set.seed(2049)
  i=sample(1:ncol(dat.both),ceiling(ncol(dat.both)*0.2),replace = F)
  reference=dat.both[,-i]
  query=dat.both[,i]
  table(reference$sex) #seurat obj
  table(query$sex) ##seurat obj
  
  # preprocess
  system.time( 
    reference <- reference %>% 
      NormalizeData() %>% 
      FindVariableFeatures(nfeatures=2000) %>% 
      ScaleData() %>% 
      RunPCA(npcs=100) %>% 
      RunUMAP(dims = 1:100)
  ) #2min
  table(reference$sex)
  #embryoFemale   embryoMale 
  #7954         5188 
  
  # Firstly, start training: Training classifiers with scPred
  ##getFeatureSpace will create a `scPred` object stored in the @misc slot. 
  system.time( reference <- getFeatureSpace(reference, "sex") ) #4s
  names(reference@misc)
  reference@misc$scPred
  #|Cell type    |    n| Features|
  #|:------------|----:|--------:|
  #|embryoFemale | 7954|      100|
  
  ## Secondly, we train the classifiers for each cell using the trainModel function. 
  ## By default, scPred will use a support vector machine with a radial kernel.
  system.time( reference<-trainModel(reference, model = "svmRadial",
                                     preProcess = c("center", "scale"),
                                     resampleMethod = "cv",
                                     number = 5) ) #8 mins
  
  tmp=reference@misc$scPred@train
  tmp #a caret output obj
  class(tmp) #list
  names(tmp)
  names(tmp$embryoFemale)
  
  head(tmp$embryoFemale$pred)
  table(tmp$embryoFemale$pred$Resample) #5-fold cv
  sum(table(reference$sex)/5) #3429 cells for test
  
  tmp$embryoFemale$ptype #PC1-100
  
  #https://stackoverflow.com/questions/17529537/example-for-svm-feature-selection-in-r
  
  
  get_probabilities(reference) %>% head()
  
  get_scpred(reference)
  
  # important (PC) features in SVM and genes with high loadings in these PCs
  # sample by PC, responses
  spmodel=get_scpred(reference)
  dim(spmodel@features$embryoFemale)  #100 PC x 3 (pc.name, pValue, pValueAdj)
  dim(spmodel@feature_loadings) #5000 HGV genes x 100 pcs
  
  dim(spmodel@cell_embeddings) #17142 cells x 100 PC embedding
  x = spmodel@cell_embeddings
  y = factor(reference$sex)
  dim(x); 
  length(y);
  
  if(F){
    system.time(
      svmProfile <- rfe(x, y,
                      sizes = c(2, 5, 10, 20, 50),
                      rfeControl = rfeControl(functions = caretFuncs,
                                              method='cv',
                                              number = 5),
                      ## pass options to train()
                      method = "svmRadial") #3.5hr
    )
    saveRDS(svmProfile,file='svmProfile_cv5.rds')
  }
  #scPred is built on top of the caret package and allows using a large set of prediction models (e.g. logistic regression, decision trees, bagging, neural networks, etc).
  #PC_1, PC_26, PC_6, PC_2, PC_14
  dim(spmodel@feature_loadings) #5000 HGV genes x 100 pcs
  gene.pc=spmodel@feature_loadings
  head(sort(gene.pc[,1],decreasing = T)) #PC1
  pc1=sort(gene.pc[,1],decreasing = T)
  grep('Sxl|roX|msl',names(pc1)) #1  43 693
  
  saveRDS(reference,'scPred_train80.rds')
}

# start predction: random sample 20% cells for prediction
reference=readRDS('scPred_train80.rds');
#query=dat.both[,sample(1:ncol(dat.both),4000,replace = F)]
query<-NormalizeData(query)
system.time(query <- scPredict(query, reference)) #12sec
table(query$sex,query$scpred_prediction)

if(F){
#https://github.com/mingwhy/bioinfo_homemade_tools/blob/main/package_101/scPred_101.R
#scPred now uses Harmony to align the query data onto the training low-dimensional space used as reference. Once the data is aligned, cells are classified using the pre-trained models.

DimPlot(query, group.by = "scpred_prediction", reduction = "scpred")
query <- RunUMAP(query, reduction = "scpred", dims = 1:30)
DimPlot(query, group.by = "scpred_prediction", label = TRUE, repel = TRUE)


FeaturePlot(query, c("scpred_B.cell", "scpred_CD4.T.cell", "scpred_CD8.T.cell", 
                     "scpred_cMono", "scpred_ncMono", "scpred_Plasma.cell", 
                     "scpred_cDC", "scpred_pDC"))
crossTab(query, "cell_type", "scpred_prediction")


# advanced options
## get classifier
get_classifiers(reference)
caret::plot.train(get_classifiers(reference)[["NK cell"]])
## use a different probability threshold
query <- scPredict(query, reference, recompute_alignment = FALSE, threshold = 0.9)

## Parallel training (See Caret parallel processing for more details)
library(doParallel)
cl <- makePSOCKcluster(2)
registerDoParallel(cl)
reference <- trainModel(reference, model = "mda", allowParallel = TRUE)

## extract the classifier from the Seurat obj
scpred <- get_scpred(reference)
query <- scPredict(query, scpred)
}

