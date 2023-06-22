
library(ggplot2)
library(gridExtra)
library(infotheo) #mutinformation() mutual info between two discrete vars
library(mpmi) #mutual infor between one continuous and one discrete
library(e1071)
library(caTools)
library(org.Dm.eg.db)
########################################################################################
## read in gene meta info
df.gene.info=data.table::fread('../external_data/2019_paper_reproduce.result/dmel_geneLength_chr.txt')
head(df.gene.info)

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

dim(expr.mat) #gene by sample 8934   49
dim(sample.meta) #sample by attributes

# gene feature selection
## expr in >5% cells
## mutual information
i=Matrix::rowSums(expr.mat>0)
sum(i<=0.05*ncol(expr.mat)) #0 gene

## read in 2011, 2015 data to subset genes
dat.2011=data.table::fread('../external_data/2011_paper_data/2011_log2.RPKM.txt')
dat.2011$gene.name #12353    24

dat.2015=data.table::fread('../external_data/2015_paper_data/2015_FPKM_normalization_6mel.txt')
dat.2015$gene.name

overlap.genes1=intersect(gene.names.2019,dat.2011$gene.name)
overlap.genes2=intersect(gene.names.2019,dat.2015$gene.name)
length(overlap.genes1) #8215
length(overlap.genes2) #5011, no Sxl in this dataset

#expr.mat.sub=expr.mat
#dim(expr.mat.sub) #gene by sample 4918   49
#dim(sample.meta) #sample by attributes

if(T){
# use `mpmi` r package to calculate mutual info
#https://www.rdocumentation.org/packages/mpmi/versions/0.43.1/topics/mmi
#mmi.pw(expr.mat[1,],sample.meta$cluster)
mi.score=mmi(t(expr.mat),as.matrix(sample.meta$cluster))
sapply(mi.score,length)
summary(mi.score$bcmi)
#plot(mi.score$bcmi,mi.score$zvalues)
}

# select top 2000 genes for SVM or use all genes
topK=1000;
expr.mat.sub=expr.mat[order(mi.score$bcmi,decreasing = T)[1:topK],]
dim(expr.mat.sub) #feature x sample
x=AnnotationDbi::select(org.Dm.eg.db,keys=rownames(expr.mat.sub),
                        keytype='FLYBASE',columns=c('SYMBOL'))
x[grep('Sxl|roX|msl',x$SYMBOL),]#all in top2000

data = data.frame(t(expr.mat.sub), y =sample.meta$cluster)
head(data)
tail(colnames(data)) #last column is sex label

# create training and valida
gene.names=colnames(data)
length(unique(gene.names))# 8935
dim(data) # 49 8935

dataset=data
dim(dataset) #49 sample x 5 feature gene
#dataset=data[,-3] #remove roX2
head(colnames(dataset))

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

if(T){
  
  set.seed(321) 
  split = sample.split(dataset$y, SplitRatio = 0.7)
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
                   cost=10,kernel = 'linear')
                   #cost=1,gamma=0.05,kernel = 'polynomial')
                   #cost=1,gamma=0.05,kernel = 'radial')
 
  classifier
  saveRDS(classifier,paste0('svm_classifier_2019_train34samples_top',topK,'.rds'))
}

classifier=readRDS(paste0('svm_classifier_2019_train34samples_top',topK,'.rds'))
y_pred = predict(classifier, newdata = test_set[-label.column])
table(test_set[,label.column],y_pred)

#https://stackoverflow.com/questions/34781495/how-to-find-important-factors-in-support-vector-machine
w<-t(classifier$coefs) %*% classifier$SV # weight vectors
w <- apply(w, 2, function(v){sqrt(sum(v^2))})  # weight
w <- sort(w, decreasing = T)
length(w) #8934
head(w)


dim(dataset) #49 x 5
top.genes=names(w)[1:25]
x=AnnotationDbi::select(org.Dm.eg.db,keys=top.genes,
                        keytype='FLYBASE',columns=c('SYMBOL'))
x$SYMBOL  
x$SYMBOL[grep('Sxl|roX|msl',ignore.case = T,x$SYMBOL)] #in top25

tmp=merge(df.gene.info,x,by.x='GENEID',by.y='FLYBASE')
table(tmp$TXCHROM) #biased to X chr

gene.names.2019=names(w)
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
#top1000 gene behaves the best
#    F10     F11     F12     U13    F14A    U14A    F14B F14B_r2    F14C F14C_r2 
#female  female  female    male    male  female  female  female  female  female 
#F14D    U14D     M10     M11     M12     M13    M14A M14A_r2    M14B M14B_r2 
#male  female  female  female  female  female  female    male    male    male 
#M14C M14C_r2 
#male    male 

##########################################
##########################################
##########################################
if(F){
  # there is a code bug in 2020 paper and this paper lack sex label. throw away this dataset
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
pred_2020 = predict(classifier, test)
table(pred_2020)

# biased to female samples
}

