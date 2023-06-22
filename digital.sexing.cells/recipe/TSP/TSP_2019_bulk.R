
library(ggplot2)
library(gridExtra)
#################################################
## read in meta information
df.meta=data.table::fread('../external_data/2019_paper_reproduce.result/sample.meta_sex.label.txt');# from 02_calTPM_PCA.R 
head(df.meta)

###################################################
## read in STAR output tab-limited read count data
files=Sys.glob('../external_data/2019_paper_reproduce.result/GSE127176_RAW/*tab')
files
raw.expr=list();
for(file in files){
  gsm=strsplit(basename(file),'_')[[1]][1]
  x=read.table(file,skip=4)
  raw.expr[[gsm]]=x
}
length(raw.expr) #54

# extract 'expressed genes', defined by >=21 counts in >=5 embryos
per.embryo.gene=lapply(raw.expr,function(x){
  x[x[,2]>=21,1]
})
sapply(per.embryo.gene,length)
expr.genes=names(which(table(unlist(per.embryo.gene))>=5))
length(expr.genes) #8983 genes
filtered.expr=lapply(raw.expr,function(x){
  x[x[,1] %in% expr.genes,]
})
sapply(filtered.expr,nrow) #8983 genes

#################################################
## TPM: read.count/gene.length
## then divided/ sum(length.corrented.count) x 10^6
df.info=data.table::fread('../external_data/2019_paper_reproduce.result/dmel_geneLength_chr.txt')
#query.genes=scan('all.gene.names.txt',what='')
query.genes=expr.genes
head(query.genes)
length(query.genes) #8983
sum(df.info$GENEID %in% query.genes) #8934
tmp=merge(filtered.expr[[1]],df.info,by.x='V1',by.y='GENEID')
length.corrected.counts=lapply(filtered.expr,function(x){
  tmp=merge(x,df.info,by.x='V1',by.y='GENEID')
  tmp$V2/tmp$exonic.gene.sizes/1000 #in kilobase
})
sapply(length.corrected.counts,length)
TPM=lapply(length.corrected.counts,function(x){
  x/sum(x)*10^6
})
log2.TPM=lapply(TPM,function(i) log(i+1,base=2))

## remove 5 PB samples
library(multiclassPairs)
log2TPM.df=Reduce(`cbind`,log2.TPM)
dim(log2TPM.df) #8934 x 54
rownames(log2TPM.df)=tmp[,1]
colnames(log2TPM.df)=names(log2.TPM)
log2TPM.df[1:3,1:3]

sum(df.meta$GSM.id==colnames(log2TPM.df)) #54
expr.mat=log2TPM.df[,df.meta$cluster!='PB']
sample.meta=df.meta[df.meta$cluster!='PB',]
sum(sample.meta$GSM.id==colnames(expr.mat)) #49

## change FBgn to gene symbol
gene.id=data.table::fread('../external_data/2019_paper_reproduce.result/validate.id_2019_paper_data.txt')
fbgn=rownames(expr.mat)
x=gene.id[match(fbgn,gene.id$FLYBASE),]
sum(x$FLYBASE==fbgn) #8934
dim(expr.mat)
rownames(expr.mat)=x$SYMBOL
expr.mat[1:3,1:3]

# due to gene name with '-' are not allowed in `multiclassPairs`
rownames(expr.mat)=gsub('-','_',rownames(expr.mat))

##########################################################################
## https://github.com/NourMarzouka/multiclassPairs/blob/master/README.md
library(multiclassPairs)

# split the data
# 60% as training data and 40% as testing data
n <- ncol(expr.mat)
set.seed(1234)
training_samples <- sample(1:n,size = n*0.6)
table(sample.meta[training_samples,]$cluster)
table(sample.meta[-training_samples,]$cluster)

train <- expr.mat[,training_samples]
test  <- expr.mat[,-training_samples]
train.meta <- sample.meta[training_samples,]
test.meta <- sample.meta[-training_samples,]

# create data object 
object <- ReadData(#Data = expr.mat,
                   Data = train,
                   Labels = train.meta$cluster,
                   verbose = FALSE)
object

filtered_genes <- filter_genes_TSP(data_object = object,
                                   filter = "one_vs_one",
                                   platform_wise = FALSE,
                                   featureNo = 1000,
                                   UpDown = TRUE,
                                   verbose = TRUE)
filtered_genes

# Let's train our model
classifier <- train_one_vs_rest_TSP(data_object = object,
                                    filtered_genes = filtered_genes,
                                    k_range = 2:500,
                                    include_pivot = FALSE,
                                    one_vs_one_scores = TRUE,
                                    platform_wise_scores = FALSE,
                                    seed = 1234,
                                    verbose = FALSE)
classifier
classifier$classifiers$female
classifier$classifiers$male

# apply on the training data
# To have the classes in output in specific order, we can use classes argument
results_train <- predict_one_vs_rest_TSP(classifier = classifier,
                                         Data = object,
                                         tolerate_missed_genes = TRUE,
                                         weighted_votes = TRUE,
                                         classes = c("female",'male'),
                                         verbose = TRUE)
# apply on the testing data
results_test <- predict_one_vs_rest_TSP(classifier = classifier,
                                        Data = test,
                                        tolerate_missed_genes = TRUE,
                                        weighted_votes = TRUE,
                                        classes = c("female",'male'),
                                        verbose = TRUE)
# get a look over the scores in the testing data
knitr::kable(head(results_train))
table(results_train$max_score)

knitr::kable(head(results_test))
table(results_test$max_score)

# Confusion Matrix and Statistics on training data
caret::confusionMatrix(data = factor(results_train$max_score, 
                                     levels = unique(object$data$Labels)),
                       reference = factor(object$data$Labels, 
                                          levels = unique(object$data$Labels)),
                       mode="everything")

# Confusion Matrix and Statistics on testing data
caret::confusionMatrix(data = factor(results_test$max_score, 
                                     levels = unique(object$data$Labels)),
                       reference = factor(test.meta$cluster,
                                          levels = unique(object$data$Labels)),
                       mode="everything")

# Visualiation
# plot for the rules and scores in the training data
pdf("TSP_out.pdf")
plot_binary_TSP(Data = object, # we are using the data object here
                classifier = classifier, 
                prediction = results_train, 
                classes =  c("female",'male'),
                margin = c(0,5,0,10),
                title = "Training data")

# plot for the rules and scores in the testing data
plot_binary_TSP(Data = test, # ExpressionSet
                ref = test.meta$cluster,
                classifier = classifier, 
                prediction = results_test, 
                classes =  c("female",'male'),
                title = "Testing data",
                margin = c(0,5,0,10))
dev.off()


