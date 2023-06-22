
library(org.Dm.eg.db)

## read in bulk 2019 svm
topK=1000;
classifier=readRDS(paste0('svm_classifier_2019_train34samples_top',topK,'.rds'))

#https://stackoverflow.com/questions/34781495/how-to-find-important-factors-in-support-vector-machine
w<-t(classifier$coefs) %*% classifier$SV # weight vectors
w <- apply(w, 2, function(v){sqrt(sum(v^2))})  # weight
w <- sort(w, decreasing = T)
length(w) #8934
head(w)

top.genes=names(w)[1:25]
x=AnnotationDbi::select(org.Dm.eg.db,keys=top.genes,
                        keytype='FLYBASE',columns=c('SYMBOL'))
x$SYMBOL  
x$SYMBOL[grep('Sxl|roX|msl',ignore.case = T,x$SYMBOL)] #in top25

set1=x

## read in bulk 2019 svm
classifier=readRDS('svm_2021sc_train80samples.rds')
dim(classifier$coefs) #2349 sample vectors x 1 (n.class-1) 
dim(classifier$SV) #2349 x 6073 original gene features
w<-t(classifier$coefs) %*% classifier$SV # weight vectors
dim(w) #1 x 6073 original gene features
w <- apply(w, 2, function(v){sqrt(sum(v^2))})  # weight
w <- sort(w, decreasing = T)
length(w) #8934
head(w,30) #top 30 contain, Sxl, msl-2, roX1, and roX2

top.genes=names(w)[1:25]
x=AnnotationDbi::select(org.Dm.eg.db,keys=top.genes,
                        keytype='FLYBASE',columns=c('SYMBOL'))
x$SYMBOL  
x$SYMBOL[grep('Sxl|roX|msl',ignore.case = T,x$SYMBOL)] #in top25
set2=x

####
sum(set1$FLYBASE %in% set2$FLYBASE) #only those 4 genes overlap
#[1] "lncRNA:roX2" "lncRNA:roX1" "Sxl"         "msl-2" 
set1$SYMBOL
set2$SYMBOL

