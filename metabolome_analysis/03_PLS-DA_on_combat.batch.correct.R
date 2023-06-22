
#BiocManager::install("mixOmics")
library(mixOmics)

# read in phenotype data
dat=readRDS('./mz.filter.dat.rds')
df.mat=dat$edata;
pheno=dat$pheno

# read in batch corrected data after `combat` function correction
df.res=readRDS('mz.filter.combat.rds')
dim(df.res) #170 mz x 60 sample

# X: nsample x nfeature
X=t(df.res)
# Y: class
#insulin sensitive: wt, wt/353
#insulin resistant: 19/74, 74/211, 211/19
pheno$insulin=rep('resistant',nrow(pheno))
pheno[pheno$genotype %in% c('wt','wt/353'),]$insulin='sensitive'

#'normal' fecundity: wt, 74/211, 211/19
#high fecundity: wt/353
#low fecundity: 19/74 
pheno$fecundity=rep('normal',nrow(pheno))
pheno[pheno$genotype=='wt/353',]$fecundity='high'
pheno[pheno$genotype=='19/74',]$fecundity='low'

table(pheno$genotype,pheno$insulin)
table(pheno$genotype,pheno$fecundity)

#############################################################
## tune the model for insulin, two groups
Y=pheno$insulin
MyResult.plsda2 <- plsda(X,Y, ncomp=10)

# select ncomp
set.seed(30) # for reproducbility in this vignette, otherwise increase nrepeat
MyPerf.plsda <- perf(MyResult.plsda2,validation="Mfold", folds=4, 
                     progressBar=FALSE, nrepeat=10) # we suggest nrepeat = 50
plot(MyPerf.plsda, col = color.mixo(5:7), sd = TRUE, legend.position = "horizontal")
MyPerf.plsda

# for eacn comp, select n.feature
list.keepX <- c(5:10,  seq(20, 100, 10))
list.keepX # to output the grid of values tested
set.seed(30) # for reproducbility in this vignette, otherwise increase nrepeat
tune.splsda.srbct <- tune.splsda(X, Y, ncomp = 3, # we suggest to push ncomp a bit more, e.g. 4
                                 validation = 'Mfold',
                                 folds = 10, dist = 'max.dist', progressBar = FALSE,
                                 measure = "BER", test.keepX = list.keepX,
                                 nrepeat = 10)   # we suggest nrepeat = 50
error <- tune.splsda.srbct$error.rate
ncomp <- tune.splsda.srbct$choice.ncomp$ncomp # optimal number of components based on t-tests on the error rate
ncomp
# 1
select.keepX <- tune.splsda.srbct$choice.keepX[1:ncomp]  # optimal number of variables to select
select.keepX
#comp1 
#5
plot(tune.splsda.srbct, col = color.jet(3))

#Based on those tuning results, we can run our final and tuned sPLS-DA model:
MyResult.splsda.final <- splsda(X, Y, ncomp = 2, keepX = select.keepX)
#MyResult.splsda.final <- splsda(X, Y, ncomp = 1, keepX = 5)
selectVar(MyResult.splsda.final,comp=1)$value
#selectVar(MyResult.splsda.final,comp=2)$value
plotIndiv(MyResult.splsda.final, ind.names = FALSE, legend=TRUE,
          ellipse = TRUE, title="sPLS-DA - final result")

#############################################################
## tune the model for fecundity
Y=pheno$fecundity
MyResult.plsda2 <- plsda(X,Y, ncomp=10)

# select ncomp
set.seed(30) # for reproducbility in this vignette, otherwise increase nrepeat
MyPerf.plsda <- perf(MyResult.plsda2, validation = "Mfold", folds = 4, 
                     progressBar = FALSE, nrepeat = 50) # we suggest nrepeat = 50

plot(MyPerf.plsda, col = color.mixo(5:7), sd = TRUE, legend.position = "horizontal")
MyPerf.plsda

# for eacn comp, select n.feature
list.keepX <- c(5:10,  seq(20, 100, 10))
list.keepX # to output the grid of values tested
set.seed(30) # for reproducbility in this vignette, otherwise increase nrepeat
tune.splsda.srbct <- tune.splsda(X, Y, ncomp = 4, # we suggest to push ncomp a bit more, e.g. 4
                                 validation = 'Mfold',
                                 folds = 10, dist = 'max.dist', progressBar = FALSE,
                                 measure = "BER", test.keepX = list.keepX,
                                 nrepeat = 10)   # we suggest nrepeat = 50
error <- tune.splsda.srbct$error.rate
ncomp <- tune.splsda.srbct$choice.ncomp$ncomp # optimal number of components based on t-tests on the error rate
ncomp

select.keepX <- tune.splsda.srbct$choice.keepX[1:ncomp]  # optimal number of variables to select
select.keepX
#comp1 comp2 comp3 
#10    60     7 
plot(tune.splsda.srbct, col = color.jet(4))

#Based on those tuning results, we can run our final and tuned sPLS-DA model:
MyResult.splsda.final <- splsda(X, Y, ncomp = 3, keepX = select.keepX)
selectVar(MyResult.splsda.final,comp=1)$value
selectVar(MyResult.splsda.final,comp=2)$value
selectVar(MyResult.splsda.final,comp=3)$value
plotIndiv(MyResult.splsda.final, ind.names = FALSE, legend=TRUE,
          ellipse = TRUE, title="sPLS-DA - final result")
plotLoadings(MyResult.splsda.final,contrib='max',method='mean')

