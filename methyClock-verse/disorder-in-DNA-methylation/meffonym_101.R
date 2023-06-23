

# https://github.com/perishky/meffonym/
library(meffonym)
# when loading library, you can see message in Rstudio
# loading /Library/Frameworks/R.framework/Versions/4.1/Resources/library/meffonym/hannum/coefs.csv  ...
# the already-constructed model can be accessed following this path
# this clock only contains 71 sites 
hannum.clock=data.table::fread('/Library/Frameworks/R.framework/Versions/4.1/Resources/library/meffonym/hannum/coefs.csv',header = T)
dim(hannum.clock) #71, site, coef

library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(minfi)
packageVersion("IlluminaHumanMethylation450kanno.ilmn12.hg19")
#[1] '0.6.0'
ann450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
head(ann450k)
dim(ann450k) #485512     33

# use ann450k to characterize the 71 site in hannum clock
x=ann450k[hannum.clock$pred.var,]
table(x$Regulatory_Feature_Group)

################################################################################################################################################
# Download 'GSE40279_series_matrix.txt.gz' from the Hannum et al. data from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE40279
source('read-gse-matrix-file.r')
hannum=read.gse.matrix.file('./GSE40279/GSE40279_series_matrix.txt.gz') #https://rdrr.io/github/perishky/meffonym/f/vignettes/age-tutorial.rmd
names(hannum) #~15min
#"dnam"    "samples"
#saveRDS(hannum$dnam,'hannum_methy.mat.rds')
#saveRDS(hannum$samples,'hannum_samples.rds')

#Extract phenotypes/exposures.
extract.characteristic <- function(ch) {
  sub("[^:]+: (.*)", "\\1", ch)
}

hannum$samples$age <- as.numeric(extract.characteristic(hannum$samples$characteristics_ch1))
hannum$samples$gender <- extract.characteristic(hannum$samples$characteristics_ch1.3)
hannum$samples$ethnicity <- extract.characteristic(hannum$samples$characteristics_ch1.4)

#Calculate scores for all available models
library(meffonym)
meffonym.models() #133 models, `models.csv` https://github.com/perishky/meffonym/ 

scores <- cbind(age=hannum$samples$age,
                sapply(meffonym.models(), function(model) {
                  meffonym.score(hannum$dnam, model)$score
                }))
dim(scores) #656 134
scores[1:3,1:3]

#Calculate correlations between scores.
cor(scores)

#Calculate differences between biological and estimated age.
apply(scores[,-1], 2, function(estimate) quantile(estimate-scores[,"age"]))


#Calculate age acceleration
#Age accelaration is obtained by adjusting age estimates for actual age.

age.scores <- scores[,c("horvath","hannum")]

acc <- apply(age.scores, 2, function(estimate) {
  residuals(lm(estimate ~ scores[,"age"]))
})

#Test age acceleration differences between the sexes. Both Hannum and Horvath scores show greater acceleration in males than females.
t(sapply(colnames(acc), function(name) {
  acc <- acc[,name]
  coef(summary(lm(acc ~ gender, data=hannum$samples)))["genderM",]
}))
##             Estimate Std. Error  t value     Pr(>|t|)
## age.horvath 1.546429  0.3908939 3.956134 8.449921e-05
## age.hannum  2.215048  0.3206417 6.908173 1.166350e-11

#Obtaining estimates identical to the online Horvath clock
#The following estimates are equivalent to that of the Horvath clock. The calculation can be time consuming, so we will obtain estimates for only the first 10.
horvath.est <- meffonym.horvath.age(hannum$dnam[,1:10]) ## 1-2 minutes

#These estimates are obtained in two steps: normalizing the methylation data to 'gold standard' used by Horvath when developing the model, and then applying the model in the normalized data.

standard <- meffonym.horvath.standard()
## step 1: normalization (1-2 minutes)
dnam.norm <- meffonym.bmiq.calibration(hannum$dnam[,1:10], standard) 
## step 2: model application
equiv.est <- meffonym.score(dnam.norm, "horvath")

#They are the same.
cor(equiv.est$score, horvath.est)
equiv.est$score - horvath.est

#Normalization in this case has a fairly minor effect on the age estimates.
cor(age.scores[1:10,"horvath"], horvath.est)

#Calculate scores for a new model
model <- meffonym.get.model("horvath")
idx <- order(abs(model$coefs), decreasing=T)[1:10]
model$coefs <- model$coefs[idx]
meffonym.add.model("age.horvath.10",
                   variables=names(model$coefs),
                   coefficients=model$coefs,
                   description="Horvath model with only the top 10 effect sizes.")

#Estimate age with that new model.
ten.est <- meffonym.score(hannum$dnam, "age.horvath.10")
#The result is obviously similar but 10 CpG sites is obviously not the same as 353.
cor(ten.est$score, scores[,"horvath"])
## [1] 0.7134123
