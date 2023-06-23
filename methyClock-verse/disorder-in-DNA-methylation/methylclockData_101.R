#https://github.com/isglobal-brge/methylclockData/tree/main
#https://bioconductor.org/packages/release/data/experiment/vignettes/methylclockData/inst/doc/methylcockData.html
library(ExperimentHub)
#BiocManager::install('methylclockData')
library(methylclockData)

# Get experimentHub records
eh <- ExperimentHub()

# Get data about methylclockData experimentHub
pData <- query(eh , "methylclockData")

# Get information rows about methylclockData
df <- mcols(pData)
df
# Retrieve data with Hobarth's clock coefficients
pData["EH6071"]
#  Hovarths CpGs to train a Bayesian Neural Network
cpgs.bn <- get_cpgs_bn()
head(cpgs.bn)
# Hannum's clock coefficients
coefHannum <- get_coefHannum()
head(coefHannum)
# Hobarth's clock coefficients
coefHorvath <- get_coefHorvath()
head(coefHorvath)
dim(coefHorvath) #intercept + 353 cpg site

#########
#https://github.com/yiluyucheng/dnaMethyAge
## Make sure 'devetools' is installed in your R
# install.packages("devtools")
#devtools::install_github("yiluyucheng/dnaMethyAge")
library('dnaMethyAge')

## prepare betas dataframe
data('subGSE174422') ## load example betas
head(betas) #cg site x samples
print(dim(betas)) ## probes in row and samples in column
# 485577 8

availableClock() ## List all supported clocks
# "HannumG2013"  "HorvathS2013" "LevineM2018"  "ZhangQ2019"   "ShirebyG2020"  "YangZ2016"    "ZhangY2017"
data('HorvathS2013')
head(coefs)

clock_name <- 'HorvathS2013'  # Select one of the supported clocks.
## Use Horvath's clock with adjusted-BMIQ normalisation (same as Horvath's paper)
horvath_age <- methyAge(betas, clock=clock_name)

print(horvath_age)
#horvath_age <- methyAge(betas, clock=clock_name, age_info=info, fit_method='Linear', do_plot=TRUE)

#https://github.com/yiluyucheng/dnaMethyAge/blob/main/R/MethyAge.R
#horvathPreprocess from: https://github.com/yiluyucheng/dnaMethyAge/blob/main/R/preprocessHorvathS2013.R
