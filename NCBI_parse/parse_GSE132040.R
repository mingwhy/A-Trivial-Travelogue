
# sample metainfo: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE132040
################################################################################################################
# homemade read in code, adapted from https://rdrr.io/github/perishky/meffonym/src/vignettes/read-gse-matrix-file.r
# used case: https://rdrr.io/github/perishky/meffonym/f/vignettes/age-tutorial.rmd
filename='./GSE132040_series_matrix.txt.gz'

dat <- readLines(filename)
str(dat) 
#chr [1:473105] "!Series_title\t\"Geno
nseries <- sum(grepl("^!Series_", dat)) #32 rows
nsamples <- sum(grepl("^!Sample_", dat)) #35 rows

ndata <- length(dat) - match("!series_matrix_table_begin", dat) - 2 #473034 CpG site

# begin read in via file connection
con <- file(filename, "r")
header <- read.table(con, sep="\t", header=F, nrows=nseries, stringsAsFactors=F)
samples <- read.table(con, sep="\t", header=F, nrows=nsamples, stringsAsFactors=F)

samples <- t(samples)
colnames(samples) <- samples[1,]
colnames(samples) <- sub("!Sample_", "", colnames(samples))
samples <- data.frame(samples[-1,], stringsAsFactors=F)
dim(samples) #947  43

rm(dat)
gc()

df=samples[,c(1,2,8,9,10,11,12,13)]
colnames(df)=c('title','geo_accession',
               'tissue_id','species','strain','tissue','age','sex')
head(df)
table(df$age,df$tissue,df$sex)
