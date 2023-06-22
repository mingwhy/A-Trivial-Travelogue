
# debug reading in GSE40279_series_matrix.txt.gz 
library("GEOquery")

# Download 'GSE40279_series_matrix.txt.gz' from the Hannum et al. data from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE40279
gse<-getGEO(filename='./GSE40279/GSE40279_series_matrix.txt.gz',GSEMatrix = TRUE,getGPL = FALSE) 
#Warning message:
# In data.table::fread(text = dat[(series_table_begin_line + 1):(series_table_end_line -  :
#  Stopped early on line 324059. Expected 657 fields but found 0. Consider fill=TRUE and comment.char=. First discarded non-empty line: <<"cg18281842"

gse
#ExpressionSet (storageMode: lockedEnvironment)
#assayData: 324057 features, 656 samples 

# I checked on line (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM989829)
# each sample should contain 473034 rows
# the unzipped GSE40279_series_matrix.txt.gz contain 473105 rows, 1~70 encode sample info, the last row is '!series_matrix_table_end'.
# 473105 - 70 - 1 = 473034 CpG sites

#getGEO source code: https://github.com/seandavi/GEOquery/blob/HEAD/R/getGEO.R
# func `parseGSEMatrix` https://github.com/seandavi/GEOquery/blob/7f416f5a0167d397c61be339bedfe35ae91e2fc1/R/parseGEO.R
# line 574
#datamat <- data.table::fread(text = dat[(series_table_begin_line + 1):(series_table_end_line -
#                                                                         1)], quote = "\"", na.strings = c("NA", "null", "NULL", "Null"))

################################################################################################################
# I further checked `GSE40279_average_beta.txt.gz`, they contain the same information as `GSE40279_series_matrix.txt.gz`.
#https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM989827

################################################################################################################
# homemade read in code, adapted from https://rdrr.io/github/perishky/meffonym/src/vignettes/read-gse-matrix-file.r
# used case: https://rdrr.io/github/perishky/meffonym/f/vignettes/age-tutorial.rmd
filename='./GSE40279/GSE40279_series_matrix.txt.gz'

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
dim(samples) #656 x 35

rm(dat)
gc()

readLines(con,1) #get rid of '!series_matrix_table_begin' row
dnam <- read.table(con, sep="\t", header=TRUE, quote="\"", dec=".", fill=TRUE,
                   na.strings = c("NA", "null", "NULL", "Null"), comment.char = "")
#dnam[nrow(dnam),]
#the last row is: `!series_matrix_table_end`
close(con)

if( dnam[nrow(dnam),1] =='!series_matrix_table_end') dnam=dnam[-nrow(dnam),]

if (ndata == 0)
  dnam <- dnam[-(1:nrow(dnam)),]

if (nrow(dnam) > 0) 
  rownames(dnam) <- dnam[,1]
dnam <- as.matrix(dnam[,-1])

rownames(samples) <- colnames(dnam)
colnames(dnam) <- samples$geo_accession

list(dnam=dnam, samples=samples)

#######################################
## source function and usage

# Download 'GSE40279_series_matrix.txt.gz' from the Hannum et al. data from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE40279
source('read-gse-matrix-file.r')
hannum=read.gse.matrix.file('./GSE40279/GSE40279_series_matrix.txt.gz') #https://rdrr.io/github/perishky/meffonym/f/vignettes/age-tutorial.rmd
names(hannum) #~15min
#"dnam"    "samples"
saveRDS(hannum$dnam,'hannum_methy.mat.rds')
saveRDS(hannum$samples,'hannum_samples.rds')



