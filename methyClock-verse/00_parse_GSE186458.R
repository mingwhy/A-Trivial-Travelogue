
############################################################################
# in Loyfer, Netanel, et al. "A DNA methylation atlas of normal human cell types." Nature (2023): 1-10.
# data can be accessed from: GSE186458
# `GSE186458_RAW.tar` is huge: 320 GB

############################################################################
# Downloading sample meta.data from GEO using GEOquery 
# http://genomicsclass.github.io/book/pages/GEOquery.html
library(GEOquery)
# Access the GEO Series Data
gse <- getGEO("GSE186458", GSEMatrix = TRUE)
show(gse)
dim(pData(gse[[1]])) #253 samples x 64 cols
head(pData(gse[[1]]))

############################################################################
# download GEO Series Data to local and process 
library("GEOquery")

# Download 'GSE186458_series_matrix.txt.gz' in Hannum et al. 
# data from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE186458
gse<-getGEO(filename='./GSE186458_series_matrix.txt.gz',GSEMatrix = TRUE,getGPL = FALSE) 

gse
#$GSE186458_series_matrix.txt.gz
#ExpressionSet (storageMode: lockedEnvironment)
#assayData: 0 features, 253 samples 

# getGEO source code: https://github.com/seandavi/GEOquery/blob/HEAD/R/getGEO.R
# func `parseGSEMatrix` https://github.com/seandavi/GEOquery/blob/7f416f5a0167d397c61be339bedfe35ae91e2fc1/R/parseGEO.R
# line 574
#datamat <- data.table::fread(text = dat[(series_table_begin_line + 1):(series_table_end_line -1)], 
#                             quote = "\"", na.strings = c("NA", "null", "NULL", "Null"))
df=pData(gse) #253 samples x 64 cols
colnames(df)
#View(df)
length(unique(df$`cell type:ch1`)) #39
df$title
tc=unique(paste(df$`tissue:ch1`,df$`cell type:ch1`,sep=':')) #79
tc[-grep('NA',tc)] #77 tc
