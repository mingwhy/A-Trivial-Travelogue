
library(zellkonverter)
library(SingleCellExperiment)
library(scater);library(scran)
library(ggplot2);library(gridExtra);
library(tidyverse);
library(zellkonverter)
library(SummarizedExperiment)
library(viridis)
library(ggpubr)


###############################################################
## use biomaRt to get gene chr info 
# https://github.com/mingwhy/bioinfo_homemade_tools/tree/main/used.cases/biomaRt_usage
if(F){ #run once
#https://stackoverflow.com/questions/13012210/find-transcription-start-sites-with-biomart
library(biomaRt)
packageVersion("biomaRt") #‘2.50.3’

biolist <- as.data.frame(listMarts())
ensembl=useMart("ensembl")
esemblist <- as.data.frame(listDatasets(ensembl))
esemblist[grep('melanogaster',esemblist$description),]
#dataset                              description  version
#56 dmelanogaster_gene_ensembl Drosophila melanogaster genes (BDGP6.32) BDGP6.32

ensembl = useDataset("dmelanogaster_gene_ensembl",mart=ensembl)
filters = listFilters(ensembl)
attributes = listAttributes(ensembl)

attributes[grep('dmel',attributes$name),]
grep('transcript',attributes$name,value=TRUE,ignore.case = TRUE)

t2g<-getBM(attributes=c('ensembl_gene_id','chromosome_name','start_position','end_position',
                        'ensembl_transcript_id','flybase_transcript_id',
                        "transcript_start", "transcript_end",
                        'transcription_start_site', "transcript_length","transcript_count"), mart = ensembl)
dim(t2g) #41209    11
head(t2g)
saveRDS(t2g,'t2g_chr.coord.rds')
}
t2g=readRDS('t2g_chr.coord.rds')

###############################################################
## read in data
sce=readH5AD('adata_headBody_S_v1.0.h5ad') # 22966 31001 
sce #SingleCellExperiment,  15992 566254 
cell.meta=colData(sce)
colnames(cell.meta)
table(cell.meta$sex,cell.meta$age)
head(cell.meta)
#table(cell.meta$afca_annotation_broad,cell.meta$afca_annotation)
length(table(cell.meta$afca_annotation)) # 163 cell types

genes=rownames(sce)
length(genes) #15992 genes

################################################################
## use AnnotationDbi and org.Dm.eg.db to get gene symbols 
# https://github.com/mingwhy/bioinfo_homemade_tools/blob/main/used.cases/biomaRt_usage/gene.id_converter_AnnotationDbi_bitr.R
library(AnnotationDbi)
library(org.Dm.eg.db)
package.version("org.Dm.eg.db") # "3.14.0"
columns(org.Dm.eg.db)
keytypes(org.Dm.eg.db)

test.inp=genes
test.out=select(org.Dm.eg.db, keys=test.inp, keytype="SYMBOL",
                #columns=c("SYMBOL","GENENAME",'FLYBASE','ALIAS','ACCNUM','ENTREZID') )
                columns=c("SYMBOL","GENENAME",'FLYBASE','ENTREZID') )
length(test.inp) #15992 genes
dim(test.out) #15992 x 2
head(test.out)
test.out$original.id=test.inp;

## keep unique mapped ones
keep1=test.out[!is.na(test.out$GENENAME),]
dim(keep1) #15937   

## save unmapped ones aned submit to flybase for more information
sum(is.na(test.out$GENENAME)) #55
write.table(test.out[is.na(test.out$GENENAME),]$SYMBOL,'genes_for_flybase.txt',quote=F,row.names = F,col.names = F)

# go to http://flybase.org/batchdownload
df=data.table::fread('FlyBase_Fields_download.txt')
dim(df) #73
df=df[!df$ANNOTATION_SYMBOL=='-',]
length(unique(df$original.id)) #54 genes could be used
df$SYMBOL
colnames(df)[1]='original.id'
#View(df[df$original.id %in% df[duplicated(df$original.id),]$original.id,])
# only 'Gtp-bp' has two different mapping result
df[df$original.id=='Gtp-bp',]
test.out[test.out$SYMBOL=='128up',] #already exists
test.out[test.out$SYMBOL=='SrpRalpha',] #not exists 
which(df$original.id=='Gtp-bp')
df=df[-26,]
df.unique=df[!duplicated(df$original.id),]
dim(df.unique) #54

sum(df.unique$SYMBOL %in% keep1$SYMBOL) #no overlap between these two
head(keep1)
head(df.unique)
sum(df.unique$FBID_KEY %in% t2g$ensembl_gene_id) #54 all exist
sum(keep1$FLYBASE %in% t2g$ensembl_gene_id) # 15937 all exist

# integrate gene symbol, flybaseId, chr, start.position, end.position
keep2=df.unique[,c('original.id','FBID_KEY','SYMBOL','NAME')]
keep1=keep1[,c(5,3,1,2)]
colnames(keep1)
colnames(keep2)=colnames(keep1)
dfc=rbind(keep1,keep2)
df.gene=merge(dfc,t2g,by.x='FLYBASE',by.y='ensembl_gene_id')
dim(df.gene) #15991
head(df.gene)
saveRDS(df.gene,'AFCA_gene.meta.rds')

