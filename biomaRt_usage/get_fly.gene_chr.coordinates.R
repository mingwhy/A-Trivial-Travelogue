library(biomaRt)
packageVersion("biomaRt") #2.48.2
  
biolist <- as.data.frame(listMarts())
ensembl=useMart("ensembl")
esemblist <- as.data.frame(listDatasets(ensembl))
#dataset                              description  version
#55 dmelanogaster_gene_ensembl Drosophila melanogaster genes (BDGP6.32) BDGP6.32


ensembl = useDataset("dmelanogaster_gene_ensembl",mart=ensembl)
filters = listFilters(ensembl)
attributes = listAttributes(ensembl)

attributes[grep('dmel',attributes$name),]

t2g<-getBM(attributes=c('ensembl_gene_id','chromosome_name','start_position','end_position'), mart = ensembl)
dim(t2g) #23932     4
head(t2g)
saveRDS(t2g,'t2g_chr.coord.rds')
