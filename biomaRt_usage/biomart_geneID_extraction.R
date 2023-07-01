
#https://stackoverflow.com/questions/48262384/increase-r-studios-connection-timeout-limit
options(timeout = max(4000000, getOption("timeout")))

#https://support.bioconductor.org/p/9139740/
#https://www.biostars.org/p/429062/
library(biomaRt)
packageVersion("biomaRt") #‘2.50.3’

biolist <- as.data.frame(listMarts()) 
biolist #biolist contain version information
#biomart                version
#1 ENSEMBL_MART_ENSEMBL      Ensembl Genes 109
#2   ENSEMBL_MART_MOUSE      Mouse strains 109
#3     ENSEMBL_MART_SNP  Ensembl Variation 109
#4 ENSEMBL_MART_FUNCGEN Ensembl Regulation 109

######################
## mouse
ensembl=useMart("ensembl") 
esemblist <- as.data.frame(listDatasets(ensembl))
esemblist[grep('musculus',esemblist$dataset),]
#dataset                              description  version
#mmusculus_gene_ensembl           Mouse genes (GRCm39)      GRCm39

mouse= useEnsembl(#version = 99, #this archived still contain dn ds values
                  version = 109,
                  biomart = 'ENSEMBL_MART_ENSEMBL', 
                  dataset = 'mmusculus_gene_ensembl')
filters = listFilters(mouse)
attributes = listAttributes(mouse)

annotLookup <- getBM(
  mart = mouse,
  attributes = c(
    'ensembl_gene_id',
    'gene_biotype',
    'external_gene_name',#'external_synonym','source_name',
    'uniprot_gn_symbol',
    'uniprot_gn_id',
    'mgi_symbol','description',
    'chromosome_name', 'start_position','end_position'),
  uniqueRows=TRUE)
head(annotLookup)

#head(subset(annotLookup, uniprot_gn_id != ''), 20)[,-4]
annotLookup2=annotLookup[annotLookup$uniprot_gn_id!='',]
dim(annotLookup) #86873
length(unique(annotLookup2$uniprot_gn_id)) #51102

saveRDS(annotLookup,'mmusculus_id_v109.rds')

######################
## rat
ensembl=useMart("ensembl") 
esemblist <- as.data.frame(listDatasets(ensembl))
esemblist[grep('norvegicus',esemblist$dataset),]
#dataset                              description  version
#rnorvegicus_gene_ensembl Rat genes (mRatBN7.2) mRatBN7.2

rat= useEnsembl(#version = 99, #this archived still contain dn ds values
  version = 109,
  biomart = 'ENSEMBL_MART_ENSEMBL', 
  dataset = 'rnorvegicus_gene_ensembl')
filters = listFilters(rat)
attributes = listAttributes(rat)
attributes[attributes$page=='feature_page',]

annotLookup <- getBM(
  mart = rat,
  attributes = c(
    'ensembl_gene_id',
    'gene_biotype',
    'external_gene_name','external_synonym','source_name',
    'uniprot_gn_symbol',
    'uniprot_gn_id',
    'mgi_symbol','description',
    'chromosome_name', 'start_position','end_position'),
  uniqueRows=TRUE)
head(annotLookup)

#head(subset(annotLookup, uniprot_gn_id != ''), 20)[,-4]
annotLookup2=annotLookup[annotLookup$uniprot_gn_id!='',]
dim(annotLookup) #50429
length(unique(annotLookup2$uniprot_gn_id)) #14983

saveRDS(annotLookup,'rnorvegicus_id_v109.rds')


######################
## human
ensembl=useMart("ensembl") 
esemblist <- as.data.frame(listDatasets(ensembl))
esemblist[grep('sapiens',esemblist$dataset),]
#dataset                              description  version
#hsapiens_gene_ensembl Human genes (GRCh38.p13) GRCh38.p13

human= useEnsembl(#version = 99, #this archived still contain dn ds values
  version = 109,
  biomart = 'ENSEMBL_MART_ENSEMBL', 
  dataset = 'hsapiens_gene_ensembl')
filters = listFilters(human)
attributes = listAttributes(human)
attributes[attributes$page=='feature_page',]

annotLookup <- getBM(
  mart = human,
  attributes = c(
    'ensembl_gene_id',
    'gene_biotype',
    'external_gene_name','external_synonym','source_name',
    'uniprot_gn_symbol',
    'uniprot_gn_id',
    #'mgi_symbol',
    'description',
    'chromosome_name', 'start_position','end_position'),
  uniqueRows=TRUE)
head(annotLookup)

#head(subset(annotLookup, uniprot_gn_id != ''), 20)[,-4]
annotLookup2=annotLookup[annotLookup$uniprot_gn_id!='',]
dim(annotLookup) #433733
length(unique(annotLookup2$uniprot_gn_id)) #71809

saveRDS(annotLookup,'hsapiens_id_v109.rds')



