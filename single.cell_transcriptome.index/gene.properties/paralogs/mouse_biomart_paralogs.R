
#Retrieving Dn Or Ds From Ensembl Using Biomart
library(biomaRt)
#Every analysis with biomaRt starts with selecting a BioMart database to use.
#https://bioconductor.org/packages/release/bioc/vignettes/biomaRt/inst/doc/accessing_ensembl.html
if(F){
  # 1. Selecting an Ensembl BioMart database and dataset
  # step1, identifying the database you need
  listEnsembl()
  ensembl <- useEnsembl(biomart = "genes")
  ensembl
  # step2, Choosing a dataset
  datasets <- listDatasets(ensembl)
  head(datasets)
  datasets[grep('musculus',datasets$dataset),]
  datasets[grep('melanogaster',datasets$dataset),]
  #Using archived versions of Ensembl
  listEnsemblArchives()
}
######################################################################################
# retrieve dn/ds between M.musculus and the rat Rattus norvegicus using their ensembly genomes
#https://www.biostars.org/p/147351/
mouse= useEnsembl(version = 99, #this archived still contain dn ds values
                  biomart = 'ENSEMBL_MART_ENSEMBL', 
                  dataset = 'mmusculus_gene_ensembl')

listAttributes(mouse) 
searchAttributes(mouse,'paralog')

## mouse and Rat (https://support.bioconductor.org/p/9141035/)
mouse_paralogs <- getBM(attributes = c('ensembl_gene_id', 
                                  'mmusculus_paralog_ensembl_gene',  
                                  'mmusculus_paralog_associated_gene_name',
                                  'mmusculus_paralog_orthology_type',
                                  'mmusculus_paralog_perc_id',
                                  'mmusculus_paralog_perc_id_r1',
                                  'mmusculus_paralog_dn',
                                  'mmusculus_paralog_ds',
                                  'mmusculus_paralog_paralogy_confidence'), 
                   mart = mouse)
head(mouse_paralogs)
dim(mouse_paralogs) #2439504     9

filter_mouse_paralogs=mouse_paralogs[mouse_paralogs$mmusculus_paralog_ensembl_gene!='',]
dim(filter_mouse_paralogs) #2405788     9
head(filter_mouse_paralogs)
table(filter_mouse_paralogs$mmusculus_paralog_orthology_type)
#gene_split          other_paralog within_species_paralog 
#64                2143962                 261762 
#https://uswest.ensembl.org/info/genome/compara/homology_types.html

filter_mouse_paralogs[filter_mouse_paralogs$ensembl_gene_id=='ENSMUSG00000064345',]
filter_mouse_paralogs[filter_mouse_paralogs$ensembl_gene_id=='ENSMUSG00000064367',]

data.table::fwrite(filter_mouse_paralogs,'filter_mouse_paralogs.txt')

