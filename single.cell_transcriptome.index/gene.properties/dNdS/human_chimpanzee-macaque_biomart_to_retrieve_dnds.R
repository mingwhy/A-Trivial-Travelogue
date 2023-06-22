
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
searchAttributes(mouse,'rnorvegicus_homolog')
searchAttributes(mouse, 'rnorvegicus_homolog_dn') #yes, there is xxx_dn
searchAttributes(mouse, 'rnorvegicus_homolog_ds') #yes, there is xxx_dn


## mouse and Rat (https://support.bioconductor.org/p/9141035/)
mouse_rat <- getBM(attributes = c('ensembl_gene_id', 
                                     'rnorvegicus_homolog_ensembl_gene',  
                                     'rnorvegicus_homolog_dn', 
                                     'rnorvegicus_homolog_ds',
                                     'rnorvegicus_homolog_orthology_type',
                                  'rnorvegicus_homolog_orthology_confidence'), 
                      mart = mouse)
head(mouse_rat)
dim(mouse_rat) #70177     5

filter_mouse_rat=mouse_rat[!is.na(mouse_rat$rnorvegicus_homolog_dn),]
dim(filter_mouse_rat) #30033     5
head(filter_mouse_rat)

data.table::fwrite(filter_mouse_rat,'mouse_rat.dnds.txt')

