
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

x=searchAttributes(mouse,'GO')
x
x$name #7 related attributes

mouse_go=getBM(attributes = c('ensembl_gene_id', 
                     x$name[1:4]), 
      mart = mouse)
dim(mouse_go)#562588  
sum(mouse_go$go_id=='') #47693
mouse_go=mouse_go[mouse_go$go_id!='',]
dim(mouse_go)#514895
length(unique(mouse_go$ensembl_gene_id)) #22260 gene IDs

saveRDS(mouse_go,'mouse_GOterms.rds') #smaller than txt
#data.table::fwrite(mouse_go,'mouse_GOterms.txt')

if(F){
#time out problem: https://www.biostars.org/p/9480321/#9513378
mouse_go = lapply(x$name,function(one.attr){
    getBM(attributes = c('ensembl_gene_id', 
                                   one.attr), 
                      mart = mouse)})
length(mouse_go)
sapply(mouse_go,dim)
head(mouse_go[[1]])
head(mouse_go[[2]])
}

