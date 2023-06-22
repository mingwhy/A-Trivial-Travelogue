#Retrieving Dn Or Ds From Ensembl Using Biomart
#Calculating dn/ds ratios

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
    #dataset                    description     version
    #18  bmusculus_gene_ensembl Blue whale genes (mBalMus1.v2) mBalMus1.v2
    #107 mmusculus_gene_ensembl           Mouse genes (GRCm39)      GRCm39
    datasets[grep('melanogaster',datasets$dataset),]
    #dataset                              description  version
    #55 dmelanogaster_gene_ensembl Drosophila melanogaster genes (BDGP6.32) BDGP6.32
    
    #Using archived versions of Ensembl
    listEnsemblArchives()
}

# as a used case show here: https://support.bioconductor.org/p/9141035/
## For the human data 
ensemblhsapiens = useEnsembl(version = 99, #this archived still contain dn ds values
                             biomart = 'ENSEMBL_MART_ENSEMBL', 
                             dataset = 'hsapiens_gene_ensembl')
ensemblhsapiens
#Object of class 'Mart':
#Using the ENSEMBL_MART_ENSEMBL BioMart database
#Using the hsapiens_gene_ensembl dataset

#https://rdrr.io/bioc/biomaRt/man/listAttributes.html
listAttributes(ensemblhsapiens) 
searchAttributes(ensemblhsapiens, 'homolog_dn') #yes, there is xxx_dn

## Input human gene list
hsapiens_GFList <- c("ENSG00000010404", "ENSG00000277796", "ENSG00000198888")

## human and Rat
hsapiens_rat <- getBM(attributes = c('ensembl_gene_id', 
                                     'rnorvegicus_homolog_ensembl_gene',  
                                     'rnorvegicus_homolog_dn', 
                                     'rnorvegicus_homolog_ds',
                                     'rnorvegicus_homolog_orthology_type'), 
                      filters = 'ensembl_gene_id', 
                      values = hsapiens_GFList, 
                      mart = ensemblhsapiens)

## human and celegans
hsapiens_celegans <- getBM(attributes = c('ensembl_gene_id', 
                                          'celegans_homolog_ensembl_gene',  
                                          'celegans_homolog_dn', 
                                          'celegans_homolog_ds',
                                          'celegans_homolog_orthology_type'), 
                           filters = 'ensembl_gene_id', 
                           values = hsapiens_GFList, 
                           mart = ensemblhsapiens)

hsapiens_rat
hsapiens_celegans

#quote"However I think it's telling that support for the dN/dS analysis has been dropped for over a year now, 
#from Ensembl 100 onwards (https://www.ensembl.info/2020/04/29/ensembl-100-has-been-released/)."


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

#####
fly= useEnsembl(version = 99, #this archived still contain dn ds values
                  biomart = 'ENSEMBL_MART_ENSEMBL', 
                  dataset = 'dmelanogaster_gene_ensembl')
x=listAttributes(fly) 
head(x,20)
head(x[grep('homolog',x$name),])

fly.genes <- getBM(attributes = c('ensembl_gene_id', 
                                  'chromosome_name',  
                                  'start_position', 
                                  'end_position',
                                  'strand',
                                  'percentage_gene_gc_content'),
                                  mart=fly)
head(fly.genes)                   
dim(fly.genes) # 17807     6
