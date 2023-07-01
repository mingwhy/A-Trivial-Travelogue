#https://www.biostars.org/p/429062/
library(biomaRt)
packageVersion("biomaRt") #‘2.50.3’

biolist <- as.data.frame(listMarts())
#               biomart                version
#1 ENSEMBL_MART_ENSEMBL      Ensembl Genes 109
#2   ENSEMBL_MART_MOUSE      Mouse strains 109
#3     ENSEMBL_MART_SNP  Ensembl Variation 109
#4 ENSEMBL_MART_FUNCGEN Ensembl Regulation 109

ensembl=useMart("ensembl")
esemblist <- as.data.frame(listDatasets(ensembl))
esemblist[grep('musculus',esemblist$dataset),]
#dataset                              description  version
#108 mmusculus_gene_ensembl           Mouse genes (GRCm39)      GRCm39

ensembl = useDataset("mmusculus_gene_ensembl",mart=ensembl)
filters = listFilters(ensembl)
attributes = listAttributes(ensembl)

annotLookup <- getBM(
  mart = ensembl,
  attributes = c(
    'ensembl_gene_id',
    'gene_biotype',
    'external_gene_name',
    'uniprot_gn_symbol',
    'uniprot_gn_id'),
  uniqueRows=TRUE)
head(annotLookup)
head(subset(annotLookup, uniprot_gn_id != ''), 20)[,-4]

annotLookup2=annotLookup[annotLookup$uniprot_gn_id!='',]
dim(annotLookup2) #51706
length(unique(annotLookup2$uniprot_gn_id)) #51102

saveRDS(annotLookup2,'mmus_id_ensembl2uniprot.rds')

####################################
## reading in gene id.mapping
id.mapping=data.table::fread('~/Documents/Data_mouse_aging_atlas/fac_20449genes_id.mapping.txt')  #readin_h5ad.R
colnames(annotLookup2)
colnames(id.mapping)
gene.meta=merge(id.mapping,annotLookup2)
dim(gene.meta) #44806    10
length(unique(gene.meta$uniprot_gn_id)) #44744
length(unique(gene.meta$mgi_symbol)) #18060, one gene -> multiple peptides

saveRDS(gene.meta,'mmus_id_ensembl2uniprot.rds')

