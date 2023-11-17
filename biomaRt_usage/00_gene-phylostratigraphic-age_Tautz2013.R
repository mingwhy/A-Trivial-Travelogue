# myTAI: https://github.com/drostlab/published_phylomaps
# Josué Barrera-Redondo, Jaruwatana Sodai Lotharukpong, Hajk-Georg Drost & Susana M. Coelho, 2023
# mouse data

########
# Download Phylostratigraphic Maps in R:
# download the Phylostratigraphic Maps
# from Neme and Tautz, 2013
if(F){
download.file( url      = "https://static-content.springer.com/esm/art%3A10.1186%2F1471-2164-14-117/MediaObjects/12864_2012_4867_MOESM1_ESM.xlsx", 
               destfile = "BMCGenomics_2013_species_PhyloMaps.xlsx" )
}
########
# Read the *.xlsx file storing the Phylostratigraphic Maps and format it for use with myTAI:
library(readxl)
MusMusculusPhyloMap.BMCGenomics <- read_excel("BMCGenomics_2013_species_PhyloMaps.xlsx", sheet = 1)
dim(MusMusculusPhyloMap.BMCGenomics) #22773
head(MusMusculusPhyloMap.BMCGenomics)
colnames(MusMusculusPhyloMap.BMCGenomics)=c('ensembl_gene_id','Phylostratum.rank','TaxID','Phylostratum')

###########################################################
## reading in gene id.mapping in mouse aging atlas data
id.mapping=data.table::fread('~/Documents/Data_mouse_aging_atlas/fac_20449genes_id.mapping.txt')  #readin_h5ad.R
head(id.mapping)
dim(id.mapping) #20449

sum(id.mapping$ensembl_gene_id %in% MusMusculusPhyloMap.BMCGenomics$ensembl_gene_id) #17614
gene.meta=merge(id.mapping,MusMusculusPhyloMap.BMCGenomics)
dim(gene.meta) #17614

saveRDS(gene.meta,file='gene-phylostratigraphic-age_Tautz2013.rds')

###########################################################
###########################################################

if(F){
  # Converting UniProtKB ids into ENSEMBL gene ids: https://support.bioconductor.org/p/9145991/#9154237
  #https://support.bioconductor.org/p/9145372/
  #BiocManager::install("biomaRt"); 
  #if error message, check to see if need to install via 
  #install.packages('httr') #https://www.r-project.org/nosvn/pandoc/httr.html
  library(biomaRt)
  mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl") # mouse
  tmp=listAttributes(mart)
  tmp[grep('uniprot',ignore.case = T,tmp$name),]
  tmp[grep('symbol',ignore.case = T,tmp$name),]
  Uniprot = getBM(
    attributes=c('ensembl_gene_id','uniprotswissprot','uniprot_gn_id',
                 'uniprot_gn_symbol'),
    #'hgnc_symbol','mgi_symbol','uniprot_gn_symbol'), 
    mart = mart)
  colnames(Uniprot) <- c("Ensembl_ID", "uniprotswissprot" ,"uniprot_gn_id",'uniprot_gn_symbol')
  dim(Uniprot) #85710
  saveRDS(Uniprot,'mouse_EnsemblID_UniProt.rds')
  
  ##
  
  Uniprot=readRDS('mouse_EnsemblID_UniProt.rds')
  dim(Mus_musculus.PhyloMap) #21985
  length(unique(Mus_musculus.PhyloMap$Uniprot_ID)) #21985
  #https://www.informatics.jax.org/marker/MGI:1919847
  #gene: Auts2, A0A087WPF7
  Uniprot[Uniprot$uniprotswissprot=='A0A087WPF7',]
  Uniprot[Uniprot$uniprot_gn_symbol=='Cyhr1',]
  Uniprot$uniprot_gn_symbol_A=toupper(Uniprot$uniprot_gn_symbol)
  
  sum(Mus_musculus.PhyloMap$Uniprot_ID %in% Uniprot$uniprotswissprot) #15402
  sum(Mus_musculus.PhyloMap$Uniprot_ID %in% Uniprot$uniprot_gn_id) #19413
  sum(Mus_musculus.PhyloMap$gene_symbol %in% Uniprot$uniprot_gn_symbol_A) #7835

}
