# myTAI: https://github.com/drostlab/published_phylomaps
# Josué Barrera-Redondo, Jaruwatana Sodai Lotharukpong, Hajk-Georg Drost & Susana M. Coelho, 2023
# mouse data

########
# Download Phylostratigraphic Maps in R:
# [Animals] from Barrera-Redondo et al., 2023
if(!file.exists('Barrera-Redondo_2023_Maps_animal.xlsx')){
  download.file( url      = "https://static-content.springer.com/esm/art%3A10.1186%2Fs13059-023-02895-z/MediaObjects/13059_2023_2895_MOESM4_ESM.xlsx", 
               destfile = "Barrera-Redondo_2023_Maps_animal.xlsx" )
}
########
# Read the *.xlsx file storing the Phylostratigraphic Maps and format it for use with myTAI:
# load package readxl
if(F){
  library(readxl)
  library(dplyr)
  #	House mouse
  Mus_musculus.data <- 
    read_excel("Barrera-Redondo_2023_Maps_animal.xlsx", sheet = "10090_gene_ages")
  head(Mus_musculus.data)
  Mus_musculus.PhyloMap <- 
    dplyr::select(
      Mus_musculus.data,
      Phylostratum = phylostratum,
      Phylostratum.rank = rank,
      GeneID = `#gene`
    ) 
  dim(Mus_musculus.PhyloMap) #21985
  
  states=unlist(lapply(strsplit(Mus_musculus.PhyloMap$GeneID,'\\|'),'[[',1))
  table(states)
  #sp    tr 
  #17046  4939 
  Mus_musculus.PhyloMap$states=states;
  
  Uniprot_ID=unlist(lapply(strsplit(Mus_musculus.PhyloMap$GeneID,'\\|'),'[[',2))
  Mus_musculus.PhyloMap$Uniprot_ID=Uniprot_ID;
  
  gene_symbol=unlist(lapply(strsplit(Mus_musculus.PhyloMap$GeneID,'\\|'),'[[',3))
  gene_symbol=gsub('_MOUSE','',gene_symbol)
  Mus_musculus.PhyloMap$gene_symbol=gene_symbol;
  saveRDS(Mus_musculus.PhyloMap,'Mus_musculus.PhyloMap.rds')
}
Mus_musculus.PhyloMap=readRDS('Mus_musculus.PhyloMap.rds')
dim(Mus_musculus.PhyloMap) #21985
summary(Mus_musculus.PhyloMap$Phylostratum.rank) 

# do id mapping online:https://www.uniprot.org/id-mapping
# UniProtKB AC/ID -> UniProtKB (idmapping_2023_11_14.tsv.gz)
# UniProtKB AC/ID -> Genome Annotation Database -> Ensembl Transcript (idmapping_2023_11_14_Ensembl-Transcript.tsv.gz)
# use this one: UniProtKB AC/ID -> Genome Annotation Database -> Ensembl  (idmapping_2023_11_14_Ensembl.tsv.gz)
# as in id.mapping file of TMS, it's `ENSMUSG0xxx`
dim(Mus_musculus.PhyloMap)
table(Mus_musculus.PhyloMap$states)
#write.table(file='protein.names.txt', Mus_musculus.PhyloMap$Uniprot_ID,row.names = F,col.names = F,quote=F)

protein.id=data.table::fread('idmapping_2023_11_14_Ensembl.tsv.gz')
dim(protein.id) #20602
head(protein.id)
colnames(protein.id)=c('Uniprot_ID','ensembl_gene_id')
Mmus.PhyloMap=merge(Mus_musculus.PhyloMap,protein.id,by.x='Uniprot_ID',by.y='Uniprot_ID')
dim(Mmus.PhyloMap) #20602

###########################################################
## reading in gene id.mapping in mouse aging atlas data
id.mapping=data.table::fread('~/Documents/Data_mouse_aging_atlas/fac_20449genes_id.mapping.txt')  #readin_h5ad.R
head(id.mapping)
dim(id.mapping) #20449
sum(id.mapping$ensembl_gene_id %in% Mmus.PhyloMap$ensembl_gene_id) #13893
Mmus.PhyloMap$ensembl_gene_id2=unlist(lapply(strsplit(Mmus.PhyloMap$ensembl_gene_id,'\\.'),'[[',1))

sum(id.mapping$ensembl_gene_id %in% Mmus.PhyloMap$ensembl_gene_id2) #16876
gene.meta=merge(id.mapping,Mmus.PhyloMap,by.x='ensembl_gene_id',by.y='ensembl_gene_id2')
dim(gene.meta) #16914

length(unique(gene.meta$ensembl_gene_id)) #16876
tmp=gene.meta[duplicated(gene.meta$ensembl_gene_id),]
keep1=gene.meta[!gene.meta$ensembl_gene_id %in% tmp$ensembl_gene_id,]

tmp=gene.meta[gene.meta$ensembl_gene_id %in% tmp$ensembl_gene_id,]
keep.gene=c()
for(gene in unique(tmp$ensembl_gene_id)){
  x=tmp[tmp$ensembl_gene_id==gene,]
  if(length(unique(x$Phylostratum.rank))==1){ keep.gene=c(keep.gene,gene)}
}
length(keep.gene) #18
keep2=tmp[tmp$ensembl_gene_id %in% keep.gene,]
keep2=keep2[!duplicated(keep2$ensembl_gene_id),]

tmp2=tmp[!tmp$ensembl_gene_id %in% keep.gene,]
# mannual check
keep3=tmp2[c(2,3,5),]

gene.meta=rbind(keep1,keep2,keep3)
dim(gene.meta)
length(unique(gene.meta$ensembl_gene_id)) #16876
saveRDS(gene.meta,file='gene-phylostratigraphic-age_Coelho2023.rds')

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
