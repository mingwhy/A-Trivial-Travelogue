###############################################################
## use biomaRt to get gene chr info 
# https://github.com/mingwhy/bioinfo_homemade_tools/tree/main/used.cases/biomaRt_usage
# add gene_biotype info: https://support.bioconductor.org/p/62441/
if(F){ #run once
  #https://stackoverflow.com/questions/13012210/find-transcription-start-sites-with-biomart
  library(biomaRt)
  packageVersion("biomaRt") #‘2.50.3’
  
  biolist <- as.data.frame(listMarts())
  ensembl=useMart("ensembl")
  esemblist <- as.data.frame(listDatasets(ensembl))
  esemblist[grep('musculus',esemblist$dataset),]
  #dataset                              description  version
  #108 mmusculus_gene_ensembl           Mouse genes (GRCm39)      GRCm39
  
  ensembl = useDataset("mmusculus_gene_ensembl",mart=ensembl)
  filters = listFilters(ensembl)
  attributes = listAttributes(ensembl)
  
  attributes[grep('mmusculus',attributes$name),]
  attributes[grep('symbol',attributes$name),]
  grep('transcript',attributes$name,value=TRUE,ignore.case = TRUE) #feature_page
  attributes[grep('ensembl_transcript_id',attributes$name),]
  grep('peptide',attributes$name,value=TRUE,ignore.case = TRUE) #feature_page
  attributes[grep('cds_start',attributes$name),] #not in feature_page
  attributes[grep('exon_chrom_start',attributes$name),] #not in feature_page
  
  #https://www.biostars.org/p/9544009/#9544241
  ensemble_data<-getBM(attributes=c('mgi_symbol','ensembl_gene_id',
                          "ensembl_transcript_id","ensembl_peptide_id",
                          "strand","gene_biotype","ensembl_exon_id"),
                          #"cds_start","cds_end","exon_chrom_start","exon_chrom_end"),
                          mart = ensembl)
  
  dim(ensemble_data) #868930     7
  head(ensemble_data)
  
  table(ensemble_data$gene_biotype)
  protein_coding_genes<-ensemble_data[ensemble_data$gene_biotype=='protein_coding',]
  dim(protein_coding_genes) #759189
  length(unique(protein_coding_genes$mgi_symbol)) #21739
  length(unique(protein_coding_genes$ensembl_peptide_id)) #66310
  
  tmp=protein_coding_genes[protein_coding_genes$ensembl_peptide_id %in% 
    protein_coding_genes[duplicated(protein_coding_genes$ensembl_peptide_id),]$ensembl_peptide_id,]
  head(tmp)
  
  #https://rdrr.io/bioc/biomaRt/man/getSequence.html
  seq = getSequence(id = protein_coding_genes$mgi_symbol[1:2], 
                    type = "mgi_symbol", 
                    seqType = "peptide", 
                    mart = ensembl)
  show(seq)
  length(unique(protein_coding_genes$mgi_symbol)) #21739
  peptide_seqs <- getSequence(id = unique(protein_coding_genes$mgi_symbol),
                              type = "mgi_symbol", 
                              seqType = "peptide", 
                              mart = ensembl)
  dim(peptide_seqs) #69799
  
  saveRDS(peptide_seqs, 'mmusculus_peptide_seqs.rds')
}

peptide_seqs=readRDS('mmusculus_peptide_seqs.rds')


