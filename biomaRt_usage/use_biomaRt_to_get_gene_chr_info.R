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
grep('transcript',attributes$name,value=TRUE,ignore.case = TRUE)

t2g<-getBM(attributes=c('mgi_symbol','ensembl_gene_id',"gene_biotype",
                        "ensembl_transcript_id","transcript_start","transcript_end",
                        #"ensembl_exon_id","exon_chrom_start","exon_chrom_end",
                        "chromosome_name","strand","transcript_biotype", 
                        'start_position','end_position',
                        'ensembl_transcript_id',"transcript_length","transcript_count"), mart = ensembl)
dim(t2g) #149482     13
head(t2g)
saveRDS(t2g,'mmusculus_t2g_chr.coord.rds')
}
t2g=readRDS('mmusculus_t2g_chr.coord.rds')
# one gene has multiple rows as muliple transcripts
df.gene.length=t2g[!duplicated(t2g$ensembl_gene_id),]
table(df.gene.length$gene_biotype)

df.gene.length=df.gene.length[df.gene.length$gene_biotype == 'protein_coding',]
dim(df.gene.length) #21959    13

table(df.gene.length$chromosome_name)
df.gene.length=df.gene.length[df.gene.length$chromosome_name %in% c(1:19,'X','Y'),]
dim(df.gene.length) #21895    13
