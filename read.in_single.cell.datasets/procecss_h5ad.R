# atlas info: https://github.com/czbiohub/tabula-muris-senis

# data source: https://figshare.com/projects/Tabula_Muris_Senis/64982
# TMS_processedFiles: Processed files (to use with scanpy)
# TMS.gene.data_final: TMS gene data (final) (code for paper: 2021-Mouse aging cell atlas analysis reveals global and cell type-specific aging signatures)
# TMS.Data.Objects: Tabula Muris Senis Data Objects


## read in h5ad file
#https://githubhot.com/repo/theislab/zellkonverter/issues/28
library(zellkonverter)
library(SummarizedExperiment)
sce1<-readH5AD('~/Documents/Data_mouse_aging_atlas/TMS_processedFiles/Kidney_facs.h5ad',use_hdf5 = TRUE)
sce2<-readH5AD('~/Documents/Data_mouse_aging_atlas/TMS.Data.Objects/Tabula Muris Senis Data Objects/tabula-muris-senis-facs-processed-official-annotations-Kidney.h5ad') #save as above
#sce<-readH5AD('~/Documents/Data_mouse_aging_atlas/TMS.gene.data_final/tabula-muris-senis-droplet-official-raw-obj.h5ad')
sce<-readH5AD('~/Documents/Data_mouse_aging_atlas/TMS.gene.data_final/tabula-muris-senis-facs-official-raw-obj.h5ad');       

# TMS.Data.Objects and TMS.gene.data_final have the same dimension
sce0<-sce[,sce$tissue=='Kidney']
dim(sce0);dim(sce1);dim(sce2)

cell.meta1=colData(sce1) #with batch information
table(cell.meta1$batch,cell.meta1$age) #batch and age are the same thing in this study

cell.meta2=colData(sce2) 
colnames(cell.meta2) #without batch information
cell.meta0=colData(sce0)
colnames(cell.meta0) #without batch information
table(cell.meta0$age)
table(cell.meta2$age)
sum(rownames(cell.meta0)==rownames(cell.meta2));dim(cell.meta2)

##
class(assay(sce))
sce

dim(assay(sce))
df.expr=assay(sce)
df.expr[1:3,1:3]
dim(df.expr) # fac: 22966 110824. droplet:10941   150.

genes=rownames(df.expr)
head(genes)

## sample meta information
cell.meta=colData(sce)
head(cell.meta)
dim(cell.meta) #fac:110824     13. droplet: 245389     12
colnames(cell.meta) #yes batch information for `TMS_processedFiles` and 

table(cell.meta$age) #6 age group: 1, 3, 18, 21, 24,30month
table(cell.meta$tissue) #16 tissue
length(table(cell.meta$cell_ontology_class)) #123
unique(paste(cell.meta$tissue,cell.meta$cell_ontology_class)) #169 unique

## more cell type annotation data (downstream analysis: https://github.com/czbiohub/tabula-muris-senis/blob/master/2_aging_signature/README.md)
meta2=data.table::fread('~/Documents/mouse_aging_atlas/TMS.gene.data_final/annotation_data/cell_ontology_class_functional_annotation.073020.tsv')
dim(meta2) #237
head(meta2)
colnames(meta2)
table(meta2$tissue); #except 'Heart_and_Aorta', 23 tissues in total, consistent with 2021 elife
table(meta2$`cell category`) # 6 cell functional classes
unique(paste(meta2$tissue,meta2$cell_ontology_class)) #237

## map gene id use biomart (https://www.biostars.org/p/301116/)
if(F){
  library(biomaRt)
  mouse= useEnsembl(version = 99, #this archived still contain dn ds values
                    biomart = 'ENSEMBL_MART_ENSEMBL', 
                    dataset = 'mmusculus_gene_ensembl')
  
  listAttributes(mouse) 
  x=listAttributes(mouse) 
  x[grep('symbol',x$name),]
  mouse_genes <- getBM(attributes = c('ensembl_gene_id', 
                                      'mgi_symbol','description',
                                      'chromosome_name',  
                                      'start_position', 
                                      'end_position'),
                       mart = mouse)
  tail(mouse_genes)
  saveRDS(mouse_genes,'mmusculus_gene_ensemblv99.rds')
}
mouse_genes=readRDS('mmusculus_gene_ensemblv99.rds')
sum('0610009B22Rik' %in% mouse_genes$mgi_symbol )

## write id mapping table
if(F){
  length(genes) #fac: 22966, droplet: 22899
  sum(genes %in% mouse_genes$mgi_symbol) #fac:20449 , droplet: 20404
  
  mouse_genes2=mouse_genes[mouse_genes$mgi_symbol %in% genes,]
  tmp=mouse_genes2[duplicated(mouse_genes2$mgi_symbol),]
  
  ok1=mouse_genes2[!mouse_genes2$mgi_symbol %in% tmp$mgi_symbol,]
  length(unique(ok1$mgi_symbol));dim(ok1) #20158 genes
  dups=mouse_genes2[mouse_genes2$mgi_symbol %in% tmp$mgi_symbol,]
  dups=dups[order(dups$mgi_symbol),]
  head(dups)
  length(unique(dups$mgi_symbol)) #291 genes
  
  dups[-grep('^\\d+$|X',dups$chromosome_name),]
  dups=dups[grep('^\\d+$|X',dups$chromosome_name),]
  length(unique(dups$mgi_symbol)) #291 genes
  
  tmp=dups[duplicated(dups$mgi_symbol),]
  ok2=dups[!dups$mgi_symbol %in% tmp$mgi_symbol,]
  length(unique(ok2$mgi_symbol));dim(ok2) #258 genes are ok now.
  dups=dups[dups$mgi_symbol %in% tmp$mgi_symbol,]
  View(dups); #choose one with longer gene description
  dim(dups)
  dups[dups$mgi_symbol=='Ndor1',]
  ok3=c()
  for(gene in unique(dups$mgi_symbol)){
    x=dups[dups$mgi_symbol==gene,]
    ok3=rbind(ok3,x[order(x$ensembl_gene_id)[1],])
  }
  dim(ok3) #33 genes
  ok=rbind(ok1,ok2,ok3)
  dim(ok) #20449 genes
  id.mapping=ok;
  sum(duplicated(id.mapping$mgi_symbol)); #0
  data.table::fwrite(id.mapping,'fac_20449genes_id.mapping.txt')
}
id.mapping=data.table::fread('fac_20449genes_id.mapping.txt')
dim(id.mapping) #20449 genes

## shown on https://github.com/czbiohub/tabula-muris-senis/blob/master/2_aging_signature/README.md
# `Differential gene expression (DGE) testing`, Per-tissue .rds format TMS data: tms_gene_data/rds_by_tissue.1e4.zip
# there are 23 `facs.normalized.XXX.rds` files, consistent with 23 tissue in 2021 elife paper
# DE analysis code: https://github.com/czbiohub/tabula-muris-senis/blob/master/2_aging_signature/job.DGE_analysis/DGE_analysis.R

##`Specifically, the gene partition statistics may be of interests to researchers tms_gene_data/result_v1/tms_gene_table/gene_stats_*.gz. The columns are:`

de=data.table::fread('./TMS.gene.data_final/result_v1/tms_gene_table/gene_stats_facs')
dim(de) #22966 genes x 987 attributes
colnames(de)[1:10]
de.sub=de[,1:10]
sum(de.sub$global) #330 GAG genes

gag=de.sub[de.sub$global,]
dim(gag) #330 x 10
head(gag)
table(gag$global.dir) #global.dir: direction of the GAG (if the gene is a GAG). This is consistent with Supp. Table. 3
#down other    up 
#190    47    9

up.genes=gag[gag$global.dir=='down',]$gene

library(clusterProfiler)
library(org.Mm.eg.db)
ego <- enrichGO(gene         = up.genes,
                 OrgDb         = org.Mm.eg.db,
                 keyType       = 'SYMBOL',
                 ont           = "BP",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.01,
                 qvalueCutoff  = 0.05)
dim(ego) 
#goplot(ego)
sum(ego@result$qvalue<0.05)
ego@result[ego@result$qvalue<0.05,]$Description
df=ego@result
df[grep('MHC',df$Description,ignore.case = T),]
df[grep('apoptosis',df$Description,ignore.case = T),]


## read in dnds data
filter_mouse_rat=data.table::fread('~/Documents/sc_transcriptome.index/dNdS/mouse_rat.dnds.txt')
head(filter_mouse_rat)
dim(filter_mouse_rat) #30033 genes
dim(id.mapping) #20449
sum(filter_mouse_rat$ensembl_gene_id %in% id.mapping$ensembl_gene_id) #20697

mouse_rat_dnds=filter_mouse_rat[!duplicated(filter_mouse_rat$ensembl_gene_id),]
dim(mouse_rat_dnds) #20048
sum(mouse_rat_dnds$ensembl_gene_id %in% id.mapping$ensembl_gene_id) #17029

gene.meta=merge(gag[,c('gene','global.dir')],id.mapping,by.x='gene',by.y='mgi_symbol')
gene.meta2=merge(gene.meta,mouse_rat_dnds,by.x='ensembl_gene_id',by.y='ensembl_gene_id')
dim(gene.meta2) #298 genes remain
gene.meta2=gene.meta2[order(gene.meta2$global.dir),]
View(gene.meta2)
