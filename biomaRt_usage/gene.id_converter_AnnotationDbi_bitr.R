
#BiocManager::install("AnnotationDbi")
#BiocManager::install("org.Dm.eg.db")
library(AnnotationDbi)
library(org.Dm.eg.db)
package.version("org.Dm.eg.db") # "3.14.0"

detach("package:clusterProfiler", unload = TRUE) #for `select` function namespace conflict
#https://bioconductor.org/packages/release/data/annotation/manuals/org.Dm.eg.db/man/org.Dm.eg.db.pdf
#Mappings were based on data provided by: Entrez Gene ftp://ftp.ncbi.nlm.nih.gov/gene/DATA With
#a date stamp from the source of: 2021-Sep13

columns(org.Dm.eg.db)
keytypes(org.Dm.eg.db)
test.inp=sample(keys(org.Dm.eg.db, keytype="SYMBOL"),6,replace = F)
test.out=select(org.Dm.eg.db, keys=test.inp, keytype="SYMBOL",
       #columns=c("SYMBOL","GENENAME",'FLYBASE','ALIAS','ACCNUM','ENTREZID') )
      columns=c("SYMBOL","GENENAME",'FLYBASE','ENTREZID') )
length(test.inp) #6 genes
dim(test.out) #6 x 2
test.out

################################################################################
## get chromosome location (prefer biomart as shown in the end of this script)
#https://davetang.org/muse/2013/12/16/bioconductor-annotation-packages/
ls("package:org.Dm.eg.db")
ls("package:org.Dm.eg.db")[grep('CHR',ls("package:org.Dm.eg.db"))]
#org.Dm.egCHRLOC is an R object that maps entrez gene identifiers to the starting position of the
#gene. The position of a gene is measured as the number of base pairs.
#The CHRLOCEND mapping is the same as the CHRLOC mapping except that it specifies the
#ending base of a gene instead of the start.
#Mappings were based on data provided by: UCSC Genome Bioinformatics (Drosophila melanogaster)
#With a date stamp from the source of: 2019-Jan15

## Bimap interface:
x <- org.Dm.egCHRLOC
dim(x) #29138
# Get the entrez gene identifiers that are mapped to chromosome locations
mapped_genes <- mappedkeys(x)
length(mapped_genes) #17129
# Convert to a list
xx <- as.list(x[test.out$ENTREZID])
length(xx) #6
x.start=xx
#x.end=xx;
x.start

###############################
## a used case: dataset S1 from
# https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.1000590#pbio.1000590.s002
# start with SYMBOL, then ALIAS, then GENENAME
# to get all mapping between input genes and FBgn IDs.
x=read.table('../Documents/single.cell_sex.differences/external_data/2011_paper_data/2011-embryo_normalized.read.count.data.txt',
             header=T,fill=T)
dim(x) #12353    79
x[1:3,1:3]
test.out=select(org.Dm.eg.db, keys=x$NAME, keytype="SYMBOL",
       columns=c("SYMBOL","GENENAME",'FLYBASE','GENENAME') )
dim(test.out) #12353     3
head(test.out)
sum(is.na(test.out$FLYBASE))#2850
sum(test.out$FLYBASE=='',na.rm=T) #0
test.out.na=test.out[is.na(test.out$FLYBASE),]
head(test.out.na)
tmp=select(org.Dm.eg.db, keys=test.out.na$SYMBOL, keytype=c("ALIAS"),
       columns=c("SYMBOL","GENENAME",'FLYBASE','GENENAME') )
sum(is.na(tmp$FLYBASE)) #63

test.out.na2=tmp[is.na(tmp$FLYBASE),]
head(test.out.na2)

AnnotationDbi::select(org.Dm.eg.db,keys='FBgn0052865',keytype='FLYBASE',
                      columns=c("SYMBOL","GENENAME",'FLYBASE','GENENAME') )

tmp=AnnotationDbi::select(org.Dm.eg.db,keys=test.out.na2$ALIAS,keytype='GENENAME',
       columns=c("SYMBOL","GENENAME",'FLYBASE','GENENAME') )
dim(test.out.na2);dim(tmp)
test.out.na3=tmp[is.na(tmp$FLYBASE),]
dim(test.out.na3) #48
head(test.out.na3)

#CG11231, http://flybase.org/reports/FBgn0039939
#This gene is referred to in FlyBase by the symbol Dmel\CG11231 (FBgn0039939).The gene previously had a gene model which is not supported by current evidence.


###############################
## use `bitr` from clusterProfiler
# BiocManager::install("clusterProfiler")
library(clusterProfiler)
query.genes = c('alphagamma-element:CR32865', 'AlstR')
bitr(query.genes, fromType="ALIAS", toType=c("ENTREZID", "FLYBASE"), OrgDb="org.Dm.eg.db")


###############################
## use gene coordinates
library(biomaRt)
packageVersion("biomaRt") #‘2.50.3’

biolist <- as.data.frame(listMarts())
ensembl=useMart("ensembl")
esemblist <- as.data.frame(listDatasets(ensembl)) #this shows the verion of reference used in biomart
#dataset                              description  version
#55 dmelanogaster_gene_ensembl Drosophila melanogaster genes (BDGP6.32) BDGP6.32


ensembl = useDataset("dmelanogaster_gene_ensembl",mart=ensembl)
filters = listFilters(ensembl)
attributes = listAttributes(ensembl)

attributes[grep('dmel',attributes$name),]

t2g<-getBM(attributes=c('ensembl_gene_id','chromosome_name','start_position','end_position'), mart = ensembl)
dim(t2g) #23932     4
head(t2g)
#saveRDS(t2g,'t2g_chr.coord.rds')
sum(t2g$ensembl_gene_id %in% test.out$FLYBASE)
t2g[t2g$ensembl_gene_id %in% test.out$FLYBASE, ]


