
library(AnnotationDbi)
library(org.Dm.eg.db)
package.version("org.Dm.eg.db") # "3.14.0"
#https://github.com/mingwhy/bioinfo_homemade_tools/blob/main/package_101/biomaRt_101/gene.id_converter_AnnotationDbi_bitr.R
columns(org.Dm.eg.db)
keytypes(org.Dm.eg.db)

###################################################
## read in STAR output tab-limited read count data
files=Sys.glob('GSE127176_RAW/*tab')
files
raw.expr=list();
for(file in files){
  gsm=strsplit(basename(file),'_')[[1]][1]
  x=read.table(file,skip=4)
  raw.expr[[gsm]]=x
}
length(raw.expr) #54

# extract 'expressed genes', defined by >=21 counts in >=5 embryos
per.embryo.gene=lapply(raw.expr,function(x){
  x[x[,2]>=21,1]
})
sapply(per.embryo.gene,length)
expr.genes=names(which(table(unlist(per.embryo.gene))>=5))
length(expr.genes) #8983 genes
filtered.expr=lapply(raw.expr,function(x){
  x[x[,1] %in% expr.genes,]
})
sapply(filtered.expr,nrow) #8983 genes

query.genes=expr.genes
##########################################################
## begin id mapping
# round 1
test.out=select(org.Dm.eg.db, keys=query.genes, keytype="FLYBASE",
                columns=c("SYMBOL","GENENAME",'FLYBASE') )
dim(test.out) #8983     3
head(test.out)
sum(is.na(test.out$FLYBASE)) #0
length(unique(test.out$SYMBOL)) #8983
test.out$query=query.genes;

keep.colnames=c("query","SYMBOL","GENENAME" ,"FLYBASE");
gene.meta=test.out[,keep.colnames]

##########################################################
## get chr location
t2g=readRDS('../t2g_chr.coord.rds')
dim(t2g) #23932 
head(t2g)
sum(gene.meta$FLYBASE %in% t2g$ensembl_gene_id) #8938 overlap
df.gene.meta=merge(gene.meta,t2g,by.x='FLYBASE',by.y='ensembl_gene_id',all.x=T)
dim(df.gene.meta) #8983
head(df.gene.meta)

data.table::fwrite(df.gene.meta,'validate.id_2019_paper_data.txt')
