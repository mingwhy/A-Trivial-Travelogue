
###############################################################
## read in data: https://flycellatlas.org/
library(zellkonverter) #https://bioconductor.org/packages/release/bioc/html/zellkonverter.html
sce=readH5AD('fat_body.h5ad') 
sce #a `SingleCellExperiment` object

library(SingleCellExperiment) #https://bioconductor.org/packages/release/bioc/html/SingleCellExperiment.html
library(SummarizedExperiment)
expr.mat=assay(sce)
class(expr.mat)
expr.mat[1:3,1:3] #gene by cell matrix

cell.meta=colData(sce) #get cell annotation data.frame
head(cell.meta)
table(cell.meta$sex,cell.meta$transf_annotation)

genes=rownames(sce)
head(genes)

## use AnnotationDbi and org.Dm.eg.db to get fly gene IDs
library(AnnotationDbi)
library(org.Dm.eg.db)
package.version("org.Dm.eg.db") # "3.14.0"

columns(org.Dm.eg.db)
keytypes(org.Dm.eg.db)
test.out=AnnotationDbi::select(org.Dm.eg.db, keys=genes, keytype="SYMBOL",
                              columns=c("SYMBOL","GENENAME",'FLYBASE','FLYBASECG','ENTREZID') )
head(test.out)
sum(is.na(test.out$FLYBASE))
test.out[is.na(test.out$FLYBASE),]

###############################################################
library(scran) #https://bioconductor.org/packages/release/bioc/html/scran.html
library(scuttle) #https://bioconductor.org/packages/release/bioc/html/scuttle.html

assayNames(sce) #'X'
assayNames(sce)<-'counts' #required by `scuttle` r package

# normalize cell-wise library size (same number of transcripts across cells)
sce<-logNormCounts(sce, log=TRUE, pseudo.count=1) #if log, the cell.lib.size range 2~3 orders
assayNames(sce) #"counts"    "logcounts"
expr.mat=assay(sce,'logcounts')
expr.mat[1:3,1:3]

genes[grep('lncRNA:roX',genes,ignore.case = T)]
tmp=expr.mat[grep('lncRNA:roX',genes,ignore.case = T),]
tmp[1:2,1:3]

library(ggplot2)
df.tmp=as.data.frame(t(tmp))
head(df.tmp)
df.tmp$sex=cell.meta$sex
df.tmp$cell.annotation=cell.meta$transf_annotation
head(df.tmp)
ggplot(df.tmp,aes(x=sex,y=`lncRNA:roX2`,col=sex))+
  facet_wrap(.~cell.annotation,scale='free')+
  geom_violin()+geom_jitter()+theme_classic()
