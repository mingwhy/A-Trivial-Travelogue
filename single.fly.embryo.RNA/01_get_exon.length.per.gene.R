#http://bioconductor.org/packages/2.10/data/annotation/html/BSgenome.Hsapiens.UCSC.hg19.html
#BiocManager::install("GenomicFeatures")
#https://www.biostars.org/p/83901/

# First, import the GTF-file that you have also used as input for htseq-count
library(GenomicFeatures)
txdb <- makeTxDbFromGFF("Drosophila_melanogaster.BDGP6.32.104.chr.gtf.gz",format="gtf")
txdb

# then collect the exons per gene id
exons.list.per.gene <- exonsBy(txdb,by="gene")
length(exons.list.per.gene) #23882
exons.list.per.gene[[1]]
# then for each gene, reduce all the exons to a set of non overlapping exons, calculate their lengths (widths) and sum then
exonic.gene.sizes <- sum(width(reduce(exons.list.per.gene)))
length(exonic.gene.sizes) #23882
df.exonic.gene.sizes=as.data.frame(exonic.gene.sizes)
head(df.exonic.gene.sizes)
df.exonic.gene.sizes$gene=rownames(df.exonic.gene.sizes)
df.exonic.gene.sizes=df.exonic.gene.sizes[,c(2,1)]

## add chromosome location info
# get chromosome location
seqlevels(txdb)
columns(txdb)
keytypes(txdb)
cols <- c("TXNAME", "TXSTRAND", "TXCHROM")
chr.info=select(txdb, keys = df.exonic.gene.sizes$gene, 
                columns=cols,
                keytype="GENEID")
head(chr.info)
dim(chr.info) # 41149 
chr.info=chr.info[!duplicated(chr.info$GENEID),]
dim(chr.info) # 23882 
sum(chr.info$GENEID==df.exonic.gene.sizes$gene) #23882
df.info=merge(chr.info,df.exonic.gene.sizes,by.x='GENEID',by.y='gene')
dim(df.info)

data.table::fwrite(df.info,"dme7_geneLength_chr.txt")

#################################################
df.info=data.table::fread('dme7_geneLength_chr.txt')
query.genes=scan('all.gene.names.txt',what='')
head(query.genes)
length(query.genes) #17481
sum(df.info$GENEID %in% query.genes) #17382

