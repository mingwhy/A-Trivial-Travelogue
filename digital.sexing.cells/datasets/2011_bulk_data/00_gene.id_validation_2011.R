
library(AnnotationDbi)
library(org.Dm.eg.db)
package.version("org.Dm.eg.db") # "3.14.0"
#https://github.com/mingwhy/bioinfo_homemade_tools/blob/main/package_101/biomaRt_101/gene.id_converter_AnnotationDbi_bitr.R
columns(org.Dm.eg.db)
keytypes(org.Dm.eg.db)

###############################
## a used case: dataset S1 from
# https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.1000590#pbio.1000590.s002
# start with SYMBOL, then ALIAS, then GENENAME
# to get all mapping between input genes and FBgn IDs.
dat=read.table('./2011_paper_data/2011-embryo_normalized.read.count.data.txt',
             header=T,fill=T)
dim(dat) #12353    79
dat[1:3,1:3]
# round 1
test.out=select(org.Dm.eg.db, keys=dat$NAME, keytype="SYMBOL",
                columns=c("SYMBOL","GENENAME",'FLYBASE') )
dim(test.out) #12353     3
head(test.out)
sum(is.na(test.out$FLYBASE))#2850
keep1=test.out[!is.na(test.out$FLYBASE),] #9503 genes
test.out.na1=test.out[is.na(test.out$FLYBASE),]# 2850
keep1$query=keep1$SYMBOL;
length(unique(keep1$query)) #9503 genes

# round 2
test.out2=select(org.Dm.eg.db, keys=test.out.na1$SYMBOL, keytype=c("ALIAS"),
           columns=c("SYMBOL","GENENAME",'FLYBASE') )
dim(test.out2) # 3037 , more than 2850
test.out.na2=test.out2[is.na(test.out2$FLYBASE),]# 63
length(unique(test.out.na2$ALIAS)) #63

test.out2=test.out2[!is.na(test.out2$FLYBASE),];
length(unique(test.out2$ALIAS)) #2787
keep2.1=test.out2[!test.out2$ALIAS %in% test.out2[duplicated(test.out2$ALIAS),]$ALIAS,] #2626
tmp=test.out2[test.out2$ALIAS %in% test.out2[duplicated(test.out2$ALIAS),]$ALIAS,] #348
head(tmp); #for ALIAS with multiple FLYBASE, choose FBgn with a bigger number
tmp$FBgnnumber=gsub('FBgn','',tmp$FLYBASE)
tmp=tmp[with(tmp,order(ALIAS,FBgnnumber)),]
keep2.2=tmp[!duplicated(tmp$ALIAS),]
dim(keep2.2) #161
keep2=rbind(keep2.1,keep2.2[,1:4]) #2787
dim(test.out.na2) #63
keep2$query=keep2$ALIAS;
length(unique(keep2$query)) #2787 genes

# round 3
test.out3=AnnotationDbi::select(org.Dm.eg.db,keys=test.out.na2$ALIAS,keytype='GENENAME',
                          columns=c("SYMBOL","GENENAME",'FLYBASE') )
dim(test.out3) #174
test.out.na3=test.out3[is.na(test.out3$FLYBASE),]
dim(test.out.na3) #48
test.out3=test.out3[!is.na(test.out3$FLYBASE),]
dim(test.out3) #15

keep3=test.out3;
keep3$query=keep3$GENENAME;
length(unique(keep3$query)) #15 genes

keep.colnames=c("query","SYMBOL","GENENAME" ,"FLYBASE");
keep=rbind(keep1[,keep.colnames],keep2[,keep.colnames],keep3[,keep.colnames]) #remove ALIAS column
dim(keep) #12305
dim(test.out.na3) #48
dim(dat) #12353

## remove duplicated FBgnxxx.
sum(duplicated(keep$FLYBASE)) #208
high.confidence=keep[!keep$FLYBASE %in% names(which(table(keep$FLYBASE)>1)),]
dim(high.confidence) #11907 genes

tmp=keep[keep$FLYBASE %in% names(which(table(keep$FLYBASE)>1)),]
tmp=tmp[order(tmp$FLYBASE),]
sum(tmp$query==tmp$SYMBOL) #96
save1=tmp[tmp$query==tmp$SYMBOL,]
length(unique(save1$query)) #96
tmp2=tmp[!tmp$query %in% save1$query,]
length(unique(tmp2$query)) #302 genes, throw away

final.keep=rbind(high.confidence,save1)
dim(final.keep) #12003
length(unique(final.keep$FLYBASE))
length(unique(final.keep$query))

gene.meta=final.keep; #12003 genes
##########################################################
## get chr location
t2g=readRDS('../t2g_chr.coord.rds')
dim(t2g) #23932 
head(t2g)
sum(gene.meta$FLYBASE %in% t2g$ensembl_gene_id) #11996 overlap
df.gene.meta=merge(gene.meta,t2g,by.x='FLYBASE',by.y='ensembl_gene_id',all.x=T)
dim(df.gene.meta) #12003 genes
head(df.gene.meta)


data.table::fwrite(df.gene.meta,'validate.id_2011_paper_data.txt')
