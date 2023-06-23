
library(zellkonverter)
library(SingleCellExperiment)
library(scater);library(scran)
library(ggplot2);library(gridExtra);
library(tidyverse);
library(zellkonverter)
library(SummarizedExperiment)
library(viridis)
library(ggpubr)

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
esemblist[grep('melanogaster',esemblist$description),]
#dataset                              description  version
#56 dmelanogaster_gene_ensembl Drosophila melanogaster genes (BDGP6.32) BDGP6.32

ensembl = useDataset("dmelanogaster_gene_ensembl",mart=ensembl)
filters = listFilters(ensembl)
attributes = listAttributes(ensembl)

attributes[grep('dmel',attributes$name),]
grep('transcript',attributes$name,value=TRUE,ignore.case = TRUE)

t2g<-getBM(attributes=c('ensembl_gene_id',"gene_biotype",
                        "ensembl_transcript_id","transcript_start","transcript_end",
                        #"ensembl_exon_id","exon_chrom_start","exon_chrom_end",
                        "chromosome_name","strand","transcript_biotype", 
                        'start_position','end_position',
                        'flybase_transcript_id',"transcript_length","transcript_count"), mart = ensembl)
dim(t2g) #41209    13
head(t2g)
saveRDS(t2g,'t2g_chr.coord.rds')
}
t2g=readRDS('t2g_chr.coord.rds')
# one gene has multiple rows as muliple transcripts
df.gene.length=t2g[!duplicated(t2g$ensembl_gene_id),]
table(df.gene.length$gene_biotype)
#miRNA                ncRNA            pre_miRNA       protein_coding           pseudogene 
#485                 2539                  262                13968                  342 
#rRNA               snoRNA                snRNA transposable_element                 tRNA 
#115                  299                   32                 5578                  312 
df.gene.length=df.gene.length[df.gene.length$gene_biotype == 'protein_coding',]
dim(df.gene.length) #13968    13

table(df.gene.length$chromosome_name)
df.gene.length=df.gene.length[df.gene.length$chromosome_name %in% c('2L','2R','3L','3R','4','X','Y'),]
dim(df.gene.length) #13950    13

##################################################################
## use AnnotationDbi and org.Dm.eg.db to get gene symbols 
df.gene=data.table::fread('~/Documents/Data_fly_FCA/fly.brain.atlas/gene.meta_brain.txt')
dim(df.gene) #12094    11
head(df.gene)
sum(df.gene$validated_id %in% df.gene.length$ensembl_gene_id) #10631


###############################################################
## read in data
sce=readRDS('~/Documents/Data_fly_FCA/fly.brain.atlas/wholebrain_filtered_valid.rds')
sce #Seurat

sum(rownames(sce) %in% df.gene$SYMBOL) #12094

cell.meta=sce@meta.data
sce <- as.SingleCellExperiment(sce)
table(cell.meta$Age)
table(cell.meta$annotation)

cell.meta$age.group='young';
cell.meta[cell.meta$Age>15,]$age.group='old'
tc=names(which(apply(table(cell.meta$age.group,cell.meta$annotation),2,function(i)sum(i>=50))==2))

cell.meta=cell.meta[cell.meta$annotation %in% tc,]
sce=sce[,sce$annotation %in% tc]

##################################################################
## undo log1p (so cells all have the same lib.size)
# and calculate gene mean expr for young and old
assayNames(sce)
mat=assay(sce,'counts')
mat[1:3,1:3]

output.file='Davie_flyBrain_res.gene.mean_log1p.rds'
#output.file='res.gene.mean_male_log1p.rds'
if(!file.exists(output.file)){
  # repeat this for each cell type separately: 5d vs 30d (contain most cells)
  all.tcs=as.character(unique(sce$annotation))
  all.tcs
  res.gene.mean<-lapply(all.tcs,function(tc){
    condition= sce$annotation==tc
    
    sub.cell.meta=cell.meta[condition,]
    
    if(nrow(sub.cell.meta)==0){return(NULL)}
    if(table(sub.cell.meta$age)[1]<20){return(NULL)} #age 5 contain less than 20 cells
    
    mat=assay(sce[,condition],'logcounts')
    
    ages=c('young','old')
    gene.mean=sapply(ages,function(age.i){
      rowMeans(mat[,sub.cell.meta$age.group==age.i])
    })
    colnames(gene.mean)=ages
    gene.mean=as.data.frame(gene.mean)
    gene.mean$cell.type=tc
    cat('cell type',tc,'is done\n')
    gene.mean
  })
  length(res.gene.mean) #163 tc
  res.gene.mean=Filter(Negate(is.null), res.gene.mean)
  saveRDS(res.gene.mean,output.file)
}

res.gene.mean=readRDS(output.file)
############################################################################
## make plots, combine with gene length information
gene.info=merge(df.gene[,c('SYMBOL','FBID_KEY')],df.gene.length,by.x='FBID_KEY',by.y='ensembl_gene_id')
dim(gene.info)

pdf('check_Davie_flyBrain_gene.mean_log1p.pdf',useDingbats = T,height = 14,width = 16)
par(mfrow=c(4,4))
for(i in 1:length(res.gene.mean)){
    test=res.gene.mean[[i]]
    colnames(test)
    
    test$SYMBOL=rownames(test) #gene id
    
    tmp2=merge(test,gene.info)
    
    
    #tmp2$gene.length=tmp2$transcript_end-tmp2$transcript_start
    tmp2$gene.length=tmp2$transcript_length;
    tmp.df=data.frame(SYMBOL=tmp2$SYMBOL,'young'=log(tmp2$young+1),'old'=log(tmp2$old+1),
                      gene.length=tmp2$gene.length)
    
    tmp.df=tmp.df[apply(tmp.df[,c('young','old')],1,sum)!=0,]
    #tmp.df=tmp.df[-grep('^mt:',tmp.df$SYMBOL),] #remove mt genes
    #tmp.df=tmp.df[-grep('RNA:',tmp.df$SYMBOL),] #remove snoRNA or lncRNA genes
    
    tmp.df=tmp.df[order(tmp.df$gene.length),]
    
    x=quantile(tmp.df$gene.length,c(0.25,0.75))
    short.genes=tmp.df[tmp.df$gene.length<x[1],]$SYMBOL
    long.genes=tmp.df[tmp.df$gene.length>x[2],]$SYMBOL
    
    tmp.df$length.group='others'
    tmp.df[tmp.df$SYMBOL %in% short.genes,]$length.group='short'
    tmp.df[tmp.df$SYMBOL %in% long.genes,]$length.group='long'
    #table(tmp.df$length.group)
    
    mycol=c('orange','blue','grey')
    tmp.df$col=mycol[factor(tmp.df$length.group,levels=c('short','long','others'))]
    tail(tmp.df)
    
    x=lm(tmp.df$old~tmp.df$young)
    x0=summary(x)
    x0
    #median(test$`5`);median(test$`30`)
    info=paste0(test$cell.type[1],'\nR^2=',round(x0$r.squared,2),',beta=',round(x0$coefficients[,1][2],2))
    
    x1=summary(lm(old~young,data=subset(tmp.df,length.group=='long') ))
    x2=summary(lm(old~young,data=subset(tmp.df,length.group=='short') ))
    info=paste0(info,'\n long.beta=',round(x1$coefficients[,1][2],2),
                '\n short.beta=',round(x2$coefficients[,1][2],2))
    plot(tmp.df$young,tmp.df$old,col=tmp.df$col,pch=16,cex=0.3,main=info,
         xlab='Mean expression in young fly (5d)',ylab='Mean expression in old fly (70d)')
    abline(a=0,b=1,lwd=1,col='black')
    
  }
dev.off()
#pick.top.bottom.200<-tmp.df[c(1:200,(nrow(tmp.df)-199):nrow(tmp.df)),]



