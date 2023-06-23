
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
# https://github.com/mingwhy/bioinfo_homemade_tools/blob/main/used.cases/biomaRt_usage/gene.id_converter_AnnotationDbi_bitr.R
df.gene=readRDS('~/Documents/Data_AgingFlyCellAtlas/AFCA_gene.meta.rds')
dim(df.gene) #15991     7
head(df.gene)
sum(df.gene$FLYBASE %in% df.gene.length$ensembl_gene_id) #13252

###############################################################
## read in data
sce=readH5AD('~/Documents/Data_AgingFlyCellAtlas/adata_headBody_S_v1.0.h5ad') # 22966 31001 
#sce=readH5AD('~/Documents/Data_AgingFlyCellAtlas/adata_head_S_v1.0.h5ad') # 22966 31001 
#sce=readH5AD('~/Documents/Data_AgingFlyCellAtlas/adata_Body_S_v1.0.h5ad') # 15992 276273
sce #SingleCellExperiment,  
assayNames(sce) #'X'
sce=sce[,sce$age %in% c(5,70)]
sce$age=droplevels(sce$age)

cell.meta=colData(sce)
colnames(cell.meta)
table(cell.meta$tissue)
#body   head 
#276273 289981 
table(cell.meta$sex,cell.meta$age)
head(cell.meta)
#table(cell.meta$afca_annotation_broad,cell.meta$afca_annotation)
length(table(cell.meta$afca_annotation)) # 163 cell types, body=73
length(table(cell.meta$afca_annotation_broad)) #17 broad cell classes
#grep('head',cell.meta$afca_annotation,value = T)

##################################################################
## undo log1p (so cells all have the same lib.size)
# and calculate gene mean expr for young and old

mat=assay(sce,'X')
mat[1:3,1:3]
test.mat=mat[,1:10]
cell.meta[1:10,]$total_counts
emat=test.mat; #use log1p directly
#emat<-(exp(test.mat)-1)
#colSums(emat) #all the same size, 1e4
#tmp=emat[,1]/10^4*cell.meta[1:10,]$total_counts[1]
#table(tmp)

#output.file='res.gene.mean_female_log1p_subsample.rds'
output.file='res.gene.mean_male_log1p_subsample.rds'

if(!file.exists(output.file)){
  # repeat this for each cell type separately: 5d vs 30d (contain most cells)
  all.tcs=as.character(unique(sce$afca_annotation))
  all.tcs
  res.gene.mean<-lapply(all.tcs,function(tc){
    condition= sce$afca_annotation==tc & sce$sex=='female';
    #condition= sce$afca_annotation==tc & sce$sex=='male';
    #condition= sce$afca_annotation==tc 
    
    mat=assay(sce[,condition],'X') #use log1p directly
    #mat=exp(assay(sce[,condition],'X'))-1 #use original raw count
    #lib.size=cell.meta[condition,]$total_counts
    #umi.mat=((exp(mat)-1)/10^4 ) * rep(lib.size, rep(nrow(mat),length(lib.size))) 
    
    sub.cell.meta=cell.meta[condition,]
    sub.cell.meta$age=droplevels(sub.cell.meta$age)
    if(nrow(sub.cell.meta)==0){return(NULL)}
    if(table(sub.cell.meta$age)[1]<20){return(NULL)} #age 5 contain less than 20 cells
    
    table(sub.cell.meta$age)
    ncell=min(table(sub.cell.meta$age))
    i.young = sample(rownames(sub.cell.meta)[sub.cell.meta$age=='5'],ncell,replace = F)
    i.old = sample(rownames(sub.cell.meta)[sub.cell.meta$age=='70'],ncell,replace = F)
    
    sub.cell.meta=sub.cell.meta[c(i.young,i.old),]
    mat=mat[,c(i.young,i.old)]
    table(sub.cell.meta$age)
    if(ncol(mat)!=nrow(sub.cell.meta)){cat('cell number does not match!!');break}
    
    ages=sort(unique(sub.cell.meta$age))
    gene.mean=sapply(ages,function(age.i){
      rowMeans(mat[,sub.cell.meta$age==age.i])
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

############################################################################
## make plots, combine with gene length information
rownames(df.gene)=df.gene$original.id
rownames(df.gene.length)=df.gene.length$ensembl_gene_id #FBgnxxx
gene.info=merge(df.gene,df.gene.length[,c(1,4,5,12,13)],by.x='FLYBASE',by.y='ensembl_gene_id')
dim(gene.info) # 13252    10

gene.info2=gene.info[-grep('Ribosomal protein',gene.info$GENENAME),]#remove ribosomal protein genes

#res.gene.mean=readRDS('res.gene.mean_female_log1p_subsample.rds');
res.gene.mean=readRDS('res.gene.mean_male_log1p_subsample.rds');

#pdf('check_5_70_gene.mean_female_log1p_subsample.pdf',useDingbats = T,height = 14,width = 16)
pdf('check_5_70_gene.mean_male_log1p_subsample.pdf',useDingbats = T,height = 14,width = 16)
par(mfrow=c(4,4))
for(i in 1:length(res.gene.mean)){
    test=res.gene.mean[[i]]
    colnames(test)
    if(sum(colnames(test) %in% c('5','70'))!=2){next}
    test$original.id=rownames(test) #gene id
    
    #tmp2=merge(test,gene.info)
    tmp2=merge(test,gene.info2)
    sum(is.na(tmp2$FLYBASE))
    
    #tmp2$gene.length=tmp2$transcript_end-tmp2$transcript_start
    tmp2$gene.length=tmp2$transcript_length;
    tmp.df=data.frame(original.id=tmp2$original.id,'young'=log(tmp2$`5`+1),'old'=log(tmp2$`70`+1),
                      gene.length=tmp2$gene.length)
    
    tmp.df=tmp.df[apply(tmp.df[,c('young','old')],1,sum)!=0,]
    #tmp.df=tmp.df[-grep('^mt:',tmp.df$original.id),] #remove mt genes
    #tmp.df=tmp.df[-grep('RNA:',tmp.df$original.id),] #remove snoRNA or lncRNA genes
    
    tmp.df=tmp.df[order(tmp.df$gene.length),]
    
    x=quantile(tmp.df$gene.length,c(0.25,0.75))
    short.genes=tmp.df[tmp.df$gene.length<x[1],]$original.id
    long.genes=tmp.df[tmp.df$gene.length>x[2],]$original.id
    
    tmp.df$length.group='others'
    tmp.df[tmp.df$original.id %in% short.genes,]$length.group='short'
    tmp.df[tmp.df$original.id %in% long.genes,]$length.group='long'
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

######################################################################
## are the same set of genes drive this pattern?

#res.gene.mean=readRDS('res.gene.mean_female_log1p.rds');
res.gene.mean=readRDS('res.gene.mean_male_log1p.rds');

all.out<-lapply(1:length(res.gene.mean),function(i){
  test=res.gene.mean[[i]]
  colnames(test)
  tc=test$cell.type[1]
  
  if(sum(colnames(test) %in% c('5','70'))!=2){next}
  test$original.id=rownames(test) #gene id
  tmp=merge(test,df.gene)
  sum(is.na(tmp$FLYBASE))
  tmp2=merge(tmp,df.gene.length,by.x='FLYBASE',by.y='ensembl_gene_id')
  
  #tmp2$gene.length=tmp2$transcript_end-tmp2$transcript_start
  tmp2$gene.length=tmp2$transcript_length;
  tmp.df=data.frame(original.id=tmp2$original.id,'young'=log(tmp2$`5`+1),'old'=log(tmp2$`70`+1),
                    gene.length=tmp2$gene.length)
  
  tmp.df=tmp.df[apply(tmp.df[,c('young','old')],1,sum)!=0,]
  #tmp.df=tmp.df[-grep('^mt:',tmp.df$original.id),] #remove mt genes
  #tmp.df=tmp.df[-grep('RNA:',tmp.df$original.id),] #remove snoRNA or lncRNA genes
  
  tmp.df=tmp.df[order(tmp.df$gene.length),]
  
  x=quantile(tmp.df$gene.length,c(0.25,0.75))
  cat(tc,x,'\n')
  short.genes=tmp.df[tmp.df$gene.length<x[1],]$original.id
  long.genes=tmp.df[tmp.df$gene.length>x[2],]$original.id
  
  tmp.df$length.group='others'
  tmp.df[tmp.df$original.id %in% short.genes,]$length.group='short'
  tmp.df[tmp.df$original.id %in% long.genes,]$length.group='long'
  
  tmp.df$cell_type=tc
  tmp.df$mean.diff=tmp.df$old-tmp.df$young
  tmp.df
})
# rank gene by age-related difference in mean gene expression
top200<-lapply(all.out,function(i){
  i=i[order(i$mean.diff),]
  i[1:200,]$original.id #over express in young
})

bottom200<-lapply(all.out,function(i){
  i=i[order(i$mean.diff,decreasing = T),]
  i[1:200,]$original.id #over express in old
})

x1=sort(table(unlist(top200))) #over express in young
x2=sort(table(unlist(bottom200))) #over express in old

Rps=gene.info[grep('Ribosomal protein',gene.info$GENENAME),] %>% arrange(transcript_length) 
nrow(Rps) #94 ribosomal proteins in fly
over.in.young<-Rps[Rps$original.id %in% names(x1),] #83
over.in.old<-Rps[Rps$original.id %in% names(x2),] #8

y=gene.info[gene.info$original.id %in% names(tail(x1,60)),] %>% arrange(transcript_length)
y=gene.info[gene.info$original.id %in% names(tail(x2,60)),] %>% arrange(transcript_length)


