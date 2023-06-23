
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
                        "transcript_length","transcript_count"), mart = ensembl)
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

###############################################################
## read in data
library(zellkonverter)
library(SingleCellExperiment)
library(scater);library(scran)
library(ggplot2);library(gridExtra)
library(tidyverse);
library(SummarizedExperiment)
library(ggpointdensity)
library(viridis)
library(patchwork); #for plot_annotation

if(!file.exists('~/Documents/aging_cell.turnover/1120_TMS_male_analysis/select.tc.h5ad')){
  inp.sce<-readH5AD('~/Documents/Data_mouse_aging_atlas/TMS.gene.data_final/tabula-muris-senis-facs-official-raw-obj.h5ad');       
  colData(inp.sce)$tissue_cell.type=paste(inp.sce$tissue,inp.sce$cell_ontology_class,sep=':')
  
  sub.sce=inp.sce[,inp.sce$sex=='male']
  x=table(sub.sce$tissue_cell.type,sub.sce$age)
  x=as.data.frame(x[,c(1,4)]) #3m and 24m
  colnames(x)=c('tc','age','ncell')
  x=x %>% spread(age,ncell)
  dim(x) #202
  data.table::fwrite(x,'ncell_per.age_per.tc_raw.txt',sep='\t',quote=F)
  
  y=x[x[,2]>=50 & x[,3]>=50, ] #both age groups contain >=50 cells
  dim(y) #72 cell states
  sum(c(y$`3m`,y$`24m`)) #47898 cells in total
  
  y[grep('Brain',y$tc),]
  y$tc=as.character(y$tc)
  unique(unlist(lapply(strsplit(y$tc,'\\:'),'[',1))) #22 unique tissues
  unique(unlist(lapply(strsplit(y$tc,'\\:'),'[',2))) #52 unique cell types
  data.table::fwrite(y,'ncell_per.age_per.tc.txt',sep='\t',quote=F)
  
  
  # #only keep tc enough cells
  #sum(y$tc %in% tc.orders) 
  #pick.cell.types=y[y$tc %in% tc.orders,]$tc 
  pick.cell.types=y$tc
  sce=sub.sce[,sub.sce$tissue_cell.type %in% pick.cell.types]
  
  unique(sce$age) #'3m','18m','21m','24m'
  sce$age=droplevels(sce$age)  
  unique(sce$age) 
  sce$age=factor(sce$age,levels=c('3m','18m','24m'))
  cell.meta=colData(sce)
  #rowData(sce)
  length(table(cell.meta$tissue_cell.type)) 
  writeH5AD(sce, 'select.tc.h5ad')
}


sce=readH5AD('~/Documents/aging_cell.turnover/1120_TMS_male_analysis/select.tc.h5ad') # 22966 31001 
sce #SingleCellExperiment,  
assayNames(sce) #'X'

cell.meta=colData(sce)
colnames(cell.meta)
table(cell.meta$tissue)

table(cell.meta$sex,cell.meta$age)
head(cell.meta)

genes=rownames(sce)
sum(genes %in% df.gene.length$mgi_symbol) #17996

##################################################################
## undo log1p (so cells all have the same lib.size)
# and calculate gene mean expr for young and old
sce=sce[,sce$age %in% c('3m','24m')]
unique(sce$age) #3m 18m 24m
assayNames(sce)<-'counts'
sce<-logNormCounts(sce)
assayNames(sce) #"counts"    "logcounts"

sce_naive=sce[rownames(sce) %in% df.gene.length$mgi_symbol,]
sce_naive # 17996 47898 
table(sce_naive$tissue)
cell.meta=colData(sce_naive)

output.file='male_log1p_organs.rds'
if(!file.exists(output.file)){
  # repeat this for each cell type separately: 5d vs 30d (contain most cells)
  all.tcs=as.character(unique(sce$tissue))
  all.tcs
  res.gene.mean<-lapply(all.tcs,function(tc){
    condition= (sce_naive$tissue==tc)
    
    sub.cell.meta=cell.meta[condition,]
    sub.cell.meta$age=droplevels(sub.cell.meta$age)
    if(nrow(sub.cell.meta)==0){return(NULL)}
    if(table(sub.cell.meta$age)[1]<20){return(NULL)} #age 5 contain less than 20 cells
    
    #mat=exp(assay(sce_naive[,condition],'logcounts'))-1
    mat=assay(sce_naive[,condition],'logcounts')
    
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
genes=rownames(sce_naive)
gene.info=df.gene.length[match(genes,df.gene.length$mgi_symbol),]
sum(gene.info$mgi_symbol==genes);length(genes)

res.gene.mean=readRDS('male_log1p_organs.rds');
length(res.gene.mean)

pdf('check_3m_24m_gene.mean_male_log1p_organs.pdf',useDingbats = T,height = 14,width = 16)
par(mfrow=c(5,5))
for(i in 1:length(res.gene.mean)){
    test=res.gene.mean[[i]]
    colnames(test)
    
    dim(test)
    sum(rownames(test)==gene.info$mgi_symbol)
    
    test$transcript_length=gene.info$transcript_length
    test$gene.length=gene.info$transcript_end-gene.info$transcript_start
    
    tmp.df=data.frame(original.id=gene.info$mgi_symbol,
                      #'young'=log(test$`3m`+1),'old'=log(test$`24m`+1),
                      #'young'=log(test$`3m`+1,base=2),'old'=log(test$`24m`+1,base=2),
                      'young'=test$`3m`,'old'=test$`24m`,
                      gene.length=test$gene.length)
    
    tmp.df=tmp.df[apply(tmp.df[,c('young','old')],1,sum)!=0,]
    #tmp.df=tmp.df[!tmp.df$young==0 | tmp.df$old==0,]
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
         #xlab='log(Mean expression in young mice (3m)+1)',ylab='log(Mean expression in old mice (24m)+1)')
         xlab='Mean expression in young mice (3m)',ylab='Mean expression in old mice (24m)')
    abline(a=0,b=1,lwd=1,col='black')
    
  }
dev.off()


######################################################################
######################################################################
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


