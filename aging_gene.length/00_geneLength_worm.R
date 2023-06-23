
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
  esemblist[grep('elegans',esemblist$dataset),]
  #dataset                              description  version
  #29 celegans_gene_ensembl Caenorhabditis elegans (PRJNA13758) genes (WBcel235) WBcel235
  
  ensembl = useDataset("celegans_gene_ensembl",mart=ensembl)
  filters = listFilters(ensembl)
  attributes = listAttributes(ensembl)
  
  attributes[grep('celegans',attributes$name),]
  attributes[grep('symbol',attributes$name),]
  
  attributes[grep('wormbase_gseqname|external_gene_name',attributes$name),]
  grep('transcript',attributes$name,value=TRUE,ignore.case = TRUE)
  
  t2g<-getBM(attributes=c('wormbase_gseqname','external_gene_name','ensembl_gene_id',"gene_biotype",
                          "ensembl_transcript_id","transcript_start","transcript_end",
                          #"ensembl_exon_id","exon_chrom_start","exon_chrom_end",
                          "chromosome_name","strand","transcript_biotype", 
                          'start_position','end_position',
                          "transcript_length","transcript_count"), mart = ensembl)
  dim(t2g) #60000    14
  head(t2g)
  saveRDS(t2g,'celegans_t2g_chr.coord.rds')
}
t2g=readRDS('celegans_t2g_chr.coord.rds')
# one gene has multiple rows as muliple transcripts
df.gene.length=t2g[!duplicated(t2g$ensembl_gene_id),]
table(df.gene.length$gene_biotype)

#df.gene.length=df.gene.length[df.gene.length$gene_biotype == 'protein_coding',]
#dim(df.gene.length) #21959    13

table(df.gene.length$chromosome_name)
df.gene.length=df.gene.length[df.gene.length$chromosome_name !='MtDNA',]
dim(df.gene.length) #46890    14

###############################################################
## read in data
library(zellkonverter)
library(SingleCellExperiment)
library(scater);library(scran)
library(ggplot2);library(gridExtra)
library(tidyverse);
library(SummarizedExperiment)

if(!file.exists('sce_worm.h5ad')){
library(reticulate)
#Sys.setenv(RETICULATE_PYTHON = "/Users/mingyang/anaconda3/envs/RCFGL/bin/python")
#py_config()
sc=import('scanpy')  #py_install('scanpy')
#ad=sc$read_h5ad('/Users/mingyang/Documents/Data_worm_aing/ad_worm_aging.h5ad') #https://scanpy.readthedocs.io/en/stable/generated/scanpy.read_h5ad.html
ad=sc$read_h5ad('/Users/mingyang/Documents/Data_worm_aing/ad_worm_aging_umi.h5ad') #https://scanpy.readthedocs.io/en/stable/generated/scanpy.read_h5ad.html
cell.meta=ad$obs
gene.meta=ad$var

ad$X$shape # 47423 20305
class(ad$X)
mat=ad$X$toarray()
class(mat) #"matrix" "array" 
dim(mat) #  47423 20305 
dim(cell.meta) #47423

# create a SingleCellExperiment
#https://www.bioconductor.org/packages/devel/bioc/vignettes/SingleCellExperiment/inst/doc/intro.html
sce <- SingleCellExperiment( list(counts=as(t(mat), "sparseMatrix")), #list(counts=mat)),
                             colData=cell.meta,
                             rowData=gene.meta);
sce
writeH5AD(sce, 'sce_worm.h5ad')
}

###########################################
## filter cell types

if(!file.exists('sce_worm_select.tc.h5ad')){
  inp.sce<-readH5AD('sce_worm.h5ad') 
  cell.meta=as.data.frame(colData(inp.sce))
  x=table(cell.meta$timepoint,cell.meta$annotate_name)
  #table(cell.meta$annotate_name)
  x[,1:3]
  # combine d1+d3, d8+d11
  
  inp.sce2=inp.sce[,cell.meta$timepoint %in% c('d1','d3','d11','d15')]
  cell.meta=as.data.frame(colData(inp.sce2))
  table(cell.meta$timepoint)
  cell.meta$age='young'
  cell.meta[cell.meta$timepoint %in% c('d11','d15'),]$age='old'
  x=table(cell.meta$age,cell.meta$annotate_name)
  pick.tc=names(which(apply(x,2,function(x) sum(x>30))==2))
  
  inp.sce2=inp.sce2[,cell.meta$annotate_name %in% pick.tc]
  table(inp.sce2$timepoint)
  age=rep('young',ncol(inp.sce2))
  age[inp.sce2$timepoint %in% c('d11','d15')]='old'
  inp.sce2$age=age;
  writeH5AD(inp.sce2, 'sce_worm_select.tc.h5ad')
}

sce=readH5AD('sce_worm_select.tc.h5ad') # 22966 31001 
sce # 20305 15004
assayNames(sce) #'counts'

cell.meta=as.data.frame(colData(sce))
colnames(cell.meta)

gene.meta=as.data.frame(rowData(sce))
sum(gene.meta$gene_ids %in% df.gene.length$wormbase_gseqname) #20248

##################################################################
## undo log1p (so cells all have the same lib.size)
# and calculate gene mean expr for young and old
table(sce$age)
sce<-logNormCounts(sce)
assayNames(sce) #"counts"    "logcounts"

sce=sce[rowData(sce)[['gene_class']]=='protein_coding',]
sce_naive=sce[ rowData(sce)[['gene_ids']] %in% df.gene.length$wormbase_gseqname,]
sce_naive #  17912 16508 
table(as.character(sce_naive$annotate_name),as.character(sce_naive$age))
cell.meta=as.data.frame(colData(sce_naive))
table(as.character(cell.meta$annotate_name),as.character(cell.meta$age))

output.file='log1p_cell.types.rds'
#output.file='log1p_cell.types_subsample.rds'

if(!file.exists(output.file)){
  # repeat this for each cell type separately: 5d vs 30d (contain most cells)
  all.tcs=as.character(unique(sce$annotate_name))
  all.tcs
  res.gene.mean<-lapply(all.tcs,function(tc){
    condition= (sce_naive$annotate_name==tc)
    #mat=exp(assay(sce_naive[,condition],'logcounts'))-1
    mat=assay(sce_naive[,condition],'logcounts')
    
    sub.cell.meta=cell.meta[condition,]
    sub.cell.meta$age=droplevels(sub.cell.meta$age)
    if(nrow(sub.cell.meta)==0){return(NULL)}
    if(table(sub.cell.meta$age)[1]<20){return(NULL)} #age 5 contain less than 20 cells
    
    #ncell=min(table(sub.cell.meta$age))
    #i.young = sample(rownames(sub.cell.meta)[sub.cell.meta$age=='young'],ncell,replace = F)
    #i.old = sample(rownames(sub.cell.meta)[sub.cell.meta$age=='old'],ncell,replace = F)
    #sub.cell.meta=sub.cell.meta[c(i.young,i.old),]
    #mat=mat[,c(i.young,i.old)]
    #table(sub.cell.meta$age)
    
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

gene.info=df.gene.length[match(rowData(sce_naive)[['gene_ids']],df.gene.length$wormbase_gseqname),]
sum(gene.info$wormbase_gseqname==rowData(sce_naive)[['gene_ids']])
dim(sce_naive)

res.gene.mean=readRDS('log1p_cell.types.rds');
#res.gene.mean=readRDS('log1p_cell.types_subsample.rds');

length(res.gene.mean)

pdf('check_gene.mean_log1p_cell.types_loglog.pdf',useDingbats = T,height = 14,width = 16)
#pdf('check_gene.mean_log1p_cell.types_loglog_subsample.pdf',useDingbats = T,height = 14,width = 16)
par(mfrow=c(5,5))
for(i in 1:length(res.gene.mean)){
    test=res.gene.mean[[i]]
    colnames(test)
    
    dim(test)
    test$transcript_length=gene.info$transcript_length
    test$gene.length=gene.info$transcript_end-gene.info$transcript_start
    
    tmp.df=data.frame(original.id=gene.info$wormbase_gseqname,
                      'young'=log(test$young+1),'old'=log(test$old+1),
                      #'young'=log(test$`3m`+1,base=2),'old'=log(test$`24m`+1,base=2),
                      #'young'=test$young,'old'=test$old,
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
         xlab='log(Mean expression in young + 1)',ylab='log(Mean expression in old + 1)')
         #xlab='Mean expression in young mice (3m)',ylab='Mean expression in old mice (24m)')
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


