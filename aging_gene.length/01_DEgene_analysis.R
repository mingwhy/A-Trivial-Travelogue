
library(zellkonverter)
library(SingleCellExperiment)
library(scater);library(scran)
library(ggplot2);library(gridExtra);
library(tidyverse);
library(zellkonverter)
library(SummarizedExperiment)
library(viridis)
library(ggpubr)
library(MAST)

###############################################################
## use biomaRt to get gene chr info 
# https://github.com/mingwhy/bioinfo_homemade_tools/tree/main/used.cases/biomaRt_usage
t2g=readRDS('t2g_chr.coord.rds')
# one gene has multiple rows as muliple transcripts
df.gene.length=t2g[!duplicated(t2g$ensembl_gene_id),]
table(df.gene.length$gene_biotype)
df.gene.length=df.gene.length[df.gene.length$gene_biotype == 'protein_coding',]
dim(df.gene.length) #13968    13

###############################################################
## use AnnotationDbi and org.Dm.eg.db to get gene symbols 
# https://github.com/mingwhy/bioinfo_homemade_tools/blob/main/used.cases/biomaRt_usage/gene.id_converter_AnnotationDbi_bitr.R
df.gene=readRDS('~/Documents/Data_AgingFlyCellAtlas/AFCA_gene.meta.rds')
dim(df.gene) #15991     7
head(df.gene)
sum(df.gene$FLYBASE %in% df.gene.length$ensembl_gene_id) #15991

###############################################################
## read in data
sce=readH5AD('~/Documents/Data_AgingFlyCellAtlas/adata_headBody_S_v1.0.h5ad') # 22966 31001 
#sce=readH5AD('~/Documents/Data_AgingFlyCellAtlas/adata_head_S_v1.0.h5ad') # 22966 31001 
#sce=readH5AD('~/Documents/Data_AgingFlyCellAtlas/adata_Body_S_v1.0.h5ad') # 15992 276273
sce #SingleCellExperiment,  
assayNames(sce) #'X'
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


##########################################################################################
## use MAST test for age effect for each gene
## `DGE_analysis.R` from https://github.com/czbiohub/tabula-muris-senis/tree/master/2_aging_signature/job.DGE_analysis
#output.file='MAST_DEtest_male.rds'
output.file='MAST_DEtest_female.rds'
if(!file.exists(output.file)){
  start_time =proc.time()
  # repeat this for each cell type separately: 5d vs 30d (contain most cells)
  all.tcs=as.character(unique(sce$afca_annotation))
  all.tcs
  
  de_res.out<-lapply(all.tcs,function(tc){
    condition= sce$afca_annotation==tc & sce$sex=='female';
    #condition= sce$afca_annotation==tc & sce$sex=='male';
    #condition= sce$afca_annotation==tc 
    sub.cell.meta=cell.meta[condition,]
    sub.cell.meta$age=droplevels(sub.cell.meta$age)
    if(nrow(sub.cell.meta)==0){return(NULL)}
    if(sum(sub.cell.meta$age=='5')<20 | sum(sub.cell.meta$age=='5')<70){return(NULL)} #age 5 contain less than 20 cells
    colnames(sub.cell.meta) #"n_genes_by_counts"
    
    sce_naive=sce[,condition]
    # Prepare sca object
    sca <- SceToSingleCellAssay(sce_naive, class = "SingleCellAssay",check_sanity = FALSE)
    colData(sca)$age_num = as.numeric(gsub('m','',colData(sca)$age))
    colData(sca)$scaled_n_genes = scale(colData(sca)$n_genes_by_counts) # n_gene (CDR)
    sca_filt = sca[rowSums(assay(sca)) != 0, ]
    
    # DGE testing (only male, so no `sex` as covariate)
    covariate = ''
    covariate = paste0(covariate, " + scaled_n_genes");
    print(paste0('covariate: ', covariate))  
    
    zlmCond <- zlm(formula = as.formula(paste0("~age_num", covariate)), sca=sca_filt)
    summaryCond <- summary(zlmCond, doLRT="age_num")
    
    # Summarize results 
    summaryDt <- summaryCond$datatable
    dt1 = summaryDt[contrast=="age_num" & component=="H", .(primerid, `Pr(>Chisq)`)]
    dt2 = summaryDt[contrast=="age_num" & component=="logFC", .(primerid, coef, z)]
    de_res = merge(dt1, dt2, by="primerid")
    colnames(de_res) <- c("gene", "age.H_p", "age.logFC", 'age.logFC_z')
    de_res$age.H_fdr <- p.adjust(de_res$age.H_p, "fdr")
    #dim(de_res);dim(sca_filt)
    #tmp=de_res[de_res$age.H_fdr>0.05,] # mean-invar genes
    #summary(tmp$age.logFC)
    
    de_res$tc=tc
    return(de_res)
    
  })
  #de_res.out=Filter(Negate(is.null), de_res.out)
  df.de_res.out=as.data.frame(Reduce(`rbind`,de_res.out))
  
  saveRDS(df.de_res.out,output.file)
  print('Finished')
  print(proc.time() - start_time) #2~3hr
  
}

df.de_res.out=readRDS(output.file)
length(unique(df.de_res.out$tc))

################################################################
## from 01_DEgene.analysis.R, DE genes length distribution per cell type
# need df.gene and df.gene.length info

#df.de_res.out=readRDS('MAST_DEtest_male.rds')
df.de_res.out=readRDS('MAST_DEtest_female.rds')
df.de_res.out=df.de_res.out[!is.na(df.de_res.out$age.logFC),]
dim(df.de_res.out)
head(df.de_res.out)
table(df.de_res.out$tc) #DE gene per cell type

tc.names=sort(unique(df.de_res.out$tc)) #39 tc
tmp=df.de_res.out %>% group_by(tc) %>% summarise(nDE=sum(age.H_fdr<0.01)/n())
tmp=tmp[order(tmp$nDE),]
summary(tmp$nDE)

df.out<-lapply(tc.names,function(tc){
  one=df.de_res.out[df.de_res.out$tc==tc,]
  one=one[one$age.H_fdr<0.01,]
  # young high exp gene
  one$DE='young'
  one[one$age.logFC>0,]$DE='old'
  #table(one$DE)
  tmp=merge(one,df.gene,by.x='gene',by.y='original.id')
  sum(is.na(tmp$FLYBASE))
  tmp2=merge(tmp,df.gene.length,by.x='FLYBASE',by.y='ensembl_gene_id')
  #tmp2$gene.length=tmp2$transcript_end-tmp2$transcript_start
  tmp2$gene.length=tmp2$transcript_length
  tmp2=tmp2[order(tmp2$gene.length),]
  tmp2$DE=factor(tmp2$DE,levels=c('young','old'));
  tmp2
  #ggplot(tmp2,aes(x=DE,y=gene.length,group=DE,col=DE))+geom_boxplot()+
  #  scale_y_log10()+theme_classic()+ggtitle(tc)+
  #  stat_compare_means(method='wilcox.test')
})
df=as.data.frame(Reduce(`rbind`,df.out))

p1<-ggplot(df,aes(x=DE,y=gene.length,group=DE,col=DE))+geom_boxplot(outlier.shape = NA)+
  facet_wrap(.~tc)+theme(strip.text = element_text(size=1.5))+
  scale_y_log10()+theme_classic()+
  stat_compare_means(method='wilcox.test',size=2)

pdf('female_DEgene_length.pdf',useDingbats = T,width = 16,height = 16)
#pdf('male_DEgene_length.pdf',useDingbats = T,width = 16,height = 16)
print(p1)
dev.off()
