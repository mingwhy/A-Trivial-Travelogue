
####################################
library(zellkonverter)
library(SingleCellExperiment)
library(SummarizedExperiment)
library(scater);library(scran)
library(ggplot2);library(gridExtra);
library(tidyverse);
library(ggpointdensity)
library(viridis)
library(patchwork); #for plot_annotation
library(lme4)
library(variancePartition)
library(ggpubr)
options(stringsAsFactors = F)
library(org.Mm.eg.db,verbose=F,quietly=T)
library(GO.db);
#https://stackoverflow.com/questions/38722202/how-do-i-change-the-number-of-decimal-places-on-axis-labels-in-ggplot2
scaleFUN <- function(x) sprintf("%.3f", x)
####################################
## get yeast - mouse orthologs
if(F){
  library(biomaRt)
  yeast= useEnsembl(#version = 99, #this archived still contain dn ds values
                    biomart = 'ENSEMBL_MART_ENSEMBL', 
                    dataset = 'scerevisiae_gene_ensembl')
  
  listAttributes(yeast) 
  searchAttributes(yeast,'mmusculus_homolog')
  #searchAttributes(mouse, 'mmusculus_homolog_dn') #yes, there is xxx_dn
  #searchAttributes(mouse, 'mmusculus_homolog_ds') #yes, there is xxx_dn
  
  yeast_mouse <- getBM(attributes = c('ensembl_gene_id', 
                                    'mmusculus_homolog_ensembl_gene',  
                                    #'mmusculus_homolog_dn', 
                                    #'mmusculus_homolog_ds',
                                    'mmusculus_homolog_associated_gene_name',
                                    'mmusculus_homolog_orthology_type',
                                    'mmusculus_homolog_orthology_confidence'), 
                     mart = yeast)
  head(yeast_mouse)
  yeast_mouse=yeast_mouse[yeast_mouse$mmusculus_homolog_ensembl_gene!='',]
  dim(yeast_mouse) #6831     5
  
  library(org.Mm.eg.db)
  yeast_mouse$mmusculus_entrezid=mapIds(org.Mm.eg.db, keys = yeast_mouse$mmusculus_homolog_ensembl_gene, keytype="ENSEMBL", column = "ENTREZID")
  
  data.table::fwrite(yeast_mouse,'yeast_mouse_orthologs.txt')
}
yeast_mouse=data.table::fread('yeast_mouse_orthologs.txt')
table(yeast_mouse$mmusculus_homolog_orthology_confidence)
yeast_mouse=yeast_mouse[yeast_mouse$mmusculus_homolog_orthology_confidence==1,]
table(yeast_mouse$mmusculus_homolog_orthology_type)

#########################################################
## get anti-aging genes from genAge database(when del, lifespan extends)
#download from: https://genomics.senescence.info/genes/models.html
#Baker's yeast
#Download baker's yeast genes: genage_models_export.tsv
#Download homologs from other model organisms: genage_models_orthologs_export.tsv
genage_models=data.table::fread('genage_models_export.tsv')
table(genage_models$Method)
table(genage_models$`Lifespan Effect`)
table(genage_models$`Longevity Influence`)

# ??? anti(del,lifespan increase) or pro(overexpr, lifespan increase) genes, fitness (del, lifespan decrease)
#del.genes<-genage_models[genage_models$Method %in% c('Deletion','depletion','Knockout','RNA interference') & genage_models$`Lifespan Effect`=='increase',]
#del.genes=genage_models[genage_models$`Longevity Influence`=='anti',]
#del.genes=genage_models[genage_models$`Longevity Influence`=='pro',]
del.genes=genage_models[genage_models$`Longevity Influence`=='fitness',]
dim(del.genes) #434
table(del.genes$`Longevity Influence`)

# ??? only keep or not keep chronological relevant genes
#del.genes=del.genes[-grep('Replicative',ignore.case = T,del.genes$`Phenotype Description`),]
#del.genes=del.genes[grep('Replicative',ignore.case = T,del.genes$`Phenotype Description`),]

# merge with ortholog data
i=(del.genes$`Ensembl ID` %in% yeast_mouse$ensembl_gene_id)
sum(i)  #98
del.genes.keep=del.genes[i,]

gene.set=yeast_mouse[yeast_mouse$ensembl_gene_id %in% del.genes$`Ensembl ID`,]
dim(gene.set) #97
length(unique(gene.set$ensembl_gene_id)) #68
length(unique(gene.set$mmusculus_homolog_ensembl_gene)) #92
table(gene.set$mmusculus_homolog_orthology_type)
#ortholog_many2many  ortholog_one2many   ortholog_one2one 
#31                 43                 23 

# ??? what type of orthologs
#gene.set=gene.set[gene.set$mmusculus_homolog_orthology_type=='ortholog_one2one',]
gene.set=gene.set[gene.set$mmusculus_homolog_orthology_confidence==1,] #100
#gene.set=gene.set[gene.set$mmusculus_homolog_orthology_confidence==1 & 
                   #gene.set$mmusculus_homolog_orthology_type=='ortholog_one2one',] #100
                   #gene.set$mmusculus_homolog_orthology_type=='ortholog_one2many',] #100
                   #gene.set$mmusculus_homolog_orthology_type=='ortholog_many2many',] #100
###########################################################
## reading in gene id.mapping in mouse aging atlas data
id.mapping=data.table::fread('~/Documents/Data_mouse_aging_atlas/fac_20449genes_id.mapping.txt')  #readin_h5ad.R
gene.meta=id.mapping

dim(gene.set) #90
sum(gene.set$mmusculus_homolog_ensembl_gene %in% gene.meta$ensembl_gene_id) #40
gene.set2=merge(gene.set,gene.meta,by.x='mmusculus_homolog_ensembl_gene',by.y='ensembl_gene_id')
dim(gene.set2) #90

####################################
## read in mouse turnover rate data
cell.lifespan=data.table::fread('../src/Dataset_S1.txt')
cell.lifespan$Ron_tc=cell.lifespan$`cell type annotation in Sender and Milo 2021`
tc.orders=cell.lifespan[order(cell.lifespan$lifespan),]$`cell type annotation in TMS`
tc.orders #39 tc

#####################################################################################
## use expression proportion (remember to change `expr.m=expr.m[,tmp$age==___ ]`)
#output_file=paste0('log1p_budding.yeast_confidence1_one2one_chronological.rds');
#output_file=paste0('log1p_budding.yeast_confidence1_one2one_replicative.rds');
#output_file=paste0('log1p_budding.yeast_confidence1_one2one.rds');
#output_file=paste0('log1p_budding.yeast_one2one.rds');
#output_file=paste0('log1p_budding.yeast_confidence1_one2many.rds');
#output_file=paste0('log1p_budding.yeast_confidence1_many2many.rds');
#output_file=paste0('log1p_budding.yeast_confidence1.rds');
output_file=paste0('log1p_budding.yeast_fitness_confidence1.rds');
#output_file=paste0('log1p_budding.yeast_pro_confidence1.rds');
#output_file=paste0('cpm_budding.yeast_confidence1_one2one.rds');
#output_file=paste0('select.tc.h5ad_log1p_budding.yeast_out.rds');

if(!file.exists(output_file)){
  
  sce=readH5AD('../0708_TMS_male_3m_24m/select.tc.h5ad') 
  unique(sce$age) #3m 18m 24m
  assayNames(sce)<-'counts'

  sce_naive=sce[,sce$age %in% c('3m','24m')]
  table(sce_naive$age)
  assayNames(sce_naive)<-'counts'
  #sce_naive<-logNormCounts(sce_naive, log=FALSE, pseudo.count=1) # use CPM, if log, the cell.lib.size range 2~3 orders
  sce_naive<-logNormCounts(sce_naive, log=TRUE, pseudo.count=1.1) #if log, the cell.lib.size range 2~3 orders
  assayNames(sce_naive)
  
  df.mouse=as.numeric();
  
  for(pick.age in c('3m','24m')){
    
    mouse_tcs_TDI.list<-lapply(tc.orders,function(tc){
      
      tmp=sce_naive[,sce_naive$tissue_cell.type==tc] #raw count data
      
      #expr.m=assay(tmp,'normcounts')
      expr.m=assay(tmp,'logcounts')
      expr.m=expr.m[,tmp$age==pick.age]
      #expr.m=expr.m[,tmp$age=='24m']
      
      n.expr.gene=Matrix::colSums(expr.m>0) #control for the number of expressed genes （https://academic.oup.com/gbe/article/12/4/300/5807614?login=true
      #expr.m=expr.m[,sce_naive$age=='24m' & n.expr.gene>=100] #keep cells which expr>=100 genes
      expr.m=expr.m[, n.expr.gene>=100] #keep cells which expr>=100 genes
      
      #used.genes=intersect(gene.meta$mgi_symbol,rownames(expr.m))
      #expr.m=expr.m[used.genes,]
      cell.expr=Matrix::colSums(expr.m)
      all.genes=rownames(expr.m)
      
      overlap.genes=intersect(all.genes,gene.set2$mgi_symbol)
      if(length(overlap.genes)==0){return(c(pick.age,0,tc,NA))}
      expr.m.tmp=expr.m[overlap.genes,,drop=F]
      go.expr=Matrix::colSums(expr.m.tmp)
      score<-mean(go.expr/cell.expr)
    
      return(c(pick.age,length(overlap.genes),tc,score))
    })
    
    tmp=as.data.frame(Reduce(`rbind`,mouse_tcs_TDI.list))
    df.mouse=rbind(df.mouse,tmp)
  }
  df.mouse=as.data.frame(df.mouse)
  colnames(df.mouse)=c('age','ngene','cell.type','score')
  saveRDS(df.mouse,output_file)
}

####################################
### see correlation
df.go.score=readRDS(output_file)
df.go.score$score=as.numeric(df.go.score$score)
df.go.score$ngene=as.numeric(df.go.score$ngene)
summary(df.go.score$ngene) 

overlap.tc=intersect(df.go.score$cell.type,cell.lifespan$`cell type annotation in TMS`)
length(overlap.tc) #39

df2=cell.lifespan[match(overlap.tc,cell.lifespan$`cell type annotation in TMS`),]

df2$duplicate='tissue-specific estimate';
df2[df2$Ron_tc %in% df2$Ron_tc[duplicated(df2$Ron_tc)],]$duplicate='non tissue-specific estimate'
table(df2$duplicate)

df1=df2[match(df.go.score$cell.type,df2$`cell type annotation in TMS`),]
sum(df1$`cell type annotation in TMS`==df.go.score$cell.type)

df3=cbind(df.go.score,df1)
colnames(df3)
#######################################################################################
## use average for non-tissue specific ones
df4=df3 %>% group_by(Ron_tc,age) %>% summarise(average=mean(score))

df5=merge(df4,df2[!duplicated(df2$Ron_tc),],all.x=TRUE)
dim(df5) #22*2, 2 age groups

df5$turnover=1-exp(-1/df5$lifespan)
cor.out=lapply(c('3m','24m'),function(pick.age){
  x=df5[df5$age==pick.age,]
  i=cor.test(x$average,x$lifespan,method='kendall',use='pairwise.complete.obs')
  c(pick.age,i$estimate,i$p.value)
})
df.cor.out<-as.data.frame(Reduce(`rbind`,cor.out))
df.cor.out[,2]=as.numeric(df.cor.out[,2])
df.cor.out[,3]=as.numeric(df.cor.out[,3])

plots=lapply(c('3m','24m'),function(pick.age){
    cor0=round(df.cor.out[df.cor.out$V1==pick.age,]$tau,4);
    pval0=round(df.cor.out[df.cor.out$V1==pick.age,3],5);
    
    ggplot(subset(df5,age==pick.age),aes(x=lifespan,y=average))+
      #geom_point(size=3,shape=16)+
      geom_point(aes(fill=duplicate),pch=21,size=5,width=6)+
      scale_x_log10()+ylab('Gene set activity score')+
      ggtitle(paste0('Age ',pick.age,'\nKendall\'s tau = ',cor0,', P.adj value=',pval0))+
      xlab('Cell lifespan (day)')+
      scale_fill_manual(name='Cell lifespan estimate',values=c('NA','black'))+
      #geom_smooth(method=lm , color="black", fill=grDevices::adjustcolor( "lightgrey", alpha.f = 0.2),
      #            se=TRUE,level=0.95) + 
      scale_shape_discrete(name='Cell lifespan estimate')+
      scale_y_continuous(labels=scaleFUN)
    #scale_color_viridis(name='Cell type annotation from Sender and Milo (2021)',option='turbo',discrete=T)
})

plots2=lapply(plots,function(i) 
  i+theme_classic(base_size = 12)+theme(legend.position = 'none',
                                        plot.margin = unit(c(1,1,1,1), "cm"),
                                        #plot.title = element_text(size = 12, face = "bold")))
                                        plot.title = element_text(size = 12)))


pdf(paste0(output_file,'.pdf'),useDingbats = T,width = 10,height = 4,pointsize=12)
#for(i in plots2){print(i)}
grid.arrange(grobs=plots2,ncol=2)
dev.off()

if(F){
plots3=lapply(c('3m','24m'),function(pick.age){
  #cor0=round(df.cor.out[df.cor.out$V1==pick.age,]$tau,4);
  #pval0=round(df.cor.out[df.cor.out$V1==pick.age,3],5);
  
  ggplot(subset(df3,age==pick.age),aes(x=lifespan,y=score))+
    geom_point(aes(shape=duplicate),size=3)+
    scale_x_log10()+ylab('Gene set activity score')+
    #ggtitle(paste0('Age ',pick.age,'\nKendall\'s tau = ',cor0,', P.adj value=',pval0))+
    xlab('Cell lifespan (day)')+
    scale_shape_manual(name='Cell lifespan estimate',values=c(1,2))+
    #geom_smooth(method=lm , color="black", fill=grDevices::adjustcolor( "lightgrey", alpha.f = 0.2),
    #            se=TRUE,level=0.95) + 
    scale_y_continuous(labels=scaleFUN)
  #scale_color_viridis(name='Cell type annotation from Sender and Milo (2021)',option='turbo',discrete=T)
})
plots4=lapply(plots3,function(i) 
  i+theme_classic(base_size = 12)+theme(legend.position = 'none',
                                        plot.margin = unit(c(1,1,1,1), "cm"),
                                        #plot.title = element_text(size = 12, face = "bold")))
                                        plot.title = element_text(size = 12)))
grid.arrange(grobs=plots4,ncol=2)
}


