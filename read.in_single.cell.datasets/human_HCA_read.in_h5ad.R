# The Tabula Sapiens: A multiple-organ, single-cell transcriptomic atlas of humans
#https://www.science.org/stoken/author-tokens/ST-495/full#sec-8
#https://tabula-sapiens-portal.ds.czbiohub.org/whereisthedata
#https://figshare.com/articles/dataset/Tabula_Sapiens_release_1_0/14267219

library(zellkonverter)
library(SingleCellExperiment)
library(scater);library(scran)
library(ggplot2);library(gridExtra);
library(tidyverse);
library(zellkonverter)
library(SummarizedExperiment)
library(ggpointdensity)
library(viridis)
library(patchwork); #for plot_annotation
library(lme4)
library(variancePartition)
library(ggpubr)

## read in data
if(F){
sce=readH5AD('14267219/TabulaSapiens.h5ad') # 
sce #58870 483152 

cell.meta=colData(sce)
saveRDS(cell.meta,'TabulaSapiens_cell.meta.rds')
}

cell.meta=readRDS('TabulaSapiens_cell.meta.rds')
colnames(cell.meta)
head(cell.meta)
length(table(cell.meta$organ_tissue)) #24 organs
table(cell.meta$compartment)
#endothelial  epithelial   germ line      immune     stromal 
#31691      104148          11      264824       82478 
length(table(cell.meta$cell_ontology_class)) #177
length(table(cell.meta$free_annotation)) #320

cell.meta$tissue_tc=paste(cell.meta$organ_tissue,cell.meta$cell_ontology_class,sep=':')
x=as.data.frame(table(cell.meta$free_annotation,cell.meta$tissue_tc))
colnames(x)=c('free_annotation','cell_ontology_class','ncell')
x=x[x$ncell!=0,]
dim(x)

tmp=cell.meta[cell.meta$cell_ontology_class=='acinar cell of salivary gland',]

tc.names=sort(unique(sce$tissue_cell.type))
tc.names #38
sce.shared=lapply(tc.names,function(tc){
  sce[,sce$tissue_cell.type==tc]
})
names(sce.shared)<-tc.names

###############################################################################################
## reading in gene id.mapping and dn.ds infomation
id.mapping=data.table::fread('~/Documents/Data_mouse_aging_atlas/fac_20449genes_id.mapping.txt')  #readin_h5ad.R
filter_mouse_rat=data.table::fread('../gene.properties/dNdS/mouse_rat.dnds.txt')
head(filter_mouse_rat)
dim(filter_mouse_rat) #30033 genes
dim(id.mapping) #20449
sum(filter_mouse_rat$ensembl_gene_id %in% id.mapping$ensembl_gene_id) #20697

x=filter_mouse_rat[filter_mouse_rat$rnorvegicus_homolog_orthology_type=='ortholog_one2one',]
anyDuplicated(x$ensembl_gene_id)
#mouse_rat_dnds=filter_mouse_rat[!duplicated(filter_mouse_rat$ensembl_gene_id),] #random choose one ortholog in rat
mouse_rat_dnds=x; #16494

head(mouse_rat_dnds)
head(id.mapping)
gene.meta=merge(mouse_rat_dnds,id.mapping,by.x='ensembl_gene_id',by.y='ensembl_gene_id')
dim(gene.meta) #14821
gene.meta$omega=gene.meta$rnorvegicus_homolog_dn/gene.meta$rnorvegicus_homolog_ds
gene.meta[is.infinite(gene.meta$omega),]$omega=NA
gene.meta[is.nan(gene.meta$omega),]$omega=NA
summary(gene.meta$omega)  #max=2.24 
gene.meta[which(gene.meta$omega>=1),] #75 genes with omega>1

#remove genes with dn/ds>1 (more likely to be under positive selection)
length(which(gene.meta$omega<=1)) #14110 genes
gene.meta=gene.meta[!is.na(gene.meta$omega),]
sum(gene.meta$omega>=1) #76 genes
#gene.meta=gene.meta[which(gene.meta$omega<1),]
dim(gene.meta) #14185 or 14109    12

summary(replicate(100,median(gene.meta[sample(1:nrow(gene.meta),500,replace = F),]$omega)))
summary(replicate(100,median(gene.meta[sample(1:nrow(gene.meta),3500,replace = F),]$omega)))
#sum(gene.meta$rnorvegicus_homolog_orthology_confidence)
#gene.meta=gene.meta[gene.meta$rnorvegicus_homolog_orthology_confidence==1,]
###############################################################
## calculate TDI per cell
sapply(sce.shared,dim)
(tc.names=names(sce.shared))

mouse_tcs_TDI.list<-lapply(tc.names,function(tc){
  sce_naive=sce.shared[[tc]]
  #x=lapply(x,function(i) i[rowData(i)$include_gene,])
  #sapply(x,dim)
  assayNames(sce_naive)<-'counts'
  
  #out<-lapply(x,function(sce_naive){
  #  summary(sizeFactors(sce_naive)) #already calculated in 02_gene_meanVar_shareGenes.R
  sce_naive <- logNormCounts(sce_naive,log = TRUE)
  assayNames(sce_naive)
  #expr.m=assay(sce_naive,'normcounts')
  expr.m=assay(sce_naive,'logcounts')
  
  overlap.genes=intersect(rownames(expr.m),gene.meta$mgi_symbol)
  cat(tc,nrow(expr.m),length(overlap.genes),'\n');
  
  expr.m=expr.m[overlap.genes,]
  n.expr.gene=Matrix::colSums(expr.m>0) #control for the number of expressed genes （https://academic.oup.com/gbe/article/12/4/300/5807614?login=true
  
  i=match(overlap.genes,gene.meta$mgi_symbol)
  gene.meta.m=gene.meta[i,]
  dim(gene.meta.m);dim(expr.m)
  sum(gene.meta.m$mgi_symbol==rownames(expr.m))
  
  expr.m.binary=expr.m;
  expr.m.binary[expr.m.binary>0]=1;
  #tmp=Matrix::rowSums(expr.m.binary)
  #tmp=tmp[order(tmp,decreasing = T)]
  #tmp1=gene.meta.m[gene.meta.m$mgi_symbol %in% names(tmp)[1:20],]
  #tmp1=tmp1[order(tmp1$omega,decreasing = T),]
  x=gene.meta.m$omega %*% as.matrix(expr.m.binary)
  #x=gene.meta.m$rnorvegicus_homolog_ds %*% as.matrix(expr.m.binary)
  index.per.cell=x/n.expr.gene;
  
  length(index.per.cell)
  cell.meta=colData(sce_naive)
  dim(cell.meta)
  cell.meta$TDI=index.per.cell;
  cell.meta$n_expr_gene=n.expr.gene
  cell.meta=as.data.frame(cell.meta)
  return(cell.meta)
})

#mouse_tcs_TDI=purrr::flatten(mouse_tcs_TDI.list)
mouse_tcs_TDI=as.data.frame(Reduce(`rbind`,mouse_tcs_TDI.list))
head(mouse_tcs_TDI)
summary(mouse_tcs_TDI$sizeFactor) #as we use subsampled UMI
summary(mouse_tcs_TDI$n_expr_gene) #44 to 4482

saveRDS(mouse_tcs_TDI,'mouse_male_binary_TDI.rds')

p0= ggplot(mouse_tcs_TDI,aes(x=n_expr_gene,y=TDI,col=age))+
  facet_wrap(.~tissue_cell.type,scale='free')+geom_jitter(size=0.2)+
  scale_x_log10()+theme_classic(base_size = 12)+
  ylab('Median dN/dS of expressed genes per cell')
pdf('test.pdf',useDingbats = T,width=16,height = 10)
print(p0);
dev.off();    

###################################################################### 
## plot
mouse_tcs_TDI=readRDS('mouse_male_binary_TDI.rds')
summary(mouse_tcs_TDI$n_expr_gene)
#df=as.data.frame(Reduce(`rbind`,mouse_tcs_TDI))
df=mouse_tcs_TDI[mouse_tcs_TDI$n_genes>=100,]
head(df)
colnames(df)

ggplot(df,aes(x=age,y=TDI))+geom_violin()+geom_jitter(size=0.2)+
  facet_wrap(.~tissue_cell.type)+
  stat_summary(fun.y=median, geom="point", size=1, color="red")+
  theme_classic()+stat_compare_means(label = "p.format",label.y = 0.25)

ggplot(df,aes(x=tissue_cell.type,y=TDI,col=tissue_cell.type))+geom_violin()+geom_jitter(size=0.2)+
  facet_wrap(.~age,nrow=3)+
  stat_summary(fun.y=median, geom="point", size=2, color="red")+
  theme_classic()

# order tc by median TDI at 3m
one.age=subset(df,age=='3m');
x=one.age %>% group_by(tissue_cell.type) %>% summarise(median.TDI=median(TDI))
x=x[order(x$median.TDI),]
x

one.age$tissue_cell.type=factor(one.age$tissue_cell.type,levels=x[order(x$median.TDI),]$tissue_cell.type)
ggplot(one.age,aes(x=tissue_cell.type,y=TDI,col=tissue_cell.type))+
  geom_jitter(size=0.2)+theme_classic()+
  #theme(axis.text.x = element_blank())+
  theme(axis.text.x = element_text(size=9,hjust=0.5,angle=45))+
  guides(color = guide_legend(override.aes = list(size = 2)))

cell.meta=df
cell.meta$tissue_cell.type=factor(cell.meta$tissue_cell.type,levels=levels(one.age$tissue_cell.type))
p1=ggplot(cell.meta,aes(x=tissue_cell.type,y=TDI,col=tissue_cell.type))+
  facet_wrap(.~age,ncol=1)+
  geom_jitter(size=0.2)+theme_classic()+
  ylab('Median dN/dS of expressed genes per cell')+
  theme(axis.text.x = element_text(size=9,hjust=1,vjust=1,angle=45),
        legend.position = 'none')+ 
  guides(color = guide_legend(override.aes = list(size = 2)))
p1

#p2=ggplot(cell.meta,aes(x=age,y=TDI,col=tissue_cell.type))+
p2=ggplot(cell.meta,aes(x=age,y=TDI))+
  facet_wrap(.~tissue_cell.type,ncol=7,scale='free')+
  geom_jitter(size=0.2)+theme_classic()+
  ylab('Median dN/dS of expressed genes per cell')+
  stat_summary(fun.y=median, geom="point", size=1, color="red")+
  theme(axis.text.x = element_text(size=9,hjust=0,vjust=0,angle=0),
        legend.position = 'none')
#guides(color = guide_legend(override.aes = list(size = 2)))
p2

pdf('mouse_male_binary_TDI.pdf',useDingbats = T,width = 18)
print(p1);
dev.off();

pdf('mouse_male_binary_TDI2.pdf',useDingbats = T,width=14,height = 16)
print(p2+stat_compare_means(label = "p.format",test='wil.cox')) 
dev.off();


