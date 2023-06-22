
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
options(stringsAsFactors = F)
library(org.Mm.eg.db,verbose=F,quietly=T)
library(GO.db);
## read in data
if(F){
  sce=readH5AD('~/Documents/aging_cell.turnover/1120_TMS_male_analysis/select.tc.h5ad') # 22966 31001 
  unique(sce$age) #3m 18m 24m
  tc.names=sort(unique(sce$tissue_cell.type))
  tc.names #38
  sce.shared=lapply(tc.names,function(tc){
    sce[,sce$tissue_cell.type==tc]
  })
  names(sce.shared)<-tc.names
}
if(T){
  sce.shared=readRDS('~/Documents/aging_cell.turnover/1120_TMS_male_analysis/sce_minCellCounts.rds')
  names(sce.shared) #38 tc
  unique(Matrix::colSums(assay(sce.shared[[1]],'counts'))) #check for cell.lib.size
}
## read in GO<->gene 
#slim.go2fb=readRDS('slim.go2SYMBOL.5_200.rds')
slim.go2fb=readRDS('Mmus_go2SYMBOL_BP.rds')

length(slim.go2fb) #1882 or 12545 (all BP GOterms)
ngene.per.go=sapply(slim.go2fb,nrow)
sum(ngene.per.go>5) # 4516 (all BP GOterms)
slim.go2fb=slim.go2fb[ngene.per.go>5]

## get GO term full name
goterms=Term(GOTERM)
length(goterms) #43851
names(slim.go2fb)[[1]]
GOTERM$"GO:0000002" #mitochondrial genome maintenance

###############################################################################################
## reading in gene id.mapping and extract MHC genes
id.mapping=data.table::fread('~/Documents/Data_mouse_aging_atlas/fac_20449genes_id.mapping.txt')  #readin_h5ad.R
HSP.genes<-id.mapping[grep('shock',id.mapping$description),]
dim(HSP.genes) #86
gene.meta=id.mapping

###############################################################################################
## calculate HSP expr per cell
sapply(sce.shared,dim)
(tc.names=names(sce.shared))

mouse_tcs_TDI.list<-lapply(tc.names,function(tc){
  sce_naive=sce.shared[[tc]]
  assayNames(sce_naive)<-'counts'
  expr.m=assay(sce_naive,'counts')
  #expr.m=expr.m[,sce_naive$age=='3m']
  
  n.expr.gene=Matrix::colSums(expr.m>0) #control for the number of expressed genes （https://academic.oup.com/gbe/article/12/4/300/5807614?login=true
  expr.m=expr.m[,sce_naive$age=='3m' & n.expr.gene>=100] #keep cells which expr>=100 genes
  
  cell.expr=Matrix::colSums(expr.m)
  all.genes=rownames(expr.m)
  
  go.median.expr<-lapply(names(slim.go2fb), function(go){
    overlap.genes=intersect(all.genes,slim.go2fb[[go]][,2])
    if(length(overlap.genes)==0){return(NA)}
    
    expr.m.tmp=expr.m[overlap.genes,,drop=F]
    go.expr=Matrix::colSums(expr.m.tmp)
    return(median(go.expr/cell.expr))  
  })
  return(unlist(go.median.expr))
})
#mouse_tcs_TDI=purrr::flatten(mouse_tcs_TDI.list)
mouse_tcs_TDI=as.data.frame(Reduce(`rbind`,mouse_tcs_TDI.list))
rownames(mouse_tcs_TDI)=tc.names
colnames(mouse_tcs_TDI)=names(slim.go2fb)

saveRDS(mouse_tcs_TDI,'mouse_male_GO.rds')

######################################################################
## read in mouse turnover rate data
#cell.lifespan=readxl::read_xlsx('~/Documents/aging_cell.turnover/2021_RonMilo_data/mouse_tc_cell.lifespan_20221117.xlsx');
cell.lifespan=readxl::read_xlsx('~/Documents/aging_cell.turnover/2021_RonMilo_data/mouse_tc_cell.lifespan_20221100.xlsx');
dim(cell.lifespan) #139
colnames(cell.lifespan);
head(cell.lifespan)

cell.lifespan=cell.lifespan[cell.lifespan$include==1,] #filter based on Kim's annotation
dim(cell.lifespan) #39 tc remain 
cell.lifespan$lifespan=cell.lifespan$lifespan.used.in.this.study

cell.lifespan[grep('-',cell.lifespan$lifespan),]
#https://onlinelibrary.wiley.com/doi/10.1111/imr.12693
cell.lifespan[cell.lifespan$lifespan=='40-70' & cell.lifespan$human_tc=='Blood:mature T cells',]$lifespan='70'; #T cell 
cell.lifespan[cell.lifespan$lifespan=='40-70' & cell.lifespan$human_tc=='Blood:mature B cells',]$lifespan='40'; #B cell 
#https://jamanetwork.com/journals/jamainternalmedicine/article-abstract/565579
cell.lifespan[cell.lifespan$lifespan=='200-400',]$lifespan='400'; #hepatocyte
cell.lifespan$lifespan=as.numeric(cell.lifespan$lifespan)
dim(cell.lifespan)  #39 tc

tmp=cell.lifespan[!duplicated(cell.lifespan$human_tc),]
dim(tmp); #unique 21 human_tc
table(tmp$species) #1 human and 20 rodents

tc.orders=cell.lifespan[order(cell.lifespan$lifespan),]$`tissue: cell.type in mouse`
tc.orders #37 tc

###################################################################### 
## plot with lifespan
mouse_tcs_TDI=readRDS('mouse_male_GO.rds')
all.go=colnames(mouse_tcs_TDI)

overlap.tc=intersect(rownames(mouse_tcs_TDI),cell.lifespan$`tissue: cell.type in mouse`)
df1=mouse_tcs_TDI[overlap.tc,]
df2=cell.lifespan[match(overlap.tc,cell.lifespan$`tissue: cell.type in mouse`),]
sum(rownames(df1)==df2$`tissue: cell.type in mouse`) #37

go.cors.values=cor(df1,df2$lifespan,method='spearman',use='pairwise.complete.obs')
length(go.cors.values) #1882 go terms
summary(go.cors.values)
df1[,which(is.na(go.cors.values))[1]]

go.cors.values=go.cors.values[!is.na(go.cors.values[,1]),]
cutoff=0.9
tmp=go.cors.values[abs(go.cors.values)>=cutoff]
length(names(tmp))
GOTERM$"GO:0043266"
sapply(names(tmp),function(i)GOTERM[[i]])
sapply(names(tmp),function(i)GOTERM[[i]]@Term)

library(GOSemSim) #https://jokergoo.github.io/simplifyEnrichment/articles/simplifyGO.html
library(simplifyEnrichment)
go_id=names(tmp)
mat = GO_similarity(go_id,ont='BP')
GO_similarity(go_id, measure = "Wang") #heatmap generated
df = simplifyGO(mat)
head(df)


#########
na.number=apply(df1,2,function(i) sum(is.na(i)))
go.cors.values=cor(df1,df2$lifespan,method='spearman',use='pairwise.complete.obs')
length(go.cors.values) #1882 go terms
summary(go.cors.values)
df1[,which(is.na(go.cors.values))[1]]

go.cors.values=go.cors.values[na.number<=5,]
go.cors.values=go.cors.values[!is.na(go.cors.values)]
summary(go.cors.values)
cutoff=0.6
sig.go.terms=go.cors.values[abs(go.cors.values)>cutoff]
sig.go.terms #24
GOTERM$"GO:0043266"
sapply(names(sig.go.terms),function(i)GOTERM[[i]])

#slim.go2fb[names(tmp)]
slim.go2fb[['GO:0071453']]

hit.genes=unique(unlist(lapply(slim.go2fb[names(sig.go.terms)],'[',2)))
length(hit.genes)
tmp=gene.meta[gene.meta$mgi_symbol %in% hit.genes,]
tmp[grep('shock',tmp$description),]


