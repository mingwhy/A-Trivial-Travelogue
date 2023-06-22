
setwd("~/Documents/sc_transcriptome.index/pseudotime_analysis/")

#BiocManager::install(c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
#                       'limma', 'S4Vectors', 'SingleCellExperiment',
#                       'SummarizedExperiment', 'batchelor', 'Matrix.utils'))

library(BiocGenerics)
library(DelayedArray)
library(DelayedMatrixStats)
library(limma)
library(S4Vectors)
library(SingleCellExperiment)
library(SummarizedExperiment)
library(batchelor)
library(Matrix.utils)
#devtools::install_github('cole-trapnell-lab/leidenbase')
library(leidenbase)
#devtools::install_github('cole-trapnell-lab/monocle3')
library(monocle3)

library(ggplot2)
library(dplyr)
library(ggpubr)
#install.packages('Seurat')
library(Seurat)

########################################################
## read in dataset
#BiocManager::install("zellkonverter")
library(zellkonverter)
library(SummarizedExperiment)

## read in data
sce<-readH5AD('~/Documents/Data_mouse_aging_atlas//TMS.gene.data_final/tabula-muris-senis-facs-official-raw-obj.h5ad');       
class(assay(sce))
sce


dim(assay(sce))
df.expr=assay(sce)
df.expr[1:10,1:10] #umi count
dim(df.expr) #22966 110824

## sample meta information
cell.meta=colData(sce)
head(cell.meta)
dim(cell.meta) #110824     13
colnames(cell.meta)


table(cell.meta$age) #6 age group: 1, 3, 18, 21, 24,30month. fac: 3,18,21,24 month
table(cell.meta$tissue) #23 tissue
length(table(cell.meta$cell_ontology_class)) #120 cell types
unique(paste(cell.meta$tissue,cell.meta$cell_ontology_class)) #207 unique #consistent with `elife-62293-supp1-v2.xlsx`
x=as.data.frame(cell.meta) %>% group_by(tissue) %>% summarise(n.cell.type=length(unique(cell_ontology_class)))
x$n.cell.type

tissue_cell.type=paste(cell.meta$tissue,cell.meta$cell_ontology_class,sep=':')
cell.meta$tissue_cell.type=tissue_cell.type

table(cell.meta[cell.meta$age=='21m',]$tissue_cell.type) #remove these cells
i=cell.meta$age!='21m'
cell.meta=cell.meta[i,]
df.expr=df.expr[,i]
cell.meta$age=factor(cell.meta$age,levels = c('3m','18m','24m'))
dim(cell.meta) # 110096     14
dim(df.expr) # 22966 110096


#######################################################################
## read in elife-62293-supp1-v2.xlsx, filter to 115 tc
library(readxl)
df=read_xlsx('~/Documents/Data_mouse_aging_atlas/SuppTables1-3/elife-62293-supp1-v2.xlsx',
             sheet='all FACS tissue-cell')
colnames(df)
dim(df) #207 tc
head(df)
unique(df$tissue) #23
unique(df$cell_ontology_class) #120
tc=paste(df$tissue,df$cell_ontology_class,sep=':')
length(unique(tc)) #207

## filter tc, contain>=20 cells in all three age groups
min.cell=20
x=df[df$`3m`>=min.cell & df$`18m`>=min.cell & df$`24m`>=min.cell,]
dim(x) #115 tc
table(x$tissue)
tc=paste(x$tissue,x$cell_ontology_class,sep=':')

df.expr=df.expr[,cell.meta$tissue_cell.type %in% tc]
cell.meta=cell.meta[cell.meta$tissue_cell.type %in% tc,]
dim(cell.meta) # 100817     14
dim(df.expr) # 22966 100817

###############################################################################################
## calculate pseudotime for each tc(https://cole-trapnell-lab.github.io/monocle3/docs/starting/)
library(grDevices);library(RColorBrewer)
library(uwot); #for UMAP
library(scran) #for clustering cells
library(slingshot, quietly = FALSE)
#https://bioconductor.org/packages/devel/bioc/vignettes/slingshot/inst/doc/vignette.html
#full quantile normalization, a well-established method which forces each cell to have the same distribution of expression values.
FQnorm <- function(counts){
  rk <- apply(counts,2,rank,ties.method='min')
  counts.sort <- apply(counts,2,sort)
  refdist <- apply(counts.sort,1,median)
  norm <- apply(rk,2,function(r){ refdist[r] })
  rownames(norm) <- rownames(counts)
  return(norm)
}

all.tcs=sort(unique(cell.meta$tissue_cell.type))

plots=list();ps.ages=list();
#for(tc in all.tcs[1:4]){
for(tc in all.tcs){
  expression_matrix=df.expr[,cell.meta$tissue_cell.type==tc]
  gene.names <- rownames(expression_matrix)
  cell.names <- colnames(expression_matrix)
  
  cell_metadata<-cell.meta[cell.meta$tissue_cell.type==tc,]
  rownames(cell_metadata)=cell.names
  
  sce <- SingleCellExperiment(assays = List(counts = expression_matrix))
  # gene filter
  geneFilter <- apply(assays(sce)$counts,1,function(x){
    sum(x >= 3) >= 10
  })
  sce <- sce[geneFilter, ]
  # quantile normalize
  assays(sce)$norm <- FQnorm(assays(sce)$counts)
  # pca
  pca <- prcomp(t(log1p(assays(sce)$norm)), scale. = FALSE)
  rd1 <- pca$x[,1:2]
  #plot(rd1, col = rgb(0,0,0,.5), pch=16, asp = 1)
  # umap
  rd2 <- uwot::umap(t(log1p(assays(sce)$norm)))
  colnames(rd2) <- c('UMAP1', 'UMAP2')
  #plot(rd2, col = rgb(0,0,0,.5), pch=16, asp = 1)
  
  #add both dimensionality reductions to the SingleCellExperiment object
  reducedDims(sce) <- SimpleList(PCA = rd1, UMAP = rd2)
  
  # cluster cells
  sce <- logNormCounts(sce)
  #clusters <- scran::clusterCells(sce, assay.type="logcounts")
  # From PCs:
  clusters <- clusterCells(sce, use.dimred="PCA")
  table(clusters)
  #plot(rd1, col = brewer.pal(9,"Set1")[clusters], pch=16, asp = 1)
  colData(sce)$GMM <- clusters
  
  # begin slingshot
  sce <- slingshot(sce, clusterLabels = 'GMM', reducedDim = 'PCA')
  summary(sce$slingPseudotime_1)
  
  colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
  plotcol <- colors[cut(sce$slingPseudotime_1, breaks=100)]
  
  dim(rd1);dim(cell_metadata)
  
  par(mfrow=c(1,2),mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
  plot(rd1, col = plotcol[as.numeric(factor(cell_metadata$age))], pch=16, asp = 1)
  #legend('topright',inset=c(-0.2,0),legend=levels(cell_metadata$age),col=plotcol[1:3],pch=16, horiz=TRUE, bty='n', cex=0.8)
  legend(0,0,legend=levels(cell_metadata$age),col=plotcol[1:3],pch=16, horiz=TRUE, cex=0.5,bty='n')
  
  plot(reducedDims(sce)$PCA, col = plotcol, pch=16, asp = 1)
  lines(SlingshotDataSet(sce), lwd=2, col='black')
  
  #############################
  #############################
  #############################
  
  x=pseudotime(cds)
  length(x) #number of cells
  sum(is.infinite(x))
  cds@colData$pseudotime=pseudotime(cds)
  
  df=cds@colData;
  p1=plot_cells(cds,
                color_cells_by = "age",
                show_trajectory_graph = FALSE,
                label_cell_groups=FALSE,
                label_leaves=FALSE,
                label_branch_points=FALSE,
                graph_label_size=1.5)+ggtitle(tc)
  p2=plot_cells(cds,
                color_cells_by = "pseudotime",
                label_cell_groups=FALSE,
                label_leaves=FALSE,
                label_branch_points=FALSE,
                graph_label_size=1.5)+ggtitle(tc)
  #Note that some of the cells are gray. This means they have infinite pseudotime, 
  #because they were not reachable from the root nodes that were picked. 
  tmp=data.frame(age=df$age,pseudotime=df$pseudotime)
  tmp$pseudotime=tmp$pseudotime+1; #avoid 0 problem
  p3=ggplot(tmp,aes(x=age,y=pseudotime,col=age))+
    geom_violin()+geom_jitter(size=0.2)+scale_y_log10()+
    theme_classic()+ggtitle(tc)
  
  plots=c(plots,list(p1,p2,p3))
  ps.ages[[tc]]=df;
}
#saveRDS(ps.ages,'pseudotime_115tc.rds')
saveRDS(ps.ages,'pseudotime_115tc_pr_nodes.rds')

#pdf('monocle3_plots.pdf',useDingbats = T,width = 12)
pdf('monocle3_plots_pr_nodes.pdf',useDingbats = T,width = 12)
for(i in seq(1,length(plots),6)){
  gridExtra::grid.arrange(grobs=plots[i:min((i+5),length(plots))],ncol=3)
}
dev.off();


#########################################################################
# https://cole-trapnell-lab.github.io/monocle3/docs/trajectories/
#It's often desirable to specify the root of the trajectory programmatically, rather than manually picking it. The function below does so by first grouping the cells according to which trajectory graph node they are nearest to. Then, it calculates what fraction of the cells at each node come from the earliest time point. Then it picks the node that is most heavily occupied by early cells and returns that as the root.
# a helper function to identify the root principal points:
get_earliest_principal_node <- function(cds, time_bin="3m"){
  cell_ids <- which(colData(cds)[, "age"] == time_bin)
  
  closest_vertex <-
    cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <-
    igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
                                                              (which.max(table(closest_vertex[cell_ids,]))))]
  
  root_pr_nodes
}
#cds <- order_cells(cds, root_pr_nodes=get_earliest_principal_node(cds))

###################################################################
###################################################################
###################################################################
## plot percentage of young, middle, old cells per ps time bin
ps.ages=readRDS('pseudotime_115tc.rds')
length(ps.ages)
names(ps.ages)
i=15
#plots=lapply(1:length(ps.ages),function(i){
tc=as.data.frame(ps.ages[[i]])
tc.name=names(ps.ages)[i]
head(tc)
summary(tc$pseudotime)
br=seq(0,max(tc$pseudotime)*1.01,len=20)
gr=cut(tc$pseudotime,breaks=br,include.lowest = T)
tc$gr=gr;
x=tc %>% group_by(gr,age) %>% summarise(n=n())
x=x %>% group_by(gr) %>% mutate(ncell=sum(n))
x$proportion=x$n/x$ncell
head(x)
ggplot(x,aes(x=gr,y=proportion,group=age,col=age))+geom_line()+theme_classic()+
  xlab('pseudotime')+theme(axis.text.x = element_blank())+ggtitle(tc.name)
#})

###################################################################
## calculate gene spearman.cor with cell pseudotime  

i=15
tc=as.data.frame(ps.ages[[i]])
tc.name=names(ps.ages)[i]

df.sub=df.expr[,cell.meta$tissue_cell.type==tc.name]
dim(df.sub)
dim(tc)
sp.cors=cor(t(as.matrix(df.sub)),tc$pseudotime,method='pearson')
length(sp.cors) #number of genes
summary(sp.cors) #there may be NA due to 0 sd for some genes
sp.cors=as.numeric(sp.cors)
names(sp.cors)=rownames(df.sub)
sp.cors=sort(sp.cors)
tail(sp.cors)


##################################################################
## reading in gene id.mapping and dn.ds infomation
id.mapping=data.table::fread('~/Documents/Data_mouse_aging_atlas/fac_20449genes_id.mapping.txt')  #readin_h5ad.R
filter_mouse_rat=data.table::fread('~/Documents/sc_transcriptome.index/gene.properties/dNdS/mouse_rat.dnds.txt')
head(filter_mouse_rat)
dim(filter_mouse_rat) #20697 genes
dim(id.mapping) #20449
sum(filter_mouse_rat$ensembl_gene_id %in% id.mapping$ensembl_gene_id) #20697

mouse_rat_dnds=filter_mouse_rat[!duplicated(filter_mouse_rat$ensembl_gene_id),]
head(mouse_rat_dnds)
head(id.mapping)
gene.meta=merge(mouse_rat_dnds,id.mapping,by.x='ensembl_gene_id',by.y='ensembl_gene_id')
dim(gene.meta) #17029
gene.meta$omega=gene.meta$rnorvegicus_homolog_dn/gene.meta$rnorvegicus_homolog_ds
gene.meta[is.infinite(gene.meta$omega),]$omega=NA
gene.meta[is.nan(gene.meta$omega),]$omega=NA
summary(gene.meta$omega) 

#remove genes with dn/ds>1 (more likely to be under positive selection)
length(which(gene.meta$omega<1)) #16192 genes
gene.meta=gene.meta[which(gene.meta$omega<1),]
dim(gene.meta) #16192    12


################################################################
## for each tc, plot gene cor.value ~ dn/ds
x=sp.cors[gene.meta$mgi_symbol]
y=gene.meta$omega
summary(x);summary(y)
summary(lm(y~x))
#plot(x,y,log='y',pch=16,cex=0.3,xlab='spearman.cor',ylab='dn/ds')
tmp=data.frame(gene=gene.meta$mgi_symbol,sp.cors=x,omega=y)
tmp=tmp[!is.na(tmp$sp.cors),]
tmp=tmp[order(tmp$sp.cors),]
ggplot(tmp,aes(x=sp.cors,y=omega))+geom_jitter(size=0.2)+
  scale_y_log10()+theme_classic()+
  xlab('spearman.cor')+ylab('dn/ds')+ggtitle(tc.name)+
  geom_smooth(method = "lm", se=FALSE, color="black", formula = y ~ x) +
  stat_cor(label.x=0,label.y = 0.8)+ #this means at 35th unit in the y axis, the r squared and p value will be shown
  stat_regline_equation(label.x=0,label.y = 1) #this means at 30th unit regresion line equation will be shown

# shuffle x and y
tmp$omega.shuffle=sample(tmp$omega,nrow(tmp),replace = F)
ggplot(tmp,aes(x=sp.cors,y=omega.shuffle))+geom_jitter(size=0.2)+
  scale_y_log10()+theme_classic()+
  xlab('spearman.cor')+ylab('dn/ds')+ggtitle(tc.name)+
  geom_smooth(method = "lm", se=FALSE, color="black", formula = y ~ x) +
  stat_cor(label.x=0,label.y = 0.8)+ #this means at 35th unit in the y axis, the r squared and p value will be shown
  stat_regline_equation(label.x=0,label.y = 1) #this means at 30th unit regresion line equation will be shown

