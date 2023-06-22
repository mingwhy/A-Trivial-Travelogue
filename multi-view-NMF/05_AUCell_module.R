library(AUCell)
library(Seurat)
library(ggplot2);library(gridExtra);
options(stringsAsFactors = F)

top=0.4;
## read in genes
input.folder='./brain_scRNA-seq/';
#input.folder='./brain_snRNA-seq_n10c0.1/';

outfile=paste0('AUCell',top,'_sc_useall_n10c0.1_T3_K25.rds');
#outfile=paste0('AUCell',top,'_sn_useall_n10c0.1_T3_K25.rds');

## read in modules 
modules=readRDS('sc_useall_n10c0.1_T3/K25_module.genes.rds');
#modules=readRDS('sn_useall_n10c0.1_T3/K17_module.genes.rds')
names(modules)=paste0('module',1:length(modules))

#modules.go=readRDS('sc_useall_n10c0.1_T3/K25_common.modules.GO.rds')
#modules.go=readRDS('sn_useall_n10c0.1_T3/K17_common.modules.GO.rds')
#head(modules.go[[3]]@result$Description)

## read in genes
genes=readRDS(paste0(input.folder,'all_common_genes.rds'))
c.genes=unlist(genes$genes)
all.genes=genes$all.genes
common.genes=genes$common.genes
#gene.names=genes$all.genes
gene.names=names(which(table(c.genes)>=length(genes$genes)/2))
gene_count=length(gene.names) #3015 for sc, 1351 for sn

## read in processed wholebrain data
if(length(grep('sn',outfile))==1){
  file="../single.cell_datasets/FCA_head/whole_head_filtered_valid.rds";
}else{
  file="../single.cell_datasets/fly.brain.atlas/wholebrain_filtered_valid.rds";
}

dat=readRDS(file);
df.meta=dat@meta.data;
#dat=NormalizeData(dat);
#df.umi=dat@assays$RNA@counts #12094 features across 56192 samples within 1 assay 
df.expr=dat@assays$RNA@data #logNormal

unique(dat@meta.data$sex)
table(dat@meta.data$annotation)

## separate female and male, select cell clusters that contain>=500cells in both sexes
min.cell=100;max.cell=Inf;
i=apply(table(dat$sex,dat$annotation),2,function(i)sum(i>=min.cell & i<max.cell))
pick.cell.clusters=names(which(i==2)) #79
pick.cell.clusters #sn,50; sc,60.

#df.expr1=df.expr[,dat$annotation %in% pick.cell.clusters]
#sparsity=sum(df.expr1==0)/nrow(df.expr1)/ncol(df.expr1)
#cat(outfile,sparsity,'\n') #sn=0.9537,sc=0.9055
#rm('df.expr1')
if(!file.exists(outfile)){
  ## effective.expressed.gene >=1umi in max(10,ncell*10%) cells
  out=list();
  for(i.sex in c('male','female')){
    for(i.cluster in pick.cell.clusters){
      mat=df.expr[,dat$sex==i.sex & dat$annotation==i.cluster]
      #gene.filter=Matrix::rowSums(mat>0) > max(10,ncol(mat)*0.1)
      #mat=as.matrix(mat[gene.filter,])
      gene.set=intersect(gene.names,rownames(mat))
      mat=mat[gene.set,]
      
      cells_rankings<-AUCell::AUCell_buildRankings(mat,nCores = 2,plotStats = F)
      #cells_AUC <- AUCell_calcAUC(modules, cells_rankings,aucMaxRank=nrow(mat)*0.05,nCores = 1)
      cells_AUC <- AUCell_calcAUC(modules, cells_rankings,aucMaxRank=nrow(mat)*top,nCores = 2)
      
      cells_AUC_matrix <- getAUC(cells_AUC)
  
      #i.cluster2=gsub('\\/','\\.',i.cluster)#there is _ and -, use . to replace /
      out[[paste0(i.sex,'__',i.cluster)]]=cells_AUC_matrix;
    }
  }
  saveRDS(out,outfile)
}

out=readRDS(outfile)
length(out) # #cell.cluster

## pick one cell cluster and plot AUCell by age
tmp=sort(sapply(out,length),decreasing = T)
(cluster.name=names(tmp[15]))#which cell cluster has most cells
pick.cluster=out[[cluster.name]]

pick.female=c('female__0','female__Ensheathing_glia','female__TmY14')
pick.male=c(gsub('female','male',pick.female));
pick=c(pick.female,pick.male);

plot_list=list();
for(cluster.name in pick[c(1,4,2,5,3,6)]){
  pick.cluster=out[[cluster.name]];
  df.meta.sub=df.meta[colnames(pick.cluster),]
  cat(cluster.name,'\n')
  print(table(df.meta.sub$Age));
  x=sapply(sort(unique(df.meta.sub$Age)),function(i){
    cells=rownames(df.meta.sub[df.meta.sub$Age==i,])
    apply(pick.cluster[,cells],1,mean)})
  x=t(x)
  rownames(x)=sort(unique(df.meta.sub$Age))
  
  x=pheatmap::pheatmap(x,cluster_rows=FALSE,
                     display_numbers = round(x,2),
                     fontsize_number=12,fontsize_row=14,
                     fontsize_column=14,
                     main=paste0('cell.cluster: ',cluster.name),
                     clustering_method='ward.D2') 
  plot_list[[cluster.name]]=x[[4]]
}
#https://www.biostars.org/p/128229/
pdf('age_aucell.pdf',useDingbats = T,height = 16)
grid.arrange(arrangeGrob(grobs= plot_list,ncol=2))
dev.off()

## plot module activity score by cell cluster
out=readRDS(outfile)
#apply(out[[3]],1,mean) #small perc ->small AUC values
#apply(out[[3]],1,median) #small perc ->small AUC values
x=sapply(1:length(out),function(i) apply(out[[i]],1,mean))
x=t(x)
rownames(x)=names(out)

if(length(grep('sn',outfile))==1){
  pdf(gsub('.rds','.pdf',outfile),useDingbats = T,height = 16)
}else{
  pdf(gsub('.rds','.pdf',outfile),useDingbats = T,height = 16)
}
print( pheatmap::pheatmap(x,clustering_method='ward.D2') )
dev.off()

#https://www.flyrnai.org/tools/single_cell/web/show_markers_table/60/70/NULL
# pick one method, whose clustering method agress the most with cell type clustering restuls
# pick some metric measuring two clustering results
# adjusted-ARI,聚类性能评估-ARI（调兰德指数)
# NMI (Normalized Mutual Information): mutual information

