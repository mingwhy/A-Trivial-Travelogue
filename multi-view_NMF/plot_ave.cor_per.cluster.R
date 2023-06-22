# plot cor transformation
#x=seq(0,1,0.001);
#y=1/(1+exp(log(9999)-2*log(9999)*x)); 
#plot(x,y)

net_transform<-function(mat){
  weight_max = max(mat);
  weight_min = 0;
  matrix_weight_1 = (mat- weight_min) / (weight_max-weight_min);
  matrix_weight_2 = 1/(1+exp(log(9999)-2*log(9999)*matrix_weight_1));    
  #matrix_weight_2[matrix_weight_1 <= 0.3] = 0;
 matrix_weight_2
}

## read in genes
genes=readRDS('brain_scRNA-seq/all_common_genes.rds')
c.genes=unlist(genes$genes)
all.genes=genes$all.genes
common.genes=genes$common.genes
#gene.names=genes$all.genes
gene.names=names(which(table(c.genes)>=length(genes$genes)/2))
gene_count=length(gene.names) #3015 for sc, 1351 for sn

#module=readRDS('sn_useall_n10c0.1_T3/K17_module.genes.rds')
module=readRDS('sc_useall_n10c0.1_T3/K25_module.genes.rds')
sapply(module,length)
files=Sys.glob('./brain_scRNA-seq/pearson_coexpr/*_pearson.rds')
cell.cluster.names=as.character(unlist(sapply(files,function(i){
  gsub('ncell_\\d+_','',basename(i),'_pearson.mat')
})))
cell.cluster.names=gsub('.rds|_pearson','',cell.cluster.names)

module.list=list();
j=1;
genes=module[[j]]
tmp=matrix(0,length(genes),length(genes))
rownames(tmp)=colnames(tmp)=genes
for(i in 1:length(files)){
  mat0=readRDS(files[[i]])
  mat=mat0
  #mat=abs(mat0)
  #mat=net_transform(mat0)
  x=intersect(genes,rownames(mat))
  tmp1=tmp;
  tmp1[x,x]=mat[x,x]
  module.list[[i]]<-tmp1;
}

MESS::cmd(module.list[[1]],module.list[[2]])
out=matrix(0,nrow=length(module.list),ncol=length(module.list))
for(i in 1:(length(module.list)-1)){
  for(j in (i+1):length(module.list)){
    out[i,j]=MESS::cmd(module.list[[i]],module.list[[j]])
    out[j,i]=out[i,j]
  }
}
colnames(out)=rownames(out)=cell.cluster.names
pdf('test.pdf',useDingbats = T,height = 16,width = 16)
pheatmap::pheatmap(out)
dev.off()

#############################
out=matrix(0,nrow=length(files),ncol=length(module))
for(i in 1:length(files)){
  
  mat=readRDS(files[[i]])
  x=intersect(gene.names,rownames(mat))
  mat=mat[x,x]
  mat2=net_transform(mat)
  for(j in 1:length(module)){
    genes=module[[j]]
    genes2=intersect(genes,rownames(mat))
    x=mat2[genes2,genes2]
    out[i,j]=mean(x[upper.tri(x)])
  }
}

colnames(out)=paste0('module',1:ncol(out))
rownames(out)=cell.cluster.names
dev.off()
pdf('test.pdf',useDingbats = T,height = 16)
print( pheatmap::pheatmap(out) )
dev.off()

