
genes=readRDS('./brain_scRNA-seq/all_common_genes.rds')
c.genes=unlist(genes$genes)
gene.names1=names(which(table(c.genes)>=length(genes$genes)/2))

genes=readRDS('./brain_snRNA-seq_n10c0.1/all_common_genes.rds')
c.genes=unlist(genes$genes)
gene.names2=names(which(table(c.genes)>=length(genes$genes)/2))

length(gene.names1); #3015
length(gene.names2); #1351
length(intersect(gene.names1, gene.names2))#1180

#m1=readRDS('sc_useall_n10c0.1_T3/K25_common.modules.genes.rds')
#m2=readRDS('sn_useall_n10c0.1_T3/K17_common.modules.genes.rds')
m1=readRDS('sc_useall_n10c0.1_T3/K25_module.genes.rds')
m2=readRDS('sn_useall_n10c0.1_T3/K17_module.genes.rds')
#m2=readRDS('sn_useall_neuron.glia_T3/K20_common.modules.genes.rds')

df.m1=data.frame(gene.symbol=unlist(m1),
  module.id=paste('module',rep(1:length(m1),sapply(m1,length))))
write.table(df.m1,'sc_genes.in.modules.txt',sep='\t',quote=F,row.names = F)
df.m2=data.frame(gene.symbol=unlist(m2),
                 module.id=paste('module',rep(1:length(m2),sapply(m2,length))))
write.table(df.m2,'sn_genes.in.modules.txt',sep='\t',quote=F,row.names = F)


m1=m2

length(m1) #4
length(m2) #4
#write.table(module[[4]],file='tmp.txt',quote=F,row.names = F,col.names = F)

gene1=unique(unlist(m1))
gene2=unique(unlist(m2))
length(gene1) #390
length(gene2) #204
length(intersect(gene1,gene2)) #113

m1.go=readRDS('sc_useall_n10c0.1_T3/K25_modules.GO.rds');
m2.go=readRDS('sn_useall_n10c0.1_T3/K17_modules.GO.rds');
#m2.go=readRDS('sn_useall_neuron.glia_T3/K20_common.modules.GO.rds');

#head(m1.go[[10]]@result$Description)

names(m1)=paste0('module',1:length(m1),',',sapply(m1,length),'g')
names(m2)=paste0('module',1:length(m2),',',sapply(m2,length),'g')
mat=matrix(0,nrow=length(m1),ncol=length(m2))
mat.count=matrix(0,nrow=length(m1),ncol=length(m2))
rownames(mat)=names(m1)
colnames(mat)=names(m2)
for(i in 1:length(m1)){
  for(j in 1:length(m2)){
    genes1=m1[[i]];
    genes2=m2[[j]]
    o=intersect(genes1,genes2);
    u=union(genes1,genes2)
    #jar=length(o)/length(u)
    jar=length(o)
    mat.count[i,j]=jar
    jar=length(o)/min(length(genes2),length(genes1))
    mat[i,j]=jar
    #mat[i,j]=length(o)
  }
}

#pdf('T3_sc_sn.pdf',useDingbats = T)
#pdf('T3_sc_sc.pdf',useDingbats = T)
pdf('T3_sn_sn.pdf',useDingbats = T,height = 4,width = 5)
print( pheatmap::pheatmap(mat, display_numbers = mat.count,
                          fontsize_number=14,fontsize_row=14,
                          #display_numbers =round(mat,3),
       main='single nuclei (column) vs single cell (row) ')
       #main='sn (column) vs sn (row) ')
       #main='sc (column) vs sc (row) ')
)
dev.off()

## check those common.modules detected in sc but not in sn
(i=which(Matrix::colSums(mat)<0.2))
(i=which(Matrix::rowSums(mat)<0.2))
(i=which(Matrix::rowSums(mat)<0.05))
m1[i]
lapply(i,function(k){
  head(m1.go[[k]]@result$Description,5)
})

sink('GO_terms.txt')
lapply(1:length(m1.go),function(k) head(m1.go[[k]]@result$Description,5))
lapply(1:length(m2.go),function(k) head(m2.go[[k]]@result$Description,5))
sink()
