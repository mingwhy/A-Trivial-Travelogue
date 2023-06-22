library(Matrix)
readMM('W_matrix.mtx')->W
readMM('H_matrix.mtx')->H
dim(W) #gene by K
dim(H) #K by cell

gene.info=data.table::fread('ndata_gene.csv',sep='\t')
dim(gene.info)
head(gene.info)

cell.info=data.table::fread('ndata_cell.csv',sep='\t')
dim(cell.info) #56875     6
head(cell.info)


## cell type info
df.cell=read.table("../single.cell_datasets/fly.brain.atlas/GSE107451_DGRP-551_w1118_WholeBrain_57k_Metadata.tsv.gz",as.is=T,header=T);
table(df.cell$annotation)
sum(df.cell$annotation=='Hsp') #668 cells annotated to 'Hsp'
# remove them due to personal email communication, they are 'stressed' cells.
dim(df.cell) #56902    31

sum(df.cell$new_barcode %in% cell.info$V1) #56875
df.cell=df.cell[df.cell$new_barcode %in% cell.info$V1,]

#############################################################
## aging vector per cell type per sex
cell.types=unique(df.cell$annotation)
sexes=unique(df.cell$sex)
out=list()
for(cell.type in cell.types){
  for(sex in sexes){
    tmp=df.cell[df.cell$annotation==cell.type & df.cell$sex==sex,]
    #table(tmp$Age )
    young=H[,cell.info$V1 %in% tmp[tmp$Age<=10,]$new_barcode,drop=F]
    old=H[,cell.info$V1 %in% tmp[tmp$Age>10,]$new_barcode,drop=F]
    if(any(ncol(young)<25,ncol(old)<25)){next}
    young.median=apply(young,1,median)
    old.median=apply(old,1,median)
    v=old.median-young.median;
    name=paste0(cell.type,'__',sex)
    out[[name]]=v
  }
}

df.out=as.data.frame(Reduce(`cbind`,out))
colnames(df.out)=names(out)
dim(df.out); #50 x 87

#BiocManager::install('lsa')
library(lsa)
df.cos=cosine(as.matrix(df.out)) #a similarity matrix
dim(df.cos)
pdf('test.pdf',useDingbats = T,height = 20,width = 20)
pheatmap::pheatmap(df.cos)
dev.off()


library(org.Dm.eg.db)
gene.df <- clusterProfiler::bitr(gene.info[W[,2]!=0,]$gene_symbols, fromType = "SYMBOL",
                                 toType = c("ENTREZID","FLYBASE","GENENAME"),
                                 OrgDb = org.Dm.eg.db)
ego <- enrichGO(gene          = gene.df$ENTREZID,
                OrgDb         = org.Dm.eg.db,
                #keyType  = 'SYMBOL',
                #ont           = category,
                ont           = "BP",
                #ont           = "MF", 
                #ont           = "CC", 
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.05,
                readable      = TRUE)
heatplot(ego)
head(ego@result$Description)
