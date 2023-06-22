# plot cor transformation
x=seq(0,1,0.001);
y=1/(1+exp(log(9999)-2*log(9999)*x)); 
plot(x,y)
#dat <- matrix(runif(10000), ncol = 100)
#image(dat, breaks = c(0.0, 0.8, 1.0), col = c("yellow", "red"))
#https://jokergoo.github.io/2020/05/06/translate-from-pheatmap-to-complexheatmap/

library(viridis)
#x=readRDS('./brain_snRNA-seq_n10c0.1/pearson_coexpr/ncell_204_female_transmedullary Y neuron TmY8_pearson.rds');
x=readRDS('./brain_scRNA-seq/pearson_coexpr/ncell_1005_male_5_pearson.rds');
#cor.mat=abs(x$cell.type)
cor.mat=abs(x);
diag(cor.mat)=NA
my.br=seq(0,1,0.01)
#my.col=viridis(length(my.br)-1)
#my.col = colorRampPalette(c("navy", "white", "firebrick3"))(length(my.br)-1)
my.col = colorRampPalette(rev(brewer.pal(n = 10, name = "RdYlBu")))(length(my.br)-1)
barplot(1:length(my.col),col=my.col)
cor.mat[1:10,1:10]
image(cor.mat[1:10,1:10], breaks = my.br, col = my.col)
#image(cor.mat, breaks = my.br, col = my.col)


#genes=readRDS('./brain_snRNA-seq_n10c0.1/all_common_genes.rds')
genes=readRDS('./brain_scRNA-seq/all_common_genes.rds')
c.genes=unlist(genes$genes)
gene.names=names(which(table(c.genes)>=length(genes$genes)/2))
length(gene.names);dim(cor.mat)
i=intersect(gene.names,rownames(cor.mat));
cor.mat=cor.mat[i,i]
gene.names=i


#module=readRDS('sn_useall_n10c0.1_T3/K17_module.genes.rds')
module=readRDS('sc_useall_n10c0.1_T3/K25_module.genes.rds')
genes=lapply(module,function(i) i[i %in% gene.names])
sum(duplicated(unlist(genes)))
genes2=gene.names[!gene.names %in% unlist(genes)]
reorder.genes=c(unlist(genes),genes2)
cor.mat2=cor.mat[reorder.genes,reorder.genes]

plot.gene.n=min(length(unlist(genes))*2,nrow(cor.mat2))
plot.mat=cor.mat2[1:plot.gene.n,1:plot.gene.n];

my.br=seq(0,max(plot.mat,na.rm=T),0.01)
#my.col=viridis(length(my.br)-1)
#my.col = colorRampPalette(c("navy", "white", "firebrick3"))(length(my.br)-1)
my.col = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 10, name = "RdYlBu")))(length(my.br)-1)

pdf('visual_sc_cor.mat.pdf',useDingbats = T)
image(plot.mat, breaks=my.br, col=my.col,
      ylab="", xlab="")
#https://stackoverflow.com/questions/15188176/adding-a-color-legend-to-an-image
#legend(grconvertX(100, "device"), grconvertY(5, "device"), as.character(my.br[-1]), fill = my.col, xpd = NA) 
dev.off();

sapply(genes,length)
