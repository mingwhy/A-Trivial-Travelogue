
#genes=readRDS('./brain_snRNA-seq_n10c0.1/all_common_genes.rds')
genes=readRDS('./brain_scRNA-seq/all_common_genes.rds')
#gene.names=genes$all.genes
c.genes=unlist(genes$genes)
gene.names=names(which(table(c.genes)>=length(genes$genes)/2))
length(gene.names);

## two feature matrix
library(R.matlab)
#strength=readMat('./sn_useall_n10c0.1_T3/strength_common.mat');
#uniformity=readMat('./sn_useall_n10c0.1_T3/uniformity_common.mat');
strength=readMat('./sc_useall_n10c0.1_T3/strength_common.mat');
uniformity=readMat('./sc_useall_n10c0.1_T3/uniformity_common.mat');
strength=strength$strength
uniformity=uniformity$uniformity
dim(strength)
colnames(strength)=rownames(strength)=gene.names
colnames(uniformity)=rownames(uniformity)=gene.names

## module info
#module=readRDS('sn_useall_n10c0.1_T3/K17_module.genes.rds')
#module.go=readRDS('sn_useall_n10c0.1_T3/K17_modules.GO.rds')
module=readRDS('sc_useall_n10c0.1_T3/k25_module.genes.rds')
module.go=readRDS('sc_useall_n10c0.1_T3/K25_modules.GO.rds')
sapply(1:length(module.go),function(i) length(module.go[[i]]@result$Description))
# some module has no GO annotation

# 41 overlap
names(which(table(c(module[[1]],module[[4]]))==2))
##########################################

library(igraph)
sapply(module,length)

#pdf('sc_4module_igraph.pdf',useDingbats = T,height = 9,width = 9)
pdf('sn_4module_igraph.pdf',useDingbats = T,height = 9,width = 9)

for(i in 1:length(module)){
  pick.genes=module[[i]]
  plot.mat=strength[pick.genes,pick.genes]
  dim(plot.mat)
  mode(plot.mat) #numeric
  net <- graph_from_adjacency_matrix(plot.mat,mode="lower",weighted = TRUE,diag = FALSE) 
  V(net)
  E(net)$weight
  
  # calcualte sum(abs(PCC)) per node
  node.weight=apply(plot.mat,1,sum)
  summary(node.weight)
  my.br=c(0,0.1,0.3,0.5,0.7,0.9,1)
  # scaled between 1 and 2
  scaled <- 1 + ((2-1) * (node.weight - min(node.weight) ) / (  max(node.weight) - min(node.weight) ) )
  legend.node.size=as.numeric(quantile(scaled,my.br));
  sizeCut=cut(node.weight,breaks=quantile(node.weight,my.br),include.lowest = T)
  sum(table(sizeCut))==length(node.weight)
  V(net)$hubness=sizeCut
  
  my.pal=rev(RColorBrewer::brewer.pal(length(levels(sizeCut)),"RdYlBu"))
  mycols=my.pal[sizeCut]
  
  ## edge color
  edge.size=E(net)$weight;
  summary(edge.size)
  edge.scaled <- 1 + ((2-1) * (edge.size - min(edge.size))/(max(edge.size) - min(edge.size) ) )
  #my.br2=c(0,0.02,0.05,0.2,0.4,0.6,0.8,1)
  my.br2=c(0,0.8,1)
  my.br2.values=as.numeric(quantile(edge.size,my.br2))
  edgeCut=cut(edge.size,breaks=my.br2.values,include.lowest = T)
  sum(table(edgeCut))==length(edge.size)
  E(net)$weighColor=edgeCut
  
  colfunc <- colorRampPalette(c("grey60", "grey30"))
  my.pal2=colfunc(length(levels(edgeCut)))
  my.pal2[1]='grey80';
  #barplot(1:length(my.pal2),col=my.pal2)
  mycols2=my.pal2[edgeCut]
  if(i==2){
    #par(mar=c(15,5,15,10))
    par(mar=c(5,4,4,2))
  }else{
    par(mar=c(5,4,4,2))
  }
  #pdf('sn_module2_phototransduction.pdf',useDingbats = T)
  plot(net,
       #vertex.size=V(net)$hubness,
       vertex.size=scaled*5,
       vertex.color = mycols,
       #layout=layout.fruchterman.reingold,
       layout=  layout_with_fr,
       #rescale=T,
       #asp=0,
       vertex.label.cex=0.5,  #2.5
       vertex.label.family="Helvetica",
       vertex.label.font=2,
       #vertex.label=t.names,
       vertex.shape="circle", 
       vertex.color="deepskyblue2",
       vertex.label.color="black",
       edge.color=mycols2,
       edge.width=0.5, #50
  )
  
  legend('topleft',legend=levels(sizeCut),pt.cex=legend.node.size,col='black',pch=21, pt.bg=my.pal)
  
  #legend('topright',legend=levels(edgeCut),col=my.pal2,lwd=3)
}
dev.off()


#https://stackoverflow.com/questions/49171958/igraph-edge-width-and-color-positive-and-negative-correlation
# edge color legend
image(5, my.br2, t(1:length(my.br2)),col=my.pal2,axes=FALSE)
axis(4,col=NA,
     at=my.br2,line=-1, #label closer to axis
     labels=round(my.br2.values,3), tck=0, cex.axis=0.6, srt=0,
     col.ticks =NULL)

