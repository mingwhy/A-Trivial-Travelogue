library(viridis)

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

names(module)=paste0('module',1:length(module))
module.go=lapply(1:length(module.go),function(i){
  if(length(module.go[[i]]@result$Description)==0) return('')
  module.go[[i]]@result$Description
})
names(module.go)=paste0('module',1:length(module))


###############################
module.size=sapply(module,length)
genes=unlist(module)
length(genes);length(unique(genes))
sum(duplicated(unlist(genes)))
plot.mat=strength[genes,genes]
#plot.mat=uniformity[genes,genes]

max(plot.mat,na.rm=T)
#my.br=seq(0,max(plot.mat,na.rm=T),0.01)
#my.br=c(0,0.02,0.05,seq(0.1,max(plot.mat,na.rm=T),0.2))
my.br=c(0,0.02,0.05,0.2,0.4,0.6,0.8,1)
#if(my.br[length(my.br)]!=max(plot.mat,na.rm=T)){my.br[length(my.br)]=max(plot.mat,na.rm=T)}
my.br
#my.col=viridis(length(my.br)-1)
#my.col = colorRampPalette(c("white","navy", "firebrick3"))(length(my.br)-1)
#my.col = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 10, name = "RdYlBu")))(length(my.br)-1)
#barplot(1:10,col=rev(RColorBrewer::brewer.pal(n = 10, name = "RdYlBu")))
#my.col = colorRampPalette((RColorBrewer::brewer.pal(n = 8, name = "YlOrRd")))(length(my.br)-1)
my.col = colorRampPalette((RColorBrewer::brewer.pal(n = 7, name = "Reds")))(length(my.br)-1)

#my.col[1]='grey60'
my.col[1]='grey90'
length(my.br);length(my.col)
#https://stat.ethz.ch/pipermail/r-help/2006-July/108765.html
axis.tick=c(1,cumsum(sapply(module,length)))
x=1:nrow(plot.mat) 
y=1:ncol(plot.mat)

pdf('T3_K17_sn_heatmap_module.pdf',useDingbats = T)
#pdf('T3_K17_sn_allGenes_heatmap_module.pdf',useDingbats = T)
#pdf('T3_K25_sc_heatmap_module.pdf',useDingbats = T)
#pdf('T3_K25_sc_allGenes_heatmap_module.pdf',useDingbats = T)
image(y,x,plot.mat, breaks=my.br,axes=FALSE, col=my.col,ylab="", xlab="")
#axis(1, at = axis.tick, labels=axis.tick,srt=45,cex.axis=0.4,tick=FALSE)
#axis(2, at = axis.tick, labels=axis.tick,srt=45,cex.axis=0.4,tick=FALSE)
for(i in 1:length(axis.tick)){
  #segments(axis.tick[i],min(y),axis.tick[i],max(y))#vertical
  #segments(min(x),axis.tick[i],max(x),axis.tick[i])#horizental
  rect(axis.tick[i],axis.tick[i],
       axis.tick[i]+module.size[i],axis.tick[i]+module.size[i])
}
box()
#abline(v=axis.tick);abline(h=axis.tick)
dev.off()

image(5, my.br, t(1:length(my.br)),col=my.col,axes=FALSE)
#axis(4,my.br,cex.axis=1)
axis(4,col=NA,
  at=my.br,line=-1, #label closer to axis
  labels=my.br, tck=0, cex.axis=0.6, srt=0,
  col.ticks =NULL)

sapply(module.go,function(i) head(i,5))


##############################################
## plot ave.module.strength VS ave.module.uni
out=lapply(module,function(i){
  x1=strength[i,i]
  x2=uniformity[i,i]
  ngene=length(i)
  return(c(mean(x1[upper.tri(x1)]),mean(x2[upper.tri(x2)]),ngene))
})
out.df=as.data.frame(Reduce(`rbind`,out))
colnames(out.df)=c('strength','uniformity','ngene')
rownames(out.df)=paste0('module',1:nrow(out.df))
out.df=out.df[order(out.df$strength,decreasing = T),]
out.df[,4]=1:nrow(out.df)
length(module) #21

cor(out.df$strength,out.df$uniformity)

sink('sn_module.strength_uni.txt');
#sink('sc_module.strength_uni.txt');
out.df
sink()


#write.table(out.df,paste0('sn_module_strength_uni.txt'),quote=F)

par(mar=c(5,5,4,3))
plot(out.df[,1],out.df[,2],pch=16,cex.lab=1.5,log='xy',
     xlab='Mean strength per module',ylab='Mean uniformity per module')
library(ggplot2)

pdf('sn_strength_uni.pdf',useDingbats = T,height =4,width =4.5)
#pdf('sc_strength_uni.pdf',useDingbats = T,height =4,width =4.5)

print(
ggplot(out.df,aes(x=strength,y=uniformity,color=factor(ngene)))+
  geom_point(size=5)+theme_classic(base_size=14)+
  scale_x_continuous(trans='log10') +
  scale_y_continuous(trans='log10')+
  #scale_color_gradientn(colours = rev(terrain.colors(7)))+
  #scale_color_viridis()+
  scale_color_brewer(name='#gene',type='div')+
  xlab('Mean strength per module')+
  ylab('Mean uniformity per module')+
  theme(legend.position = "top")
)
dev.off();


out.df=out.df[order(out.df[,2],decreasing = T),]
out.df[1:2,]
head(module.go[[out.df[1,3]]])
module[[out.df[1,3]]]

head(module.go[[out.df[2,3]]])
module[[out.df[2,3]]]

head(module.go[[out.df[3,3]]])
module[[out.df[3,3]]]

###############################
## reorder, get module cluster
if(F){
  tmp=diag(0,nrow=length(module))
  for(i in 1:(length(module)-1)){
    for(j in (i+1):length(module)){
      m1=head(module.go[[i]],5);
      m2=head(module.go[[j]],5);
      over=sum(m1 %in% m2)
      tmp[i,j]=tmp[j,i]=1-over/5
    }
  }
  out=hclust(as.dist(tmp))
  i=out$order
  
  module.go<-module.go[i]
  module<-module[i]
}
