
library(ggplot2);library(gridExtra)
library(dplyr);library(ggpubr)
library(Seurat);
library(grDevices);library(RColorBrewer)
library(transport)
opt.dist<-function(x,y,n.repeat=300,perc=0.8){
  n1=nrow(x);n2=nrow(y);
  n.sample.cell=min(min(n1,n2)*perc,100)
  #set.seed(123)
  out=lapply(1:n.repeat,function(i){
    sx <- pp(x[sample(1:n1,n.sample.cell,replace = F),])
    sy <- pp(y[sample(1:n2,n.sample.cell,replace = F),])
    wasserstein(sx,sy,p=1)
    #wasserstein(sx,sy,p=2)
  })
  return(unlist(out))
}

###############################################################################################
#nmf.out=readRDS('nmf_marrow.heart.rds')
#cell.meta=readRDS('nmf.cellmeta_marrow.heart.rds')
nmf.out=readRDS('nmf_115tc.rds')
cell.meta=readRDS('nmf.cellmeta_115tc.rds')
length(nmf.out);length(cell.meta)

all.OT.distances=list();
for(tc in names(nmf.out)){
  nmf_model=nmf.out[[tc]]
  cell_metadata=cell.meta[[tc]]
  
  nmf.h=t(nmf_model@h)
  colnames(nmf.h) <- paste0("NMF_", 1:choose.k)
  #sum(rownames(nmf.h)==rownames(cell_metadata));dim(nmf.h)
  
  n1=sum(cell_metadata$binary.age=='young')
  n2=sum(cell_metadata$binary.age=='old')
  if(n1<100 | n2<100){next}
  #expr.pc=expression_matrix
  x=nmf.h[cell_metadata$binary.age=='young',]
  y=nmf.h[cell_metadata$binary.age=='old',]
  
  out.xx=opt.dist(x,x,n.repeat=300,perc=0.8)
  out.yy=opt.dist(y,y,n.repeat=300,perc=0.8)
  out.xy=opt.dist(x,y,n.repeat=300,perc=0.8)
  
  #boxplot(out.xx,out.yy,out.xy)
  df=data.frame(group=rep(c('young-young','old-old','young-old'),each=300),OT.dist=c(out.xx,out.yy,out.xy))
  all.OT.distances[[tc]]=df
  cat('tc',tc,'is done\n')
}
length(all.OT.distances) #73
#saveRDS(all.OT.distances,'nmf_marrow.heart_OT.dist.rds')
#saveRDS(all.OT.distances,'nmf_115tc_OT.dist.rds')
saveRDS(all.OT.distances,'nmf_73tc_OT.dist.rds')


#all.OT.distances=readRDS('nmf_marrow.heart_OT.dist.rds')
#all.OT.distances=readRDS('nmf_115tc_OT.dist.rds')
all.OT.distances=readRDS('nmf_73tc_OT.dist.rds')
df=Reduce(`rbind`,all.OT.distances)
df$tc=rep(names(all.OT.distances),sapply(all.OT.distances,nrow))
head(df)

plotcol <- brewer.pal(8,'Dark2')
#pdf('nmf_marrow.heart_OT.dist.pdf',useDingbats = T,width = 24,height=14)
#pdf('nmf_115tc_OT.dist.pdf',useDingbats = T,width = 24,height=14)
pdf('nmf_73tc_OT.dist.pdf',useDingbats = T,width = 24,height=14)
ggplot(df,aes(x=group,y=OT.dist,col=group))+
  scale_color_manual(values=plotcol)+
  facet_wrap(.~tc)+geom_violin()+geom_jitter(size=0.5)+theme_classic()
dev.off()

tc=df[grep('Liver:endo',df$tc),]$tc[1]
ggplot(df[grep('Liver:endo',df$tc),],aes(x=group,y=OT.dist,col=group))+
  scale_color_manual(values=plotcol)+ggtitle(tc)+
  geom_violin()+geom_jitter(size=0.5)+theme_classic()


