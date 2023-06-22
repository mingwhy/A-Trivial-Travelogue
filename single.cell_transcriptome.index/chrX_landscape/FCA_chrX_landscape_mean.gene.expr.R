
library(tidyverse)
library(Seurat)
library(gridExtra)
#for each tc, female-mean, male-mean
#50  gene on chrX as a window, 25 gene as a step.
#mean(~50 gene female-mean), mean(~50 gene male-mean)
##########################################################
## read in gene meta info for FCA
gene.meta=data.table::fread('~/Documents/Data_fly_FCA/gene.meta.txt')
dim(gene.meta) #16276     11
gene.meta=as.data.frame(gene.meta)
rownames(gene.meta)=gene.meta$current_symbol
dim(gene.meta) #16276    11

gene.meta$LOCATION_MAX=as.numeric(gene.meta$LOCATION_MAX)
gene.meta$LOCATION_MIN=as.numeric(gene.meta$LOCATION_MIN)
gene.meta[is.na(gene.meta$LOCATION_MAX),] #one gene, no information for chr arm

#####################################################################
## add chr position per gene for landscape
all.tcs=readRDS('expr.genes_cellpop.rds') #`FCA_find.expr.genes_cellpop.R`
length(all.tcs) # per tc per sex
names(all.tcs)
tc.names=gsub(';female$|;male$','',names(all.tcs))
pick.tc.names=names(which(table(tc.names)==2)) #235

pdf('test.pdf',useDingbats = T,height = 5,width = 16)
for(tc.name in pick.tc.names){
  female=all.tcs[[which(paste0(tc.name,';female')== names(all.tcs))]]
  male=all.tcs[[which(paste0(tc.name,';male')== names(all.tcs))]]
  
  dim(female);dim(male)
  x=female[,c('SYMBOL','mean','LOCATION_ARM','LOCATION_MIN')]
  y=male[,c('SYMBOL','mean','LOCATION_ARM','LOCATION_MIN')]
  colnames(x)[2]='female_mean'
  colnames(y)[2]='male_mean'
  xy=merge(x,y)
  dim(xy) #commonly expr genes in both F and M
  xy$FMratio=xy$female_mean/xy$male_mean
  xy$LOCATION_MIN=as.numeric(xy$LOCATION_MIN)
 
  
  print( ggplot(xy,aes(x=LOCATION_MIN,y=FMratio))+
            geom_point(size=0.2,pch=21)+#geom_line()+
            facet_wrap(.~LOCATION_ARM,scale='free',nrow=1)+
            theme(axis.text.x = element_blank())+
            scale_y_continuous(limits = c(min(xy$FMratio), max(xy$FMratio)))+
            #geom_point(size=0.2,pch=21, position = position_jitterdodge())+
            theme_classic()+theme(legend.position="top")+ggtitle(tc.name)
  )
}
dev.off()
  

    p2= ggplot(df.windows,aes(x=window.index,y=ngene,col=sex,group=sex))+
             geom_point(size=0.2,pch=21)+geom_line()+
             theme(legend.position="top")+
             #geom_point(size=0.2,pch=21, position = position_jitterdodge())+
             theme_classic()+theme(legend.position="top")+ggtitle(paste0(tissue,'\n',tc))
    
    df.windows.ratios=df.windows %>% group_by(window.index) %>% 
        summarise(ratio = gene.mean[sex=="female"]/gene.mean[sex=="male"])
    #max(df.windows.ratios$ratio)
    
    p3= ggplot(df.windows.ratios,aes(x=window.index,y=log2(ratio)))+
      geom_point(size=0.2,pch=21)+geom_line()+
      #geom_point(size=0.2,pch=21, position = position_jitterdodge())+
      theme_classic()+ggtitle(paste0(tissue,'\n',tc))
  
    print( grid.arrange(p1,p2,p3,nrow=3) )
  }
  dev.off()
}



##################################################################
## plot chrX gene expr for all cells
colnames(df.gene.table)
df.gene.table$LOCATION_MIN=as.numeric(df.gene.table$LOCATION_MIN)
df.gene.table$LOCATION_MAX=as.numeric(df.gene.table$LOCATION_MAX)
dim(df.gene.table)

dat.mat=umi.mat;
#dat.mat=sce.count;
#dat.mat=log(umi.mat+1) #raw count
#dat.mat=log10.sce.1.mat #sce normalized count
dim(dat.mat)
dat.mat=dat.mat[df.gene.table$current_symbol,]

chrA.mat=dat.mat[!df.gene.table$LOCATION_ARM %in% c('X','Y'),]
per.cell.A.median=apply(chrA.mat,2,function(x){ median(x[x!=0]) })

chrX.mat=dat.mat[df.gene.table$LOCATION_ARM=='X',]
dim(chrX.mat) #1679

relative.toA.chrX.mat=sapply(1:ncol(chrX.mat),function(i){ chrX.mat[,i]/per.cell.A.median[i]})

chrX.gene.coords=df.gene.table[df.gene.table$LOCATION_AR=='X',] #1679 X genes
sum(chrX.gene.coords$current_symbol==rownames(relative.toA.chrX.mat)) #1679 X genes

# order along X chr
chrX.gene.coords=chrX.gene.coords[order(chrX.gene.coords$LOCATION_MIN),]
relative.toA.chrX.mat=relative.toA.chrX.mat[chrX.gene.coords$current_symbol,]
sum(chrX.gene.coords$current_symbol == rownames(relative.toA.chrX.mat)) #1679 X genes
#plot(1:nrow(chrX.gene.coords),chrX.gene.coords$LOCATION_MIN)

##
plot.mat=relative.toA.chrX.mat
tmp=Matrix::rowSums(plot.mat)
summary(tmp) #min 4

(y.max=(max(plot.mat[plot.mat!=0]))*1.2)
(y.min=(min(plot.mat[plot.mat!=0]))*1.2)
ymin=0;

pdf('plot_chrX_landscape.pdf',width = 16)
par(mfrow=c(3,1))
plot.ncell=c(20,100,1297);
#plot.ncell=c(20,100);
for(ncell in plot.ncell){
  
  plot(chrX.gene.coords$LOCATION_MIN,
       plot.mat[,1],
       #ylab='log10(normalized.count)',
       ylab='X.genes/median.A.per.cell',
       ylim=c(y.min,y.max),
       type='l',xlab='X chromosome coordinates',
       pch=16,col=rgb(0,0,0,0.2),cex=0.3,main=paste0(ncell,'cells') )
  for(i in 2:ncell){
    points(chrX.gene.coords$LOCATION_MIN,
           log10(ncell[,i]),type='l',
           pch=16,col=rgb(0,0,0,0.2),cex=0.3)
  }
  points(chrX.gene.coords$LOCATION_MIN,rep(y.min,length(chrX.gene.coords$LOCATION_MIN)),
         pch=16,col='red',cex=0.5)

}

dev.off()


