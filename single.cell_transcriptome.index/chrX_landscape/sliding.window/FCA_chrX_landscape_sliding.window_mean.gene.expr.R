
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


# chrX: set up sliding window (#gene per window)
chrX=gene.meta[gene.meta$LOCATION_ARM=='X',]
chrX=chrX[order(chrX$LOCATION_MIN),]
min(chrX$LOCATION_MIN)
max(chrX$LOCATION_MAX)

size=50; step=25; #50 genes, step = 25genes
(ngene=nrow(chrX)); #2488
window.start=seq(1,ngene,25)
window.end=window.start+49;
length(window.start);length(window.end)
window.start=window.start[1:99]
window.end=window.end[1:99]
if(window.end[length(window.end)]>ngene){
  window.end[length(window.end)]=ngene
}
chrX.window.start=window.start
chrX.window.end=window.end


# chr2L: set up sliding window (#gene per window)
chr2L=gene.meta[gene.meta$LOCATION_ARM=='2L',]
chr2L=chr2L[order(chr2L$LOCATION_MIN),]
size=50; step=25; #50 genes, step = 25genes
(ngene=nrow(chr2L)); # 3155
chr2L.window.start=seq(1,ngene,25)
chr2L.window.end=chr2L.window.start+49;
length(chr2L.window.start);length(chr2L.window.end)
chr2L.window.start=chr2L.window.start[1:126]
chr2L.window.end=chr2L.window.end[1:126]

chr2L.window.end[length(chr2L.window.end)]=ngene

chr2L.window.start
chr2L.window.end

#####################################################################
## read in FCA data, average gene expr per tissue
folders=folders=Sys.glob('~/Documents/Data_fly_FCA/*FCA*')
out.path='gene.mean.var_per.tc';
dir.create(out.path)

for(folder in folders){
  tissue=strsplit(basename(folder),'\\_')[[1]][2]
  input.file=Sys.glob(paste0(folder,'/*valid.rds'))
  cat('begin',tissue,input.file,'\n')
  
  ## read in data
  #tissue='gut'
  #input.file = "../../single.cell_datasets/FCA_gut/whole_gut_filtered_valid.rds"
  sc = readRDS(input.file)
  dim(sc)
  
  # antenna and body has female, male, mix, three class labels
  # only keep cell types with both female and male cells profiled
  sc=sc[,sc$sex!='mix' & sc$annotation!='unannotated']
  sc$tc_sex=paste(sc$annotation,sc$sex,sep=';')
  
  sc <- NormalizeData(sc, normalization.method = "LogNormalize", scale.factor = 1e6,verbose = F) #ln(CPM+1)
  cat(tissue,'tc_sex #',length(unique(sc$tc_sex)),'\n');
  
  # mean gene expr across cells per sex per tc
  out=lapply(unique(sc$tc_sex),function(i){
    x=sc[,sc$tc_sex==i];
    x.mean=Matrix::rowMeans(x@assays$RNA@data)
    x.median=apply(x@assays$RNA@data,1,median)
    x.var=apply(x@assays$RNA@data,1,var)
    x=data.frame(SYMBOL=names(x.mean),n.cell=ncol(x),
                 mean=x.mean,median=x.median,var=x.var,tc_sex=i)
    merge(x,gene.meta)
  })
  length(out)
  
  #ggplot(out[[1]],aes(x=LOCATION_ARM,y=var^0.5/mean))+geom_violin()+geom_jitter(size=0.2)+theme_classic()
  out.file=paste0(out.path,'/',tissue,'_gene.var.rds')
  saveRDS(out,out.file)
  cat(tissue,'is done\n')
}

#####################################################################
## add chr position per gene for sliding window
files=Sys.glob('gene.mean.var_per.tc/*.rds')
(files=files[c(7,17,4)]) #head, fat.body, wing

pick.chr=chr2L
window.start=chr2L.window.start
window.end=chr2L.window.end
for(file in files){
  tissue=gsub('_gene.var.rds','',basename(file))
  x=readRDS(file)
  x=as.data.frame(Reduce(`rbind`,x))
  x=x[x$LOCATION_ARM %in% c('X','2L','2R','3L','3R'),]
  x$tc=unlist(lapply(strsplit(x$tc_sex,';'),'[',1))
  x$sex=unlist(lapply(strsplit(x$tc_sex,';'),'[',2))
  #table(x$sex);#table(x$tc);
  
  #pdf(paste0('FCA_chrX_',tissue,'.pdf'),useDingbats = T,width=16,height = 12)
  pdf(paste0('FCA_chr2L_',tissue,'.pdf'),useDingbats = T,width=16,height = 12)
  
  for(tc in unique(x$tc)){
    tc.gene=x[x$tc==tc,]
    if(length(unique(tc.gene$sex))<2){next}
    if(sum(tc.gene$n.cell<20)>0){next} #per tc_sex, >=20 cells in computing gene.mean.expr
    windows=lapply(1:length(window.start),function(i){
      genes=pick.chr[window.start[i]:window.end[i],]$SYMBOL
      window1=tc.gene[tc.gene$SYMBOL %in% genes,]
      window1 %>% filter(mean>0) %>% group_by(sex) %>% 
        summarise(gene.mean=mean(mean),ngene=n()) %>% mutate(window.index=i)
    })
    df.windows=as.data.frame(Reduce(`rbind`,windows))
    
    p1<- ggplot(df.windows,aes(x=window.index,y=gene.mean,col=sex,group=sex))+
             geom_point(size=0.2,pch=21)+geom_line()+
             #geom_point(size=0.2,pch=21, position = position_jitterdodge())+
             theme_classic()+theme(legend.position="top")+ggtitle(paste0(tissue,'\n',tc))
    
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


