
library(tidyverse)
library(Seurat)

#####################################################################
## read in gene dn/ds info
dnds=as.data.frame(data.table::fread('~/Documents/sc_transcriptome.index/gene.properties/dNdS/Dmel_Dsim_analysis_results_flydivas_v1.2.txt'))
head(dnds)
dnds[duplicated(dnds$id),]

dim(dnds);#10765
str(dnds)
summary(dnds$omega) 
sum(dnds$omega>2,na.rm=T) #16
dnds=dnds[dnds$omega<=2,]
dim(dnds)#10749    10
rownames(dnds)=dnds$id;

#####################################################################
## read in FCA data
## check for sex ncell in each dataset
if(F){
  folders=Sys.glob('~/Documents/Data_fly_FCA/*FCA*')
  for(folder in folders){
    tissue=strsplit(basename(folder),'\\_')[[1]][2]
    input.file=Sys.glob(paste0(folder,'/*valid.rds'))
    cat('begin',tissue,input.file,'\n')
    sc = readRDS(input.file)
    print(table(sc$sex))
  }
}
## body and antenna has: female, male, and mix, three class labels

folders=folders=Sys.glob('~/Documents/Data_fly_FCA/*FCA*')
#folders=folders[c(1,7)] #antenna and body
#folders=folders[c(2,9)]  #gut and fat.body

for(folder in folders){
  tissue=strsplit(basename(folder),'\\_')[[1]][2]
  input.file=Sys.glob(paste0(folder,'/*valid.rds'))
  cat('begin',tissue,input.file,'\n')
  
  ## read in data
  #tissue='gut'
  #input.file = "../../single.cell_datasets/FCA_gut/whole_gut_filtered_valid.rds"
  sc = readRDS(input.file)
  dim(sc)
  
  #table(sc$annotation,sc$sex)
  if(length(table(sc$sex))==1){
    cat('skip: one sex present',tissue,input.file,'\n');
    next
  }
  # antenna and body has female, male, mix, three class labels
  # only keep cell types with both female and male cells profiled
  sc=sc[,sc$sex!='mix' & sc$annotation!='unannotated']
  meta.info=sc@meta.data
  cell.type.names=lapply(unique(sc$annotation),function(i){
    if(length(unique(meta.info[meta.info$annotation==i,]$sex))>1){return(i)}
  })
  cell.type.names=Filter(Negate(is.null), cell.type.names)
  sc=sc[,sc$annotation %in% cell.type.names]
  
  out.file=paste0('fca_',tissue,'_TDI.rds');
  if(!file.exists(out.file)){
    sc=NormalizeData(sc)
    expr.mat=sc@assays$RNA@data
    
    #sum(rownames(expr.mat) %in% dnds$symbol)
    overlap.genes=intersect(rownames(expr.mat),dnds$symbol)
    expr.mat=expr.mat[overlap.genes,]
    tmp=dnds[match(overlap.genes,dnds$symbol),]
    #sum(tmp$symbol==overlap.genes)
    
    x1=as.numeric(tmp$omega %*% expr.mat)
    #sum(tmp$omega * expr.mat[,1]);x1[1]
    x2=as.numeric(Matrix::colSums(expr.mat))
    index.per.cell=x1/x2;
    cell.meta=sc@meta.data;
    cell.meta$TDI=index.per.cell;
    if(F){
      cell.meta$tc_sex=paste(cell.meta$annotation,cell.meta$sex)
      ggplot(cell.meta,aes(x=annotation,y=TDI,col=sex,group=tc_sex))+
        geom_violin()+#geom_jitter(size=0.2)+
        ggtitle(tissue)+
        #stat_summary(fun.y=median, geom="point", size=2, color="red")+
        theme_classic()+
        theme(axis.text.x = element_text(angle=45,vjust=1,hjust=1))
    }
   
    saveRDS(cell.meta,out.file)
  }
  cat('finished',tissue,input.file,'\n')
}


##############################################
## plot
files=Sys.glob('fca_*_TDI.rds')
pdf('FCA_TDI.pdf',useDingbats = T,width = 16 )
for(file in files){
  cell.meta=readRDS(file)
  tissue=strsplit(file,'\\_')[[1]][2]
  
  cell.meta$tc_sex=paste(cell.meta$annotation,cell.meta$sex)
  print( ggplot(cell.meta,aes(x=annotation,y=TDI,col=sex,group=tc_sex))+
    geom_violin()+#geom_jitter(size=0.2)+
    ggtitle(tissue)+
    #stat_summary(fun.y=median, geom="point", size=2, color="red")+
    theme_classic()+
    theme(axis.text.x = element_text(angle=45,vjust=1,hjust=1))
  )
}
dev.off()

