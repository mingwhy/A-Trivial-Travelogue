
library(ggplot2);library(gridExtra)

# from mz to gene, or from gene to mz from `00_Fruitfly-GEM_metabolite_neighbors.R`
mz_gene_list=readRDS('mz_gene_list.rds')
names(mz_gene_list)
########################################################
## read in gene symbol mapping
gene.meta=readRDS('~/Documents/Data_AgingFlyCellAtlas/AFCA_gene.id.rds')
#gene.meta=readRDS('~/Documents/Data_AgingFlyCellAtlas/AFCA_gene.meta.rds')
gene.meta=gene.meta[gene.meta$FLYBASE!='',]
head(gene.meta)
dup.names=names(which(table(gene.meta$original.id)>1))

##########################################################
## read in mz measurement data for trajectory comparison
load('~/Documents/aging_metabolism/batch_block_corrected_mzdata')
ls() #'dat' 'mzs'
head(dat)
dim(dat) #181  93
dat[1:3,1:8]
sort(unique(dat$AgeNum)) #4 10 24 45 69 80 days
colnames(dat)[1:8]
sample.info=metabolome.dat[,1:6]
metabolome.dat=dat;

load('~/Documents/aging_metabolism/age-associated.metabolites.for.Ming')
ls() #"cmat",age-independent covariance between metabolites 
# "for.Ming",  mz ~ age betas 
mz.age.betas=for.Ming
dim(mz.age.betas) #86

######################################################################
## calculate mean gene expr score per cell type per age
sce=readH5AD('sce_filteredBy_ncell100.h5ad')

assayNames(sce)='counts'
sce=logNormCounts(sce)
assayNames(sce) # "counts"    "logcounts"

out.file='bootstrap_log1p_expr.rds'
unique(sce$tc)
bootstrap=100;
ncell=80;

all.mz.out=list();
for(mz in names(mz_gene_list)){
  
  # obs metabolome data
  mz.name=mz.age.betas[mz.age.betas$KEGGid==mz,]$mz
  obs.values=metabolome.dat[, c(1:6,which(colnames(metabolome.dat)==mz.name))]
  p.obs<-ggplot(obs.values,aes(x=AgeNum,y=obs.values[,7]))+geom_jitter()+
    theme_classic()+ggtitle(paste0('mz:',mz.name))+ylab('abundance')+
    scale_x_continuous(breaks=obs.values$AgeNum)+
    stat_summary(
      geom = "point",fun.y = "median",
      col = "black",size = 2,
      shape = 24,fill = "red"
    )+theme(legend.position = 'none')# +scale_y_log10()
  #p.obs
  
  # extract this mz associated gene expr from scRNA-seq
  df.g=mz_gene_list[[mz]]
  genes=c(df.g$from,df.g$to)
  genes=gsub("^C[0-9]{5}", "", genes)
  genes=genes[genes!='']  
  genes=unlist(lapply(genes,function(i){unlist(strsplit(i,'or'))}))
  genes=gsub('^\\s+|\\s+$','',genes)
  genes=unique(genes)
  
  #sce.sub=sce[intersect(genes,rownames(sce)),,drop=F]
  sce.sub=sce[intersect(genes,rownames(sce)),]
  
  one.mz.out=list()
  for(gene.name in rownames(sce.sub)){
    out<-lapply(as.character(unique(sce.sub$tc)),function(tc){
      sce.tc=sce.sub[gene.name,sce.sub$tc==tc]
      expr.m=assay(sce.tc,'logcounts')
      #expr.m[expr.m==0]=NA #remove zero expr gene, transform NA here too slow
      #class(expr.m) # "dgCMatrix"
      #table(sce.tc$age)
      
      x=lapply(unique(sce.tc$age),function(age){
        sce.tc.age=expr.m[,sce.tc$age==age,drop=F]
        expr.values<-sapply(1:bootstrap,function(i){
          tmp=sce.tc.age[,sample(1:ncol(sce.tc.age),ncell,replace = F),drop=F]
          #if(TRUE){ #remove 0, then compute mean
            tmp1=tmp[rowSums(tmp)!=0,]
            if(nrow(tmp1)==0){return(0)} #sampled values are all 0
            gene.na=tabulate(tmp1@i + 1L, nrow(tmp1)) ### nnz per row,number of non-zeros
            exprs=rowMeans_drop0(tmp)
            names(exprs)=rownames(tmp1)
            exprs
          #}
          #rowMeans(tmp)
        })
        data.frame(age=age,expr.values=expr.values)
      })
      df=as.data.frame(Reduce(`rbind`,x))
      df$cell.type=tc
      cat('cell.type',tc,'is done\n')
      return(df)
    })
      
    df.out=as.data.frame(Reduce(`rbind`,out))
    df.out$gene=gene.name;
    one.mz.out[[gene.name]]=df.out
  }
  df.one.mz=as.data.frame(Reduce(`rbind`,one.mz.out))
  all.mz.out[[mz]]=df.one.mz
}
saveRDS(all.mz.out,out.file);

names(all.mz.out)
###################################################################  
## per mz, per cell.type, select the connected genes with the max cor with metabolome
## record the chosen gene, and record observed P value
all.mz.out=readRDS('bootstrap_log1p_expr.rds')

  p=ggplot(df.out,aes(x=age,y=expr.values))+
    facet_wrap(.~cell.type, scale='free')+
    #facet_wrap(.~title,scale='free')+
    geom_point()+theme_classic()+
    ggtitle(paste0('mz:',mz.name))+
    stat_summary(
      geom = "point",fun.y = "median",
      col = "black",size = 2,
      shape = 24,fill = "red"
    )+theme(legend.position = 'none') #+scale_y_log10()
  p

