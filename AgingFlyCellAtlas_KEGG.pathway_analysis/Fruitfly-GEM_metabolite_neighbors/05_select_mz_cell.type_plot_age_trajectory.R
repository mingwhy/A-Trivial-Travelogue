library(ggplot2);
library(gridExtra);
##########################################################
## read in mz measurement data for trajectory comparison
load('~/Documents/aging_metabolism/batch_block_corrected_mzdata')
ls() #'dat' 'mzs'
head(dat)
dim(dat) #181  93
dat[1:3,1:8]
sort(unique(dat$AgeNum)) #4 10 24 45 69 80 days
colnames(dat)[1:8]
metabolome.dat=dat;
sample.info=metabolome.dat[,1:6]

load('~/Documents/aging_metabolism/age-associated.metabolites.for.Ming')
ls() #"cmat",age-independent covariance between metabolites 
# "for.Ming",  mz ~ age betas 
mz.age.betas=for.Ming
dim(mz.age.betas) #86
############################################################
df.mz.out=readRDS('metabolite_cell.type_Pvalue.rds')
#df.mz.out=readRDS('metabolite_cell.type_Pvalue_remove0.rds')
head(df.mz.out)

sum(df.mz.out$Pval<0.01,na.rm=T) #34 for remove0
df.mz.out.sig=df.mz.out[df.mz.out$Pval<0.01 & df.mz.out$cor.coeff>0,]
dim(df.mz.out.sig) #16 for remove0, or 27 otherwise
ranked.mz<-sort(table(df.mz.out.sig$mz))
length(ranked.mz) #12 unique mz
df.mz.out.sig$mz=factor(df.mz.out.sig$mz,levels=names(ranked.mz))
df.mz.out.sig=df.mz.out.sig[order(df.mz.out.sig$mz),]
unique(df.mz.out.sig$mz)
ranked.mz

######################################################################
# plot observed metabolome trajectory VS gene expression trajectory
# code copied from `01_metabolite_gene.neighbors_maxR.R`


all.mz.out=readRDS('bootstrap_log1p_expr.rds')
#all.mz.out=readRDS('bootstrap_log1p_expr_remove0.rds')
names(all.mz.out) #53 mz

obs.result=readRDS('obs_metabolome_sc_trajectory.rds')
#obs.result=readRDS('obs_metabolome_sc_trajectory_remove0.rds')

pdf('sig_metabolome_sc_trajectory.pdf',height = 3.5,width = 8)
#pdf('sig_metabolome_sc_trajectory_remove0.pdf',height = 3.5,width = 8)
for(i in 1:nrow(df.mz.out.sig)){
  row=df.mz.out.sig[i,]
  mz=as.character(row$mz)
  cat(mz,'\n')
  
  ## obs metabolome data
  mz.name=mz.age.betas[mz.age.betas$KEGGid==mz,]$mz
  df.met=data.frame(cell.type='metabolome',age=metabolome.dat$AgeNum,
                    mz.name=mz.name,value=metabolome.dat[,mz.name])
  fit.metabolome=loess(value ~ age, data=df.met, span=1)
  newd <- data.frame(age=0:90)
  newd$pred <- predict(fit.metabolome, newd)
  newd$dat='mz'
  
  p.obs<-ggplot(df.met,aes(x=age,y=value))+
    geom_jitter(size=1)+geom_violin(aes(group=age),fill=NA)+
    theme_classic()+ggtitle(paste0(mz,', ',mz.name))+ylab('abundance')+
    scale_x_continuous(breaks=df.met$age)+
    stat_summary(
      geom = "point",fun.y = "median",
      col = "black",size = 2,
      shape = 24,fill = "red"
    )+theme(legend.position = 'none',
            plot.title = element_text(size=10))+ # +scale_y_log10()
    geom_line(data=newd,aes(x=age,y=pred),color='blue')
  #p.obs
  
  ## obs per cell type per associated gene expr value
  ## for each mz, record cell.type, gene, obs.r
  df.one.mz=all.mz.out[[mz]]
  df.one.mz=df.one.mz[df.one.mz$cell.type==row$cell.type,]
  
  x=obs.result[[mz]]
  picked.gene=x[x$cell.type==row$cell.type,]$gene
  df.one.mz<-df.one.mz[df.one.mz$gene==picked.gene,]
  
  df1=df.one.mz
  df1$age=as.numeric(as.character(df1$age))
  fit.sc=loess(expr.values ~ age, data=df1, span=1)
  #if(is.null(fit.sc) || is.nan(fit.sc$residuals[1])){return(NA)}
  newd2 <- data.frame(age=0:90)
  newd2$pred <- predict(fit.sc, newd2)
  newd2$dat='sc'
    
  df1$age=as.numeric(as.character(df1$age))
  p=ggplot(df1,aes(x=age,y=expr.values))+
      facet_wrap(.~cell.type, scale='free')+#geom_violin(aes(group=age),fill=NA)+
      #facet_wrap(.~title,scale='free')+
      ylab(paste0('mz:',mz.name))+
      geom_jitter(size=1)+theme_classic()+#+scale_y_log10()
      stat_summary(
        geom = "point",fun.y = "median",
        col = "black",size = 2,
        shape = 24,fill = "red"
      )+theme(legend.position = 'none') +
      ggtitle(paste0('Pearson r = ',round(row[[3]],3),', gene ',picked.gene,
                     ', P = ',round(row$Pval,5)))
  p.cell.type<-p+geom_vline(xintercept = as.numeric(as.character(unique(df.met$age))),linetype = "dashed")+
      theme(axis.title.y = element_text(size=5),
            plot.title = element_text(size=10))+
      geom_line(data=newd2,aes(x=age,y=pred),color='blue')+
      scale_x_continuous(breaks=sort(unique(df1$age,df.met$age)))

  grid.arrange(p.obs,p.cell.type,ncol=2)
}
dev.off()



