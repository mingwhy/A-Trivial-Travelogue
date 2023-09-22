
library(ggplot2);library(gridExtra)
library(tidyverse)

###################################################################  
## per mz, per cell.type, select the connected genes with the max cor with metabolome
## record the chosen gene, and record observed P value
all.mz.out=readRDS('bootstrap_log1p_expr.rds')
names(all.mz.out) #53 mz, gene expression

obs=readRDS('obs_metabolome_sc_trajectory.rds')
length(obs) # n.mz
head(obs[[1]]) #cell.type, selected.best.fitted gene

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

##########################################################
## for each metabolite, permute age labels among the 6 age groups, 
# which maintained the relationships between genotypes within each age group, 

factorial(6) #720
n.rep=100; #rep100, 5min, rep200，10min

mz.permu.out=list()

pdf('permu_metabolome_head_body_out.pdf',useDingbats = T,height = 9,width = 14)
for(mz in names(all.mz.out)){
  #for(mz in names(all.mz.out)[1:2]){
  cat(mz,'\n')
  
  ## obs metabolome data
  mz.name=mz.age.betas[mz.age.betas$KEGGid==mz,]$mz
  df.met=data.frame(cell.type='metabolome',age=metabolome.dat$AgeNum,
                    mz.name=mz.name,value=metabolome.dat[,mz.name])
  fit.metabolome=loess(value ~ age, data=df.met, span=1)
  newd <- data.frame(age=0:90)
  newd$pred <- predict(fit.metabolome, newd)
  newd$dat='mz'
  
  ## obs per cell type per associated gene expr value
  ## for each mz, record cell.type, gene, obs.r
  df.one.mz=all.mz.out[[mz]] #sample multiple times per age group per cell.type
  #tmp=expand.grid(unique(df.one.mz$cell.type),unique(df.one.mz$gene))
  obs.gene.cor=obs[[mz]] #observed correlation values across cell types
  
  permu.out<-lapply(1:n.rep,function(irep){
    df.met.pseudo=df.met;
    x1=unique(as.character(df.met.pseudo$age));
    x2=sample(x1)
    names(x2)=x1
    df.met.pseudo$new.age=x2[as.character(df.met.pseudo$age)]
    #table(df.met.pseudo$age,df.met.pseudo$new.age)
    df.met.pseudo$age=as.numeric(as.character( df.met.pseudo$new.age))
    
    fit.metabolome.pseudo=loess(value ~ age, data=df.met.pseudo, span=1)
    newd.pseudo <- data.frame(age=0:90)
    newd.pseudo$pred <- predict(fit.metabolome.pseudo, newd.pseudo)
    newd.pseudo$dat='mz'
    
    tmp=obs.gene.cor;
    cor.values=apply(tmp,1,function(row){
      df1=subset(df.one.mz, cell.type==row[[1]] & gene==row[[2]])
      df1$age=as.numeric(as.character(df1$age))
      fit.sc=loess(expr.values ~ age, data=df1, span=1)
      if(is.null(fit.sc) || is.nan(fit.sc$residuals[1])){return(NA)}
      
      newd2 <- data.frame(age=0:90)
      newd2$pred <- predict(fit.sc, newd2)
      newd2$dat='sc'
      #newdf=rbind(newd,newd2)
      #ggplot(newdf,aes(x=age,y=pred,color=dat))+geom_point()+facet_wrap(.~dat,scale='free_y')+theme_classic(base_size = 15)
      if(sd(newd2$pred,na.rm=T)<1e-22 | sd(newd$pred,na.rm=T)<1e-22){return(NA)}
      cor.value=cor(newd2$pred, newd.pseudo$pred,use='complete')
      return(cor.value)
    })
    tmp$cor.coeff=cor.values
    colnames(tmp)=c('cell.type','gene','cor.coeff')
    tmp$irep=irep;
    tmp
  })
  
  df.permu.out=as.data.frame(Reduce(`rbind`,permu.out))
  mz.permu.out[[mz]]<-df.permu.out;
  
  #calculate P value across cell types
  par(mfrow=c(5,6))
  for(cell.type in unique(obs.gene.cor$cell.type)){
    cat(cell.type,'\n')
    
    x=obs.gene.cor[obs.gene.cor$cell.type==cell.type,]
    
    permu.cor.values=df.permu.out[df.permu.out$cell.type==cell.type,]$cor.coeff
    
    
    pval=min((sum(permu.cor.values>x$cor.coeff,na.rm=T)+1)/(length(permu.cor.values)+1),
             (sum(permu.cor.values<x$cor.coeff,na.rm=T)+1)/(length(permu.cor.values)+1));
    #pval
    
    hist(permu.cor.values,
         main=paste0(cell.type,'\nobs.cor=',round(x$cor.coeff,3),
                     #'\none tail P=',formatC(pval, format = "e", digits = 2)),
                     '\none tail P=',round(pval,3)),
         #https://stackoverflow.com/questions/39623636/forcing-r-output-to-be-scientific-notation-with-at-most-two-decimals
         cex.main=1,
         breaks=12,#xlab='Pearson\'s r with metabolome data')
         xlab=paste0(mz,', ',mz.name),
         xlim=c(min(permu.cor.values,x$cor.coeff,na.rm=T),max(permu.cor.values,x$cor.coeff,na.rm=T)))
    abline(v=x$cor.coeff,col='darkred',lwd=2)
  }
}

dev.off()

saveRDS(mz.permu.out,'shuffle_metabolome_out.rds')



