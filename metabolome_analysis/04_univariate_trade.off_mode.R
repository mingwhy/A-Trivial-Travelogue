
options(stringsAsFactors = F)
library(ggplot2)
library(gridExtra)
library(car)

plot_mz<-function(df.inp,mz.name){
  mz=df.inp[rownames(df.inp)==mz.name,]
  tmp=data.frame(mz=as.numeric(mz),age=factor(pheno$age),mode=factor(pheno$mode),batch=factor(pheno$batch));
  tmp$group=factor(paste(tmp$mode,tmp$age))
  tmp$group=factor(tmp$group,levels=group.levels)
  tmp$mode=factor(tmp$mode,levels=mode.levels)
  
  plots=list();
  plots[[1]]=ggplot(tmp,aes(x=mode,y=mz,col=age,group=group,shape=batch))+
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(size=3,position = position_jitterdodge())+
    theme_bw(base_size = 14)+ggtitle(mz.name)
  
  plots[[2]]=ggplot(tmp,aes(x=age,y=mz,col=mode,group=group,shape=batch))+
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(size=3,position = position_jitterdodge())+
    theme_bw(base_size = 14)+ggtitle(mz.name)

  do.call(grid.arrange,c(plots,ncol=2))
}

mz.out.file='tradeoff_candidate.mz.csv'

##################################################
dat=readRDS('./mz.filter.dat.rds')
df.mat=dat$edata;
pheno=dat$pheno

## log transform
#summary(as.numeric(unlist(df.mat))) #no zero, log directly
df.mat.log=log(df.mat)

## normalize by sample
df.scaled=scale(df.mat.log,center=T,scale=T)
#apply(df.scaled,2,mean) #confirm
#apply(df.scaled,2,sd) #confirm

## mode setting
if(T){
# trade-off
#wt is wildtype
#wt/353 is the ‘Mode 2’, no trade-off
#the other three are the 'Mode 1', extend lifespan by reducing costs of reproduction.
pheno$mode='mode.1';
pheno[pheno$genotype=='wt/353',]$mode='mode.2'
pheno[pheno$genotype=='wt',]$mode='wt'
#table(pheno$mode)
#mode.1 mode.2     wt 
#36     12     12 
}
if(F){
#insulin sensitive: wt, wt/353
#insulin resistant: 19/74, 74/211, 211/19
pheno$mode=rep('insulin.resistant',nrow(pheno))
pheno[pheno$genotype %in% c('wt','wt/353'),]$mode='insulin.sensitive'
}

if(F){
#'normal' fecundity: wt, 74/211, 211/19
#high fecundity: wt/353
#low fecundity: 19/74 
pheno$mode=rep('normal.fecundity',nrow(pheno))
pheno[pheno$genotype=='wt/353',]$mode='high.fecundity'
pheno[pheno$genotype=='19/74',]$mode='low.fecundity'
}
if(F){
#genotype
pheno$mode=pheno$genotype
}


group.levels=sort(unique(paste(pheno$mode,pheno$age)))
mode.levels=sort(unique(pheno$mode))
mode.levels=c('wt','mode.1','mode.2') #re-level the 3 modes

# select two metabolites
if(F){
pick.mz=c('Sarcosine','Reduced Glutathione');
df.tmp=df.scaled[pick.mz,]
plots=lapply(1:nrow(df.tmp),function(i){
    mz=df.tmp[i,]
    mz.name=rownames(df.tmp)[i]
    tmp=data.frame(mz=as.numeric(mz),age=factor(pheno$age),mode=factor(pheno$mode),batch=factor(pheno$batch));
    tmp$group=factor(paste(tmp$mode,tmp$age))
    tmp$group=factor(tmp$group,levels=group.levels)
    tmp$mode=factor(tmp$mode,levels=mode.levels)
    
    ggplot(tmp,aes(x=mode,y=mz,col=age,group=group,shape=batch))+
      geom_boxplot(outlier.shape = NA) +
      geom_jitter(size=3,position = position_jitterdodge())+
      theme_bw(base_size = 14)+ggtitle(mz.name)
  })
  do.call(grid.arrange,c(plots,ncol=2))

  plots=lapply(1:nrow(df.tmp),function(i){
    mz=df.tmp[i,]
    mz.name=rownames(df.tmp)[i]
    tmp=data.frame(mz=as.numeric(mz),age=factor(pheno$age),mode=factor(pheno$mode),batch=factor(pheno$batch));
    tmp$group=factor(paste(tmp$mode,tmp$age))
    tmp$group=factor(tmp$group,levels=group.levels)
    tmp$mode=factor(tmp$mode,levels=mode.levels)
   ggplot(tmp,aes(x=age,y=mz,col=mode,group=group,shape=batch))+
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(size=3,position = position_jitterdodge())+
    theme_bw(base_size = 14)+ggtitle(mz.name)
  })
  do.call(grid.arrange,c(plots,ncol=2))
}  


## --------------------------------------------------------------------------------------------------------------------
df.res=readRDS('mz.filter.combat.rds')
#df.tmp=as.data.frame(df.res);
#df.tmp$metabolite=rownames(df.tmp)
#data.table::fwrite(df.tmp,file='batch.effect.corrected_InR_metaboliomics.csv',row.names = T,col.names = T)
#data.table::fwrite(pheno,file='pheno.txt',row.names = F,col.names = T)
  
### use genotype and age as predictor
out=lapply(1:nrow(df.res),function(i){
  mz=df.res[i,]
  mz.name=rownames(df.res)[i] # fit model using df.res with age and mode as predictor, including interaction term
  tmp=data.frame(mz=as.numeric(mz),age=factor(pheno$age),mode=factor(pheno$mode));
  tmp$mode=factor(tmp$mode,levels=mode.levels) #re-level 'mode'
  x=lm(mz~age*mode,tmp )
  x
})
names(out)=rownames(df.res)

# use Anova to get p values for each factor of each mz
#Anova(out[[1]])
sig.out=lapply(out,function(x){ #make a table
  x1=Anova(x)
  #x1=summary(x) #this gives you coef, coef.p.value for each term
  pvalue=x1$`Pr(>F)`
  #terms=rownames(x1)
  pvalue
})

df.sig.out=Reduce(`rbind`,sig.out)
x1=Anova(out[[1]])
colnames(df.sig.out)=rownames(x1)
rownames(df.sig.out)=rownames(df.res)
df.sig.out=df.sig.out[,-4] #remove Residuals column
#df.sig.out.adj=data.frame(matrix(p.adjust(df.sig.out,method = 'BH'),ncol=ncol(df.sig.out)))
df.sig.out.adj=data.frame(apply(df.sig.out,2,function(i) p.adjust(i,method='BH')))
#df.sig.out.adj=data.frame(apply(df.sig.out,2,function(i) p.adjust(i,method='fdr')))
rownames(df.sig.out.adj)=rownames(df.sig.out)
colnames(df.sig.out.adj)=colnames(df.sig.out)


## --------------------------------------------------------------------------------------------------------------------
var.explain=sapply(out,function(i) {x=summary(i);x$adj.r.squared})
hist(var.explain,main='variance explained')


## --------------------------------------------------------------------------------------------------------------------
df.sig.out.adj[c(rownames(df.sig.out.adj)[1:3],pick.mz),]


## --------------------------------------------------------------------------------------------------------------------
par(mfrow=c(1,3))
hist(df.sig.out.adj$age,main='p.adj for Age')
hist(df.sig.out.adj$mode,main='p.adj for Mode')
hist(df.sig.out.adj$`age:mode`,main='p.adj for Age:Mode')
p.value.cutoff=0.01;


## --------------------------------------------------------------------------------------------------------------------
df.sig.out.adj2=df.sig.out.adj;
apply(df.sig.out.adj2,2,function(i){sum(i<p.value.cutoff)})
#age     mode age:mode 
#76      127        4 
#rownames(df.sig.out.adj[df.sig.out.adj$mode<p.value.cutoff,])
#df.sig.out.adj['Sarcosine',]

## --------------------------------------------------------------------------------------------------------------------
length(out) #170 regression result
nrow(df.sig.out.adj2) #170 mz's p.adj

x=out[[1]]
summary(x) #age15 is the ref level

# for those metabolites which has a sig 'age' factor
#  decrease or increase with age
test.mz=rownames(df.sig.out.adj[df.sig.out.adj$age<0.01,])
length(test.mz) #76

pool=c('b0','age','mode','mode','age:mode','age:mode')
beta=t(sapply(test.mz,function(mz.name){
  x=out[[mz.name]]
  r=summary(x)$coefficients
  r[2,]
  #tmp=r[c(2,5,6),]
  #tmp[which(tmp[,4]<0.05)[1],]
}))
sum(beta[,4]<0.05) #only 44 
sum(beta[,1]>0) #39 increase with age
sum(beta[,1]<0) #37 decrease with age     
sum(beta[,1]>0 & beta[,4]<0.05)  #21 increase with age
sum(beta[,1]<0 & beta[,4]<0.05)  #23 decrease with age

## look at mode in detail
test.mz=rownames(df.sig.out.adj[df.sig.out.adj$mode<0.01,])
length(test.mz) #127
beta=t(sapply(test.mz,function(mz.name){
  x=out[[mz.name]]
  r=summary(x)$coefficients
  #r[2,]
  tmp=r[c(2,3),]
  tmp1=tmp[order(tmp[,4]),]
  tmp1[1,]
}))
sum(beta[,4]<0.05) #101 
# extract mz with sig mode2 and non-sig mode1
beta=lapply(test.mz,function(mz.name){
  x=out[[mz.name]]
  r=summary(x)$coefficients
  #r[2,]
  tmp=r[c(3,4),] #mode1 and mode2
  #if(tmp[1,4]>0.05 & tmp[2,4]<0.05){
  if((tmp[1,4]<0.05 | tmp[2,4]<0.05) & tmp[1,1]*tmp[2,1]<0){ #one sig and opposite direction
  #if(tmp[1,1]*tmp[2,1]<0){ 
    tmp[2,]
  }
})
names(beta)=test.mz
length(beta) #127
beta=Filter(Negate(is.null), beta) #remove NULL
length(beta) #39
sort(names(beta))
df.beta=data.frame(Reduce(`rbind`,beta))
df.beta=cbind(names(beta),df.beta)
colnames(df.beta)[1]='metabolite'
data.table::fwrite(df.beta,mz.out.file)

## plot some mz as examples
pick=grep('Methionine|glutath|Sarcosine',df.beta$metabolite,ignore.case = T)
df.beta[pick,]
pdf('test.pdf',useDingbats = T,width = 9,height = 5)
for(mz.name in df.beta[pick,]$metabolite){
  plot_mz(df.inp=df.res,mz.name=mz.name)
}
dev.off()

nrow(df.beta) #39
pdf('sig.mz_39.pdf',useDingbats = T,width = 9,height = 5)
for(mz.name in df.beta$metabolite){
  plot_mz(df.inp=df.res,mz.name=mz.name)
}
dev.off()


### boxplots for all 170 metabolites
test.mz=rownames(df.sig.out.adj)
length(test.mz) #170

pdf('all.mz_170.pdf',useDingbats = T,width = 9,height = 5)
for(mz.name in test.mz){
  plot_mz(df.inp=df.res,mz.name=mz.name)
}
dev.off()

