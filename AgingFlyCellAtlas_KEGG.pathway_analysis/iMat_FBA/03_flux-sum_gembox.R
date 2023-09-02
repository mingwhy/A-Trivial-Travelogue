
library(ggplot2);library(gridExtra)
library(tidyverse);library(Matrix)

##############################################
## load FruitMat GEM
library(R.matlab) #in matlab `load(Fruitfly-GEM.mat)`
gem=readMat('~/Documents/bioinfo_software/GEM_metabolic.models/Fruitfly-GEM-main/model/Fruitfly-GEM.mat')
# run GEM_convertOldStyleModel.m to get `rules` in GEM fields, used in exprs2fluxes.R
rules=readMat('rules.mat')
length(rules)

ifly=gem[[1]][,,1]
length(ifly)
lapply(ifly,dim)
names(ifly)
dim(ifly$S) # 8135mz X 11898 fluxes (how many overlapped with metabolome assay)
#genes=unlist(ifly$genes)
#length(genes) #1753 genes
ifly$rules=rules$rules
length(ifly$subSystems)
ifly$subSystems=unlist(ifly$subSystems)

ifly$genes=unlist(ifly$genes)
length(ifly$genes)
length(ifly$rules)
ifly$ub=ifly$ub[,1]
ifly$lb=ifly$lb[,1]
ifly$rowub=rep(0,nrow(ifly$S)) #used in iMat, check with human data `data("recon1")`, all 0
ifly$rowlb=rep(0,nrow(ifly$S))

##
#not do this: ifly$rxnNames=unlist(ifly$rxnNames) #unlist would drop NULL automatically
ifly$rxnNames=as.character(ifly$rxnNames)
length(unlist(ifly$rxns)) #11898
ifly$rxns=unlist(ifly$rxns)
ifly$mets=unlist(ifly$mets) #all members non-empty, checked length
ifly$metNames=unlist(ifly$metNames)

#some quick access test
length(table(ifly$subSystems)) #149
ifly$subSystems[grep('hist',ifly$subSystems)]
ifly$subSystems[grep('hist',ifly$subSystems,ignore.case = T)]

col.idx=grep('hist',ifly$subSystems,ignore.case = T);
tmp=ifly$S[,grep('hist',ifly$subSystems,ignore.case = T)]
row.idx=which(Matrix::rowSums(tmp)!=0)

ifly$rxns[col.idx]
ifly$rxnNames[col.idx]
ifly$metFormulas[col.idx]
ifly$grRules[col.idx]
ifly$metNames[row.idx]
tmp=ifly$S[row.idx,col.idx]
rowSums(tmp)
rowSums(ifly$S[row.idx,])
tmp=as.matrix(tmp)
rownames(tmp)=ifly$metNames[row.idx]
colnames(tmp)=ifly$rxnNames[col.idx]
tmp

#############################################
# id Cross references
#kegg id C01672: https://metabolicatlas.org/explore/Fruitfly-GEM/gem-browser/metabolite/MAM02123e
#model folder in: https://github.com/SysBioChalmers/Fruitfly-GEM
mets.mapping=data.table::fread('~/Documents/bioinfo_software/GEM_metabolic.models/Fruitfly-GEM-main/model/metabolites.tsv')
head(mets.mapping)
mets.mapping[mets.mapping$metKEGGID=='C00388',]#Cytosol,Extracellular #ifly$compNames
tmp=ifly$S[ifly$mets %in% mets.mapping[mets.mapping$metKEGGID=='C00388',]$mets,]
col.index=which( colSums(abs(tmp))!=0 ) 
tmp[,col.index] #all rxn which involve `C00388`, not necessariliy histine metabolism
ifly$rxns[col.index] #"MAR04428" "MAR00619"
#https://metabolicatlas.org/explore/Fruitfly-GEM/gem-browser/metabolite/MAM02124c
#https://metabolicatlas.org/explore/Fruitfly-GEM/gem-browser/metabolite/MAM02124e
ifly$subSystems[col.index]

rxns.mapping=data.table::fread('~/Documents/bioinfo_software/GEM_metabolic.models/Fruitfly-GEM-main/model/reactions.tsv')
rxns.mapping[rxns.mapping$rxns=="MAR01442",]

##############################################
## iMAT
#https://github.com/mingwhy/A-Trivial-Travelogue/tree/main/metabolic_flux/gembox-ming
(files=Sys.glob('./my.utils/R/*R')) #download `ImNotaGit/my.utils` from github
# https://github.com/ImNotaGit/my.utils
# comment out `rmote::xxx` in `utils.R` then source files
#install.packages("httpgd")
lapply(files, source)

nc <- 2L # number of CPUs/cores to use for parallelization
library(gembox)
model=ifly;

##############################################
## iMAT result
res=readRDS('head_body_flux_vec.rds')
df.info=as.data.frame(res$df.info)

datf=res$df.flux.vec
dim(datf) #11898   108
ncol(datf)==nrow(df.info)
max(datf) #1000
min(datf) #-1000

length(abs(ifly$S) %*% abs(datf [,1])) #mz length
mz.fluxsum=abs(ifly$S) %*% abs(datf)
dim(mz.fluxsum) #8135   108 cell.type_per.age
max(mz.fluxsum) # 533992.8
min(mz.fluxsum) #0
which(as.matrix(mz.fluxsum)==max(mz.fluxsum),arr.ind=T)
#3409  97
ifly$metNames[3409] #"H2O" 

cell.type=paste(df.info$tissue,df.info$cell.type,sep=';')
age=df.info$age
rownames(mz.fluxsum)=ifly$mets
#df.mz.fluxsum=as.data.frame(as.matrix(t(mz.fluxsum)))
#colnames(df.mz.fluxsum)=ifly$mets
#df.mz.fluxsum$cell.type=cell.type
#df.mz.fluxsum$age=age

##########################################################
## read in mz measurement data for trajectory comparison
load('../aging_metabolic.models/batch_block_corrected_mzdata')
ls() #'dat' 'mzs'
head(dat)
dim(dat) #181  93
dat[1:3,1:8]
sort(unique(dat$AgeNum)) #4 10 24 45 69 80 days
colnames(dat)[1:8]


load('../aging_metabolic.models/age-associated.metabolites.for.Ming')
ls() #"cmat",age-independent covariance between metabolites 
# "for.Ming",  mz ~ age betas 
mz.age.betas=for.Ming
dim(mz.age.betas) #86

share.mz=intersect(mz.age.betas$KEGGid, mets.mapping$metKEGGID)
share.mz=share.mz[share.mz!='']
length(share.mz) #61

#as.character(ifly$comps) 
#as.character(ifly$compNames)
##########################
# how many mz show 0 flux across all cell types

keep.mz<-lapply(share.mz,function(mz){
  cat(mz,'\n')
  ids=mets.mapping[mets.mapping$metKEGGID==mz,]$mets
  df=as.data.frame(t(as.matrix(mz.fluxsum[ids,,drop=F])))
  colnames(df)=ids;
  df$age=age; df$cell.type=cell.type
  df1=reshape::melt(df,id.vars=c('cell.type','age'))
  colnames(df1)[3]='mz.name'
  tmp=df1 %>% group_by(cell.type,age) %>% summarise(value=sum(value))
  if(sum(tmp$value)==0){return(NA)}
  else{return(mz)}
})
length(keep.mz) #61
keep.mz=unlist(keep.mz)
keep.mz=keep.mz[!is.na(keep.mz)]
length(keep.mz) #55

#########################

res.cor<-lapply(keep.mz,function(mz){
  cat(mz,'\n')
  ids=mets.mapping[mets.mapping$metKEGGID==mz,]$mets
  df=as.data.frame(t(as.matrix(mz.fluxsum[ids,,drop=F])))
  colnames(df)=ids;
  df$age=age;df$cell.type=cell.type
  df1=reshape::melt(df,id.vars=c('cell.type','age'))
  colnames(df1)[3]='mz.name'
  df1$mz.name=substring(df1$mz.name,1,8)
  df11=df1 %>% group_by(cell.type,age,mz.name) %>% summarise(value=sum(value))
  df11$age=as.numeric(as.character(df11$age))
  
  # add metabolome data for plot
  mz.name=mz.age.betas[mz.age.betas$KEGGid==mz,]$mz
  #df2=dat[,c('AgeNum',mz.name),drop=F]
  df2=data.frame(cell.type='metabolome',age=dat$AgeNum,
                 mz.name=mz.name,value=dat[,mz.name])
  fit.metabolome=loess(value ~ age, data=df2, span=1)
  newd <- data.frame(age=0:90)
  newd$pred <- predict(fit.metabolome, newd)
  newd$dat='mz'
  
  # order cell types by trajectory similarity 
  cor.values<-lapply(unique(df11$cell.type),function(i){
    fit.sc=loess(value ~ age, data=df11[df11$cell.type==i,], span=1)
    if(is.null(fit.sc) || is.nan(fit.sc$residuals[1])){return(NULL)}
    newd2 <- data.frame(age=0:90)
    newd2$pred <- predict(fit.sc, newd2)
    newd2$dat='sc'
    newdf=rbind(newd,newd2)
    #ggplot(newdf,aes(x=age,y=pred,color=dat))+geom_point()+facet_wrap(.~dat,scale='free_y')+theme_classic(base_size = 15)
    cor.value=cor(newd2$pred,newd$pred,use='complete')
    return(c(i,cor.value))
  })
  
  cor.values2<-as.data.frame(Reduce(`rbind`,cor.values))
  cor.values2=cor.values2[order(cor.values2[,2],decreasing = T),]
  colnames(cor.values2)=c('cell.type','cor')
  cor.values2$cor=as.numeric(cor.values2$cor)
  
  #colnames(df11);colnames(df2)
  dfc=rbind(df11,df2)
  dfc=merge(dfc,cor.values2,all.x=T)
  dfc$cell.type=factor(dfc$cell.type,levels=c('metabolome',cor.values2$cell.type))
  dfc$age=as.numeric(as.character(dfc$age))
  
  dfc$mz.kegg.id=mz
  return(dfc)
  
})
length(res.cor) #55
head(res.cor[[1]])
table(res.cor[[1]]$cell.type) #`metabolome`
names(res.cor)=keep.mz
saveRDS(res.cor, 'obs_mz_age.cor.rds')

##############################
plots<-lapply(keep.mz,function(mz){
  cat(mz,'\n')
  ids=mets.mapping[mets.mapping$metKEGGID==mz,]$mets
  df=as.data.frame(t(as.matrix(mz.fluxsum[ids,,drop=F])))
  colnames(df)=ids;
  df$age=age;df$cell.type=cell.type
  df1=reshape::melt(df,id.vars=c('cell.type','age'))
  colnames(df1)[3]='mz.name'
  df1$mz.name=substring(df1$mz.name,1,8)
  df11=df1 %>% group_by(cell.type,age,mz.name) %>% summarise(value=sum(value))
  df11$age=as.numeric(as.character(df11$age))
  
  # add metabolome data for plot
  mz.name=mz.age.betas[mz.age.betas$KEGGid==mz,]$mz
  #df2=dat[,c('AgeNum',mz.name),drop=F]
  df2=data.frame(cell.type='metabolome',age=dat$AgeNum,
             mz.name=mz.name,value=dat[,mz.name])
  fit.metabolome=loess(value ~ age, data=df2, span=1)
  newd <- data.frame(age=0:90)
  newd$pred <- predict(fit.metabolome, newd)
  newd$dat='mz'
  
  # order cell types by trajectory similarity 
  cor.values<-lapply(unique(df11$cell.type),function(i){
    fit.sc=loess(value ~ age, data=df11[df11$cell.type==i,], span=1)
    if(is.null(fit.sc) || is.nan(fit.sc$residuals[1])){return(NULL)}
    newd2 <- data.frame(age=0:90)
    newd2$pred <- predict(fit.sc, newd2)
    newd2$dat='sc'
    newdf=rbind(newd,newd2)
    #ggplot(newdf,aes(x=age,y=pred,color=dat))+geom_point()+facet_wrap(.~dat,scale='free_y')+theme_classic(base_size = 15)
    cor.value=cor(newd2$pred,newd$pred,use='complete')
    return(c(i,cor.value))
  })
  
  cor.values2<-as.data.frame(Reduce(`rbind`,cor.values))
  cor.values2=cor.values2[order(cor.values2[,2],decreasing = T),]
  colnames(cor.values2)=c('cell.type','cor')
  cor.values2$cor=as.numeric(cor.values2$cor)
  
  #colnames(df11);colnames(df2)
  dfc=rbind(df11,df2)
  dfc=merge(dfc,cor.values2,all.x=T)
  dfc$cell.type=factor(dfc$cell.type,levels=c('metabolome',cor.values2$cell.type))
  dfc$age=as.numeric(as.character(dfc$age))
  
  dfc$title=paste0(dfc$cell.type,'\n',round(dfc$cor,3))
  dfc$tmp=as.integer(dfc$cell.type)
  dfc=dfc[order(dfc$tmp),]
  dfc$title=factor(dfc$title,levels=unique(dfc$title))
  
  p=ggplot(dfc,aes(x=age,y=value))+
    #facet_wrap(.~cell.type,scale='free')+
    facet_wrap(.~title,scale='free')+
    geom_point()+theme_classic()+
    ggtitle(paste0('mz:',mz.name))+
    stat_summary(
      geom = "point",fun.y = "median",
      col = "black",size = 2,
      shape = 24,fill = "red"
    )+theme(legend.position = 'none')
  p
  
})

pdf('head_body_out.pdf',useDingbats = T,height = 9,width = 14)
#pdf('head_body_out_mzCompartment.pdf',useDingbats = T,height = 9,width = 14)
for(i in plots){print(i)}
dev.off()


