
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

##########################################################
## read in mz measurement data for trajectory comparison
load('~/Documents/aging_metabolism/batch_block_corrected_mzdata')
ls() #'dat' 'mzs'
head(dat)
dim(dat) #181  93
dat[1:3,1:8]
sort(unique(dat$AgeNum)) #4 10 24 45 69 80 days
colnames(dat)[1:8]
metabolome_age=dat;

load('~/Documents/aging_metabolism/age-associated.metabolites.for.Ming')
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
# how share.mz in-degree, out-degree, beta from metabolome_age

keep.mz<-lapply(share.mz,function(mz){
  cat(mz,'\n')
  ids=mets.mapping[mets.mapping$metKEGGID==mz,]$mets
  
  row.index=which(ifly$mets %in% ids)
  tmp=ifly$S[row.index,,drop=F]
  x=apply(tmp,1,function(x)c(sum(x>0,na.rm=T),sum(x<0,na.rm=T)))
  #2 x metz, two rows are: in-degree(>0, product), out-degree(<0, comsuption)
  rownames(x)=c('kin','kout')
  colnames(x)<-ids
  return(x)
})
names(keep.mz)=share.mz
keep.mz[[1]]

tmp=lapply(keep.mz,function(x){
  sum(x[,grep('c$',colnames(x))])
  #sum(x[1,grep('c$',colnames(x))]) #only kin
  #sum(x[2,grep('c$',colnames(x))]) #only kout
})
tmp=lapply(keep.mz,sum)
mz.degree=as.data.frame(unlist(tmp))
colnames(mz.degree)='degree'
mz.degree$KEGGid=rownames(mz.degree)
df=merge(mz.degree,mz.age.betas,all.x=T)
dim(df) #61
df=df[order(df$degree,decreasing = T),]
summary(df$degree)
cor.test(df$degree,df$beta)
cor.test(df$degree,abs(df$beta))
par(mfrow=c(2,1))
plot(df$degree,df$beta,log='x')
plot(df$degree,abs(df$beta),log='x')
#########################
