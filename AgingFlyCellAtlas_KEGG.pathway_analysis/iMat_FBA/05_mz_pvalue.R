
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

############################################
## iMAT obs result
obs.res=readRDS('obs_mz_age.cor.rds')
obs.kegg.id<-names(obs.res) #kegg.id
length(obs.kegg.id) #55

##############################################
## iMAT shuffled result
#(files=Sys.glob('imat_out_shuffle/*.rds'))
(files=Sys.glob('imat_out_shuffle2/*.rds'))
#files=Sys.glob('imat_out_shuffle_2tc_rep20/*.rds')
if(length(grep('_P.rds',files))!=0){ files=files[-grep('_P.rds',files)] }
files

# one cell type at one time
for(file in files){
  simu.out=readRDS(file) #one cell type
  out.file=paste0(file,'_P.rds')
  dim(simu.out) # rxns X rep_cell.type_age
  
  #transform #flux x condition into #mz.flux-sum x condition  view
  mz.fluxsum=abs(ifly$S) %*% abs(simu.out)
  dim(mz.fluxsum) #8135mz X  80 rep_cell.type_age
  
  x=colnames(simu.out)
  x=strsplit(x,'_')  
  rep=unlist(lapply(x,'[[',1))
  cell.type=unlist(lapply(x,'[[',3))
  age=unlist(lapply(x,'[[',2))
  age.num=as.numeric(gsub('age','',age))
  rownames(mz.fluxsum)=ifly$mets

  all.mz.res.cor<-lapply(obs.kegg.id,function(mz){
    cat(mz,'\n')
    ids=mets.mapping[mets.mapping$metKEGGID==mz,]$mets
    df=as.data.frame(t(as.matrix(mz.fluxsum[ids,,drop=F])))
    colnames(df)=ids;
    
    df$age=age.num; df$cell.type=cell.type
    df$rep=rep;
    
    rep.cor.values<-lapply(unique(rep),function(irep){
      tmp=df[df$rep==irep,];
      df1=reshape::melt(tmp[,-ncol(tmp)],id.vars=c('cell.type','age'))
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
        #return(c(i,cor.value))
        cor.value
      })
      unlist(cor.values)
    })
    df.rep.cor.values=data.frame(rep=unique(rep),cor=unlist(rep.cor.values))
    df.rep.cor.values$mz.kegg.id=mz;
    df.rep.cor.values
  })
  length(all.mz.res.cor) # number of mz in obs
  df.all.mz.res.cor=Reduce(`rbind`,all.mz.res.cor)
  saveRDS(df.all.mz.res.cor, out.file)
}
############
out.file='imat_out_shuffle_2tc_rep20/shuffle_head;adult brain perineurial glial cell.rds_P.rds'
#out.file='imat_out_shuffle_2tc_rep20/shuffle_head;gamma Kenyon cell.rds_P.rds'

df.all.mz.res.cor=readRDS(out.file)
simu.dist=df.all.mz.res.cor[df.all.mz.res.cor$mz.kegg.id=='C00388',]
head(simu.dist)

tmp=obs.res$C00388
#tmp[tmp$cell.type=='head;adult brain perineurial glial cell',]
tmp=tmp[tmp$cell.type==tc,]
pval=min(sum(tmp$cor[1]>simu.dist$cor)/length(simu.dist$cor),
         sum(tmp$cor[1]<simu.dist$cor)/length(simu.dist$cor));
hist(simu.dist$cor,main=paste0(tc,'\none tail P=',round(pval,2)),
     breaks=12,xlab='Pearson\'s r with metabolome data')
abline(v=tmp$cor[1],col='darkred',lwd=2)
######################

length(names(obs.res)) #55

out.file='imat_out_shuffle/shuffle_head;gamma Kenyon cell.rds_P.rds'
out.file='imat_out_shuffle2/shuffle_head;gamma Kenyon cell.rds_P.rds'

(files1=Sys.glob('imat_out_shuffle/shuffle_*P.rds'))
(files2=Sys.glob('imat_out_shuffle/shuffle_*P.rds'))
df.simu=c();
for(i in 1:length(files1)){
  out.file1=files1[i]; 
  (tc1=strsplit(basename(out.file1),'shuffle_|.rds')[[1]][[2]]) 
  out.file2=files2[i];
  (tc2=strsplit(basename(out.file2),'shuffle_|.rds')[[1]][[2]])
  if(tc1!=tc2){ cat('tc name not match',tc,'\n');break}
  
  df1=readRDS(out.file1)
  df2=readRDS(out.file2)
  df=rbind(df1,df2)
  df$tc=tc1;
  df.simu=rbind(df.simu,df)
}
nrow(df.simu)
length(table(df.simu$mz.kegg.id)) #55
x=df.simu %>% group_by(mz.kegg.id,tc) %>% summarise(n=n())
table(x$n)


pdf('simu_head_body_out.pdf',useDingbats = T,height = 9,width = 14)
for(mz in names(obs.res)){
  cat('mz',mz,'\n')
  mz.name=mz.age.betas[mz.age.betas$KEGGid==mz,]$mz
  simu.dist=df.simu[df.simu$mz.kegg.id==mz,]
  head(simu.dist)
  
  tmp=obs.res[[mz]]
  tmp=tmp[!duplicated(tmp$cell.type),]
  #dim(tmp)
  tmp=tmp[tmp$cell.type!='metabolome',]
  tmp=tmp[order(tmp$cor,decreasing = T),]
  
  par(mfrow=c(5,6))
  for(tc in tmp$cell.type){
    cat('tc',tc,'\n')
    tmp1=tmp[tmp$cell.type==tc,]
    tmp2=simu.dist[simu.dist$tc==tc,]
    pval=min((sum(tmp1$cor[1]>tmp2$cor,na.rm=T)+1)/length(tmp2$cor),
             (sum(tmp1$cor[1]<tmp2$cor,na.rm=T)+1)/length(tmp2$cor));
    #pval
    if(sum(is.na(tmp2$cor))==length(tmp2$cor)){next}
    hist(tmp2$cor,main=paste0(tc,'\none tail P=',round(pval,2)),
         breaks=12,#xlab='Pearson\'s r with metabolome data')
         xlab=paste0(mz,', ',mz.name),
         xlim=c(min(tmp2$cor,tmp1$cor,na.rm=T),max(tmp2$cor,tmp1$cor,na.rm=T)))
    abline(v=tmp1$cor[1],col='darkred',lwd=2)
  }
}
dev.off()


