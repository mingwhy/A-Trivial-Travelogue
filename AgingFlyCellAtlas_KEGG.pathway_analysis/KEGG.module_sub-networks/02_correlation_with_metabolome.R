
library(ggplot2);library(gridExtra)
library(tidyverse);library(Matrix)
#https://www.preprints.org/manuscript/202101.0280/v1
source('propagationFuncs.R')

load("./files/all_producers.RData")
length(all_producers) #1559, each one <=> one metabolite
names(all_producers)
load("./files/rxn2gene.RData")
dim(rxn2gene) #df,6420 obs. of  4 variables, reaction and its associated enzymes
head(rxn2gene) #3rd column: Entrez 
length(unique(rxn2gene$Reaction)) #1901

# end product of each (sub)module
mz.produced=lapply(1:length(all_producers),function(i){
  all_producers[[i]]$isProduced
})
table(unlist(lapply(mz.produced,length))) #1559,each sub-graph only produce one met
mz.produced=unlist(mz.produced)
tail(sort(table(mz.produced))) #there are metabolites which are end product of multiple module
length(unique(mz.produced)) #1270
which(mz.produced=='cpd:C00388')
all_producers[[718]] #sub-network which produce Histamine

##########################################################
## read in mz measurement data for trajectory comparison
load('../batch_block_corrected_mzdata')
ls() #'dat' 'mzs'
head(dat)
dim(dat) #181  93
dat[1:3,1:8]
sort(unique(dat$AgeNum)) #4 10 24 45 69 80 days
colnames(dat)[1:8]
metabolome_ages=dat;

load('~/Documents/aging_metabolism/age-associated.metabolites.for.Ming')
ls() #"cmat",age-independent covariance between metabolites 
# "for.Ming",  mz ~ age betas 
mz.age.betas=for.Ming
dim(mz.age.betas) #86
head(mz.age.betas)
tmp=paste0('cpd:', mz.age.betas$KEGGid)
sum(tmp %in% mz.produced) #54 overlapped
mz.age.betas[tmp %in% mz.produced,]

##################################################
endNode.res=readRDS('endNode.res')
names(endNode.res)
out=lapply(1:length(endNode.res),function(i){
  x=names(endNode.res[[i]])
  mz.name=names(endNode.res)[i]
  labels=strsplit(x,';')
  tissues=unlist(lapply(labels,'[',1))
  cell.types=unlist(lapply(labels,'[',2))
  ages=unlist(lapply(labels,'[',3))
  tc=paste(tissues,cell.types,sep=';')
  df=data.frame(mz.name=mz.name,cell.type=tc,age=ages,
                raw.value=endNode.res[[i]])
  #ggplot(df,aes(x=age,y=value))+facet_wrap(.~cell.type,scale='free')+
  #  geom_point()+theme_classic()+ggtitle(paste0('mz:',mz.name))
  df
})
df.all=as.data.frame(Reduce(`rbind`,out))
head(df.all)
df.all$mz=sapply(strsplit(df.all$mz.name,'-'),'[',1)
#df.all[df.all$cell.type=='body;adult ventral nervous system' & df.all$mz=='cpd:C00388',]

# when there are metabolites as product from multiple sub-networks, sum them up
df.end=df.all %>% group_by(cell.type,age,mz) %>% summarise(value=sum(raw.value))
head(df.end)
#df.end[df.end$cell.type=='body;adult ventral nervous system' & df.end$mz=='cpd:C00388',]
#df.end[df.end$cell.type=="body;adult oenocyte" & df.end$mz=='cpd:C00388',]

##################################################################
df.end$KEGGid=gsub('cpd:','',df.end$mz)
keep.mz=unique(df.end$KEGGid)

res.cor<-lapply(keep.mz,function(mz){
  cat(mz,'\n')
  mz.name=mz.age.betas[mz.age.betas$KEGGid==mz,]$mz
  df1=as.data.frame(df.end[df.end$KEGGid==mz,c('cell.type','age','value','KEGGid')])
  
  # add metabolome data for plot
  mz.name=mz.age.betas[mz.age.betas$KEGGid==mz,]$mz
  #df2=metabolome_ages[,c('AgeNum',mz.name),drop=F]
  df2=data.frame(cell.type='metabolome',age=metabolome_ages$AgeNum,
                 mz.name=mz.name,value=metabolome_ages[,mz.name])
  fit.metabolome=loess(value ~ age, data=df2, span=1)
  newd <- data.frame(age=0:90)
  newd$pred <- predict(fit.metabolome, newd)
  newd$dat='mz'
  
  # order cell types by trajectory similarity 
  cor.values<-lapply(unique(df1$cell.type),function(i){
    fit.sc=loess(value ~ age, data=df1[df1$cell.type==i,], span=1)
    if(is.null(fit.sc) || is.nan(fit.sc$residuals[1])){return(c(i,NA))}
    newd2 <- data.frame(age=0:90)
    newd2$pred <- predict(fit.sc, newd2)
    newd2$dat='sc'
    #newdf=rbind(newd,newd2)
    #ggplot(newdf,aes(x=age,y=pred,color=dat))+geom_point()+facet_wrap(.~dat,scale='free_y')+theme_classic(base_size = 15)
    if(sd(newd2$pred,na.rm=T)<1e-22 | sd(newd$pred,na.rm=T)<1e-22){return(c(i,NA))}
    cor.value=cor(newd2$pred,newd$pred,use='complete')
    return(c(i,cor.value))
  })
  
  cor.values2<-as.data.frame(Reduce(`rbind`,cor.values))
  cor.values2=cor.values2[order(cor.values2[,2],decreasing = T),]
  colnames(cor.values2)=c('cell.type','cor')
  cor.values2$cor=as.numeric(cor.values2$cor)
  
  #colnames(df1);colnames(df2)
  dfc=rbind(df1[,c(1,2,3)],df2[,c(1,2,4)])
  dfc$KEGGid=mz
  dfc$mz.name=mz.name
  dfc=merge(dfc,cor.values2,all.x=T)
  dfc$cell.type=factor(dfc$cell.type,levels=c('metabolome',cor.values2$cell.type))
  dfc$age=as.numeric(as.character(dfc$age))
  
  return(dfc)
  
})
length(res.cor) #54
head(res.cor[[1]])
table(res.cor[[1]]$cell.type) #`metabolome`
names(res.cor)=keep.mz
saveRDS(res.cor, 'obs_mz_age.cor.rds')

sapply(res.cor,length)
##############################
plots<-lapply(keep.mz,function(mz){
  cat(mz,'\n')
  mz.name=mz.age.betas[mz.age.betas$KEGGid==mz,]$mz
  df1=as.data.frame(df.end[df.end$KEGGid==mz,c('cell.type','age','value','KEGGid')])
  
  # add metabolome data for plot
  mz.name=mz.age.betas[mz.age.betas$KEGGid==mz,]$mz
  #df2=metabolome_ages[,c('AgeNum',mz.name),drop=F]
  df2=data.frame(cell.type='metabolome',age=metabolome_ages$AgeNum,
                 mz.name=mz.name,value=metabolome_ages[,mz.name])
  fit.metabolome=loess(value ~ age, data=df2, span=1)
  newd <- data.frame(age=0:90)
  newd$pred <- predict(fit.metabolome, newd)
  newd$dat='mz'
  
  # order cell types by trajectory similarity 
  cor.values<-lapply(unique(df1$cell.type),function(i){
    fit.sc=loess(value ~ age, data=df1[df1$cell.type==i,], span=1)
    if(is.null(fit.sc) || is.nan(fit.sc$residuals[1])){return(c(i,NA))}
    newd2 <- data.frame(age=0:90)
    newd2$pred <- predict(fit.sc, newd2)
    newd2$dat='sc'
    newdf=rbind(newd,newd2)
    #ggplot(newdf,aes(x=age,y=pred,color=dat))+geom_point()+facet_wrap(.~dat,scale='free_y')+theme_classic(base_size = 15)
    if(sd(newd2$pred,na.rm=T)<1e-22 | sd(newd$pred,na.rm=T)<1e-22){return(c(i,NA))}
    cor.value=cor(newd2$pred,newd$pred,use='complete')
    return(c(i,cor.value))
  })
  
  cor.values2<-as.data.frame(Reduce(`rbind`,cor.values))
  cor.values2=cor.values2[order(cor.values2[,2],decreasing = T),]
  colnames(cor.values2)=c('cell.type','cor')
  cor.values2$cor=as.numeric(cor.values2$cor)
  
  #colnames(df1);colnames(df2)
  dfc=rbind(df1[,c(1,2,3)],df2[,c(1,2,4)])
  dfc$KEGGid=mz
  dfc$mz.name=mz.name
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
  p+scale_y_log10()
  
})

pdf('head_body_out.pdf',useDingbats = T,height = 9,width = 14)
for(i in plots){print(i)}
dev.off()


