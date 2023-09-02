
library(ggplot2);library(gridExtra);
library(tidyverse);
library(viridis)
library(ggpubr)
##############################################
## load FruitMat GEM
library(R.matlab) #in matlab `load(Fruitfly-GEM.mat)`
gem=readMat('./Fruitfly-GEM-main/model/Fruitfly-GEM.mat')
ifly=gem[[1]][,,1]
names(ifly)
dim(ifly$S) # 8135mz X 11898 fluxes (how many overlapped with metabolome assay)
genes=unlist(ifly$genes)
length(genes) #1753 genes

##############################################
load('age-associated.metabolites.for.Ming')
ls()
for.Ming$mz
metabolome=for.Ming %>% arrange(desc(beta),FDR) 
head(metabolome)
dim(metabolome) #86

# id Cross references
#kegg id C01672: https://metabolicatlas.org/explore/Fruitfly-GEM/gem-browser/metabolite/MAM02123e
#model folder in: https://github.com/SysBioChalmers/Fruitfly-GEM
all.reactions=unlist(ifly$rxnNames)
all.reactions[grep("\\bglucose\\b",all.reactions)]
all.mets=unlist(ifly$mets)
mets.mapping=data.table::fread('./Fruitfly-GEM-main/model/metabolites.tsv')
head(mets.mapping)
sum(metabolome$KEGGid %in% mets.mapping$metKEGGID) #63

share.mz=merge(metabolome[metabolome$KEGGid!='',],mets.mapping[,c('mets','metKEGGID')],by.x='KEGGid',by.y='metKEGGID')
unique(share.mz$KEGGid)

# subset Sto matrix
Smat=ifly$S[mets.mapping$metKEGGID %in% metabolome$KEGGid,]
dim(Smat) #4635 11898
sum(colSums(Smat)!=0) #2608 reactions are assocaited with those metabolies
##############################################
## after iMat is done
res=R.matlab::readMat("results_rep1/iMat_fly.body.mat") 
#res=R.matlab::readMat("results_rep1/iMat_fly.head.mat") 
experiment=res$experiment[,,1]
names(experiment)
experiment$name
cell.type_per.age=unlist(experiment$conditions)
f=experiment$fluxes.exp.all
fmat=lapply(f,'[[',1)
sapply(fmat,dim)
fmat=as.data.frame(Reduce(`cbind`,fmat))
summary(fmat[,2])
apply(fmat,2,summary)
colnames(fmat)=cell.type_per.age

# remove `unannotated` label
fmat=fmat[,-grep('unannotated', cell.type_per.age)]
cell.type_per.age=cell.type_per.age[-grep('unannotated', cell.type_per.age)]

dim(fmat) #11898 reactions x cell.type.age 
rownames(fmat)=unlist(ifly$rxns)
length(cell.type_per.age)

grp.cell.type=unlist(lapply(strsplit(cell.type_per.age,';'),'[[',2))
grp.age=unlist(lapply(strsplit(cell.type_per.age,';'),'[[',3))

##
share.mz[share.mz$mz=='HISTAMINE',]
mz.id='C00388'
mets=share.mz[share.mz$KEGGid==mz.id,]$mets

Si=ifly$S[all.mets %in% mets,]
dim(Si)
all.reactions[colSums(Si)!=0]
fmati=fmat[colSums(Si)!=0,]

plots<-lapply(unique(grp.cell.type),function(i){
  x=fmati[,grp.cell.type==i]
  age.labels=grp.age[grp.cell.type==i]
  colnames(x)=age.labels
  x$rxn=rownames(x)
  df=reshape2::melt(x)
  colnames(df)[2:3]=c('age','flux')
  df$age=factor(df$age,levels=sort(as.numeric(as.character(unique(df$age)))))
  ggplot(df,aes(x=age,y=flux,group=rxn,col=rxn))+geom_jitter()+ #geom_point()+geom_line()
    theme_classic()+ggtitle(paste(mz.id,i,sep='\n'))
})
length(plots)
grid.arrange(grobs=plots,ncol=5)

