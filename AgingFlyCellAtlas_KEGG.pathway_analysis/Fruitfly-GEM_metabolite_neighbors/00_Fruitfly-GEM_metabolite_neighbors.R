
library(ggplot2);library(gridExtra)

##############################################
## load FruitMat GEM
library(R.matlab) #in matlab `load(Fruitfly-GEM.mat)`
gem=readMat('~/Documents/bioinfo_software/GEM_metabolic.models/Fruitfly-GEM-main/model/Fruitfly-GEM.mat')

ifly=gem[[1]][,,1]

length(ifly$subSystems)==ncol(ifly$S)
ifly$subSystems=unlist(ifly$subSystems)
length(table(ifly$subSystems)) #149

#not do this: ifly$rxnNames=unlist(ifly$rxnNames) #unlist would drop NULL automatically
ifly$rxnNames=as.character(ifly$rxnNames)
length(unlist(ifly$rxns))==ncol(ifly$S) #11898
ifly$rxns=unlist(ifly$rxns)

length(unlist(ifly$mets))==nrow(ifly$S)
ifly$mets=unlist(ifly$mets) #all members non-empty, checked length
ifly$metNames=unlist(ifly$metNames)

length(unlist(ifly$grRules))==ncol(ifly$S) #FALSE
ifly$grRules=as.character(ifly$grRules)

# compartment short names and full names
as.character(ifly$comps)
as.character(ifly$compNames)
#c:Cytosol, m:Mitochondria, e:Extracellular

x=substr(ifly$mets,9,9)
table(x)
y=ifly$subSystems[x=='e']
table(y)
## test extracting one metabolite
ifly$subSystems[grep('hist',ifly$subSystems,ignore.case = T)]

col.idx=grep('hist',ifly$subSystems,ignore.case = T);
tmp=ifly$S[,grep('hist',ifly$subSystems,ignore.case = T)]
row.idx=which(Matrix::rowSums(tmp)!=0)

ifly$rxns[col.idx]
ifly$rxnNames[col.idx]
#ifly$metFormulas[col.idx]
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

load('~/Documents/aging_metabolism/age-associated.metabolites.for.Ming')
ls() #"cmat",age-independent covariance between metabolites 
# "for.Ming",  mz ~ age betas 
mz.age.betas=for.Ming
dim(mz.age.betas) #86

share.mz=intersect(mz.age.betas$KEGGid, mets.mapping$metKEGGID)
share.mz=share.mz[share.mz!='']
length(share.mz) #61

mets.mapping[mets.mapping$metKEGGID=='C00388',] #Histamine
mz_gene_list=list()

pdf('metabolite_genes_graph.pdf')
for(mz in share.mz){
  #mz='C00388';
  cat(mz,'\n')
  # might be multiple ids due to different compartments
  ids=mets.mapping[mets.mapping$metKEGGID==mz,]$mets
  ids
  #mz.age.betas[mz.age.betas$KEGGid==mz,]
  
  bait=match(ids,ifly$mets)
  tmp=ifly$S[bait,,drop=F]
  col.idx=which(Matrix::colSums(abs(tmp))!=0) #rxn target mz is involved
  #ifly$subSystems[col.idx] #there might be  "Exchange/demand reactions" and "Transport reactions"
  x=ifly$grRules[col.idx]
  col.idx=col.idx[!duplicated(x)] #rxn via the same set of genes, only keep one
  x=ifly$grRules[col.idx]
  col.idx=col.idx[x !="list(character(0))"]
  if(length(col.idx)==0){cat('Genes-associated-rxns not found\n');next}
  enzyme.genes=ifly$grRules[col.idx]
  
  tmp=ifly$S[,col.idx,drop=F]
  row.idx=which(Matrix::rowSums(abs(tmp))!=0) #metabolite these traget rxn requires
  tmp=ifly$S[row.idx,col.idx,drop=F]
  
  tmp=as.matrix(tmp)
  #rownames(tmp)=ifly$metNames[row.idx]
  rownames(tmp)=ifly$mets[row.idx]
  colnames(tmp)=ifly$rxns[col.idx]
  
  #combine ids to one metabolite and ignore stoi number, only keep 1, -1,0
  enzyme.gene.list=lapply(enzyme.genes,function(x){
    substr(x,7, nchar(x)-2)
  })
  genes.with.expr.info=intersect(unlist(enzyme.gene.list) , gene.meta$original.id)
  #cat('reaction with avaiable enzymes ',ncol(tmp),', gene ',length(genes.with.expr.info),'\n')
  if(length(genes.with.expr.info)==0){cat('Genes-associated-rxns not found\n');next}
  #target.mz.row.ids<-match(ids,rownames(tmp))
  
  i=match(rownames(tmp),mets.mapping$mets)
  #mets.mapping[i,]$mets==rownames(tmp)
  mets=unique(mets.mapping[i,]$metKEGGID)
  tmp1=lapply(mets[mets!=''],function(met){
    Matrix::colSums(tmp[mets.mapping[i,]$metKEGGID==met,,drop=F])
  })
  tmp2=as.matrix(Reduce(`rbind`,tmp1))
  rownames(tmp2)=mets[mets!='']
  
  #if it's transport reactions, tmp1 would have column sum 0, eg, C01081
  tmp2=tmp2[,colSums(abs(tmp2))!=0,drop=F]
  tmp3=tmp2[,tmp2[which(rownames(tmp2)==mz),,drop=F]!=0,drop=F]
  cat('reaction with avaiable enzymes ',ncol(tmp3),', gene ',length(genes.with.expr.info),'\n')
  
  # create igraph object
  all.df=lapply(1:ncol(tmp3),function(i){
    #cat(i,'\n')
    df=c()
    x=tmp3[,i]
    from=names(x[x<0])
    if(length(from)!=0){
      df=rbind(df,data.frame(from=from,to=enzyme.gene.list[[i]]))
    }
    to=names(x[x>0])
    if(length(to)!=0){
      df=rbind(df,data.frame(from=enzyme.gene.list[[i]],to=to))
    }
    df$reaction=colnames(tmp3)[i]
    df
  })
  all.df2=as.data.frame(Reduce(`rbind`,all.df))
  mz_gene_list[[mz]]=all.df2
  
  df.g=igraph::graph.data.frame(all.df2,directed=TRUE)
  v.names=igraph::V(df.g)$name
  v.shapes=ifelse(v.names %in% rownames(tmp3), 'circle','square')
  v.color=ifelse(v.names ==mz, 'orange','lightblue')
  df.g$shapes=v.shapes;
  df.g$color=v.color;
  plot(df.g, vertex.shape=df.g$shapes,vertex.size=10,  #edge.color="darkgreen", edge.label=links$value,
          vertex.label.font=1, vertex.color = df.g$color,
       edge.arrow.size=0.5,edge.width=1,
         vertex.label.cex = 1 , layout=igraph::layout_with_gem,#layout=igraph::layout.fruchterman.reingold,
        #main=paste0(mz,', ',colnames(tmp3)[i]))
       main=paste0(mz))
    
}
dev.off()

length(mz_gene_list) #53
saveRDS(mz_gene_list,'mz_gene_list.rds')
