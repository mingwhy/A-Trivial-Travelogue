
library(igraph);library(tidyverse);
library(ggplot2);library(gridExtra)
#https://www.preprints.org/manuscript/202101.0280/v1
#https://github.com/babelomics/Metabolica
source('../propagationFuncs.R')

load("../files/all_producers.RData")
length(all_producers) #1559, each one <=> one metabolite
names(all_producers)

load("../files/rxn2gene.RData")
dim(rxn2gene) #df,6420 obs. of  4 variables, reaction and its associated enzymes
head(rxn2gene) #3rd column: Entrez 
length(unique(rxn2gene$Reaction)) #1901

mz.produced=lapply(1:length(all_producers),function(i) all_producers[[i]]$isProduced )
table(unlist(lapply(mz.produced,length))) #each sub-graph only produce one met
mz.produced=unlist(mz.produced)
(i=which(mz.produced=='cpd:C00388')) #Histamine, https://www.genome.jp/entry/C00388
all_producers[[i]] #sub-network which produce Histamine
length(unique(mz.produced)) #1270
length(mz.produced) #1559
# one metabolite might be generated via multiple modules

mz.start=lapply(1:length(all_producers),function(i) all_producers[[i]]$initialNodes )
table(unlist(lapply(mz.start,length))) #each sub-graph can have multiple start mz.start
############################################################
# read in fly-human orthologs (00_process_h5ad_AFCA.R)
orthologs=readRDS('~/Documents/Data_AgingFlyCellAtlas/human_fly_ortholog_entrezid.rds')
colnames(orthologs)
sum(orthologs$human_entrezid=='')
sum(is.na(orthologs$human_entrezid)) #1558
orthologs=orthologs[!is.na(orthologs$human_entrezid),]
dim(orthologs) # 24085    10
gene.order=unique(orthologs$human_entrezid)
length(gene.order)#11310
##################################################
#overlap with metabolome measured data
load('~/Documents/aging_metabolism/batch_block_corrected_mzdata')
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

overlap.mz=mz.age.betas[tmp %in% mz.produced,]$KEGGid
sum(overlap.mz=='C00388') #exist 
###################################################
# pathview: https://bioconductor.org/packages/release/bioc/html/pathview.html
library(pathview)
library(igraph)
load("../files/all_producers.RData")
length(all_producers) #1559, each one <=> one metabolite
names(all_producers)

i='C00388'
i='C00135'
#for(i in overlap.mz){
pathway.index=grep(paste0('cpd:',i),names(all_producers))
names(all_producers)[pathway.index]

if(length(pathway.index)==0){next}
for(j in pathway.index){
  
  my_producer <- all_producers[[j]] # i can 1 to length(all_producers)
  pathway.id=gsub('hsa|_SIF','',my_producer$cluster)
  cpd=gsub('cpd:','',my_producer$isProduced)
  #g<-my_producer$clusteRigraph
  g<-my_producer$produceRigraph
  node.names=V(g) #E(g)
  node.names=names(node.names)
  node.names=node.names[grep('cpd',node.names)]
  node.names
  inp.cpd=rep(1,length(node.names))
  names(inp.cpd)=gsub('cpd:','',node.names)
  my_producer$cluster #"hsa00340_SIF"
  
  # figure: kegg.native 
  r=tryCatch(
    pv.out <- pathview(cpd.data = inp.cpd, 
                       pathway.id =pathway.id, species = "hsa", 
                       #out.suffix = "hsa00340_C00388.cpd");
                       out.suffix = cpd),
    error=function(c)'error')
  if(r=='error'){next }
  str(pv.out)
  names(pv.out)
  head(pv.out$plot.data.gene)
  head(pv.out$plot.data.cpd)
  # figure: graphviz
  pv.out <- pathview(cpd.data = inp.cpd, 
                     pathway.id =pathway.id, species = "hsa", 
                     kegg.native = F,
                     out.suffix = cpd);
}
#}  

########################################################
## replace fly gene symbols to human gene ENTREZID
dataset.name='fly.female'
dat=readRDS('../log1p_female_gene.mean.expr.rds')
length(dat) #27
length(dat[[1]]) #4
names(dat[[1]])
# tissue;cell.type;age, 4 age groups
head(dat[[1]][[1]])

sapply(dat[[2]],length)
head(names(dat[[2]][[1]]))
expr.mat=matrix(0,nrow=length(orthologs$human_entrezid),ncol=27*4)
rownames(expr.mat)=orthologs$human_entrezid
col.names.vec=c();
z=0;
for(i in 1:length(dat)){
  for(j in 1:length(dat[[i]])){
    z=z+1;
    expr.vec=dat[[i]][[j]]
    expr.vec=expr.vec[names(expr.vec) %in% orthologs$original.id]
    #length(expr.vec)
    tmp=orthologs[match(names(expr.vec),orthologs$original.id),]
    #sum(tmp$original.id==names(expr.vec))
    df.tmp=data.frame(gene.id=tmp$human_entrezid,expr=expr.vec[tmp$original.id])
    #sum(table(df.tmp$gene.id)>1)
    df.tmp2=df.tmp %>% group_by(gene.id) %>% summarize(expr.value=mean(expr))
    #length(unique(df.tmp2$gene.id))==nrow(df.tmp2)
    expr.mat[df.tmp2$gene.id,z]=df.tmp2$expr.value
    #colnames(expr.mat)[z]=names(dat[[i]])[[j]]
    col.names.vec=c(col.names.vec,names(dat[[i]])[[j]])
  }
}
dim(expr.mat)
colnames(expr.mat)=col.names.vec
expr.mat[1:3,1:3]
###################################################
## look at enzyme mapped to this sub-network
## nodeVals only need to compute once
nodeVals <- do.calc.rxnvals(all_producers,rxn2gene,expr.mat)
dim(nodeVals); #3056 108
length(unique(unlist(rownames(nodeVals)))) #3056, #node=(rxn+cpd)
#nodeVals[1:10,1:2]
#table(Matrix::rowSums(nodeVals))

# pick one metabolite
sum(overlap.mz=='C00388') #exist 
kegg.id='C00388'; #histamine
kegg.id='C00135' #histidine
#kegg.id='C01181';
pathway.index=grep(paste0('cpd:',kegg.id),names(all_producers))
names(all_producers)[pathway.index]
length(pathway.index)

for(j in pathway.index){
  my_producer <- all_producers[[j]] # i can 1 to length(all_producers)
  pathway.id=gsub('hsa|_SIF','',my_producer$cluster)
  cpd=gsub('cpd:','',my_producer$isProduced)
  
  nodes=c(V(my_producer$produceRigraph)) #extract nodes
  g<-my_producer$produceRigraph
  node.names=V(g) #E(g)
  node.names=names(node.names)
  rxn.node.names=node.names[grep('rn',node.names)]
  rxn.node.names
  
  mapped.human.genes=rxn2gene[rxn2gene$Reaction %in% rxn.node.names,]
  mapped.human.genes
  
  dim(expr.mat) #gene by tissue;cell.type;age
  gene.index=intersect(mapped.human.genes$Entrez , rownames(expr.mat)) # 3 out of 4 genes mapped
  sc.expr=expr.mat[gene.index,]
  
  ## data.frame, gene expr 
  df.sc=reshape2::melt(sc.expr)
  colnames(df.sc)=c('gene','tc_age','value')
  
  tmp=strsplit(as.character(df.sc$tc_age),';')
  tissues=sapply(tmp,'[[',1)
  cell.types=sapply(tmp,'[[',2)
  ages=sapply(tmp,'[[',3)
  df.sc$age=ages;
  df.sc$age=factor(df.sc$age,levels=sort(unique(as.numeric(ages))))
  df.sc$tc=paste(tissues,cell.types,sep=';')
  
  ## data.frame, substrate node values from metabolica
  node.values=nodeVals[rxn.node.names,,drop=F]
  
  df.node=reshape2::melt(node.values)
  colnames(df.node)=c('reaction','tc_age','value')
  tmp=strsplit(as.character(df.node$tc_age),';')
  tissues=sapply(tmp,'[[',1)
  cell.types=sapply(tmp,'[[',2)
  ages=sapply(tmp,'[[',3)
  df.node$age=ages;
  df.node$age=factor(df.node$age,levels=sort(unique(as.numeric(ages))))
  df.node$tc=paste(tissues,cell.types,sep=';')
  
  ## data.frame, metabolite product values
    trial<-tryCatch(
      res <- metabolica(nodes.vals=nodeVals, subgraph = my_producer$produceRigraph, 
                        ininodes = my_producer$initialNodes, 
                        endnode = my_producer$isProduced, algR = "minprod", algM = "nsqrt", 
                        maxnum = my_producer$maxnum*1.5),
      error=function(c) 'error')
    if(trial=='error'){return(NA)}
    res$node.signal
    
  df.product=data.frame(tc_age=names(res$node.signal),value=res$node.signal)
  tmp=strsplit(as.character(df.product$tc_age),';')
  tissues=sapply(tmp,'[[',1)
  cell.types=sapply(tmp,'[[',2)
  ages=sapply(tmp,'[[',3)
  df.product$age=ages;
  df.product$age=factor(df.product$age,levels=sort(unique(as.numeric(ages))))
  df.product$tc=paste(tissues,cell.types,sep=';')
  
  # data.frame, add metabolome data for plot
  mz.name=mz.age.betas[mz.age.betas$KEGGid==kegg.id,]$mz
  #df2=dat[,c('AgeNum',mz.name),drop=F]
  df.metabolome=data.frame(cell.type='metabolome',age=metabolome_ages$AgeNum,
                 mz.name=mz.name,value=metabolome_ages[,mz.name])
  fit.metabolome=loess(value ~ age, data=df.metabolome, span=1)
  newd <- data.frame(age=0:90)
  newd$pred <- predict(fit.metabolome, newd)
  newd$dat='mz'
  
  ##############
  ## ready for plot, select one cell type
  
  #pick.tc='head;photoreceptor cell R7';
  #pick.tc='head;adult optic chiasma glial cell';
  #pick.tc='body;adult oenocyte';
  
  pdf(paste0('hsa',pathway.id,'_',kegg.id,'.pdf'),useDingbats = T,height = 8,width = 10)
  for(pick.tc in unique(df.product$tc)){
  
    df.tmp1=df.sc[df.sc$tc==pick.tc,]
    p1=ggplot(df.tmp1,aes(x=age,y=value))+
      facet_wrap(.~gene,scale='free')+ geom_point()+theme_classic()+
      ggtitle(pick.tc)+stat_summary(
        geom = "point",fun.y = "median",
        col = "black",size = 2,
        shape = 24,fill = "red"
      )#+theme(legend.position = 'none')
    
    df.tmp2=df.node[df.node$tc==pick.tc,]
    p2=ggplot(df.tmp2,aes(x=age,y=value))+
      facet_wrap(.~reaction,scale='free')+ geom_point()+theme_classic()+
      ggtitle(pick.tc)+stat_summary(
        geom = "point",fun.y = "median",
        col = "black",size = 2,
        shape = 24,fill = "red"
      )#+theme(legend.position = 'none')
    
    df.tmp3=df.product[df.product$tc==pick.tc,]
    p3=ggplot(df.tmp3,aes(x=factor(age,levels=sort(df.tmp3$age)),y=value))+
      geom_point()+theme_classic()+
      ggtitle(pick.tc)+scale_y_log10()+xlab('age')+
      #scale_y_continuous(trans=scales::pseudo_log_trans(base = 10))+
      stat_summary(
        geom = "point",fun.y = "median",
        col = "black",size = 2,
        shape = 24,fill = "red"
      )+theme(legend.position = 'none')
    p3
    
    df.product$age=as.numeric(as.character(df.product$age))
    fit.sc=loess(value ~ age, data=df.product[df.product$tc==pick.tc,], span=1)
    if(is.null(fit.sc) || is.nan(fit.sc$residuals[1])){return(NULL)}
    newd2 <- data.frame(age=0:90)
    newd2$pred <- predict(fit.sc, newd2)
    newd2$dat='sc'
    
    if(sd(newd2$pred,na.rm=T)<1e-22 | sd(newd$pred,na.rm=T)<1e-22){
      cat('SD==0, not for calculating cor\n')
      cor.value='NA';
    }else{
      cor.value=cor(newd2$pred,newd$pred,use='complete')
      cor.value=round(cor.value, 3)
    }
    #df.metabolome$age=factor(df.metabolome$age,levels=sort(unique(df.metabolome$age)))
    p4=ggplot(df.metabolome,aes(x=age,y=value))+
      geom_point()+theme_classic()+
      ggtitle(paste0('mz:',mz.name,', Pearson\'r:',cor.value))+
      stat_summary(
        geom = "point",fun.y = "median",
        col = "black",size = 2,
        shape = 24,fill = "red"
      )+#theme(legend.position = 'none')
      geom_vline(xintercept = as.numeric(as.character(unique(df.product$age))),linetype = "dashed")
    p4
    print(grid.arrange(p1,p2,p3,p4,ncol=2))
    
  }
  dev.off()
}



