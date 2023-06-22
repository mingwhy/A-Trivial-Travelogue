
## for a give cell classification, calculate the information content of a GO term
options(stringsAsFactors = F)

library(Seurat)
#library(iterators)
#library(doParallel)
#numCores <- parallel::detectCores() # Requires library(parallel)
#registerDoParallel(2)
source('src_normalized.mutual.information.R')


## extract a GO term and its assocaited genes
library(org.Dm.eg.db,verbose=F,quietly=T)
library(GO.db);
packageVersion('org.Dm.eg.db') #"3.13.0"
packageVersion('GO.db') #"3.13.0"

all.fly.genes<-keys(org.Dm.eg.db,"FLYBASE")
length(all.fly.genes) #25097 genes
length(unique(all.fly.genes)) #25097 genes
fb2go=select(org.Dm.eg.db,keys=all.fly.genes,keytype = 'FLYBASE',
             columns = c('GO'))
pick.go=c('GO:0065004','GO:1902275')
goterms=Term(GOTERM)
goterms[[pick.go[[1]]]] #"protein-DNA complex assembly"
goterms[[pick.go[[2]]]] #"regulation of chromatin organization"

go2fb=lapply(pick.go,function(i) fb2go[which(fb2go$GO==i),]$FLYBASE)
go2fb
# change FBgn to symbol
go2symbol=lapply(go2fb,function(i){
  AnnotationDbi::select(org.Dm.eg.db,keys=i,keytype='FLYBASE',c('SYMBOL'))
})
names(go2symbol)=pick.go
go2symbol      

## read in single cell data
dat.test=readRDS('test.rds')
dat.test #12878 features across 985 samples
dat.test<-NormalizeData(dat.test,normalization.method = "LogNormalize")

df.test=dat.test@assays$RNA@data
#df.test=dat.test@assays$RNA@counts

# log2(CPM+1) transformation
i=Matrix::rowSums(df.test) #remove genes with all 0
df.test=df.test[!(i==0),]
#df.test = log(df.test+1,base=2);

# get sex label
table(dat.test$sex)
x=as.numeric(factor(dat.test$sex)) #male=2,female=1
label.test=data.frame(name=colnames(df.test),label=x)
if(sum(label.test$name!=colnames(df.test))==ncol(df.test)){
  cat('warning: cell number of label and df discrete matrix do not match!\n')
}

df.test=as.matrix(df.test)
genes=rownames(df.test)
n.rep=10;
mi.out=list()

## begin GO term information content calculation
for(go.id in names(go2symbol)){
  go.genes=go2symbol[[go.id]]$SYMBOL
  go.genes2=go.genes[go.genes %in% genes]
  if(length(go.genes2)<=2){mi.out[[go.id]]=NA;next}
  
  mat=df.test[go.genes2,]
  #dist.x=as.matrix(dist(t(mat),method='manhattan'))
  dist.x=as.matrix(dist(t(mat)))
  #sum(attr(x,'Labels')==colnames(mat))
  #sum(colnames(x)==colnames(mat))
  #ncol(mat);sum(label.test$name==colnames(mat))
  # select one female and male cell 
  out=list();
  for(i.rep in 1:n.rep){
    tmp=matrix(0,ncol=2,nrow=2)
    colnames(tmp)=rownames(tmp)=c('1','2')
    i=sample(which(label.test$label==1),1) #female
    j=sample(which(label.test$label==2),1) #male
    dist.x1=dist.x[-c(i,j),c(i,j)]
    true.label=label.test[-c(i,j),2]
    assign=apply(dist.x1,1,function(k){
      ifelse(k[1]==k[2],sample(1:2,1),which.min(k))
    })
    x=table(assign,true.label)    
    if(length(rownames(x))==1){tmp[rownames(x),]=x
    }else{tmp=x}          
    out[[i.rep]]=tmp
  }
  x=lapply(out,NMI)
  mi=mean(unlist(x),na.rm=T)
  #out1=Reduce(`+`,out)
  #mi=NMI(out1)
  
  mi.out[[go.id]]=mi;
}
cat('cell.type',cell.type,'is done\n')
mi.out



