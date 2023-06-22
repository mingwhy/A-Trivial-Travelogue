
#source('cal_bigScale2.R') 
#source('cal_CDI.R')
#source('cal_BaCo.R')
library(scLink)
library(Seurat)
library(scWGCNA)

cal_cor.mat<-function(method='pearson',m1,m2){

  if(sum(method %in% c('pearson','propr.phi_s',
                       'propr.rho_p','baco',
                       'sclink','bigscale','CDI','pseudocell'))==0){
    cat('method not matched\n')
    break
  }
  if(method=='pearson'){
    #cor.mat1= Rfast::cora(t(m1)) #cor.mat
    #cor.mat2= Rfast::cora(t(m2)) #cor.mat
    cor.mat1= coop::tpcor(m1) #faster
    cor.mat2= coop::tpcor(m2) #faster
  }else if(method=='propr.phi_s'){
    cor.mat1 = -1.0 * propr::phis(t(m1))@matrix
    cor.mat2 = -1.0 * propr::phis(t(m2))@matrix
  }else if(method=='propr.rho_p'){
    cor.mat1 = propr::perb(t(m1))@matrix
    cor.mat2 = propr::perb(t(m2))@matrix
  }else if(method=='baco'){
    ## baco
    cor.mat1= BaCo(m1)
    cor.mat2= BaCo(m2)
  }else if(method=='sclink'){
    ## sclink
    genes=rownames(m1)
    count.norm1 = sclink_norm(t(m1), scale.factor = 1e6, filter.genes = FALSE, gene.names = genes)
    count.norm2 = sclink_norm(t(m2), scale.factor = 1e6, filter.genes = FALSE, gene.names = genes)
    cor.mat1 = sclink_cor(expr = count.norm1, ncores = 2)
    cor.mat2 = sclink_cor(expr = count.norm2, ncores = 2)
  }else if(method=='bigscale'){
    ## bigscale 
    gene.names=rownames(m1)
    net.out1=compute.network(expr.data = m1,gene.names = gene.names,clustering='direct',speed.preset = "fast")
    net.out2=compute.network(expr.data = m2,gene.names = gene.names,clustering='direct',speed.preset = "fast")
    #net.out1=compute.network(expr.data = m1,gene.names = gene.names,clustering='recursive',speed.preset = "fast")
    #net.out2=compute.network(expr.data = m2,gene.names = gene.names,clustering='recursive',speed.preset = "fast")
    cor.mat1=net.out1$Dp.mat
    cor.mat2=net.out2$Dp.mat
  }else if(method=='CDI'){
    cor.mat1=CDI(m1)
    cor.mat2=CDI(m2)
  }else if(method=='pseudocell'){    
     mydata1 <- CreateSeuratObject(count=m1,project = "group1")
     mydata2 <- CreateSeuratObject(count=m2,project = "group2")
     mat.out=list();
     for(i in 1:2){
        if(i==1) cell=mydata1 else cell=mydata2
        cell <- NormalizeData(cell, normalization.method = "LogNormalize", scale.factor = 10000)
        cell <- FindVariableFeatures(cell, selection.method = "vst", nfeatures = 200)
        all.genes <- rownames(cell)
        cell <- ScaleData(cell, features = all.genes)
        cell <- RunPCA(cell, features = VariableFeatures(object = cell))
        # now cell@reductions has $pca 
        
        cell  = RunUMAP(cell,dims=1:20)
        # generate pseudo cells
        ps.cell=construct_metacells(cell, k=10,reduction='umap')
        #ps.cell@assays$RNA@counts[1:3,1:3]
        #ps.cell@assays$RNA@data[1:3,1:3]
        
        #mat.out[[i]]=Rfast::cora(t(as.matrix(ps.cell@assays$RNA@counts)))    
        mat.out[[i]]=coop::tpcor(as.matrix(ps.cell@assays$RNA@counts))
    }
    cor.mat1=mat.out[[1]];cor.mat2=mat.out[[2]]
  }
  return(list(cor.mat1,cor.mat2))
}

cal_agreement<-function(topk=10000,cor.mat1,cor.mat2){
  # top 50%,40,30,20,10,5,1, edge.overlap/n.edge
  x1=cor.mat1[upper.tri(cor.mat1)]
  x2=cor.mat2[upper.tri(cor.mat2)]
  right1=sort(x1,decreasing = T)[topk]
  right2=sort(x2,decreasing = T)[topk]
  left1=sort(x1,decreasing = F)[topk]
  left2=sort(x2,decreasing = F)[topk]
 
  tmp1=matrix(0,nrow=nrow(cor.mat1),ncol=ncol(cor.mat1))
  tmp2=tmp1
  tmp1[cor.mat1>right1]=1
  tmp2[cor.mat2>right2]=1
  diag(tmp1)=0;diag(tmp2)=0
  #sum(tmp1==1);sum(tmp2==1)
  tmp3=tmp1+tmp2
  aggrement1=sum(tmp3==2)/sum(tmp1==1)
  
  tmp1=matrix(0,nrow=nrow(cor.mat1),ncol=ncol(cor.mat1))
  tmp2=tmp1
  tmp1[cor.mat1<left1]=1
  tmp2[cor.mat2<left2]=1
  diag(tmp1)=0;diag(tmp2)=0
  #sum(tmp1==1);sum(tmp2==1)
  tmp3=tmp1+tmp2
  aggrement2=sum(tmp3==2)/sum(tmp1==1)
  return(c(aggrement1,aggrement2))
}

splitData <- function(df, propShared=0.5){
  sampleCount <- ncol(df)
  colIdx <- sample.int(sampleCount)
  sharedSampleCount <- round(propShared*sampleCount)
  sharedColsIdx <- 1:sharedSampleCount
  specificSampleCount <- floor(0.5*(sampleCount-sharedSampleCount))

  s1=sharedSampleCount+1;
  e1=sharedSampleCount+1+specificSampleCount;
  s2=e1+1;
  e2=min(s2+specificSampleCount,sampleCount);
  df1 <- df[,colIdx[c(sharedColsIdx, s1:e1) ]]
  df2 <- df[,colIdx[c(sharedColsIdx, s2:e2)]]
  dfList <- list(df1, df2)
  return(dfList)
}

## if you want to use scran to normalize data
library(scran)
scran_norm<-function(expr.mat){
  sce <- SingleCellExperiment(list(counts=expr.mat),
                              #colData=DataFrame(cell.type=cell.type),
                              rowData=DataFrame(gene=rownames(expr.mat)) )
  #sce
  clusters <- quickCluster(sce)
  sce <- computeSumFactors(sce, clusters=clusters)
  #summary(sizeFactors(sce))
  sce <- logNormCounts(sce)
  log.sce=sce@assays@data$logcounts
  #dim(expr.mat);dim(log.sce);
  return(log.sce)
}


