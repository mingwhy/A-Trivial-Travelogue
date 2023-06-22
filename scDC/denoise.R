
# this function takes a expression UMI count matrix and a confidence score cutoff (0~1)
# after fitting a mixture model to each gene separately, it assigns a confidence score to each matrix element.
# based on the input confidence cutoff value, a denoised matrix is generate.
# this function returns a list,
# denoised matrix: denoised.mat=input.mat,
# fitted_paras: mixture model fitted parameters for each gene
# droprate: droprate per gene per cell
# the mixture model was developed in: Li, Wei Vivian, and Yanzeng Li. "sclink: Inferring sparse gene co-expression networks from single-cell expression data." Genomics, proteomics & bioinformatics 19.3 (2021): 475-492.

# other input parameter: 
# gene.meta, a data.frame, one gene per row, gene symbol, chr.poi, etc
# ncores: the number of CPU cores used in mclapply function from `get_mix_parameters.R`.
# plot: if you want to save some middel result plots and they would be saved in the plot.fig file.

# usage:
# x=readRDS('~/Documents/Data_Jay_fly_development/SupplementaryOnlineMaterial/seurat_objects/pred_windows/NNv1/18_20_finished_processing.rds')
# table(x$seurat_clusters) # cell# per cluster
# umi.mat=x@assays$RNA@counts
# input.mat=umi.mat[,x$seurat_clusters=='17']
# gene.meta=data.table::fread('validated_17831genes.txt') 
# start_time <- Sys.time()
# res_out=denoise(input.mat=input.mat,conf.cutoff=0.999,gene.meta=gene.meta,ncores=1)
# end_time <- Sys.time()
# end_time-start_time # ncores=1, 360 cells takes 48secs
# names(res_out);
# denoise.mat=res_out$denoised.mat


library(parallel) #for mclapply function used in `get_mix_parameters.R`
source('simu_gamma_norm.R')
source('get_mix_parameters.R')
source('src_plot_one_cell.R')

denoise<-function(input.mat=input.mat,conf.cutoff=0.999,gene.meta=gene.meta,
                  ncores=1,plot=TRUE,plot.fig='test.pdf'){
  
  #remove non-expressed genes in the input matrix
  input.mat=input.mat[Matrix::rowSums(input.mat)!=0,] 
  sce <- SingleCellExperiment(list(counts=input.mat))
  
  # use computeSumFactors(https://rdrr.io/bioc/scran/man/computeSumFactors.html) to do cell-wise normalization 
  clusters <- quickCluster(sce)
  sce <- computeSumFactors(sce, clusters=clusters)
  sce.count <- normalizeCounts(sce, log=FALSE)
  
  log10.sce.1.mat=log10(sce.count+1)
  log10.sce.1.01.mat=log10(sce.count+1.01)
  
  ## fit mixture model to each gene
  # get_mix_parameters() return rate as the proportion of gamma in the mixture.
  log.mat=log10.sce.1.01.mat; #gene by cell matrix
  pa = get_mix_parameters(log.mat, ncores = ncores) #using log.mat for fit
  #dim(pat) #one row per gene
  rownames(pa)=rownames(log.mat)
  
  # plot two cell examples (one cell one row, left:obs, right:simu based on estiamted parameter values)
  if(plot){
    pdf(plot.fig,useDingbats = T)
    tmp=pa[!is.na(pa[,1]),]
    pick.gene=rownames(tmp[sample(1:nrow(tmp),2,replace = F),])
    gene.para=pa[rownames(pa) %in% pick.gene,] #estiamted para values for gene.i
    
    # simulate distribution based on estiamted parameters
    # as in r_gamma_normal(), pi is the proportion of Normal in the mixture.
    # get_mix_parameters() return rate as the proportion of gamma in the mixture.
    plot.mat=log10.sce.1.mat
    
    par(mfrow=c(2,2))
    for(i in 1:nrow(gene.para)){
      simu.dat=r_gamma_normal(n=nrow(log.mat),pi=gene.para[i,1],
                              mu=gene.para[i,4],sigma2=gene.para[i,5]^2,alpha=gene.para[i,2],beta=gene.para[i,3])
      
      x=plot.mat[pick.gene[i],]
      hist.data = hist(x, plot=F)
      hist.data$counts = log10(hist.data$counts+1)
      plot(hist.data, ylab='log10(Frequency)',main=paste0('gene ',pick.gene[i],'\nobs.data, log10(cpm+1)'),
           xlab='obs.data')
      
      x=simu.dat
      hist.data = hist(x, plot=F)
      hist.data$counts = log10(hist.data$counts+1)
      plot(hist.data, ylab='log10(Frequency)',
           main=paste0('gene ',pick.gene[i],'\nsimu.data, rate=',round(gene.para[i,1],2)));
    }
  }
  
  ## calculate cell-wise dropout rate for each gene
  I = nrow(log.mat) #number of genes
  J= ncol(log.mat) #number of cells
  # calculate_weight<-function (x, paramt): https://github.com/Vivianstats/scLink/blob/master/R/get_mix_parameters.R
  # pz1 = paramt[1] * dgamma(x, shape = paramt[2], rate = paramt[3])
  # pz2 = (1 - paramt[1]) * dnorm(x, mean = paramt[4], sd = paramt[5])
  # pz = pz1/(pz1 + pz2)
  
  droprate = sapply(1:I, function(i) {
    if(is.na(pa[i,1])) return(rep(0,J))
    wt = calculate_weight(as.numeric(log.mat[i, ]),as.numeric(pa[i, ]) )
    return(wt[, 1])
  })
  droprate=t(droprate) #gene by cell
  colnames(droprate)=colnames(log.mat)
  rownames(droprate)=rownames(log.mat)
  
  # transform droprate into confidence score
  confidence = 1-droprate;
  if(plot){
    tmp=apply(confidence,1,function(i) sum(i>0.1 & i<0.9)) #per gene
    tmp1=order(tmp,decreasing = T)
    pick.gene=tmp1[1:2]
    
    tmp=apply(confidence,2,function(i) sum(i>0.1 & i<0.9)) #per gene
    tmp1=order(tmp,decreasing = T)
    pick.cell=tmp1[1:2]
      
    plot.mat=log10.sce.1.mat; #cell by gene
    par(mfrow=c(2,2))
    plot(jitter(plot.mat[,pick.cell[1]]),confidence[,pick.cell[1]],pch=16,col=rgb(0,0,1,alpha=0.5),
         xlab='log10(cpm+1)',ylab='measure confidence',main='cell 1',cex.lab=1.5)
    plot(jitter(plot.mat[,pick.cell[2]]),confidence[,pick.cell[2]],pch=16,col=rgb(0,0,1,alpha=0.5),
         xlab='log10(cpm+1)',ylab='measure confidence',main='cell 2',cex.lab=1.5)
    
    plot(jitter(plot.mat[pick.gene[[1]],]),confidence[pick.gene[[1]],],pch=16,col=rgb(0,0,1,alpha=0.5),
         xlab='log10(cpm+1)',ylab='measure confidence',main='gene 1',cex.lab=1.5)
    plot(jitter(plot.mat[pick.gene[[2]],]),confidence[pick.gene[[2]],],pch=16,col=rgb(0,0,1,alpha=0.5),
         xlab='log10(cpm+1)',ylab='measure confidence',main='gene 2',cex.lab=1.5)
    
    cutoffs=c(0,0.6,0.9,0.999)  
    plots=list();iplot=0
    for(test.cutoff in cutoffs){
      i.cell=plot.mat[,pick.cell[[1]]]
      f=confidence[,pick.cell[[1]]]
      keep=1*(f>test.cutoff)
      i.cell.keep=i.cell[keep==1]
      df.gene.tmp=gene.meta[gene.meta$SYMBOL %in% names(i.cell.keep),]
      xs=i.cell.keep[which(df.gene.tmp$chromosome_name == 'X')]
      as=i.cell.keep[df.gene.tmp$chromosome_name %in% c('2L','2R','3L','3R','4')]
      iplot=iplot+1
      plots[[iplot]]=plot_one_cell_hist(i,xs,as,test.cutoff) #from src_plot_one_cell.R
    }
  
    grid.arrange(grobs=plots,ncol=2,nrow=2)
    dev.off()
  }
  
  # apply confience cutoff to generate final output matrix
  keepind = 1* (confidence > conf.cutoff) #1,keep
  keepind[which(input.mat==0,arr.ind = T)]=0
  cat('number of 0 in the input matrix ',  sum(input.mat==0,na.rm=T),'\n')
  cat('denoise remove',sum(keepind==0,na.rm=T),'elements\n')
  input.mat[keepind==0]=NA
  return(list(denoised.mat=input.mat, fitted_paras=pa, droprate=droprate))
}
  
  
