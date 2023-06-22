library(R.matlab)
source('pipeline_codes/src_featureNets.R')
source('pipeline_codes/src_SNMF.R')
source('pipeline_codes/src_SNMFforView.R')
source('pipeline_codes/src_multiViewNMF.R')
source('pipeline_codes/src_moduleNodesSelection.R')
source('pipeline_codes/src_setSimilarity.R')
source('pipeline_codes/src_significantModules.R')
source('pipeline_codes/src_significantModules_diff.size.R')
source('pipeline_codes/src_fly.gene_GOenrich.R')

#for(num.of.latent.factors in c(15,17,20,25,30)){
#for(num.of.latent.factors in c(17,20,25,30)){ 
for(num.of.latent.factors in c(25)){  
  library(R.matlab)
  code.files=Sys.glob('pipeline_codes/src*R')
  for(i in code.files) source(i)
  
  numCores=12;
  permu_times=200;
  #num.of.latent.factors=25;# sc=17,sn=25
  T_cutoff=3;
  input.folder='./brain_scRNA-seq/';
  out.folder=paste0('./sc_useall_n10c0.1_T',T_cutoff,'/');
  #input.folder='./brain_snRNA-seq_n10c0.1/';
  #out.folder=paste0('./sn_useall_n10c0.1_T',T_cutoff,'/');
  
  if(!dir.exists(out.folder)) dir.create(out.folder)
  
  ## read in genes
  genes=readRDS(paste0(input.folder,'all_common_genes.rds'))
  c.genes=unlist(genes$genes)
  all.genes=genes$all.genes
  common.genes=genes$common.genes
  #gene.names=genes$all.genes
  gene.names=names(which(table(c.genes)>=length(genes$genes)/2))
  gene_count=length(gene.names) #3015 for sc, 1351 for sn
  cat('use',gene_count,'genes in feature matrix construction\n');
  ## generate 2 feature matrix
  feature1=paste0(out.folder,'/strength_common.mat');
  feature2=paste0(out.folder,'/uniformity_common.mat');
  
  if(any(!file.exists(feature1),!file.exists(feature2))){
    
    files=Sys.glob(paste0(input.folder,'/pearson_coexpr/*rds'));
    cell.cluster.names=as.character(unlist(sapply(files,function(i){
      gsub('ncell_\\d+_','',basename(i),'_pearson.mat')
    })))
    cat('read in nets',length(files),'\n')
    
    min.cell=100;max.cell=Inf;
    nets=list();
    for(file in files){
      mat=abs(readRDS(file))
      ncell=as.numeric(strsplit(basename(file),'ncell\\_|\\_')[[1]][[2]])
      if( any(ncell>max.cell , ncell<min.cell) ){next}
      diag(mat)=0;
      x=intersect(gene.names,rownames(mat))
      nets[[basename(file)]]=mat[x,x];
    }
    
    network_count=length(nets); #number of input gene co-expression nets 
    
    nets.transform=list()
    for(k in 1:network_count){
       weight_max = max(nets[[k]]);
       weight_min = 0;
       matrix_weight_1 = (nets[[k]] - weight_min) / (weight_max-weight_min);
       matrix_weight_2 = 1/(1+exp(log(9999)-2*log(9999)*matrix_weight_1));    
       #matrix_weight_2[matrix_weight_1 <= 0.3] = 0;
       nets.transform[[k]]=matrix_weight_2
    }
    rm('nets')
    gc()
    
    # compute strength matrix
    all.weigh=matrix(0,nrow=gene_count,ncol=gene_count)
    rownames(all.weigh)=colnames(all.weigh)=gene.names
    for(i in 1:length(nets.transform)){
      net=nets.transform[[i]];
      all.weigh[rownames(net),colnames(net)]=all.weigh[rownames(net),colnames(net)]+net
    }
    strength=all.weigh/network_count
    diag(strength)=0;
    writeMat(feature1,strength=strength)
    cat('feature: strength is done\n')
    rm('all.weigh');rm('strength');
    gc()
    
    # compute uniformity matrix
    ref.mat=matrix(0,nrow=gene_count,ncol=gene_count)
    rownames(ref.mat)=colnames(ref.mat)=gene.names
    uni.nets=array(0,dim=c(network_count,gene_count,gene_count)) #number of matrix, nrow, ncol
    for(i in 1:network_count){
      net=nets.transform[[i]];
      tmp.mat=ref.mat;
      tmp.mat[rownames(net),colnames(net)]=net;
      uni.nets[i,,]=tmp.mat
    }
    max.weigh <- apply(uni.nets, c(2,3), max)
    dim(uni.nets);dim(max.weigh)
    #max.weigh[1:3,1:3]
    
    C.nets=list();
    for(i in 1:network_count){
      mat=uni.nets[i,,] #all pos values
      C.nets[[i]]=1-(mat/max.weigh)^2 
      #C.nets[[i]]=1-(mat/max.weigh) 
    }
    
    all.C=Reduce(`+`,C.nets)
    uni=1 - all.C/(network_count-1) #simlify to xi/max-1 / n-1, as xi/max>=1 as there is at least one i whose value == max
    diag(uni)=0;
    sum(is.na(uni))
    sum(is.nan(uni))
    sum(is.infinite(uni))
    uni[is.nan(uni)]=0;
    uni[is.infinite(uni)]=0;
    if(min(uni)<0){ stop('entries in uniformity <0\n')}
    uniformity=uni
   
    writeMat(feature2,uniformity=uniformity)
    cat('feature: uniformity is done\n')
    rm('nets.transform');rm('uni.nets');
    rm('C.nets');rm('all.C');
    gc();
  }
  
  #################################################
  # go to matlab to get Hc matrix by multi-view NMF
  # ConMod_Matlab.m
  #################################################
  
  #################################################
  ## or stay in R, slowerly
  Hc.outfile=paste0(out.folder,'K',num.of.latent.factors,'_Hc.mat');
  if(!file.exists(Hc.outfile)){
    if(!any(ls()=='strength')){
      x=readMat(feature1)
      strength=x$strength
    }
    if(!any(ls()=='uniformity')){
      x=readMat(feature2)
      uniformity=x$uniformity
    }
      
    #cat('similarity:',MESS::cmd(uniformity,strength),'\n'); #0.3
    
    start=Sys.time()
    #X=list(strength=strength,information=information)
    X=list(strength=strength,uniformity=uniformity)
    lambda = c(0.01, 0.01);
    #X=list(information=information,participate=participate)
    #lambda = c(0.01, 0.01);
    maxIter=50;
    K=num.of.latent.factors;
    Hc.out=multiViewNMF( X, num.of.latent.factors, lambda, maxIter)
    names(Hc.out)
    end=Sys.time()
    print(end-start)
    writeMat(Hc.outfile,Hc.mat=Hc.out$Hc)
  }
  
  #################################################
  if(!any(ls()=='Hc.out')) Hc.out=readMat(Hc.outfile)
  
  Hc=Hc.out$Hc
  dim(Hc) #ngene by K
  sum(Hc<0) #0, non-negative matrix
  
  
  xita = T_cutoff; #cutoff for node seletion in a module
  modules_final = moduleNodesSelection( Hc, xita ); #module merge in this step
  length(modules_final) #sn=17,
  sapply(modules_final,length)
  tmp=setSimilarity(modules_final)
  cat('max.similarity',max(tmp),'\n')

  modules_final_genes=lapply(modules_final,function(i) 
    gene.names[unlist(i)])
  sapply(modules_final_genes,length)
  cat('#module',length(modules_final_genes),
  '#gene',length(unique(unlist(modules_final_genes))),'\n'); #533
  
  saveRDS(modules_final_genes,
          paste0(out.folder,'K',num.of.latent.factors,'_module.genes.rds'))
  gc()
  
  ## Module validation
  if(!any(ls()=='nets')){
    files=Sys.glob(paste0(input.folder,'/pearson_coexpr/*rds'));
    cell.cluster.names=as.character(unlist(sapply(files,function(i){
      gsub('ncell_\\d+_','',basename(i),'_pearson.mat')
    })))
    cat('read in nets',length(files),'\n')
    
    min.cell=100;max.cell=Inf;
    nets=list();
    for(file in files){
      mat=abs(readRDS(file))
      ncell=as.numeric(strsplit(basename(file),'ncell\\_|\\_')[[1]][[2]])
      if( any(ncell>max.cell , ncell<min.cell) ){next}
      diag(mat)=0;
      x=intersect(gene.names,rownames(mat))
      nets[[basename(file)]]=mat[x,x];
    }
  }
  sapply(nets,dim)
  num_Nodes=nrow(Hc)
  
  #module.in.nets=significantModules_diff.size(modules_final_genes, nets,num_Nodes,permu_times=permu_times);
  #system.time({
    module.in.nets=significantModules_diff.size(modules=modules_final_genes, 
            multiNetworks=nets,
            permu_times=permu_times,numCores=numCores)
  #})
  
  saveRDS(module.in.nets,
          paste0(out.folder,'K',num.of.latent.factors,'_module.fdr.rds'))
  
  names(module.in.nets)
  dim(module.in.nets$FDR) #100 module by 41 nets
  
  if(FALSE){
  #if(max(module.in.nets$FDR)>0){
  colnames(module.in.nets$FDR)=cell.cluster.names
  rownames(module.in.nets$FDR)=paste0('module',1:nrow(module.in.nets$FDR))
  if(length(grep('sc',input.folder))>0){
    pdf('sc_module_by_cell.type.pdf',useDingbats = T,height = 16)
  }else{
    pdf('sn_module_by_cell.type.pdf',useDingbats = T,height = 16)
  }
  print(pheatmap::pheatmap(t(module.in.nets$FDR)))
  dev.off()
  }
  
  i=apply(module.in.nets$FDR,1,function(i) sum(i<0.05))
  i;
  M=ncol(module.in.nets$FDR);
  cat(sum(i==M),'common modules detected\n')
  
  if(sum(i==M)>0){
    common.modules.genes=modules_final_genes[i==M]
    sapply(common.modules.genes,length)
    common.modules.genes[[1]]
    
    saveRDS(common.modules.genes,
            paste0(out.folder,'K',num.of.latent.factors,'_common.modules.genes.rds'))
  }
  
  #################################################
  ## GO enrich
  #GO.out=lapply(1:length(common.modules.genes),function(i){
  GO.out=lapply(1:length(modules_final_genes),function(i){
    genes=modules_final_genes[[i]]
    x=GOenrich(genes)
    x
  })
  #sapply(GO.out,nrow)
  #head(GO.out[[1]]$Description)
  for(i in 1:length(GO.out)){
    cat('module',i,'\n')
    print(head(GO.out[[i]]$Description))
  }
  saveRDS(GO.out,
          paste0(out.folder,'K',num.of.latent.factors,'_modules.GO.rds'))
  
  i=apply(module.in.nets$FDR,1,function(i) sum(i<0.05))
  M=ncol(module.in.nets$FDR);
  if(sum(i==M)>0){
    common.GO.out=GO.out[i==M]
    saveRDS(common.GO.out,
          paste0(out.folder,'K',num.of.latent.factors,'_common.modules.GO.rds'))
  }
  rm(list=setdiff(ls(), "nets"))
}


#############################################

