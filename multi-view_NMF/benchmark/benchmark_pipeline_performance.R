
## generate simu nets
#library(R.matlab)
source('src_syn_dataset_common.R')
source('src_evaluation.R')
source('src_syn_dataset_overlap.R')

(files=Sys.glob('../brain_core_scran/pipeline_codes/*'))
for(file in files){source(file)}
numCores=12;

nrep=20; #7hrs
for(alpha.i in c(0.1,0.3,0.5,0.7,0.9)){
#for(alpha.i in c(0,0.01,0.05,seq(0.1,1,0.1))){
  T_cutoff=1.5; #1.5
  #dir.create(paste0('common_T',T_cutoff))
  #out.file=paste0('common_T',T_cutoff,'/common_alpha_',alpha.i,'.txt')
  
  dir.create(paste0('overlap_T',T_cutoff))
  out.file=paste0('overlap_T',T_cutoff,'/overlap_alpha_',alpha.i,'.txt')
  
  rep.out=list();rep.i=0; 
  while(length(rep.out)<=nrep){
    #out=syn_dataset_common(alpha=alpha.i);K=5;
    out=syn_dataset_overlap(alpha=alpha.i);K=2;
    names(out)
    sapply(out$dataset,dim)
    
    #out=syn_dataset_overlap(alpha=0.6) 
    
    if(F){
      my.br=c(0,seq(0.1,1,0.2));
      my.col = colorRampPalette((RColorBrewer::brewer.pal(n = 8, name = "YlOrRd")))(length(my.br)-1)
      my.col[1]='white'
      barplot(1:length(my.col),col=my.col)
      x=1:nrow(out$dataset[[1]]) 
      y=1:ncol(out$dataset[[1]])
      #pdf('syn_common.nets.pdf',useDingbats = T)
      pdf('syn_over.nets.pdf',useDingbats = T)
      sapply(1:length(out$dataset),function(i)
        image(y,x,out$dataset[[i]],breaks=my.br,axes=FALSE, col=my.col,ylab="", xlab=""))
      dev.off()
    }
    
    realLabels=out$realLabels
    labels_specific=out$labels_specific
    dim(labels_specific)
    
    simu.nets=out$dataset
    
    ## detect modules, strength and uniformity
    multiNetworks=simu.nets;
    network_count=length(multiNetworks) #number of simulated matrix
    gene_count=nrow(multiNetworks[[1]])
    
    
    # data transform
    if(T){ #yes transform
      multiNetworks.transform=list()
      for(k in 1:network_count){
        weight_max = max(multiNetworks[[k]]);
        weight_min = 0;
        matrix_weight_1 = (multiNetworks[[k]] - weight_min) / (weight_max-weight_min);
        matrix_weight_2 = 1/(1+exp(log(9999)-2*log(9999)*matrix_weight_1));
        
        #matrix_weight_2[matrix_weight_1 <= 0.3] = 0;
        multiNetworks.transform[[k]]=matrix_weight_2
      }
    }
    
    #multiNetworks.transform=multiNetworks;
    # compute strength matrix
    all.weigh=Reduce(`+`,multiNetworks.transform)
    strength=all.weigh/network_count
    diag(strength)=0;
    
    # compute uniformity matrix
    tmp=array(0,dim=c(network_count,gene_count,gene_count)) #number of matrix, nrow, ncol
    for(i in 1:network_count){
      tmp[i,,]=multiNetworks.transform[[i]]
    }
    max.weigh <- apply(tmp, c(2,3), max)
    dim(max.weigh)
    max.weigh[1:3,1:3]
    
    C.nets=list();
    for(i in 1:network_count){
      mat=multiNetworks.transform[[i]] #all pos values
      C.nets[[i]]=1-(mat/max.weigh)^2 #no need square
      #C.nets[[i]]=1-(mat/max.weigh)
    }
    
    all.C=Reduce(`+`,C.nets)
    uni=1 - all.C/(network_count-1)
    diag(uni)=0;
    uni[is.infinite(uni)]=0;
    uni[is.nan(uni)]=0;
    if(min(uni)<0){ stop('entries in uniformity <0\n')}    
    uniformity=uni
    
    #library(MESS)
    #cat('similarity:',MESS::cmd(uniformity,strength),'\n'); #0.3
    
    
    ## mvNMF
        
    start=Sys.time()
    X=list(strength=strength,uniformity=uniformity)
    #X=list(strength=strength,information=information) #working
    #X=list(strength=strength,participate=participate) #working
    #X=list(participate=participate,information=information) #not working at all
    lambda = c(0.01, 0.01);
    maxIter=50;
    #K=5;
    out=multiViewNMF( X, K, lambda, maxIter)
    names(out)
    end=Sys.time()
    print(end-start)
    
    ## select node to become modules
    Hc=out$Hc
    dim(Hc)
    sum(Hc<0) #0, non-negative matrix
    
    xita = T_cutoff; #cutoff for node seletion in a module
    #modules_final = moduleNodesSelection( Hc, xita ); 
    r=tryCatch(
      modules_final <- moduleNodesSelection( Hc, xita ),
      error=function(c) 'error')    
    if(r[1]=='error'){cat('no gene pass xita,redo this rep\n');next }

    length(modules_final)
    sapply(modules_final,length)
    
    ## evaluate
    num_Nodes=nrow(Hc)
    if(length(modules_final)==0){next}
    out1=evaluation(modules_final, realLabels, num_Nodes)
    out1
    
    ## see what's their significance in each simulated net
    module.in.nets=significantModules(modules_final, simu.nets,
      permu_times=100,numCores=numCores);
    
    names(module.in.nets)
    out2=apply(module.in.nets$FDR,1,function(i) sum(i<0.05)) ##modules detected in each net
    names(out2)=rep('n.net',length(out2))
    out2
    
    x=data.frame('name'=c(names(out1),names(out2)),
               'value'=c(as.numeric(out1),as.numeric(out2)))
    x=rbind(x,c('alpha',alpha.i))
    rep.i=rep.i+1;
    x$rep=rep.i
    rep.out[[rep.i]]=x
    #write.table(x,file=paste0('alpha_',alpha.i,'.txt'),quote=F,row.names=F)
  }
  rep.df=Reduce(`rbind`,rep.out)
  write.table(rep.df,file=out.file,sep='\t',quote=F,row.names=F)
}

###########################################
if(F){
library(tidyverse)
library(gridExtra)
  
plot.list=list();i=0;
dir.name=Sys.glob('noise*')
for(x.dir in dir.name){
  g.mean=strsplit(x.dir,'_')[[1]][2]
  g.var=strsplit(x.dir,'_')[[1]][3]
  for(y in c('common_T1.5','overlap_T1.5')){
    files=Sys.glob(paste0(x.dir, '/',y,'/*txt'))
    #files=Sys.glob(paste0(x.dir, '/common_T1.5/*txt'))
    #files=Sys.glob('noise_0.25_0.1/overlap_T1.5/*txt')
    out=lapply(files,function(file){
      df=read.table(file,header=T)
      alpha=df[df$name=='alpha',]$value[1]
      df=df[!df$name %in% c('alpha','n.net'),]
      x=df %>% group_by(name) %>% summarise(n=n(),mean=mean(value),sd=sd(value))
      x$alpha=alpha;
      x
    })
    
    if(length(grep('common',y))==1){
      type='common'
    }else{
      type='overlapping'
    }
    plot.title=paste0('Gaussian noise with mean ',g.mean,', var ',g.var,'\n',
                      type,' module simulation scheme');
    df.x=as.data.frame(Reduce(`rbind`,out))
    df.x$name=factor(df.x$name,levels=c('TPR','FPR','Accuracy','MCC'))
    i=i+1;
    plot.list[[i]]<- ggplot(subset(df.x,name!='Accuracy'),aes(x=alpha,y=mean))+geom_point()+
      facet_wrap(.~name,ncol=4)+theme_classic()+
        xlab('Alpha')+ylab('')+
      scale_x_continuous(breaks=seq(0.1,0.9,0.2),
                        labels=seq(0.1,0.9,0.2))+
      geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=0.05,
                    position=position_dodge(0))+
      #ggtitle(dirname(files[1]))
        ggtitle(plot.title)
      
    }
  
}

pdf('simulated.network_plot.pdf',useDingbats = T,height = 6,width = 8)
grid.arrange(arrangeGrob(grobs= plot.list,ncol=2))
dev.off()

tiff("simulated.network_plot.tif", width = 500, height =800,
     units = 'in', res = 300)
grid.arrange(arrangeGrob(grobs= plot.list,ncol=1))
dev.off()

}
