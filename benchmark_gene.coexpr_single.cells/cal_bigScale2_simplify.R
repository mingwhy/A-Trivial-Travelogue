
#' Gene regulatory network
#'
#' Infers the gene regulatory network from single cell data
#'
#' @param expr.data matrix of expression counts. Works also with sparse matrices of the \pkg{Matrix} package.
#' @param gene.names character of gene names, now it supports Gene Symbols or Ensembl, Mouse and Human.
#' @param clustering type of clustering and correlations computed to infer the network.
#' \itemize{
#'  \item {\bold{recursive}} Best quality at the expenses of computational time. If the dataset is larger than 10-15K cells and is highly heterogeneous this can lead to very long computational times (24-48 hours depending of the hardware).
#'  \item {\bold{direct}} Best trade-off between quality and computational time. If you want to get a quick output not much dissimilar from the top quality of \bold{recursive} one use this option. Can handle quickly also large datasets (>15-20K cells in 30m-2hours depending on hardware)
#'  \item {\bold{normal}} To be used if the correlations (the output value \bold{cutoff.p}) detected with either \bold{direct} or \bold{recursive} are too low. At the moment, bigSCale displays a warning if the correlation curoff is lower than 0.8 and suggests to eithe use \bold{normal} clustering or increase the input parameter \bold{quantile.p}
#' }
#' @param quantile.p To select how many correlatoins are retained. By default \eqn{quantile.p=0.9} meaning that only Pearson Correlatins>0.9 are retained to create the network. If the number is higher than 1 then it is interpreted as a quantile, meaning tha only the first \eqn{1 - quantile.p} correlations are used to create the edges of the network. 
#' If \eqn{0<quantile.p<1} ten it is interprested directly as the minimum Pearson correlation. If the networ is too sparse(dense) decrease(increase) \eqn{quantile.p}
#' @param speed.preset Used only if  \code{clustering='recursive'} . It regulates the speed vs. accuracy of the Zscores calculations. To have a better network quality it is reccomended to use the default \bold{slow}.
##' \itemize{
#'   \item {\bold{slow}} {Highly reccomended, the best network quality but the slowest computational time.} 
#'   \item {\bold{normal}} {A balance between network quality and computational time. }
#'   \item {\bold{fast}} {Fastest computational time, worste network quality.}
#' }
#' @param previous.output previous output of \code{compute.network()} can be passed as input to evaluate networks with a different quantile.p without re-running the code. Check the online tutorial at https://github.com/iaconogi/bigSCale2.
#' 
#' @return  A list with the following items:
#' \itemize{
#' \item {\bold{centrality}} {Main output: a Data-frame with the network centrality (Degree,Betweenness,Closeness,PAGErank) for each gene(node) of the network}
#' \item {\bold{graph}} {The regulatory network in iGraph object}
#' \item {\bold{correlations}} {All pairwise correlations between genes. The correlation is an average between \emph{Pearson} and \emph{Spearman}. Note that it is stored in single precision format (to save memory space) using the package \pkg{float32}.To make any operation or plot on the correlations first transform it to the standard double precisione by running \code{correlations=dbl(correlations)} }
#' \item {\bold{cutoff.p}} {The adptive cutoff used to select significant correlations}
#' \item {\bold{tot.scores}} {The Z-scores over which the correlations are computed. The visually check the correlation between to genes \emph{i} and \emph{j} run  \code{plot(tot.scores[,i],tot.scores[,j])} }
#' \item {\bold{clusters}} {The clusters in which the cells have been partitioned}
#' \item {\bold{model}} {Bigscale numerical model of the noise}
#' }
#'
#'
#' @examples
#' out=compute.network(expr.data,gene.names)
#'
#' @export

  
compute.network = function (expr.data,gene.names,modality='pca',model=NA,clustering='recursive',quantile.p=0.9,speed.preset='slow',previous.output=NA){

  if (is.na(previous.output))  
    {
    
    if (ncol(expr.data)>20000)
      warning('It seems you are running compute.network on a kind of large dataset, it it failed for memory issues you should try the compute.network for large datasets')
    expr.data=as.matrix(expr.data)
    
    
    # Part 1) Initial processing of the dataset  ************************************************************
    print('Pre-processing) Removing null rows ')
    exp.genes=which(Rfast::rowsums(expr.data)>0)
    if ((nrow(expr.data)-length(exp.genes))>0)
      print(sprintf("Discarding %g genes with all zero values",nrow(expr.data)-length(exp.genes)))
    expr.data=expr.data[exp.genes,]
    gc()
    gene.names=gene.names[exp.genes]
    tot.cells=Rfast::rowsums(expr.data>0)
    
    print('PASSAGE 1) Setting the size factors ....')
    lib.size = Rfast::colsums(expr.data)
      
      
    
    if (is.list(model))
    {
    print('Using edges given as input')
    edges=model$edges
    }
    else
    {
      print('PASSAGE 2) Setting the bins for the expression data ....')
      edges=generate.edges(expr.data)
    }
      
    print('PASSAGE 3) Storing in the single cell object the Normalized data ....')
    avg.library.size=mean(lib.size)
    for (k in 1:ncol(expr.data)) expr.data[,k]=expr.data[,k]/lib.size[k]*avg.library.size
    expr.data.norm=expr.data
    rm(expr.data)
    expr.data.norm=Matrix::Matrix(expr.data.norm)
    gc()
    
    if (is.list(model))
      {
      print('Using Model given in the input')
      model=model$model
      }
    else
      if (!(modality=='pca' & clustering=='direct'))
      {  
      print('PASSAGE 4) Computing the numerical model (can take from a few minutes to 30 mins) ....')
      model=fit.model(expr.data.norm,edges,lib.size)$N_pct
      gc()
      }
      

      
    print('PASSAGE 5) Clustering ...')
    if (clustering=="direct" | clustering=="recursive")
      {
      if (clustering=="direct")
        mycl=bigscale.recursive.clustering(expr.data.norm = expr.data.norm,model = model,edges = edges,lib.size = lib.size,fragment=TRUE,modality=modality)$mycl

      if (clustering=="recursive")
        mycl=bigscale.recursive.clustering(expr.data.norm = expr.data.norm,model = model,edges = edges,lib.size = lib.size,modality=modality)$mycl
      }
    else
      {
      ODgenes=calculate.ODgenes(expr.data.norm)
      dummy=as.matrix(ODgenes[[1]])
      ODgenes=which(dummy[,1]==1)
      print('Computing distances ...')
      D=compute.distances(expr.norm = expr.data.norm,N_pct = model,edges = edges,driving.genes = ODgenes,lib.size = lib.size)
      mycl=bigscale.cluster(D,plot.clusters = FALSE,method.treshold = 0.5,verbose=FALSE)$clusters
      }
    tot.clusters=max(mycl)
    
    
    
    
    #filtering the number of genes
    #pass.cutoff=which(tot.cells>(max(15,ncol(expr.data.norm)*0.005)))
    pass.cutoff=which(tot.cells>0) #ming, use this line, keep all genes without removing, 20210925  
    gene.names=gene.names[pass.cutoff]
    
    
    if (clustering=="direct")
        {
        print(sprintf('Assembling cluster average expression for %g genes expressed in at least %g cells',length(pass.cutoff),max(15,ncol(expr.data.norm)*0.005)))
        tot.scores=matrix(0,length(pass.cutoff),tot.clusters)
        for (k in 1:(tot.clusters))
            tot.scores[,k]=Rfast::rowmeans(as.matrix(expr.data.norm[pass.cutoff,which(mycl==k)]))
        tot.scores=log2(tot.scores+1)
        }
      else
        {
        print(sprintf('Calculating Zscores for %g genes expressed in at least %g cells',length(pass.cutoff),max(15,ncol(expr.data.norm)*0.005)))  
        Mscores=calculate.marker.scores(expr.norm = expr.data.norm[pass.cutoff,], clusters=mycl, N_pct=model, edges = edges, lib.size = lib.size, speed.preset = speed.preset)$Mscores
    
        rm(expr.data.norm)
        gc()
        
        tot.scores=matrix(0,length(Mscores[[1,2]]),(tot.clusters*(tot.clusters-1)/2))
        
        # assembling the matrix of Mscores
        counts=1
        for (k in 1:(tot.clusters-1))
          for (h in (k+1):tot.clusters)
          {
            tot.scores[,counts]=Mscores[[k,h]]
            counts=counts+1 
          }
        }
    tot.scores=t(tot.scores)
    }
    
    else
    {
      print('It appears you want to tweak previously created networks with a different quantile.p, proceeding ....')
      tot.scores=previous.output$tot.scores
      gene.names=colnames(previous.output$correlations)
      mycl=previous.output$clusters
      model=previous.output$model
      rm(previous.output)
      gc()
    }
    
    




  o=gc()
  #print(o)
  print('Calculating Pearson ...')
  rm(list=setdiff(setdiff(ls(), lsf.str()), c('tot.scores','quantile.p','model','gene.names','tot.scores','mycl')))
  o=gc()
  #print(o)

  Dp=Rfast::cora(tot.scores)
  o=gc()
  #print(o)
    if(F){ #ming
    if (quantile.p>1)
      {
      print('Calculating correlation quantile associated with your quantile.p ...')
      cutoff.p=quantile(Dp,quantile.p/100,na.rm=TRUE)
      }
    else
      {
      print('Detacted quantile.p<1, using it directly as correlation treshold ...')
      cutoff.p=quantile.p
      }

    print(sprintf('Using %f as cutoff for pearson correlation',cutoff.p))


    Dp=float::fl(Dp)
    gc()
  }

  #print('Calculating Spearman ...') #ming, 2021-09-01
  #Ds=cor(x = tot.scores,method = "spearman") #ming, 2021-09-01
  if(F){ #ming
    Ds=float::fl(Ds)
    gc()

    print('Calculating the significant links ...')
    rm(list=setdiff(setdiff(ls(), lsf.str()), c("Dp","Ds","cutoff.p",'gene.names','tot.scores','mycl','model')))  
    gc()
    network=((Dp>cutoff.p & Ds>0.7) | (Dp<(-cutoff.p) & Ds<(-0.7)))
    diag(network)=FALSE
    network[is.na(network)]=FALSE
    degree=Rfast::rowsums(network)
    to.include=which(degree>0)

    G=igraph::graph_from_adjacency_matrix(adjmatrix = network[to.include,to.include],mode = 'undirected')
    G=igraph::set_vertex_attr(graph = G,name = "name", value = gene.names[to.include])
    rm(network)
    gc()


    print('Calculating the final score ...')
    Df=(Ds+Dp)/float::fl(2)
    rm(Dp)
    rm(Ds)
    gc()
    rownames(Df)=gene.names
    colnames(Df)=gene.names


    print(sprintf('Inferred the raw regulatory network: %g nodes and %g edges (ratio E/N)=%f',length(igraph::V(G)),length(igraph::E(G)),length(igraph::E(G))/length(igraph::V(G))))

    G=polish.graph(G)
    comp=igraph::components(G)
    small.comp=which(comp$csize<0.01*sum(comp$csize))
    to.remove=which(is.element(comp$membership,small.comp))
    G=igraph::delete_vertices(G, to.remove)
    comp=igraph::components(G)
    cat('\n')
    print(sprintf('Final network after GO filtering: %g nodes and %g edges (ratio E/N)=%f and %g components',length(igraph::V(G)),length(igraph::E(G)),length(igraph::E(G))/length(igraph::V(G)), comp$no))


    print('Computing the centralities')
    Betweenness=igraph::betweenness(graph = G,directed=FALSE,normalized = TRUE)
    Degree=igraph::degree(graph = G)
    PAGErank=igraph::page_rank(graph = G,directed = FALSE)$vector
    Closeness=igraph::closeness(graph = G,normalized = TRUE)

    G=igraph::set.vertex.attribute(graph = G,name = 'bigSCale.Degree',value = Degree)
    G=igraph::set.vertex.attribute(graph = G,name = 'bigSCale.PAGErank',value = PAGErank)
    G=igraph::set.vertex.attribute(graph = G,name = 'bigSCale.Closeness',value = Closeness)
    G=igraph::set.vertex.attribute(graph = G,name = 'bigSCale.Betweenness',value = Betweenness)

    if (cutoff.p<0.7)
      warning('bigSCale: the cutoff for the correlations seems very low. You should either increase the parameter quantile.p or select clustering=normal (you need to run the whole code again in both options,sorry!). For more information check the quick guide online')
  }

  rownames(Dp)=gene.names;colnames(Dp)=gene.names #ming
  #rownames(Ds)=gene.names;colnames(Ds)=gene.names ##ming, 2021-09-01
  #return(list(graph=G,correlations=Df,tot.scores=tot.scores,clusters=mycl,centrality=as.data.frame(cbind(Degree,Betweenness,Closeness,PAGErank)),cutoff.p=cutoff.p,model=model))
  return(list(#graph=G,
    Dp.mat=Dp,#Ds.mat=Ds,
    tot.scores=tot.scores,clusters=mycl,
    #centrality=as.data.frame(cbind(Degree,Betweenness,Closeness,PAGErank)),cutoff.p=cutoff.p,
    model=model)); #ming

}


bigscale.recursive.clustering = function (expr.data.norm,model,edges,lib.size,
          fragment=FALSE,create.distances=FALSE,modality){
    
  gc()
    
  num.samples=ncol(expr.data.norm)
    
  if (fragment==FALSE)
    {
    # Adjusting max_group_size according to cell number
    if (num.samples<1000) dim.cutoff=10
    if (num.samples>=1000 & num.samples<5000) dim.cutoff=50
    if (num.samples>=5000 & num.samples<10000) dim.cutoff=100
    if (num.samples>=10000) dim.cutoff=150
    }
  else
  {
    if (fragment==TRUE)
        {
        if (num.samples<1000) dim.cutoff=10
        else dim.cutoff=50
        }
    else
      {
      dim.cutoff=fragment*ncol(expr.data.norm)
      #cat(sprintf('\nClustering to groups of at most %g cells (%g %%)',dim.cutoff,fragment*100))
      }
      
  }
      

  print(sprintf('Clustering cells down to groups of approximately %g-%g cells',dim.cutoff,dim.cutoff*5))   
  #dim.cutoff=ncol(expr.data.norm)*min.group.size  

  mycl=rep(1,ncol(expr.data.norm))
  tot.recursive=1
  #current.cutting=40
  unclusterable=rep(0,length(mycl))

  if (create.distances==TRUE) D.final=matrix(0,ncol(expr.data.norm),ncol(expr.data.norm))

  while(1){

    #cat(sprintf('\nRecursive clustering, beginning round %g ....',tot.recursive))
    action.taken=0
    mycl.new=rep(0,length(mycl))
    
    
    for (k in 1:max(mycl))
    {
      #print(sprintf('Checking cluster %g/%g, %g cells',k,max(mycl),length(which(mycl==k))))
      if (length(which(mycl==k))>dim.cutoff & sum(unclusterable[which(mycl==k)])==0 ) # then it must be re-clustered
      {
        #print('Computing Overdispersed genes ...')
        ODgenes=calculate.ODgenes(expr.data.norm[,which(mycl==k)],verbose = FALSE,min_ODscore = 2)
        
        if (is.null(ODgenes))
          {
          D=NA
          unclusterable[which(mycl==k)]=1   
          }
        else
          {
          dummy=as.matrix(ODgenes[[1]])
          ODgenes=which(dummy[,1]==1)
          #print('Computing distances ...')
          D=compute.distances(expr.norm = expr.data.norm[,which(mycl==k)],N_pct = model,edges = edges,driving.genes = ODgenes,lib.size = lib.size[which(mycl==k)],modality=modality)
          temp.clusters=bigscale.cluster(D,plot.clusters = FALSE,clustering.method = 'low.granularity',granularity.size=dim.cutoff,verbose=FALSE)$clusters #cut.depth=current.cutting,method.treshold = 0.2
          if (max(temp.clusters)>1) 
            action.taken=1
          else
            unclusterable[which(mycl==k)]=1  
          }
        
        
        if (create.distances==TRUE)
          {
          D.final[which(mycl==k),which(mycl==k)]=D # assigning distances
          max.inter=compute.max.inter(D,temp.clusters)
          D.final[which(mycl==k),which(mycl!=k)]=D.final[which(mycl==k),which(mycl!=k)]+max.inter
          D.final[which(mycl!=k),which(mycl==k)]=D.final[which(mycl!=k),which(mycl==k)]+max.inter
          }
        #print(sprintf('Partitioned cluster %g/%g in %g sub clusters ...',k,max(mycl),max(temp.clusters)))
        if (is.null(ODgenes)){mycl.new[which(mycl==k)]=1+max(mycl.new)} #ming, 2021-09-01, check out '2021-09-01-debug-bigscale2.R'
        else{ mycl.new[which(mycl==k)]=temp.clusters+max(mycl.new) } #ming, 2021-09-01, check out '2021-09-01-debug-bigscale2.R'
         
        
      }
      else
      {
        #print(sprintf('Cluster %g/%g is of size %g and cannot be further partitioned',k,max(mycl),length(which(mycl==k))))
        mycl.new[which(mycl==k)]=1+max(mycl.new)
      }
    }
    
    tot.recursive=tot.recursive+1  
    
    #cat(sprintf('\nRecursive clustering, after round %g obtained %g clusters',tot.recursive,max(mycl.new)))
    # if (modality=='jaccard' & tot.recursive==5)
    # {
    #   cat(sprintf('\n\nYou are analyzing ATAC-seq data, I stop the clustering to avoid creating too many clusters. If you want to customize the analysis with more/less rounds of recursive contact the developer at gio.iacono.work@gmail.com'))
    #   break
    # }
    #current.cutting=current.cutting+10
    if (action.taken==0) 
      break
    else
      mycl=mycl.new 
    
    
    
    #if(tot.recursive==5) break
  }

  # Preventing small clusters to be used, with a ugly trick

  bad.cells=c()
  bad.clusters=c()
  indexes=list()
  for (k in 1:max(mycl))
      {
      indexes[[k]]=which(mycl==k)
      if (length(which(mycl==k))<5)
        {
        bad.clusters=c(bad.clusters,k)
        bad.cells=c(bad.cells,which(mycl==k))
        }
      }
  if (length(bad.clusters)>0)
    {
    indexes=indexes[-bad.clusters]
    mycl=rep(0,length(mycl))
    for (k in 1:length(indexes))
      mycl[indexes[[k]]]=k
    warning(sprintf('Hidden %g clusters of %g cells total to be used: they are too small',length(bad.clusters),length(bad.cells)))
  }

  if (create.distances==FALSE)
    D.final=NA
  return(list(mycl=mycl,D=D.final))
}


compute.network.model = function (expr.data)
{
  sce = SingleCellExperiment(assays = list(counts = expr.data))  
  sce = preProcess(sce)
  sce = storeNormalized(sce)
  sce=setModel(sce)
  return(list(model=sce@int_metadata$model,edges=sce@int_metadata$edges))
}
  

generate.edges<-function(expr.data){
  
  edges=c(0,0.0000001)
  

  # if (pipeline=='atac')
  #   {
  #   step=0.75
  #   for (k in 1:40)
  #     edges[length(edges)+1]=edges[length(edges)]+step
  #   edges[length(edges)+1]=Inf
  #   return(edges)
  #   }
  if (ncol(expr.data)>6000)
  {
  print('Subsetting dataset...')
  ix=sample(c(1:ncol(expr.data)),5000)
  expr.data=expr.data[,ix]    
  }
  print('Creating edges...')
  percentage.exp=(sum(expr.data<10 & expr.data>0)/sum(expr.data>0))
  # Testing how many nonzeros values are below 70, if number is large than we have UMIs like distribution
  if (percentage.exp>0.9)
      { 
      modality='UMIs'
      num.edges=80
      period=10
      print(sprintf('%.1f %% of elements < 10 counts, therefore Using a UMIs compatible binning',percentage.exp*100))
      }#'UMIs'
    else
    { 
      modality='reads'
      num.edges=80
      period=5
      print(sprintf('%.1f %% of elements < 10 counts, therefore Using a reads compatible binning',percentage.exp*100))
    }
  
  if (modality=='UMIs')
      {
      step=0.75
      period.reset=period
      for (k in 3:num.edges)
        {
        edges[k]=edges[k-1]+step
        period=period-1
        if (period==0) 
          {
          period=period.reset
          step=step*2
          }
        }
      }
  else
      {
        count.periods=0
        step=5
        period.reset=period
        for (k in 3:num.edges)
        {
          edges[k]=edges[k-1]+step
          period=period-1
          if (period==0) 
          {
            period=period.reset
            step=step+10
          }
        }
      }
  edges[length(edges)+1]=Inf
  return(edges)
}

calculate.ODgenes = function(expr.norm,min_ODscore=2.33,verbose=FALSE,use.exp=c(0,1)) {

  if (ncol(expr.norm)<15000)
    downsample.od.genes = c(1:ncol(expr.norm))
  else
    downsample.od.genes = sample(ncol(expr.norm),15000)
  
  if (class(expr.norm)=='big.matrix')
    expr.norm=bigmemory::as.matrix(expr.norm[,downsample.od.genes])
  else
    expr.norm=as.matrix(expr.norm[,downsample.od.genes])
  
  #start.time <- Sys.time() 
  
    
  num.samples=ncol(expr.norm) 
  num.genes=nrow(expr.norm) 
  min.cells=max( 15,  round(0.002*length(expr.norm[1,]))) 
  skwed.cells=max( 5,  round(0.002*length(expr.norm[1,]))) 
  
  if (verbose)
    print(sprintf('Analyzing %g cells for ODgenes, min_ODscore=%.2f',ncol(expr.norm),min_ODscore))
  
  
  # Discarding skewed genes
  if (verbose)
    print('Discarding skewed genes')
  
  if (max(expr.norm>1))
      {
      expr.row.sorted=Rfast::rowSort(expr.norm) #MEMORY ALERT with Rfast::sort_mat
      a=Rfast::rowmeans(expr.row.sorted[,(num.samples-skwed.cells):num.samples])
      La=log2(a)
      B=Rfast::rowVars(expr.row.sorted[,(num.samples-skwed.cells):num.samples], suma = NULL, std = TRUE)/a
      rm(expr.row.sorted)
      gc()
      f=smooth.spline(x=La[La>0],y=B[La>0],df = 12) # smoothing sline approximatinf the La/B relationship
      yy=approx(f$x,f$y,La)$y # smoothing spline calculate on each exact value of La
    
      skewed=rep(0,length(yy))
      skewed_genes=which(B/yy>4)
      skewed[skewed_genes]=1
       df=as.data.frame(t(rbind(La,B,yy,skewed)))
      df$skewed=as.factor(df$skewed)
      g1=ggplot2::ggplot(df) + ggplot2::ggtitle('Discarding skewed genes')  + ggplot2::geom_point(ggplot2::aes(x=La, y=B,color=skewed),alpha=0.5)  + ggplot2::scale_color_manual(breaks = c("0", "1"),values=c("gray", "red")) + ggplot2::geom_line(ggplot2::aes(x=La, y=yy),colour="black")
      }
    else
      {
      skewed_genes=c()
      g1=plot(1,1)
      }
 
  
  
  
  
  # Fitting OverDispersed genes
  okay=which( Rfast::rowsums(expr.norm>0)>min.cells )
  if (verbose)
    print(sprintf('Using %g genes detected in at least >%g cells',length(okay),min.cells))

  okay=setdiff(okay,skewed_genes)
  if (verbose)
    print(sprintf('Further reducing to %g geni after discarding skewed genes', length(okay)))

  if (length(okay)<200)
    {
    #print('Returning no highly variable genes, too few genes, too much noise!')
    return(NULL)
    }
  
  # STEP1: local fit of expression and standard deviation
  expr.norm=expr.norm[okay,]
  #expr.norm.odgenes<<-expr.norm
  gc()
  a=Rfast::rowmeans(expr.norm)
  La=log2(a)
  B=Rfast::rowVars(expr.norm,std = TRUE)/a
  f=smooth.spline(x=La,y=B,df = 8) # smoothing sline approximatinf the La/B relationship
  yy=approx(f$x,f$y,La)$y # smoothing spline calculate on each exact value of La
  
  df=as.data.frame(t(rbind(La,B,yy)))
  g2=ggplot2::ggplot(df) + ggplot2::ggtitle('STEP1: local fit of expression and standard deviation')  + ggplot2::geom_point(ggplot2::aes(x=La, y=B),colour="brown",alpha=0.5)  + ggplot2::geom_line(ggplot2::aes(x=La, y=yy),colour="black") +
    ggplot2::xlab('Log2(expression)') + ggplot2::ylab('Standard Deviation')
  
  
  
  # STEP2: moving standard deviation of fit1
  B_corr1=B-yy
  sLa=sort(La, index.return = TRUE)
  sd_width=100
  movSD=zoo::rollapply(B_corr1[sLa$ix], width = sd_width,FUN=trim_sd, fill = NA)
  f=smooth.spline(x=sLa$x[!is.na(movSD)],y=movSD[!is.na(movSD)],df = 16)
  yy=approx(f$x,f$y,La)$y
  
  pos=max(which(!is.na(yy[sLa$ix])))
  yy[sLa$ix[pos:length(yy)]]=yy[sLa$ix[pos]]

  
  # movSD2=runSD(B_corr1[sLa$ix],sd_width)
  # f2=smooth.spline(x=sLa$x[!is.na(movSD)],y=movSD2[!is.na(movSD)],df = 16)
  # yy2=approx(f$x,f$y,La)$y
  
  
  ODscore=B_corr1/yy
  od_genes=which(ODscore>min_ODscore)

  
  od_genes=intersect(od_genes,which( La>=quantile(La,use.exp[1]) &  La<=quantile(La,use.exp[2]) ))
      
  ODgenes=rep(0,length(La))
  ODgenes[od_genes]=1
  
  df=as.data.frame(t(rbind(La,B_corr1,yy,ODgenes,ODscore)))
  df$ODgenes=as.factor(df$ODgenes)
  
  g3 = ggplot2::ggplot(df) + ggplot2::ggtitle('Determining overdispersed genes (OGgenes)')  + 
    ggplot2::geom_point(ggplot2::aes(x=La, y=B_corr1,color=ODgenes),alpha=0.5)  + ggplot2::scale_color_manual(breaks = c("0", "1"),values=c("gray", "red")) + ggplot2::geom_line(ggplot2::aes(x=La, y=yy),colour="black") +
  ggplot2::xlab('Log2(expression)') + ggplot2::ylab('Adjusted standard deviation ')
  
  # B_corr1.out=B_corr1
  # yy.out<<-yy
  # La.out<<-La
  
  g4=ggplot2::ggplot(df) + ggplot2::ggtitle('Determining overdispersed genes (OGgenes)')  + ggplot2::geom_point(ggplot2::aes(x=La, y=ODscore,color=ODgenes),alpha=0.5)  + ggplot2::scale_color_manual(breaks = c("0", "1"),values=c("gray", "red"))
  
  # generatnig a data.frame for the export
  dummy=matrix(0,num.genes,2)
  #print(dim(dummy))
  #print(length(okay))
  #print(dim(t(rbind(ODgenes,ODscore))))
  dummy[okay,]=t(rbind(ODgenes,ODscore))
  DFout=as.data.frame(dummy)
  colnames(DFout)=c('ODgenes','ODscore')
  
  
  if (verbose)
    print(sprintf('Determined  %g overdispersed genes',length(od_genes)))
  
  #end.time <- Sys.time()
  #print(end.time - start.time)
  
  gc()
  return(list(DFout,g1,g2,g3,g4))
  
  
}

trim_sd<-function(x){
  x=x[!is.na(x)]
  
  tmp=quantile(abs(x),c(0.05,0.95))
  out=sd(x[x>tmp[1] & x<tmp[2]])
}


compute.distances = function (expr.norm, N_pct , edges, driving.genes , genes.discarded,lib.size,modality='pca',pca.components=25){
  
  
  #print(sprintf('Proceeding to calculated cell-cell distances with %s modality',modality))
  
  if (modality=='jaccard')
      {
      print("Calculating Jaccard distances ...")
      D=as.matrix(jaccard_dist_text2vec_04(x = Matrix::Matrix( t(as.matrix(expr.norm[driving.genes,]>0)), sparse = T  )))
      return(D)
      }
 
  if (modality=='pca')
  {
    
    
    #print(sprintf('Using %g PCA components for %g genes and %g cells',pca.components,length(driving.genes),ncol(expr.norm)))
    if(max(expr.norm)==1)
      {
      print('Detecting ATAC-seq data...')
      dummy=svd(expr.norm[driving.genes,],0,pca.components)  
      }
      else
      dummy=svd(log10(expr.norm[driving.genes,]+1),0,pca.components)
    
    for (k in 1:ncol(dummy$v))
      dummy$v[,k]=dummy$v[,k]*dummy$d[k]
    #print('Computing distance from PCA data...')
    D=dist(dummy$v,method = 'euclidean')
    return(D) # is a distance object
  }  
  
  # normalize expression data for library size without scaling to the overall average depth
  if (class(expr.norm)=='big.matrix')
    expr.driving.norm=bigmemory::as.matrix(expr.norm[driving.genes,])/mean(lib.size)
  else
    expr.driving.norm=as.matrix(expr.norm[driving.genes,])/mean(lib.size)
  gc()
  
  if (modality=='correlation')
    {
    print("Calculating normalized-transformed matrix ...")
    expr.norm.transformed = transform.matrix(expr.driving.norm , 2 )
    gc()
    print("Calculating Pearson correlations ...")
    D=1-Rfast::cora(expr.norm.transformed)
    rm(expr.norm.transformed)
    gc()
    return(D)
  }
  
  
  
  
  if (!hasArg(genes.discarded)) genes.discarded =c()

  
  log.scores=get.log.scores(N_pct)
  
  # Consumes several Gb of memory this step!
  # Vector is a trick to increase speed in the next C++ part

  indA.size=1000000

  critical.value=max(lib.size)*max(expr.driving.norm)
  if (critical.value>indA.size/10)
    {
    indA.size=indA.size*50
      warning(sprintf('Critical value very high (%g): Increased memory usage!!',critical.value))
    if (critical.value>indA.size/10)
      stop(sprintf('Critical value way too high (%g): Stopping the analysis!!',critical.value))
    }
  
  #print('Proceding to allocate large vector')
  vector=c(0:indA.size)/10 # increase if code blocks, It can assign a gene exprssion level up to 10000000
  gc()
  ind.A=as.integer(cut(vector,edges,include.lowest = TRUE))
  rm(vector)
  gc()
  ind.A=ind.A-1
  gc()
  
  

  
  if (length(genes.discarded)>0)
  
    print('Detect that you want to remove genes')
  # vector.weights=ones(num_genes,1)
  # vector_weights(remove_genes)=0;
  # tic; D_ls=SC_dist_MEX_double_rv(total_data_norm,log_scores,ind_A,sum_ex,0,vector_weights); t=toc
  else
    {
    start.time=Sys.time() 
    D=C_compute_distances(expr.driving.norm,log.scores,ind.A,lib.size)
    #print("Time elapsed to calculate distances")
    #print(Sys.time()-start.time)
    }

  #D=Dist(Rfast::squareform(D), method = "euclidean", square = FALSE,vector=FALSE)
  #D=(1-fastCor(Rfast::squareform(D)))
  #diag(D)=0
  D=Rfast::squareform(D)
  
  gc()
  
  return(D)  

}

bigscale.cluster = function (D,plot.clusters=FALSE,cut.depth=NA,method.treshold=0.5,clustering.method='high.granularity',granularity.size,verbose=FALSE){

  gc()
    
    if (class(D)=='big.matrix')
      D=bigmemory::as.matrix(D)
    if (class(D)=='dgCMatrix')
      error('Error of the stupid biSCale2 programmer!')
    
  ht=hclust(as.dist(D),method='ward.D')
  ht$height=round(ht$height,6)


  if (is.na(cut.depth) & clustering.method=='high.granularity')
    {
    # Trying a variation of the elbow method to automatically set the cluster number
    result=rep(0,100)
    for (k in 1:100)
      {
      mycl <- cutree(ht, h=max(ht$height)*k/100)
      result[k]=max(mycl)
      }
    movAVG=zoo::rollapply(data = -diff(result),10,mean, fill = NA)

    cut.depth=max(which(movAVG>method.treshold)) #1

    if (is.infinite(cut.depth)) cut.depth=50
    
    }

  if (is.na(cut.depth) & clustering.method=='low.granularity')
  {

    
    result=c()
    progressive.depth=c(90, 80, 70, 60, 50, 40, 30, 20, 15, 10)
    
    for (k in 1:length(progressive.depth))
    {
      mycl <- cutree(ht, h=max(ht$height)*progressive.depth[k]/100)
      result[k]=max(mycl)
    }
    #print(result)
    result= ( diff(result)>0 )
    
    #print(result)
    
    detected.positions=c()
    for (k in 1:(length(result)-2))
      if (result[k]==TRUE & result[k+1]==TRUE)
        detected.positions=c(detected.positions,k)
        
    
    if (sum(result)==0)
      cut.depth=30
    else  
    {
      if (length(detected.positions)==0) 
        cut.depth=progressive.depth[min(which(result==TRUE))]
      else
        cut.depth=progressive.depth[min(detected.positions)]
    }

    
    
    #print(sprintf('Cut depth before checking granularity: %g percent',cut.depth))
    if (cut.depth>=30 & length(ht$order)<granularity.size) cut.depth=100 # do not cluster at all if it starts fragemnting so high alreay
    if (cut.depth>=40 & length(ht$order)>=granularity.size & length(ht$order)<(granularity.size*2)) cut.depth=100 # do not cluster at all if it starts fragemnting so high alreay
    if (cut.depth>=50 & length(ht$order)>=(granularity.size*2) & length(ht$order)<(granularity.size*3)) cut.depth=100 # do not cluster at all if it starts fragemnting so high alreay
    if (cut.depth>=60 & length(ht$order)>=(granularity.size*3) & length(ht$order)<(granularity.size*4)) cut.depth=100 # do not cluster at all if it starts fragemnting so high alreay
    if (cut.depth>=70 & length(ht$order)>=(granularity.size*4) & length(ht$order)<(granularity.size*5)) cut.depth=100 # do not cluster at all if it starts fragemnting so high alreay
    if (cut.depth>=80 & length(ht$order)>=(granularity.size*5) & length(ht$order)<(granularity.size*6)) cut.depth=100 # do not cluster at all if it starts fragemnting so high alreay
  }


  if (cut.depth<100)
    mycl = cutree(ht, h=max(ht$height)*cut.depth/100)
  else
    mycl = rep(1,length(ht$order)) 

  if (verbose)
  print(sprintf('Automatically cutting the tree at %g percent, with %g resulting clusters',cut.depth,max(mycl)))

  #if (length(unique(mycl))>1 & cut.depth==100) stop('Error in the clustering, wake up!')

  # Fixing the ordering of the clusters
  vicious.order=unique(mycl[ht$order])
  clusters=rep(0,length(mycl))
  for (k in 1:max(mycl))
    clusters[which(mycl==k)]=which(vicious.order==k)

  if (plot.clusters) # plotting the dendrogram
    {
    #print('We are here')
    d=as.dendrogram(ht)
    tot.clusters=max(clusters)
    palette=set.quantitative.palette(tot.clusters)
    d = dendextend::color_branches(d, k = tot.clusters,col = palette)  
    plot(d)
    }

  gc()
  return(list(clusters=clusters,ht=ht))

  }


fit.model = function(expr.norm,edges,lib.size,plot.pre.clusters=TRUE){
  # for clustering='recursive'
  # Performing down-sampling, model does not require more than 5000 cells
  if (ncol(expr.norm)>5000)
    {
    selected=sample(1:ncol(expr.norm),5000)
    if (class(expr.norm)=='big.matrix')
      expr.norm=bigmemory::as.matrix(expr.norm[,selected])
    else
      expr.norm=as.matrix(expr.norm[,selected])
    }
  else
    if (class(expr.norm)=='big.matrix')
      expr.norm=bigmemory::as.matrix(expr.norm)
    else
      expr.norm=as.matrix(expr.norm)
  
  print(dim(expr.norm))
  gc()

  # if (pipeline=='atac')
  #   cutoff=0
  # else
    cutoff=1
  #Remove genes with 1 or 0 UMIs/reads in the whole dataset.
  genes.low.exp =which(Rfast::rowsums(expr.norm>0)<=cutoff)
  if (length(genes.low.exp))
    {
    print(sprintf('I remove %g genes not expressed enough', length(genes.low.exp)))
    expr.norm=expr.norm[-genes.low.exp,]
    }
  
  num.genes=nrow(expr.norm)
  num.samples=ncol(expr.norm)
  
  
  # if (pipeline=='atac')
  # {
  #   print("Calculating distances for ATAC-seq pre-clustering...")
  #   D=as.matrix(jaccard_dist_text2vec_04(x = Matrix::Matrix(t(expr.norm>0))))
  #   gc()
  #   }  
  # else
  #   {
    print("Calculating normalized-transformed matrix ...")
    expr.norm.transformed = transform.matrix(expr.norm , 2 )
    gc()
    print("Calculating Pearson correlations ...")
    #D=fastCor(expr.norm) #improvable with BLAS, or cor(expr.norm,method = 'pearson')
    D=Rfast::cora(expr.norm.transformed)
    rm(expr.norm.transformed)
    gc()
    
    
  
  D=as.dist(1-D)
  print("Clustering  ...")
  ht=hclust(D,method='ward.D')
  
  
  if (plot.pre.clusters) # plotting the dendrogram
    plot(as.dendrogram(ht))
  
  
  # Adjusting max_group_size according to cell number
  if (num.samples<1250) max.group.size=0.10
  if (num.samples>=1250 & num.samples<=2000) max.group.size=0.08
  if (num.samples>=2000 & num.samples<=3000) max.group.size=0.06
  if (num.samples>=3000 & num.samples<=4000) max.group.size=0.05
  if (num.samples>=4000) max.group.size=0.04
  
  # Calculating optimal cutting depth for the pre-clustering
  print('Calculating optimal cut of dendrogram for pre-clustering')
  MAX_CUT_LEVEL=0.9#0.6
  cut.level=MAX_CUT_LEVEL
  while (TRUE)
    {
    while (TRUE)
      {
      if (cut.level<0) break
      T = cutree(ht, h=max(ht$height)*cut.level)
      if ( max(table(as.factor(T)))<(max.group.size*num.samples) & length(unique(T))<(0.05*num.samples) ) break
      cut.level=cut.level-0.01
      }
    if (cut.level>0 | abs(cut.level-0)<1e05) break
    else
      {
      max.group.size=max.group.size+0.01;
      cut.level=MAX_CUT_LEVEL;
      }
  }
  
   #g.dendro = fviz_dend(ht, h=max(ht$height)*cut.level)
   mycl <- cutree(ht, h=max(ht$height)*cut.level)
   
   if (plot.pre.clusters) # plotting the dendrogram
   {
     print('We are here')
     d=as.dendrogram(ht)
     tot.clusters=max(mycl)
     palette=set.quantitative.palette(tot.clusters)
     d = dendextend::color_branches(d, k = tot.clusters,col = palette)
     plot(d)
   }
   
   
   print(sprintf('Pre-clustering: cutting the tree at %.2f %%: %g pre-clusters of median(mean) size %g (%g)',cut.level*100,max(mycl),median(as.numeric(table(mycl))),mean(as.numeric(table(mycl))) ))
   # For loop over the clusters

   pb <- progress::progress_bar$new(format = "Analyzing cells [:bar] :current/:total (:percent) eta: :eta", total = length(mycl))
   
   N=matrix(0,length(edges)-1,length(edges)-1)
   for (i in 1:max(mycl))
      {
       pos=which(mycl==i)
      if (length(pos)>3) # just to aviod bugs
        {
        #print(sprintf("Pre-cluster %g/%g",i,max(mycl)))
        N=N+enumerate(expr.norm[,pos],edges,lib.size)
        }
       pb$tick(length(pos))
      }
     
   
   # Calculating percentiles
   N_pct=matrix(NA,length(edges)-1,length(edges)-1)
   for (k in 1:nrow(N))
     for (j in 1:nrow(N))
       N_pct[k,j] = sum(N[k,j:ncol(N)])/sum(N[k,])
  
  N_pct[is.na(N_pct)]=1 # fix 0/0
  print(sprintf("Computed Numerical Model. Enumerated a total of %g cases",sum(N)))
  gc()
  

  
  return(list(N_pct=N_pct,N=N))
  
}


transform.matrix<-function(expr.norm,case){
  ## for clustering='recursive'
  print('Computing transformed matrix ...')
  print('Converting sparse to full Matrix ...')
  
  expr.norm=as.matrix(expr.norm)
  
  
  if (case==2)# model=2. Log(x+1), 
    expr.norm=log2(expr.norm+1)
  
  if (case==4)# capped to 95% expression
    {
    print('Capping expression gene by gene ...') 
    for (k in 1:nrow(expr.norm))
      expr.norm[k,]=cap.expression(expr.norm[k,])
    }
  gc()
  print('Normalizing expression gene by gene ...')  
  # each row (gene) normalized between [0:1]  
  
  for (k in 1:nrow(expr.norm))
    expr.norm[k,]=shift.values(expr.norm[k,],0,1)
  #t(apply(expr.norm, 1,shift.values,A=0,B=1))
  
  
  gc()
  return(expr.norm)
  
}

shift.values<-function(x,A,B){
  ## for clustering='recursive'
  a=min(x)
  b=max(x)
  x_out=(x-a)/(b-a)*(B-A)+A
}

set.quantitative.palette = function(tot.el){
  ## for clustering='recursive'
  gc()
  if (tot.el<=11)
      {
      palette=c('#ffe119', '#4363d8', '#f58231', '#e6beff', '#800000', '#000075', '#a9a9a9', '#000000','#FF0000','#fffac8','#f032e6')
      palette=palette[1:tot.el]
      #https://sashat.me/2017/01/11/list-of-20-simple-distinct-colors/
      }
    else
    {
      #palette=c('#e6194B', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#42d4f4', '#f032e6', '#bfef45', '#fabebe', '#469990', '#e6beff', '#9A6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#a9a9a9', '#ffffff', '#000000')
      #if (tot.el>length(palette))
      #  stop('You have more than 22 clusters and I cannot find enough colors for them. Contact the developer at gio.iacono.work@gmail.com to fix this issue')
      #palette=palette[1:tot.el]
      #https://sashat.me/2017/01/11/list-of-20-simple-distinct-colors/
      color = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
      palette=sample(color, tot.el, replace = FALSE)
    }

    
  
    gc()                    
return(palette)
  
}

enumerate = function(expr.norm,edges,lib.size) {
  ### for clustering='recursive'
  gc()
  num.genes=nrow(expr.norm)
  num.samples=ncol(expr.norm) 
  min.cells = min ( round(num.samples/40) , 10 )

  
  # normalize expression data for library size without scaling to the overall average depth
  expr.norm=expr.norm/mean(lib.size)
  gc()
  
  # Removing lowly expressed genes
  genes.low.exp =which(Rfast::rowsums(expr.norm>0)<=min.cells)
  #print(sprintf('Enumeration of combinations, removing  %g genes not expressed enough', length(genes.low.exp)))
  expr.norm=expr.norm[-genes.low.exp,]
  
  # Initiating variable N
  #print(sprintf('Calculating %g couples over %g cells',(num.samples*num.samples - num.samples)/2,num.samples))
  N=matrix(0,length(edges)-1,length(edges)-1)
  
  # Enumerating  combinations
  for (i in 1:num.samples)
    for (j in i:num.samples)
      if (i!=j)
        {
        factor = mean ( lib.size[c(i,j)] )
        xgroup=cut(expr.norm[,i]*factor,edges,include.lowest = TRUE)
        ygroup=cut(expr.norm[,j]*factor,edges,include.lowest = TRUE)
        xy = table(xgroup, ygroup)
        N=N+xy
      }
  gc()
  return(N)
}


calculate.marker.scores = function (expr.norm, clusters, N_pct, edges, lib.size, speed.preset,cap.ones=F){
  ### for clustering='recursive'
  gc()  
  tot.clusters=max(clusters)
  
  # Initializing the 2D list (a matrix list)
  Mscores = matrix(vector('list',tot.clusters*tot.clusters),tot.clusters,tot.clusters)
  Fscores = matrix(vector('list',tot.clusters*tot.clusters),tot.clusters,tot.clusters)
    
 
  to.calculate=tot.clusters*(tot.clusters-1)/2
  computed=0
  # Rolling trought the pairwise cluster comparisons
  
  pb <- progress::progress_bar$new(format = "Running Diff.Expr. [:bar] :current/:total (:percent) eta: :eta", total = tot.clusters*(tot.clusters-1)/2)
  pb$tick(0)

  if (cap.ones==T)
    expr.norm[expr.norm>0]=1
  
  for (j in 1:(tot.clusters-1))
    for (k in (j+1):tot.clusters)
    {
      out=bigscale.DE(expr.norm = expr.norm, N_pct = N_pct, edges = edges, lib.size = lib.size, group1 = which(clusters==j),group2 = which(clusters==k),speed.preset = speed.preset)
      Mscores[[j,k]] = out[,1] # DEscores.final
      Fscores[[j,k]] = out[,2] # Log2(Fold change)
      Mscores[[k,j]] = -out[,1] # DEscores.final
      Fscores[[k,j]] = -out[,2] # Log2(Fold change)
      computed=computed+1
      pb$tick(1)
    }

  out=list(Mscores=Mscores,Fscores=Fscores)
  gc()
  return(out)
  
}

bigscale.DE = function (expr.norm, N_pct, edges, lib.size, group1,group2,speed.preset='slow',plot.graphic=FALSE){
  ## for clustering='recursive'
  #print('Starting DE')  
  gc.out=gc()


  num.genes.initial=nrow(expr.norm)

  if (class(expr.norm)=='big.matrix')
    expr.norm=bigmemory::as.matrix(expr.norm[,c(group1,group2)])
  else
    expr.norm=as.matrix(expr.norm[,c(group1,group2)])


  tot.lib.size=lib.size
  lib.size=lib.size[c(group1,group2)]
  group1=c(1:length(group1))
  group2=c((length(group1)+1):(length(group1)+length(group2)))

  #remove the genes with zero reads
  #genes.expr =which(Rfast::rowsums(expr.norm[,c(group1,group2)])>0)
  genes.expr =which(Rfast::rowsums(expr.norm)>0)
  #print(sprintf('I remove %g genes not expressed enough', (nrow(expr.norm)-length(genes.expr))))
  expr.norm=expr.norm[genes.expr,]
  gc()
  num.genes=nrow(expr.norm)
  num.samples=ncol(expr.norm)
  tot.el=nrow(N_pct)  


  # normalize expression data for library size without scaling to the overall average depth
  if (speed.preset!='fast')
    expr.norm=expr.norm/mean(tot.lib.size)



  # saving the FOLD changes
  F.change=log2(Rfast::rowmeans(expr.norm[,group2])/Rfast::rowmeans(expr.norm[,group1]))
  dummy=matrix(0,num.genes.initial,1)
  dummy[genes.expr]=F.change
  F.change=dummy
  rm(dummy)

  # if ( (length(group1)+length(group2))==1774 )
  # {
  #   # expr.norm.out<<-expr.norm
  #   # group2.out<<-group2
  #   # group1.out<<-group1
  #   # num.genes.initial.out<<-num.genes.initial
  #   F.change.out<<-F.change
  #   #print('Fold change of positin 40764')
  #   #print(F.change[40764])
  # }
  #   


  # Calculating Wilcoxon
  dummy=rep(0,num.genes)
  DE.scores.wc=rep(0,num.genes.initial)
  info.groups=c(rep(TRUE,length(group1)),rep(FALSE,length(group2)))
  input.matrix=expr.norm[,c(group1,group2)]


  #print(sprintf('Launching the Wilcoxon, %g VS %g cells',length(group1),length(group2)))
  dummy=BioQC::wmwTest(t(input.matrix), info.groups , valType = 'p.two.sided' )
  # fixing the zeroes in wilcoxon
  zeroes=which(dummy==0)
  if (length(zeroes)>0)
    dummy[zeroes]=min(dummy[-zeroes])
  #print('Computed Wilcoxon')
  dummy[dummy>0.5]=0.5
  DE.scores.wc[genes.expr]=-qnorm(dummy)
  DE.scores.wc=abs(DE.scores.wc)*sign(F.change)
  rm(dummy,input.matrix)

  if (speed.preset=='fast')
    return(cbind(DE.scores.wc,F.change))

  problem.size=mean(length(group1),length(group2))

  if (speed.preset=='normal')
    {
    max.size=2000
    if (problem.size<=200) 
      max.rep =  5
    else if (problem.size>200 & problem.size<=500)
      max.rep =  3
    else if (problem.size>500 & problem.size<=1000)
      max.rep =  2
    else if (problem.size>1000)
     max.rep =  1
    }
    
    
  if (speed.preset=='slow')
    {
    max.size=5000
    if (problem.size<=1000) 
      max.rep =  4
    else if (problem.size>1000 & problem.size<=2500)
      max.rep =  3
    else if (problem.size>2500 & problem.size<=4000)
      max.rep =  2
    else if (problem.size>4000)
      max.rep =  1
    } 


  # Downsampling according to MAX-SIZE
  if (length(group1)>max.size)
  {
    group1=sample(group1)
    group1=group1[1:max.size]
  }
  if (length(group2)>max.size)
  {
    group2=sample(group2)
    group2=group2[1:max.size]
  } 


  # Making N_pct symmetric
  for (k in 1:tot.el)
    for (h in 1:tot.el)
      if (k>h) N_pct[k,h]=N_pct[h,k]


  # 3) Calculating log_scores and removing infinite values
  log.scores=-log10(N_pct)
  for (k in 1:tot.el)
  {
    max.val=max(log.scores[k,is.finite(log.scores[k,])])
    infinite.el=which(is.infinite(log.scores[k,]))
    log.scores[k,infinite.el]=max.val
    log.scores[infinite.el,k]=max.val
  }

  # zeroing the diagonal and setting pval<dereg zone
  diag(log.scores)=0
  for (k in 1:tot.el)
    for (h in 1:tot.el)
      if (k>h) log.scores[k,h]=-log.scores[h,k]

  # Vector is a trick to increase speed in the next C++ part
  indA.size=1000000
  critical.value=max(lib.size)*max(expr.norm)

  if (critical.value>indA.size/10)
  {
    indA.size=indA.size*50
    warning(sprintf('Critical value very high (%g): Increased memory usage!!',critical.value))
    if (critical.value>indA.size/10)
      stop(sprintf('Critical value way too high (%g): Stopping the analysis!!',critical.value))
  }

  #print('Proceding to allocate large vector')
  vector=c(0:indA.size)/10 # increase if code blocks, It can assign a gene exprssion level up to 10000000
  gc()
  ind.A=as.integer(cut(vector,edges,include.lowest = TRUE))
  rm(vector)
  gc()
  ind.A=ind.A-1
  gc()



  # Calculating scores of real DE
  #print('Rcpp computing real DE')
  results.DE.real=DE_Rcpp(expr.norm[,group1],expr.norm[,group2], log.scores,  ind.A, lib.size[group1], lib.size[group2])
  gc()
  DE.scores.real=results.DE.real[1:num.genes]
  DE.counts.real=results.DE.real[(num.genes+1):length(results.DE.real)]



  # Calculating scores of random permutations
  #print('Rcpp computing randomly reshuffled DEs')
  idx=c(group1,group2)
  Total.scores.fake=c()
  Total.counts.fake=c()
  for (repetitions in 1:max.rep)
    {
    #print(sprintf('Random repetition %d',repetitions))
    idx=sample(idx)
    fake.A=idx[1:length(group1)]
    fake.B=idx[(length(group1)+1):length(idx)]
    results.random.DE=DE_Rcpp(expr.norm[,fake.A],expr.norm[,fake.B], log.scores,  ind.A, lib.size[fake.A], lib.size[fake.B])
    gc()
    Total.scores.fake=c(Total.scores.fake,results.random.DE[1:num.genes])
    Total.counts.fake=c(Total.counts.fake, results.random.DE[(num.genes+1):length(results.random.DE)])
    gc()
    }

    
  # Fitting the random DEs: Moving standard deviation DE_scores against DE_counts
  #print('Fitting the random DEs')
  sa=sort(Total.counts.fake, index.return = TRUE)
  sd_width=200
  movSD=zoo::rollapply(Total.scores.fake[sa$ix], width = sd_width,FUN=sd, fill = NA)
  last.value.y=movSD[(length(movSD)-sd_width/2)] # From now a quick fix to the two ends of movSD
  before.last.value.y=movSD[(length(movSD)-sd_width)]
  last.value.x=sa$x[(length(movSD)-sd_width/2)]
  before.last.value.x=sa$x[(length(movSD)-sd_width)]
  estimated.slope=(last.value.y-before.last.value.y)/(last.value.x-before.last.value.x)
  estimated.slope=max(0,estimated.slope)
  estimated.increase=estimated.slope*(sa$x[length(sa$x)]-last.value.x)
  movSD[length(movSD)]=last.value.y+estimated.increase
  movSD[1:sd_width/2]=min(movSD[(sd_width/2+1):sd_width])
  f=smooth.spline(x=sa$x[!is.na(movSD)],y=movSD[!is.na(movSD)],df = 32)


  yy=approx(f$x,f$y,DE.counts.real)$y

    # f.out<<-f
    #yy.out<<-yy
    #DE.counts.real.out<<-DE.counts.real
    # x.out<<-sa$x[!is.na(movSD)]
    # y.out<<-movSD[!is.na(movSD)]


  if (any(is.na(yy))) # fixing the NA in the yy, when they coincide with the highest DE.counts.real
    {
      null.fits=which(is.na(yy))
      for (scroll.null.fits in 1:length(null.fits))
      {
        result=unique(DE.counts.real[-null.fits]<DE.counts.real[null.fits[scroll.null.fits]])
        if (length(result)>1) stop('Fix this bug, programmer!')
        if (result=='FALSE') 
          yy[null.fits[scroll.null.fits]]=min(yy[-null.fits])
        else
          yy[null.fits[scroll.null.fits]]=max(yy[-null.fits])
      }
    }
   

  # I want at least 5 non zeros cell in one group if the other is full empty
  treshold=min( 5*max(length(group1),length(group2)),(length(group1)*length(group2))/2) # values fitted below the 400 comparisons (eg.20cell*20cells) are considered noisy and replace with the lower value of next interval
  if (sum(DE.counts.real>treshold)==0)
    error('Clusters too small for DE')
  yy[DE.counts.real<=treshold]=min(yy[DE.counts.real>treshold])


  # Setting to zero the scores of all gene with less than 5 cells.
  #DE.scores.real.out<<-DE.scores.real
  #DE.scores.wc.out<<-DE.scores.wc
  #yy.out<<-yy

  DE.scores=DE.scores.real/yy
  ## debugging


  # saving the DE scores
  DE.scores=rep(0,num.genes.initial)
  DE.scores[genes.expr]=DE.scores.real/yy
  DE.scores=abs(DE.scores)*sign(F.change)

  #DE.scores.out<<-DE.scores


  if (plot.graphic)
    {
    yy2=approx(f$x,f$y,Total.counts.fake)$y
    df=as.data.frame(t(rbind(DE.counts.real,DE.scores.real,yy,Total.scores.fake,Total.counts.fake,yy2)))
    g1=ggplot2::ggplot(df) + ggplot2::ggtitle('Fitting over random resampling')  + ggplot2::geom_point(ggplot2::aes(x=Total.counts.fake, y=Total.scores.fake,color='gray'),alpha=0.5)  + ggplot2::geom_line(ggplot2::aes(x=Total.counts.fake, y=yy2),colour="black") +
      ggplot2::xlab('DE_counts') + ggplot2::ylab('DE_score')
    print(g1)
    g2=ggplot2::ggplot(df) + ggplot2::ggtitle('DE results')  + ggplot2::geom_point(ggplot2::aes(x=DE.counts.real, y=DE.scores.real,color='gray'),alpha=0.5)  + ggplot2::geom_line(ggplot2::aes(x=DE.counts.real, y=yy),colour="black") + ggplot2::xlab('DE_counts') + ggplot2::ylab('DE_score')
    print(g2)
    }


       

  # Merging Wilcoxon and bigscale pvalues to increase DE quality

  factor1=max(DE.scores.wc[!is.infinite(DE.scores.wc)])/max(DE.scores[!is.na(DE.scores)])
  factor2=min(DE.scores.wc[!is.infinite(DE.scores.wc)])/min(DE.scores[!is.na(DE.scores)])
  factor=mean(c(factor1,factor2))
  # print(max(DE.scores.wc[!is.infinite(DE.scores.wc)]))
  # print(max(DE.scores[!is.na(DE.scores)]))
  # print(factor1)
  # print(min(DE.scores.wc[!is.infinite(DE.scores.wc)]))
  # print(min(DE.scores[!is.na(DE.scores)]))
  # print(factor2)
  if (is.infinite(factor) | factor<0)
  {
    stop('Problems with the factor')
  }
  #print(sprintf('Factor1 %.2f, factor2 %.2f, average factor %.2f',factor1,factor2,factor))

  DE.scores.wc=DE.scores.wc/factor
  DE.scores.final=sqrt(DE.scores.wc^2+DE.scores^2)*sign(F.change)

  #print(Sys.time()-start)

  gc()

  return(cbind(DE.scores.final,F.change))  
}

