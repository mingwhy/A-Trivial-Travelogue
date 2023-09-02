#############################################################
#         Metabolica: A systems biology approach			#
#           to prioritize METABOLItes in CAncer.			#         
#															#        
#         cankutcubuk [at] {gmail} [dot] {com}				#       
#                     2018-2020								#      
#              http://clinbioinfosspa.es					#       
#                @ FPS, Sevilla , Spain						#    
############################################################# 


#' Metabolica: Integrates gene expression into metabolic networks to estimate metabolite abundances. 
#'
#' @importFrom igraph get.edgelist incident
#' @param nodes.vals 
#' You can use do.calc.rxnvals to calculate node values. nodes.vals <- do.calc.rxnvals(all_producers,rxn2gene,combat.vals)
#' @param subgraph precomputed subpathways of KEGG metabolic networks. Each subpathway produces a metabolite. 
#' Precomputed for hsa :/files/all_producers.RData
#' @param ininodes initial nodes, the metabolites of subgraph where the flux propagation starts.
#' @param endnode end node, the metabolite which is produced by the subgraph. Flux propagation ends at this node.
#' @param maxnum the enough number of iterations to visit all the nodes of subpathway and convergence the loops (includes reversible reactions). Default=5000.
#' To reduce the computation time and find the maximum number of the required iterations you can use grid.search option.
#' Analysis of 1500 samples using 20 threads takes 2 days when we use default parameters 
#' @param tol a threshold to determine convergenced loops. Default=1.0e-10
#' Flux propagation stops when the difference of flux value before and after visiting the nodes in the loop does not change more than tol value.
#' @param response.tol is an obsolete argument which will removed in the next release.
#' @param algR integrates the values of arriving fluxes and the reaction node value. In other words estimates the amount of substrate+enzyme complex. Default="minprod".
#' @param algM integrates amount of the products (metabolites). Default="nsqrt".
#' @param grid.search used to search optimum maxnum value
#' 
#' @examples
#' combat.vals is gene expression matrix (batch corrected, normalized and rescaled between 0-1). 
#' load("/files/combat.vals.LUAD.RData")
#' load("/files/all_producers.RData")
#' load("/files/rxn2gene.RData")
#' my_producer <- all_producers[[i]] # i can 1 to length(all_producers)
#' nodeVals <- do.calc.rxnvals(all_producers,rxn2gene,combat.vals)
#' res <- metabolica(nodes.vals=nodeVals, subgraph = my_producer$produceRigraph, ininodes = my_producer$initialNodes, endnode = my_producer$isProduced, algR = "minprod", algM = "nsqrt", maxnum = my_producer$maxnum*1.5)
#' 
#' @return Returns list with node.signal and sig.dif
#' 
#' @export
#'


metabolica <- path.value <- function(
		nodes.vals, 
		subgraph, 
		ininodes, 
		endnode, 
		maxnum = 5000, 
		tol = 1.0e-10, 
		response.tol = 0,
        algR=c("minprod","maxprod","maxmin","minmin"), 
        algM=c("nsqrt", "sum", "mean"), 
        grid.search=F){
			
  # Remove response.tol
  # tol = 1.0e-20 starts to return NaNs
  suppressMessages(library(igraph))
  algR <- match.arg(algR)
  algM <- match.arg(algM)
  
  # Initialize lists
  ready <- ininodes
  processed <- list()
  notAllowed <- list()
  actnodeVex <- c()
  
  # Initialize node values
  node.signal <- nodes.vals
  # node.signal[,] <- NA # this line has been inactivated on 09 May 2018. Otherwise it causes an error with cpd:C00022-hsa00010_SIF cpd:C00103 (actnode) <-> R00959 (NA)
  endnode.signal.dif <- 10
  num <- 0
  reached_last <- F
  
  while(length(ready) > 0 && num <= maxnum && endnode.signal.dif[length(endnode.signal.dif)]>tol){
    num <- num + 1
    actnode <- ready[[1]]
    old.signal <- node.signal[actnode,]
    actnodeVex <- c(actnodeVex,actnode)
    
    node.signal[actnode,] <- compute.node.flux(actnode, nodes.vals, node.signal, subgraph, response.tol, algR, algM)
    # Updates compound node vals
    if(length(grep("rn:",actnode,invert = T))>0){  nodes.vals[actnode,] <- node.signal[actnode,]}
    
    # Transmit flux
    nxt <- do.get.nextnodes(subgraph,actnode,notAllowed)
    nextnodes <- nxt$nextnodes
    notAllowed <- nxt$notAllowed
    
    dif <- old.signal - node.signal[actnode,]
    
    if(actnode==endnode){
      reached_last <- T
      if(!all(is.na(dif))){
        endnode.signal.dif <- c(endnode.signal.dif, sqrt(sum(dif^2,na.rm = T)))
      }
      ready <- ready[ready!=endnode]
    }
    if(all(is.na(old.signal)) || endnode.signal.dif[length(endnode.signal.dif)] > tol){
      
      if(actnode!=endnode){ # If the flux reaches endnode, it does not initialize the propagation from there again. while node.signal[,] <- NA and since we did "trim outgoing edges of endnode": we do not need this if statement.
        ready <- unique(c(ready, nextnodes))
      }
      ready <- ready[ready!=actnode] # The flux is under propagation through several branches of the network (e.g. multifurcations). The fluxes which reach to the endnode must stop. But the rest needs to continue.
    }
  }
  if(reached_last==F){
    endnode.signal.dif <- NA
  }
  if(grid.search){
    # Get minimum number for argument maxnum.  %95 (arbitary) percent of nodes should be visited
    if(length(unique(actnodeVex))==length(V(subgraph)$name)){complete <- T}else{complete <- F}
    if(endnode.signal.dif[length(endnode.signal.dif)]==tol){convergence <- T}else{convergence <- F}
    not_visited <- V(subgraph)$name[which(!V(subgraph)$name %in% unique(actnodeVex))]
    return(list("maxnum"=num,"notVisited"=not_visited,"iscompleted"=complete,"isconvergence"=convergence))
  }else{
    # return(list("node.signal"=node.signal[endnode,], "sig.dif"=endnode.signal.dif,"node.vals"=nodes.vals))
    return(list("node.signal"=node.signal[endnode,], "sig.dif"=endnode.signal.dif))
  }
}


#' @param actnode
#' @param nodes.vals 
#' @param node.signal
#' @param subgraph
#' @param response.tol = 0
#' @param algR
#' @param algM

compute.node.flux <- function(actnode, nodes.vals, node.signal, subgraph, response.tol = 0, algR, algM){
  
  node.val <- nodes.vals[actnode,]
  incis <- incident(subgraph, actnode, mode="in")
  
  if(length(incis)==0){
    signal <- rep(1, length(node.val))
  }else{
    # Get flux - From the direction of substrate to reaction
    prevs <- get.edgelist(subgraph)[incis,1]
    input_signals <- node.signal[prevs,,drop=F]
    nas <- is.na(input_signals[,1])
    prevs <- prevs[!nas]
    incis <- incis[!nas]
    input_signals <- input_signals[!nas,,drop=F]
    #typeincis <- E(subgraph)$relation[incis]
    
    # Flux percentages for multifurcations -  Flux is divided at the percentage of their outgoing node concentration - From the direction of substrate to reaction
    if(length(grep("rn:",actnode))>0){
      percents <- fluxPercentages(subgraph,prevs,nodes.vals,actnode,mymode = "out")
      for(p in rownames(input_signals)){
        input_signals[p,] <-  input_signals[p,] * percents[[p]]
      }
    }
    signal <- do.compute.flux(actnode, input_signals, node.val, algR, algM)
    # If the signal is too low, it will not be propagated
    # if(sum(nas) == 0 && signal < response.tol){
    #   signal <- rep(0,length(node.val))
    # }
  }
  return(signal)
}


# here node.signal or nodes.vals
fluxPercentages <- function(subgraph,prevs,nodes.vals,actnode,mymode="out"){
  prevsoutsCoef <- list()
  for(i in 1:length(prevs)){
    p <- prevs[i]
    pout <- get.edgelist(subgraph)[incident(subgraph, p, mode = mymode),2]
    if(length(pout)>1){
      poutpercent <- apply(nodes.vals[pout,],2, function(x) x/sum(x))
      # !!! if all outgoing reactions are close (if they are all zero) 0/0=NA.
      # So we replace NAs with zero and switch off this branch of the circuit. 
      # subgraph=cpd:C00468-hsa00140_SIF, actnode=rn:R04854
      poutpercent[is.na(poutpercent)] <- 0
      prevsoutsCoef[[p]] <- poutpercent[actnode,]
    }else{
      poutpercent <- rep(1,ncol(nodes.vals))
      prevsoutsCoef[[p]] <- poutpercent
    }
  }
  return(prevsoutsCoef)
}


prettyifelse <- function(test,una,olaotra){
  if(test){
    return(una)
  } else {
    return(olaotra)
  }
}


do.compute.flux <- function(actnode, input_signals, node.val, algR, algM){
  
  if(length(grep("rn:",actnode))>0){
    # There is no stoichiometry. For that reason, our generalized rule is 1 unit enzyme needs 1 unit metabolite.
    # if(enzyme > substrate){ product <- substrate }else{ product <- enzyme }
    # If a node has more than one indegree than we use the maximum or minimum one.
    
    if(algR=="maxmin"){
      signal <- prettyifelse(nrow(input_signals)>0, apply(input_signals, 2, max), node.val)
      signal <- cbind(node.val,signal)
      signal <- apply(signal,1,min)
    }
    if(algR=="minmin"){
      signal <- prettyifelse(nrow(input_signals)>0, apply(input_signals, 2, min), node.val)
      signal <- cbind(node.val,signal)
      signal <- apply(signal,1,min)
    }
    if(algR=="maxprod"){
      signal <- prettyifelse(nrow(input_signals)>0, apply(input_signals, 2, max), node.val)
      signal <- cbind(node.val,signal)
      signal <- apply(signal,1,prod)
    }
    if(algR=="minprod"){
      signal <- prettyifelse(nrow(input_signals)>0, apply(input_signals, 2, min), node.val)
      signal <- cbind(node.val,signal)
      signal <- apply(signal,1,prod)
    }
    
  }else{
    # sqrt is to suppress the effect of drastic compound concentration changes 
    # This needs node.val between [0-1] interval. And it states a dynamic system.
    if(algM=="nsqrt"){
      n <- nrow(input_signals) + 1
      # If the reaction is closed, the amount of the product which contributes to metabolite node will be zero...
      input_signals[input_signals==0] <- min(input_signals[input_signals!=0]) * 0.1
      signal <- apply(input_signals, 2, prod)
      signal <- do.nth.root(x = (signal * node.val), n)
    }else if(algM=="sum"){
      signal <- apply(input_signals, 2, sum)
      signal <- (signal + node.val)
      if(any(signal>1)){warning("saturation occurs"); signal[signal>1] <- 1}
    }else if(algM=="mean"){
      n <- nrow(input_signals) + 1
      signal <- apply(input_signals, 2, sum)
      signal <- (signal + node.val)/n
    }else{print("algR: nsqrt, sum or mean")}
  }
  return(signal)
}


do.nth.root <- function(x, n){
  return(exp(log(x)/n))
}

#' do.calc.rxnvals: Calculates the reaction node values using GPR like rules and gene expression.
#'
#' @param all_producers precomputed subpathways of KEGG metabolic networks. Each subpathway produces a metabolite. 
#' Precomputed for hsa :/files/all_producers.RData
#' @param rxn2gene data.frame to map between reaction, enzyme, entrez gene id and hgnc gene symbols
#' for hsa :/files/rxn2gene.RData
#' @param combat.vals gene expression matrix (batch corrected, normalized and rescaled between 0-1). 
#' Columns are the samples and rows are the genes. Gene IDs (rownames) need to be in Entrez/NCBI ID
#' @return Reaction node values

do.calc.rxnvals <- function(all_producers,rxn2gene,combat.vals){
  
  RxnNodesList  <- c()
  CpdNodesList <- c()
  
  for(l in names(all_producers)){
    my_producer <- all_producers[[l]]
    nodes <- unique(c(my_producer$procuderSubpath$V1, my_producer$procuderSubpath$V3))
    rxnodes <- nodes[grep("rn:",nodes)]
    cpdnodes <- nodes[grep("rn:",nodes,invert = T)]
    RxnNodesList <- c(RxnNodesList,rxnodes)
    CpdNodesList <- c(CpdNodesList,cpdnodes)
  }
  
  CpdNodesList <- unique(CpdNodesList)
  RxnNodesList <- unique(RxnNodesList)
  
  rxn2gene <- rxn2gene[rxn2gene$Entrez!="",]
  combat.vals <- combat.vals[which(rownames(combat.vals) %in% unique(rxn2gene$Entrez)),]
  combat.vals[combat.vals<0]<-0
  combat.vals <- t(apply(combat.vals, 1, function(x) (x-min(x,na.rm=T))/(max(x,na.rm=T)-min(x,na.rm=T))))
  
  nodeVals <- mat.or.vec(nr = length(RxnNodesList)+length(CpdNodesList),nc = ncol(combat.vals))
  rownames(nodeVals) <- c(RxnNodesList,CpdNodesList)
  colnames(nodeVals) <- colnames(combat.vals)
  
  count <- 1
  for(rn in rownames(nodeVals)){
    rxns <- strsplit(rn, split = " ")[[1]]
    isTRUE <- length(rxns)>1 & length(grep("rn",rxns))>1
    gexCL <- list()
    for(r in rxns){
      gexC <- c()
      genes <- unique(rxn2gene$Entrez[rxn2gene$Reaction==r])
      genes <- genes[which(genes %in% rownames(combat.vals))]
      gex <- combat.vals[genes,, drop=F]
      if(nrow(gex)>1){
        # Before, it was an arbitary value (quantile 0.9).
        gexC <- apply(gex,2, function(x){quantile(x,0.5, na.rm = TRUE)})
      }else if(nrow(gex)==1){
        gexC <- gex
      }else if(nrow(gex)==0){
        gexC <- rep(NA, ncol(nodeVals))
      }
      if(isTRUE){
        gexCL[[r]] <- gexC
      }
    }
    if(isTRUE){ 
      gexCL <- do.call(gexCL, what = "rbind")
      gexC <- apply(gexCL,2, function(x) max(x,na.rm = T))
    }else{
      gexC <- gexC
    }
    nodeVals[rn,] <- gexC
    count <- count+1
    if((count %% 50) == 0){ print(count)}
  }
  nodeVals[is.na(nodeVals)] <- 1
  return(nodeVals)
}


do.simplify.results <- function(res, algR, algM){
  # res is the result object of metabolica
  
  if(names(res)[length(res)]=="parameters"){
    print(res$parameters)
    res <- res[-length(res)]
  }
  res <- t(sapply(res, function(x) x$node.signal))
  # isNovarMets <- which(apply(res,1, function(x) sd(x)!=0))
  # res <- res[isNovarMets,]
  mets <- sapply(rownames(res), function(x) strsplit(x,split = "-")[[1]][1])
  mtab <- table(mets)  
  dups <- names(mtab)[mtab > 1]
  
  sumMet <- sapply(dups, function(x){ 
    input_signals  <- res[grep(x, rownames(res)),]
    if(algM!="sum"){ node.val <- rep(1, ncol(input_signals)) }else{ node.val <- rep(0, ncol(input_signals))}
    sM  <- do.compute.flux(x, input_signals, node.val, algR=algR, algM)
    sM
  })
  
  sumMet <- t(sumMet)
  res <- rbind(res[-grep(paste0(dups,collapse = "|"),rownames(res)),],sumMet)
  
  # print info
  isAlsoKnown<- rownames(res)[grep(" ",rownames(res))]
  if(length(isAlsoKnown)>0){
    message("is also known:")
    message(paste0(isAlsoKnown,collapse = "\n"))
  }
  cpds <- rownames(res)
  drugs <- grep("dr:", cpds)
  glycans <- grep("gl:", cpds)
  compounds <- grep("cpd:",cpds[-c(drugs, glycans)])
  message(paste0(length(drugs)," Drug(s), ", length(glycans), " Glycan(s), ", length(compounds), " Compound(s)"))
  message(paste0("There is no other type of metabolite: ",length(cpds)==length(c(drugs, glycans, compounds))))
  rownames(res) <- sapply(rownames(res), function(x) strsplit(x,split = "-")[[1]][1])
  return(res)
}


do.fold.change <- function(simplify.results,designmat,class1, class2){
  # value of do.simplify.results
  # and design matrix col1=name, col2=class
  foldvec <- c()
  for(i in 1:nrow(simplify.results)){
    x <- simplify.results[i,]
    mT <- mean(x[designmat$class==class1])
    mN <- mean(x[designmat$class==class2])
    if(is.na(mT) || is.na(mN)){
      fold <- NA
    }else{
      if(mT==0 & mN>0){
        fold <- "DOWN"
      }else if(mT>0 & mN==0){
        fold <- "UP"
      }else{
        f <- mT/mN
        fold <- ifelse(f>1, "UP", "DOWN")
      }
    }
    foldvec <- c(foldvec,fold)
  }
  names(foldvec) <- rownames(simplify.results)
  return(foldvec)
}


get.metabolitename <- function(metabolicaIDs,compound_info){
  #names(all_producers)
  #load("/Metabolica/Annotation/LiteratureSearch/KEGG_compound_info.RData")
  metaboliteNames <- sapply(metabolicaIDs, function(x) compound_info[[strsplit(x, "-")[[1]][1]]][1])
  metaboliteNames[sapply(metaboliteNames, is.null)] <- "No name"
  return(unlist(metaboliteNames))
}


get.pathname <- function(metabolicaIDs, pathwayList){
  #names(all_producers)
  #pathwayList <- read.delim("/Metabolica/Annotation/keggpathway.txt",sep = "\t",header = F,stringsAsFactors = F)
  pathwayList$V1 <- gsub("path:map","hsa",pathwayList$V1)
  pathIDs <- sapply(metabolicaIDs, function(x) strsplit(x, "-")[[1]][2])
  pathNames <- pathwayList[match(pathIDs, pathwayList$V1),2]
  return(pathNames)
}


do.get.nextnodes <- function(subgraph,actnode,notAllowed){
  # A <-> Reaction <-> B shows reversible reaction.
  # A and B are metabolites.
  # If the flux goes from A to B it can also go from B to A but not Reaction to A...
  # Here we control the over-propagation from reaction node to metabolite node ...
  # Alternetively, we can build a network with weighted edges where the edges are representing reactions and nodes are only metabolites....
  # In this case, we will not need this function
  nextnodes <- get.edgelist(subgraph)[incident(subgraph, actnode, mode="out"),2]
  
  if(length(grep("rn:",nextnodes))>0){
    notAllowed[nextnodes] <- actnode
  }
  if(length(grep("rn:",actnode))>0){
    if(!is.null(notAllowed[[actnode]])){
      nextnodes <- nextnodes[nextnodes!=notAllowed[[actnode]]]
      notAllowed[[actnode]] <- NULL
    }
  }
  return(list("nextnodes"=nextnodes, "notAllowed"=notAllowed))
}


#' @importFrom foreach foreach
#' @importFrom doMC registerDoMC

do.calcMaxNum <- function(all_producers, nodeVals, samplesizes, rep, outDir,verbose){
  # samplesizes =  a vector of different sample sizes. samplesizes <- c(20,40,100,200,400)
  # rep is number of repeats,. rep <- 10
  
  if(any(ncol(nodeVals)<samplesizes)){ print("sample sizes can not be bigger than total sample size in nodeVals matrix"); stop()}
  
  # different sample sizes
  for(s in samplesizes){ # sample size had no affect on maxnum value. Tested on sample sizes of 20,40,100,200,400
    pstart <- Sys.time()
    myMaxnum <- list()
    for(j in 1:rep){
      maxnumList <- foreach(i=1:length(all_producers), .combine = c) %dopar% {
        if(verbose){print(names(all_producers)[i])}
        my_producer <- all_producers[[i]] 
        mysamples <- sample(x = 1:ncol(nodeVals),s)
        start_time <- Sys.time()
        maxnum <- path.value(nodes.vals = nodeVals[,mysamples], subgraph = my_producer$produceRigraph ,ininodes = my_producer$initialNodes, endnode = my_producer$isProduced,grid.search = T)
        end_time <- Sys.time()
        maxnum$time <- end_time - start_time
        return(list("maxnum"=maxnum))
      }
      names(maxnumList) <- names(all_producers)
      myMaxnum[[paste0(s,"_samples")]][[paste0("rep",j)]] <- maxnumList
      if(verbose){print(paste0("############ rep ", j," finished ############"))}
    }
    estart <- Sys.time()
    if(verbose){
      print(paste0("################## saved sample size ", s," ################## "))
      print(round(estart - pstart))
    }
    save(myMaxnum,file = paste0(outDir,s,"_maxnum.RData"))
  }
}


do.updateMaxNumMat <- function(maxnumfiles, all_producers){
  maxmat <- matrix(nrow = length(all_producers) ,ncol = length(maxnumfiles))
  minmat <- maxmat
  dimnames(minmat) <- dimnames(maxmat) <- list(names(all_producers), paste0(unlist(regmatches(basename(maxnumfiles), gregexpr("[[:digit:]]+", basename(maxnumfiles)))),"_samples"))
  
  load(maxnumfiles)
  m1 <- unlist(regmatches(basename(maxnumfiles), gregexpr("[[:digit:]]+", basename(maxnumfiles))))
  m1 <- paste0(m1,"_samples")
  for(p in names(all_producers)){
    repmax <- sapply(myMaxnum[[m1]], function(x) x[[p]]$maxnum)
    repmax <- as.numeric(as.character(repmax))
    # there are some NULLs returned from calculateMaxNum.R, it requires double check
    # if(length(unique(repmax))>1){
    #   print(p)
    #   print(unique(unlist(repmax)))
    # }
    maxmat[p,m1] <- max(repmax,na.rm = T)
    minmat[p,m1] <- min(repmax,na.rm = T)
  }
  #print(paste0(m1,"...finished"))
  return(list("maxmat"=maxmat,"minmat"=minmat))
}
