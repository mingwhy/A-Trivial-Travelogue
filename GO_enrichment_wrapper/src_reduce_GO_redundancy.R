
########################################
#https://rdrr.io/bioc/clusterProfiler/man/simplify-methods.html
#References:
#https://github.com/YuLab-SMU/clusterProfiler/issues/162
# two must included column in res: res$ID and res$p.adjust
#avaiable similarity measures:https://www.bioconductor.org/packages/devel/bioc/manuals/GOSim/man/GOSim.pdf 
#measure One of "Resnik", "Lin", "Rel", "Jiang" "TCSS" and "Wang" methods.

# usage
#res$ID=res$GO
#res$p.adjust=res$FDR
#MmGO <- GOSemSim::godata('org.Mm.eg.db', ont="BP") #https://yulab-smu.top/biomedical-knowledge-mining-book/GOSemSim.html

# by=p.adjust, use min(p.adjust) to select retained GO terms
#simplify_internal(res,cutoff=0.7, by="p.adjust", select_fun=min, measure="Rel", ontology='BP',semData=MmGO)

# by=ngene, use max(ngene) to select retained GO terms
#res.simple=simplify_internal(res,cutoff=0.1, by="ngene", select_fun=max, measure="Wang", ontology='BP',semData=MmGO)

# by=abs(rho), use max(rho) to select retained GO terms
#res$abs.rho=abs(res$Spearman.rho)
#res.simple=simplify_internal(res,cutoff=0.25, by="abs.rho", select_fun=max,measure="Wang", ontology='BP',semData=MmGO)

simplify_internal <- function(res, cutoff=0.7, by="p.adjust", select_fun=min, 
                              measure="Wang", ontology='BP', semData) {
  if (missing(semData) || is.null(semData)) {
    if (measure == "Wang") {
      semData <- GOSemSim::godata(ont = ontology)
    } else {
      stop("godata should be provided for IC-based methods...")
    }
  } else {
    if (ontology != semData@ont) {
      msg <- paste("semData is for", semData@ont, "ontology, while enrichment result is for", ontology)
      stop(msg)
    }
  }
  
  sim <- GOSemSim::mgoSim(res$ID, res$ID,
                          semData = semData,
                          measure=measure,
                          combine=NULL)
  
  ## to satisfy codetools for calling gather
  go1 <- go2 <- similarity <- NULL
  
  sim.df <- as.data.frame(sim)
  sim.df$go1 <- row.names(sim.df)
  sim.df <- tidyr::gather(sim.df, go2, similarity, -go1)
  
  sim.df <- sim.df[!is.na(sim.df$similarity),]
  
  ## feature 'by' is attached to 'go1'
  sim.df <- merge(sim.df, res[, c("ID", by)], by.x="go1", by.y="ID")
  sim.df$go2 <- as.character(sim.df$go2)
  
  ID <- res$ID
  
  GO_to_remove <- character()
  for (i in seq_along(ID)) {
    ii <- which(sim.df$go2 == ID[i] & sim.df$similarity > cutoff) #similarity larger than cutoff
    ## if length(ii) == 1, then go1 == go2
    if (length(ii) < 2)
      next
    
    sim_subset <- sim.df[ii,]
    
    jj <- which(sim_subset[, by] == select_fun(sim_subset[, by])) #select one with minimal p.value, 
    #jj may contain multiple GO ids if they are share the same minimal p.value
    
    ## sim.df <- sim.df[-ii[-jj]]
    GO_to_remove <- c(GO_to_remove, sim_subset$go1[-jj]) %>% unique
  }
  GO_to_remove
  
  res[!res$ID %in% GO_to_remove, ]
}

########################################
