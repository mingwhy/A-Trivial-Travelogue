
calc_cumulative_information_of_set<-function(df_discrete,my_genes,df_labels){
  H_naive=entropy(df_labels)/log(2)
  df_temp= t(df_discrete[my_genes,,drop=F])
  
  x=apply(df_temp,1,function(i) paste(i,collapse = ''))
  names(x)=NULL
  
  H = 0
  # Find unique combinations of expression levels
  for(status in unique(x)){
    cell.names=rownames(df_temp[x==status,,drop=F]) ## Get cell names having this unique combination of expression levels
    labels_cond=df_labels[cell.names] ## Get class labels of cells
    
    # Calculate entropy of classification (conditional on expression levels)
    H_cond = entropy(labels_cond)/log(2)
    weight = length(cell.names) /nrow(df_temp) # Weight by fraction of cells in this set
    H = H+ H_cond*weight
  }
  
  I = H_naive - H
  return(I)
}


find_nonredundant_gene_set<-function(df_discrete, genes,
                               df_labels,df_info,
                               H_naive, N_constrain=length(genes),
                               cumulative_information_cutoff=0.999,
                                 verbose=FALSE){
    
  cis = c()
  nonredundant_genes = c()
  current_ci = 0
  current_relative_ci = 0.0
  
  # Rank genes by mutual information
  remaining_genes=names(sort(df_info[genes],decreasing = T))
  
  # begin loop, select gene in non-redundant set, then recalculate my_info_gains removing that selected gene
  i = 0
  while(1){
    i = i+1;
    
    my_info_gains<- foreach(i=1:length(remaining_genes),.combine="c") %dopar% {
      my_genes = c(nonredundant_genes,remaining_genes[i])
      my_ci = calc_cumulative_information_of_set(df_discrete, my_genes, df_labels)
      my_info_gain = my_ci - current_ci;
      my_info_gain}
    
    # too slow
    if(F){
      my_info_gains=c()
      for(gene in remaining_genes){
        my_genes = c(nonredundant_genes,gene)
        my_ci = calc_cumulative_information_of_set(df_discrete, my_genes, df_labels)
        my_info_gain = my_ci - current_ci
        my_info_gains=c(my_info_gains,my_info_gain)
      }
    }    
    names(my_info_gains)=remaining_genes
  
    # Sort genes by information gain
    my_info_gains=sort(my_info_gains,decreasing = TRUE) #there may be more than 1 top genes (equal info_gain)
  
    # Take best gene name
    hit = names(my_info_gains[1])
    
    nonredundant_genes=c(nonredundant_genes,hit)
    remaining_genes=remaining_genes[remaining_genes!=hit] #sequentially remove genes
    
    current_ci = current_ci + my_info_gains[hit]
    current_relative_ci = current_ci / H_naive
    
    if(verbose){
      print(hit)
      print(current_relative_ci)
    }
    if(length(remaining_genes) == 0 || current_relative_ci > cumulative_information_cutoff){
      cat('length of remaining_genes=',length(remaining_genes),'\n');
      cat('current_relative_ci=',current_relative_ci,'\n');
      return(nonredundant_genes)
    }
  }
}


calc_cumulative_informations<-function(df_discrete, genes, df_labels, N=5){
  # Calculate total information for each gene set defined by iteratively adding a single gene
  # starting from the beginning of the list genes
  # Calculates total information of gene sets starting from top of df
  cis = c();
  for(i in 1:N){
    my_ci = calc_cumulative_information_of_set(df_discrete, genes[1:i], df_labels)
    cis=c(cis,my_ci)
  }
  df_result = data.frame(symbol=genes,cumulative_mutual_information=cis)
  return(df_result)
}
