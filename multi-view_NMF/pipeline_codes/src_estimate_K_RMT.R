#' Title Estimate the number of components
#' an Random Matrix theory based component number estimation
#'
#' @param data.m a matrix object which need to be estimated its number of components
#' @param default boolean variable to determine whether restrict the maximum number of decomposed components
#' @param svd.num an integer value to indicate maximum number of decomposed components
#'
#' @return a integer value which indicating the estimated number of ica components
#' @export
#' @importFrom coop tpcor
#' @importFrom rARPACK svds
#'
estimate_K <- function(data.m){    
    require(coop) #faster than Rfast
    require(rARPACK)
    data.m <- data.m[rowSums(data.m)>0,] #gene by cell

    M <- apply(data.m,2,function(X){ (X - mean(X))/sqrt(var(X))}); #cell-wise scale    

    sigma2 <- var(as.vector(M));
    Q <- nrow(data.m)/ncol(data.m);
    lambdaMAX <- sigma2*(1+1/Q + 2*sqrt(1/Q));
    lambdaMIN <- sigma2*(1+1/Q - 2*sqrt(1/Q));
    
    C <- coop::tpcor(t(M))
    #C <- Rfast::cora(M)
     
    if(ncol(data.m)>1000){
        eigen.o <- svds(C,k=1000)$d
    }else{
        eigen.o <- svd(C)$d
    }

    intdim <- length(which(eigen.o > lambdaMAX))
    
    if(intdim==50){
        cat('warning! the number of significant component is too large')
    }

    rm(M,data.m)

    #return
    intdim
}
