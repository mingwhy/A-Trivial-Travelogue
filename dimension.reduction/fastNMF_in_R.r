#https://cran.r-project.org/web/packages/RcppML/vignettes/RcppML.html
sessionInfo()
#R version 4.0.2 (2020-06-22) or R version 4.1.3 (2022-03-10)
#Platform: x86_64-apple-darwin17.0 (64-bit)
#Running under: OS X  12.1
if(F){
  devtools::install_github("zdebruine/RcppSparse")  #https://github.com/zdebruine/RcppSparse
  devtools::install_github("zdebruine/RcppML") # #https://github.com/zdebruine/RcppML
}

#https://github.com/zdebruine/RcppML
library(RcppML)
library(Matrix)
packageVersion("RcppML") # 0.5.4
?RcppML::nmf

#############################
# basic NMF
A =rsparsematrix(1000, 100, 0.1)
dim(A) #1000 100
model <- nmf(A, k = 10)
model
dim(model@w) #1000 10
model@d #10 eigenvalues
dim(model@h) #10 100

# compare rank-2 NMF to second left vector in an SVD
data(iris)
A <- as(as.matrix(iris[,1:4]), "dgCMatrix")
nmf_model <- nmf(A, 2, tol = 1e-5)
nmf_model

bipartitioning_vector <- apply(nmf_model$w, 1, diff)
second_left_svd_vector <- base::svd(A, 2)$u[,2]
abs(cor(bipartitioning_vector, second_left_svd_vector))

# compare rank-1 NMF with first singular vector in an SVD
abs(cor(nmf(A, 1)$w[,1], base::svd(A, 2)$u[,1]))

# symmetric NMF
A <- crossprod(rsparsematrix(100, 100, 0.02))
dim(A)
model <- nmf(A, 10, tol = 1e-5, maxit = 1000)
plot(model$w, t(model$h))
# see package vignette for more examples

#############################
#https://cran.r-project.org/web/packages/RcppML/vignettes/RcppML.html
A <- rsparsematrix(100, 100, 0.1)
model <- RcppML::nmf(A, 10, verbose = F)

w <- model$w
d <- model$d
h <- model$h
model_tolerance <- tail(model$tol, 1)

#Mean squared error of a factorization can be calculated for a given model using the RcppML::mse function:
#RcppML::mse(A_sym, model$w, model$d, model$h)
#above not working, refer to https://cran.r-project.org/web/packages/RcppML/RcppML.pdf
#c_mse <- mse(A, model$w, model$d, model$h)
R_mse <- mean((A - model$w %*% Diagonal(x = model$d) %*% model$h)^2)
#all.equal(c_mse, R_mse)
R_mse
mse.out=sapply(1:10, function(i){
  #R_mse <- mean((A - model$w %*% Diagonal(x = model$d) %*% model$h)^2)
  R_mse <- mean((A - model$w[,1:i] %*% Diagonal(x = model$d[1:i]) %*% model$h[1:i,])^2)
  R_mse
})
plot(1:10,mse.out)

