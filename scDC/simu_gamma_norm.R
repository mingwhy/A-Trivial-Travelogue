
#simulation scheme: param.est.pnb.R function from https://github.com/anlingUA/scDoc 
#f(y) = pi*Gamma(y;alpha,beta) +(1-pi)*Normal(y;mu,sigma)
r_gamma_normal <- function(n, pi, mu, sigma2, alpha, beta) {
  if (length(n) > 1) n <- length(n)
  u <- runif(n) 
  y <- apply(as.matrix(u), 1, function(x) 
    ifelse(x<=pi, rgamma(1, shape=alpha, rate=beta),
           rnorm(1, mu,sd=sqrt(sigma2)) ) )
  #y1 <- rnorm(floor(n*pi), mu, sd=sqrt(sigma2))
  #y2 <- rgamma(n-length(y1), shape=alpha, rate=beta)
  #y=c(y1,y2)
  return(y)
}

# usage
#y1=r_gamma_normal(n=2000,pi=0.3,mu=10,sigma2=1,alpha=1,beta=1)
#y2=r_gamma_normal(n=2000,pi=0.1,mu=6,sigma2=2,alpha=0.5,beta=1)
#count=t(data.frame(y1=y1,y2=y2)) 
#dim(count) #two genes, each gene has 2000 cell info

