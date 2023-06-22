#sample data
a = 1:10;
b = as.matrix(rbind(a,a^2,a^3,a^4,a^5))

# Creating functions
#Vn = @(Aij, Bij) (1/(n*(n-3))).*(sum(sum(Aij.*Bij)) - (n/(n-2))*diag(Aij)'*diag(Bij));
#Rn = @(Aij, Bij) Vn(Aij, Bij)./sqrt(Vn(Aij, Aij).*Vn(Bij, Bij));

#x1=matrix(1:16,ncol=4,byrow =T)
#x2=matrix(rep(1:4,each=4),ncol=4,byrow=T)

Vn = function(x1,x2,n) {1/(n*(n-3)) * (sum(apply(x1*x2,2,sum)) -(n/(n-2))*diag(x1) %*% diag(x2)) }
Rn = function(x1,x2,n) { Vn(x1,x2,n) / sqrt( Vn(x1,x1,n)*Vn(x2,x2,n) ) }

gcl = function(data,num_divisions){
  num_genes=nrow(data)
  n=ncol(data)
  gcl_output=as.numeric();
  
  for(i in 1:num_divisions){
    
    # Random divisions
    random_genes = sample(1:num_genes,num_genes,replace=F);
    X1 = t(data[random_genes[1:floor(num_genes/2)],])
    X2 = t(data[random_genes[(floor(num_genes/2)+1):num_genes],])
    
    # Aij matrices
    d1 = as.matrix(dist(X1,method='euclidean',upper=T,diag=T))
    m1 = apply(d1,2,mean);
    M1 = mean(d1);
    Aij1 = d1 - m1 %o% rep(1,n);
    Aij1 = Aij1-rep(1,n) %o% m1;
    Aij1 = Aij1 + M1;
    Aij1 = Aij1 - d1/n;
    diag(Aij1) = m1 - M1;
    Aij1 = (n/(n-1))*Aij1;
    
    d2 = as.matrix(dist(X2,method='euclidean',upper=T,diag=T))
    m2 = apply(d2,2,mean);
    M2 = mean(d2);
    Aij2 = d2 - m2 %o% rep(1,n);
    Aij2 = Aij2-rep(1,n) %o% m2;
    Aij2 = Aij2 + M2;
    Aij2 = Aij2 - d2/n;
    diag(Aij2) = m2 - M2;
    Aij2 = (n/(n-1))*Aij2;
    
    # Calculating bcdcorr
    gcl_output[i] = Rn(Aij2, Aij1,n);
  }
  gcl_output
}

gcl(b,5)
