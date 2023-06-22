#Co-dependency index score
#https://academic.oup.com/nar/article/49/18/e104/6324613#304882741

#dim(m1)
#cor.mat=matrix(0,nrow=nrow(m1),ncol=nrow(m1))
#cor.mat1=CDI(m1)
#for(i in 1:(nrow(m1)-1)){
#  for(j in (i+1):nrow(m1)){
#    cor.mat[i,j]=CDI_2gene(m1[i,],m1[j,])
#  }
#}

CDI_2gene<-function(i,j){
  i=i>0;
  j=j>0
  ncell=length(i)
  pi=sum(i)*sum(j)/ncell^2
  obs=sum(i+j==2) #joint occur
  #pbinom(obs,ncell,pi)
  #sum(dbinom(obs:ncell,ncell,pi))
  p.tail=pbinom(obs-1,ncell,pi,lower.tail = F) 
  cdi=-1*log10(p.tail)
  return(cdi)
}

CDI<-function(mat){
  bm1=mat;
  bm1[mat>0]=1
  cov1=bm1 %*% t(bm1)
  ncell=ncol(mat)
  ngene=nrow(mat)

  gene.freq=Matrix::rowSums(bm1)/ncell
  gene.freq.mat=gene.freq %*% t(gene.freq)
  #dim(gene.freq.mat)

  pi=gene.freq.mat[upper.tri(gene.freq.mat)]
  #sum(cov1==ncol(m1))

  obs=cov1[upper.tri(cov1)]
  p.tail=pbinom(obs-1,ncell,pi,lower.tail = F) 
  cdi=-1*log10(p.tail)
  #summary(cdi)

  cdi.mat=diag(ngene)
  cdi.mat[upper.tri(cdi.mat)]<-cdi
  cdi.mat=t(cdi.mat)
  cdi.mat[upper.tri(cdi.mat)]<-cdi
  #isSymmetric(cdi.mat)
  return(cdi.mat)
}
