
filter_by_frac_xa_ratio<-function(cutoff,x_expr,a_expr){
  #cutoff=0.1; #expr.gene cutoff
  (frac1=sum(x_expr>cutoff)/length(x_expr))
  (frac2=sum(a_expr>cutoff)/length(a_expr))
  select.frac=ifelse(frac1>frac2,frac2,frac1)
  x_expr=sort(x_expr,decreasing = T);
  a_expr=sort(a_expr,decreasing = T);
  x_expr.inp=x_expr[1:floor(select.frac*length(x_expr))]
  a_expr.inp=a_expr[1:floor(select.frac*length(a_expr))]
  #c(median(x_expr.inp)/median(a_expr.inp),length(x_expr.inp),length(a_expr.inp))
  c(mean(x_expr.inp)/mean(a_expr.inp),length(x_expr.inp),length(a_expr.inp))
  #mean(x_expr.inp)/mean(a_expr.inp)
  #ks.test(x_expr.inp,a_expr.inp)
}

filter_by_expr_xa_ratio<-function(cutoff,x_expr,a_expr){
  #cutoff=0.1; #expr.gene cutoff
  x_expr.inp=x_expr[x_expr>cutoff]
  a_expr.inp=a_expr[a_expr>cutoff]
  #c(median(x_expr.inp)/median(a_expr.inp),length(x_expr.inp),length(a_expr.inp))
  c(mean(x_expr.inp)/mean(a_expr.inp),length(x_expr.inp),length(a_expr.inp))
  #mean(x_expr.inp)/mean(a_expr.inp)
  #ks.test(x_expr.inp,a_expr.inp)
}