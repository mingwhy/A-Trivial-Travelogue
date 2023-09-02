#two matrix with colname
#ub: col,cell.type_age, row, reactions
#lb

(files=Sys.glob('imat_out/*rds'))
df.info=c();
df.ub=c();
df.lb=c();

for(file in files){
  x=strsplit(file,';')[[1]]
  cell.type=x[[3]];
  age=x[[4]];
  age=gsub('.rds','',age)
  #cat(cell.type,',', age,'\n')
  df.info=rbind(df.info,c(cell.type,age))
  
  x=readRDS(file)
  names(x)
  ub=x$result.model$ub
  lb=x$result.model$lb
  unique(ub);unique(lb)
  df.ub=cbind(df.ub,ub)
  df.lb=cbind(df.lb,lb)
}
colnames(df.info)=c('cell.type','age')
data.table::fwrite(df.info, 'body_cell.type_age.txt')

library(R.matlab)
writeMat('body_lb_ub.mat', imat_ub=df.ub,imat_lb=df.lb)
#load('body_lb_ub.mat') in matlab
#size(imat_ub) 

