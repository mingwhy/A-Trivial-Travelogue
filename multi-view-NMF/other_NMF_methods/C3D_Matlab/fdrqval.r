#! /usr/bin/Rscript
# written by Xiaolin Xiao, 2013 

getwd()
a=getwd()
library(fdrtool)
fin=paste(a,"vx.txt",sep="/")
vin=read.table(fin)
vin=vin[,1]
fdr.out=fdrtool(vin,statistic="normal",plot=FALSE)
qval_all=fdr.out$qval
fout_qall=paste(a,"vx_outqval_all.txt",sep="/")
write.table(qval_all,fout_qall)

