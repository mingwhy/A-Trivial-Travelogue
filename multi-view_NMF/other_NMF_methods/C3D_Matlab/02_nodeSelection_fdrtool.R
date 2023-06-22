#df=readRDS('../brain_scRNA-seq/all_common_genes.rds')
#genes=df$common.genes
#length(genes)
cut=0.01;#qvalu<0.001
min.size=5; #size>min.size
library(R.matlab) 
#x=readMat('sc_brain_58/vx_lambda.mat');
#x=readMat('sn_brain_100/vx_lambda.mat');
#x=readMat('sc_male_29/vx_lambda.mat');
#x=readMat('sc_female_29/vx_lambda.mat');
#x=readMat('sn_female_50/vx_lambda.mat');
#x=readMat('sn_male_50/vx_lambda.mat');
#x=readMat('sn_male_50_rep2/vx_lambda.mat');
#x=readMat('sn_1kcells/vx_lambda.mat');
x=readMat('sn_500cells//vx_lambda.mat');
#x=readMat('sc_500cells//vx_lambda.mat');
V=x$V
dim(V)

candidate.modules=list()
for(i in 1:ncol(V)){
        r<-tryCatch(
                {fdr.out=fdrtool::fdrtool(V[,i],statistic="normal",plot=FALSE)},
                error=function(c)'error')
        if(r=='error'){ 
                cat('failed')
                next
        }
        
        #plot(V[,i],fdr.out$qval)
        qval_all=fdr.out$qval
        #hist(qval_all)
        if(sum(qval_all<cut)<min.size){next}
        candidate.modules[[as.character(i)]]=which(qval_all<cut)
}

length(candidate.modules)
sapply(candidate.modules,length)

id=rep(1:length(candidate.modules),
       sapply(candidate.modules,length))
df=data.frame(id=id,gene=unlist(candidate.modules))
#writeMat("sc_male_candi.modules.mat",module.df=df)
#writeMat("sc_female_candi.modules.mat",module.df=df)
#writeMat("sn_1kcells_modules.mat",module.df=df)
writeMat("sn_500cells_modules.mat",module.df=df)
#writeMat("sc_500cells_modules",module.df=df)

