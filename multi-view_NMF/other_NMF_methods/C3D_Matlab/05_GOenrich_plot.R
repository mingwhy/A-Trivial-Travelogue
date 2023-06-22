
#genes=readRDS('brain_scRNA-seq/all_common_genes.rds')
genes=readRDS('../brain_snRNA-seq/all_common_genes.rds')

names(genes)
common.genes=genes$common.genes
length(common.genes)

library(R.matlab)
#x=readMat('ConMod_SP_female_male_common.mat') #strength, participate
#x=readMat('ConMod_IP_female_male_common.mat') #participate, information
#x=readMat('brain_snRNA-seq/ConMod_SI_female_male_common.mat') #strength,information
#x=readMat('brain_snRNA-seq/ConMod_SI_female_common.mat') #strength,information
#x=readMat('brain_snRNA-seq/ConMod_SI_male_common.mat') #strength,information

#x=readMat('./C3D_Matlab/sc_C3D_out.mat') #strength,information
#x=readMat('./C3D_Matlab/sn_C3D_out.mat') #strength,information
#x=readMat('./sn_1kcells_modules.fdr.mat')
x=readMat('./sn_500cells_modules.fdr.mat')
#x=readMat('./C3D_Matlab/sc_500cells_modules.fdr.mat')
#x=readMat('./C3D_Matlab/sn_C3D_500.mat') 
#x=readMat('./C3D_Matlab/sc_C3D_female_out.mat') #strength,information
#x=readMat('./C3D_Matlab/sc_C3D_male_out.mat') #strength,information
#x=readMat('./C3D_Matlab/sn_C3D_female_out.mat') #strength,information
#x=readMat('./C3D_Matlab/sn_C3D_male_out.mat') #strength,information
#x=readMat('./C3D_Matlab/sn_C3D_male_out_rep2.mat') #strength,information

names(x)
cell.types=unlist(x$filenames)
dim(x$FDR2) #module by cell.type
#dim(x$Hc)
length(x$modulesFinal)
x$modulesFinal[[1]]
sapply(x$modulesFinal,function(i)nrow(i[[1]]))

dim(x$pvalues.modulePerNet) #module by cell.types


dim(x$FDR2) #module by cell.types
i=apply(x$FDR2,1,function(i) sum(i<0.05))
max(i)
sum(i==ncol(x$FDR2)) #12 module common to all
sum(i>=ncol(x$FDR2)*0.9) #12 module common to all
sum(i>=ncol(x$FDR2)*0.8) #12 module common to all
sum(i!=0) #28


#common.modules=x$modulesFinal[i==ncol(x$FDR2)]
#common.modules=x$modulesFinal[i>=1 ]
common.modules=x$modulesFinal[i>=ncol(x$FDR2)*0.8 ]
length(common.modules)

#common.modules[[1]]
common.modules.genes=lapply(common.modules,function(i) 
  common.genes[unlist(i)])
sapply(common.modules.genes,length)

#saveRDS(common.modules.genes,'C3D_sc_modules.rds')
#saveRDS(common.modules.genes,'C3D_sc_female_modules.rds')
#saveRDS(common.modules.genes,'C3D_sc_male_modules.rds')
#saveRDS(common.modules.genes,'C3D_sn_female_modules.rds')
#saveRDS(common.modules.genes,'brain_snRNA-seq/female_male_common.rds')
#saveRDS(common.modules.genes,'brain_snRNA-seq/female_common.rds')
#saveRDS(common.modules.genes,'brain_snRNA-seq/male_common.rds')

i=apply(x$FDR2,1,function(i) sum(i<0.05))
tmp=x$FDR2[i>0,]
tmp1=tmp;
tmp[tmp1<0.05]=1;
tmp[tmp1>=0.05]=0
dim(tmp)
colnames(tmp)=unlist(x$filenames)
pdf('test.pdf',useDingbats = T)
print(pheatmap::pheatmap(tmp, fontsize = 6))
dev.off()


## GO enrich
source('../src_fly.gene_GOenrich.R')
GO.out=lapply(common.modules.genes,function(genes){
  GOenrich(genes)
})
sapply(GO.out,nrow)
head(GO.out[[1]]$Description)
for(i in 1:length(GO.out)){
  cat('module',i,'\n')
  print(head(GO.out[[i]]$Description))
}

saveRDS(GO.out,'sn_500cells_moduleGO.rds')


