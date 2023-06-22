
source('pipeline_codes/src_significantModules_diff.size.R')

input.folder='./brain_scRNA-seq/';
genes=readRDS(paste0(input.folder,'all_common_genes.rds'))
c.genes=unlist(genes$genes)
all.genes=genes$all.genes
common.genes=genes$common.genes
#gene.names=genes$all.genes
gene.names=names(which(table(c.genes)>=length(genes$genes)/2))
gene_count=length(gene.names) #3015 for sc, 1351 for sn


## test common.modules.genes on fake modules
out.folder='./sc_useall_n10c0.1_T3/';

common.modules.genes=readRDS(paste0(out.folder,'K25_common.modules.genes.rds'))

fake=lapply(common.modules.genes,function(i){
  sample(gene.names,length(i),replace = F)})


## read in pearson.mat
files=Sys.glob(paste0(input.folder,'/pearson_coexpr/*rds'));
files=files[1:2]
cell.cluster.names=as.character(unlist(sapply(files,function(i){
  gsub('ncell_\\d+_','',basename(i),'_pearson.mat')
})))
cat('read in nets',length(files),'\n')


min.cell=100;max.cell=Inf;
nets=list();
for(file in files){
  mat=abs(readRDS(file))
  ncell=as.numeric(strsplit(basename(file),'ncell\\_|\\_')[[1]][[2]])
  if( any(ncell>max.cell , ncell<min.cell) ){next}
  diag(mat)=0;
  x=intersect(gene.names,rownames(mat))
  nets[[basename(file)]]=mat[x,x];
}
sapply(nets,dim)


permu_times=100;
fake.module.in.nets=significantModules_diff.size(modules=fake, 
                                            multiNetworks=nets,
                                           permu_times=permu_times,numCores=2)

i=apply(fake.module.in.nets$FDR,1,function(i) sum(i<0.05))
i;

saveRDS(fake.module.in.nets,
        paste0(out.folder,'K25_negative.control_test.on.nets.rds'))


