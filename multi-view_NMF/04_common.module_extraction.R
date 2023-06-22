
input.folder='brain_snRNA-seq_n10c0.1/'
#input.folder='brain_scRNA-seq//'
files=Sys.glob(paste0(input.folder,'/pearson_coexpr/*rds'));
cell.cluster.names=as.character(unlist(sapply(files,function(i){
  gsub('ncell_\\d+_','',basename(i),'_pearson.mat')
})))

x=readRDS('./sn_useall_n10c0.1_T3/K17_module.fdr.rds')
#x=readRDS('./sc_useall_n10c0.1_T3/K25_module.fdr.rds')
module.fdr=x$FDR
dim(module.fdr) #7module by 108clusters
colnames(module.fdr)=cell.cluster.names
rownames(module.fdr)=paste0('module',1:nrow(module.fdr))
apply(module.fdr,1,function(i) sum(i<0.05))
ncol(module.fdr)*0.8
#log.module.fdr=-log10(module.fdr)
module.fdr1=module.fdr
module.fdr1[module.fdr<0.05]=0
module.fdr1[module.fdr1!=0]=1

#pdf("sc_pvalue.pdf",useDingbats = T,height = 16)
pdf("sn_pvalue.pdf",useDingbats = T,height = 16)
pheatmap::pheatmap(t(module.fdr1),main='sn')
dev.off()

##
module=readRDS('./sc_useall_n10c0.1_T3/K25_module.genes.rds')
module.go=readRDS('./sc_useall_n10c0.1_T3/K25_modules.GO.rds')
sink('T3_sc_GO_terms.txt');
module
lapply(1:length(module.go),function(k) head(module.go[[k]]@result$Description,10))
sink()

module=readRDS('./sn_useall_n10c0.1_T3/K17_module.genes.rds');
module.go=readRDS('./sn_useall_n10c0.1_T3/K17_modules.GO.rds')
sink('T3_sn_GO_terms.txt');
module
lapply(1:length(module.go),function(k) head(module.go[[k]]@result$Description,10))
sink()

