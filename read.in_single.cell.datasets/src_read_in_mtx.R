
read_mtx<-function(mtx,cells,genes){
  
  # Load Data
  ## mtx sparse matrix
  #df.expr=readMM("GSE135810/GSE135810_RAW/GSM4030593_matrix_sample_1.mtx.gz")
  data=as.matrix(readMM(mtx))
  #dim(df.expr) # 17492   630
  #df.expr[1:3,1:3]
  
  ## gene ID 
  #df.gene=read.table("GSE135810/GSE135810_RAW/GSM4030593_genes_sample_1.tsv.gz");
  feature=read.table(genes)
  #dim(df.gene) #17492 genes
  #colnames(df.gene)=c("flybase",'symbol')
  
  ## cell barcode
  #barcode=read.table("GSE135810/GSE135810_RAW/GSM4030593_barcodes_sample_1.tsv")
  barcode=read.table(cells)
  
  ## check dimension
  dim(data)
  dim(feature)
  dim(barcode)
  
  colnames(data) <- barcode[,1]
  rownames(data) <- feature[,2]
  data <- as(data, Class = "dgCMatrix")
  return(data)
}
