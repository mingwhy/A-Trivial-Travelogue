
setwd("~/Downloads/metabolic-network/");
load('XQmzdataforMing')
ls() #"gmeans" "panel" 
head(panel)
all.mz=panel$KEGG.ID[panel$KEGG.ID!='']
length(all.mz)  #81

setwd("~/Downloads/metabolic-network/nicepath/nicepath")

if(FALSE){  #run once
  #81 mz, 5hr
  library(reticulate)
  Sys.setenv(RETICULATE_PYTHON = "/Users/mingyang/anaconda3/envs/network/bin/python")
  py_config()
  sys=import('sys')
  sys$path
  # make sure you have 'decibel_module_gcl.py' in the working env
  sys$path=c(sys$path,'/Users/mingyang/Downloads/metabolic-network/nicepath/nicepath/') #sys.path.append('./decibel/module/')
  sys$path 
  
  #mz1='C01197'; #target
  #mz2='C00021'; #precursor
  for(mz1 in all.mz){
    for(mz2 in all.mz){
      if(mz1==mz2){next}
      
      #output.file='Pathways_C00378_C02494.txt'
      output.file=paste0('Pathways_',mz1,'_',mz2,'.txt')
      if(file.exists(output.file)){next}
      
      cmd <- paste("perl", "../input/KEGG/modify_parameters.pl", mz1,mz2,
                   '../input/KEGG/parameters_template.txt > ../input/KEGG/parameters.txt')
      system(cmd)
      system('python main.py KEGG')
      
    }
  }

}
# check for 'symmetry'
file1='../output/KEGG/Pathways_C10164_C00019.txt';
file2='../output/KEGG/Pathways_C00019_C10164.txt';
df1=as.data.frame(data.table::fread(file1))
df2=as.data.frame(data.table::fread(file2))
df1[10,];df2[10,]
select.col=c(1,3,4,5,7,8,12,13);
View(df1[,select.col])

##########################################################
# change 'KEGG' to 'KEGG_81mz'
setwd("~/Downloads/metabolic-network/nicepath/output/KEGG_81mz/");
# organize files into a R list
files=Sys.glob('Pathways_C*_C*.txt')
length(files) #4528 files
length(all.mz)*(length(all.mz)-1) #6480
choose(length(all.mz),2) #3240

paths.between.mz.pairs=list()
for(file in files){
  a=unlist(strsplit(gsub('.txt','',file),'\\_'))
  mz1=a[2]
  mz2=a[3]
  df=as.data.frame(data.table::fread(file))
  select.col=c(1,3,4,5,7,8,12,13);
  df=df[,select.col]
  pair.name=gsub('.txt','',file)
  paths.between.mz.pairs[[pair.name]]<-df
}
length(paths.between.mz.pairs) #4528
sort(sapply(paths.between.mz.pairs,nrow))
table(sapply(paths.between.mz.pairs,nrow))
paths.between.mz.pairs[['Pathways_C00666_C00388']]
saveRDS(paths.between.mz.pairs,'paths.between.mz.pairs.rds')

paths.between.mz.pairs=readRDS('paths.between.mz.pairs.rds')
names(paths.between.mz.pairs) #all mz pairs with at least one pathway identified result
paths.between.mz.pairs[['Pathways_C00666_C00388']]
paths.between.mz.pairs[['Pathways_C00388_C00666']] #they are the same, both are included in the R object

###########################################
if(F){
#for(i in 1:(length(all.mz)-1)){
#  mz1=all.mz[i]
#  for(j in (i+1):length(all.mz)){
#    mz2=all.mz[j]
for(mz1 in all.mz){
  for(mz2 in all.mz){
    if(mz1==mz2){next}
    #output.file='Pathways_C00378_C02494.txt'
    output.file=paste0('Pathways_',mz1,'_',mz2,'.txt')
    file.exists(output.file)
    if(!file.exists(output.file)){next}
    df=as.data.frame(data.table::fread(output.file))
    select.col=c(1,3,4,5,7,8,12,13);
    df=df[,select.col]
    pair.name=gsub('.txt','',output.file)
    paths.between.mz.pairs[[pair.name]]<-df
  }
}
length(paths.between.mz.pairs) #1753
}
