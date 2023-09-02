#two matrix with colname
#ub: col,cell.type_age, row, reactions
#lb

(files=Sys.glob('imat_out/*rds'))
df.info=c();
df.flux.vec=c();

#x=data.table::fread('log1p_female_head.csv')
#y=data.table::fread('log1p_female_body.csv')
#pick.tc=c(colnames(x)[-1],colnames(y)[-1])
#length(pick.tc) #108/4=27 
#pick.tc=gsub('\\/','\\-',pick.tc) #filename adjustment

for(file in files){
  x=strsplit(file,';')[[1]]
  #tmp=gsub('.rds','',paste(x[-1],collapse = ';'))
  #if(sum(tmp %in% pick.tc)==0){next}
  
  tissue=x[[2]]
  cell.type=x[[3]];
  age=x[[4]];
  age=gsub('.rds','',age)
  #cat(tissue,cell.type,',', age,'\n')
  df.info=rbind(df.info,c(tissue,cell.type,age))
  
  x=readRDS(file)
  names(x);
  # imat.R in gembox
  # the mode 1 of imat (the fva-like approach): solve the iMAT MILP under the forced inactivaton/activation of each rxn, determine (de)activated rxns based on the resulting objective values
  # return a list(solver.out, flux.int.imat), solver.out is a data.table of the optimal iMAT objectives for all rxns, flux.int.imat is a vector in the order of the model rxns, 
  # with values 0/9/1/-1 representing a rxn being inactive, activity level not enforced, active in the forward direction, and active in the backward direction as determined by iMAT
  
  #length(x$imat.model$fluxes.int) # 11898, 1  0 -1
  length(x$imat.model$fluxes.int.imat) # 11898, 9  1  0 -1
  #length(x$imat.model$solver.out[[1]]$xopt) #14495
  flux.vec=x$imat.model$solver.out[[1]]$xopt[1:length(x$imat.model$fluxes.int.imat)]
  #summary(flux.vec[x$imat.model$fluxes.int.imat==0]) # (-0.1,0.1)
  #summary(flux.vec[x$imat.model$fluxes.int.imat==-1]) #(-1000,-1)
  #summary(flux.vec[x$imat.model$fluxes.int.imat==1]) #(1,1000)
  #summary(flux.vec[x$imat.model$fluxes.int.imat==9]) #(-1000,1000)
  df.flux.vec=cbind(df.flux.vec,flux.vec)
}
colnames(df.info)=c('tissue','cell.type','age')
saveRDS(list(df.info=df.info,df.flux.vec=df.flux.vec),
             'head_body_flux_vec.rds')
data.table::fwrite(df.info, 'head_body_cell.type_age.txt')

#library(R.matlab)
#writeMat('head_body_flux_vec.mat', df.flux.vec=df.flux.vec)
#load('body_lb_ub.mat') in matlab
#size(imat_ub) 

