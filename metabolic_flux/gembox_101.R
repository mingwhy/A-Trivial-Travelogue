# https://github.com/ruppinlab/gembox
# need to install gurobi (solver algorithm) (or Rcplex2: https://github.com/ruppinlab/Rcplex2)
# library(gurobi)
# check out install_gurobi.R
#devtools::install_github("ruppinlab/gembox")
library(gembox)
## lee et al. BMC Syst Biol 2012: https://github.com/ruppinlab/gembox/blob/master/R/exprs.R

#README.md from https://github.com/ruppinlab/covid_metabolism
# 1, download data
#$ cd covid_metabolism-main/data
#$pip3 install gdown #for use gdown in ./dnload.sh
#$ ./dnload.sh
# covid_metabolism_files.tar.gz (270MB), liao.RDS(3.6G), chua.RDS(1.7G) files.

# process data `collect.validation.data.R`
#BiocManager::install('GSA')    
#devtools::install_github("ImNotaGit/my.utils")
#library(my.utils)
(files=Sys.glob('../my.utils_R/*R')) #download `ImNotaGit/my.utils` from github
# comment out `rmote::xxx` then source files
#install.packages("httpgd")
lapply(files, source)
