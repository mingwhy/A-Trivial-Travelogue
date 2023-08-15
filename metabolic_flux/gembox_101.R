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
(files=Sys.glob('../my.utils/R/*R')) #download `ImNotaGit/my.utils` from github
# https://github.com/ImNotaGit/my.utils
# comment out `rmote::xxx` in `utils.R`, install 'httpgd', then source files.
#install.packages("httpgd")
lapply(files, source)

# an example from Genome-scale metabolic modeling reveals SARS-CoV-2-induced metabolic changes and antiviral targets
# https://www.embopress.org/doi/full/10.15252/msb.202110260
# full code repo: https://github.com/ruppinlab/covid_metabolism
# download required data on google drive 
# `Alternatively, one can manually download the data files from this Google Drive link then decompress them into the data folder.`
# you'd have `covid_metabolism_files.tar.gz`.
# there is `de.and.gsea.res.RData` in folder `expression`.
# go to `covid_metabolism/GEM` folder, run `prepare.data.R` and `imat.and.mta.R` scripts.


