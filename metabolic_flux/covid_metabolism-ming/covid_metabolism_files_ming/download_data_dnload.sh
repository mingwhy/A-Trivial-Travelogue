#!/bin/bash

#install gdown: pip3 install gdown: https://pypi.org/project/gdown/
# download required data, part 1 from a shared Google Drive link, the downloaded file name will be covid_metabolism_files.tar.gz
gdown https://drive.google.com/uc\?id\=1bVPCQlDR3G8TTMx09jkN3IuGwzq9n9hi
# decompress
tar xzvf covid_metabolism_files.tar.gz -C ../

# download required data, part 2: the two scRNA-seq datasets
# Liao et al.
wget http://cells.ucsc.edu/covid19-balf/nCoV.rds -O liao.RDS
# Chua et al.
wget https://ndownloader.figshare.com/files/22927382 -O chua.RDS

# github: https://github.com/ruppinlab/covid_metabolism
# download `covid_metabolism_files.tar.gz` from google drive link:
# https://drive.google.com/file/d/1bVPCQlDR3G8TTMx09jkN3IuGwzq9n9hi/view

