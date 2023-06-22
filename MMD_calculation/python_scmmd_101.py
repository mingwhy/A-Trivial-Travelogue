################################################################################
# https://github.com/calico/scmmd
# scMMD -- Comparing cell populations using the maximum mean discrepancy (MMD)

git clone https://github.com/josipd/torch-two-sample.git
cd torch-two-sample
pip install .

pip install scmmd

################################################################################
## debug 
# in ipython
import scmmd

error message: OMP: Error #15: Initializing libomp.dylib, but found libomp.dylib already initialized.
# https://stackoverflow.com/questions/53014306/error-15-initializing-libiomp5-dylib-but-found-libiomp5-dylib-already-initial

## debug 
pip install configargparse

## try and works
import os
os.environ['KMP_DUPLICATE_LIB_OK']='True'
import scmmd
import anndata

# load an anndata object
adata = anndata.read_h5ad(
    'kang_2017_stim_pbmc.h5ad',
)

distances, p_values = scmmd.compute_mmd_contrast(
    adata=adata, # [Cells, Genes] object
    representation='X_pca', # representation to use, "X" or key in `adata.obsm`.
    groupby='cell', # a categorical grouping variable in `adata.obs`
    contrast='stim', # a binary contrast in `adata.obs`
    n_iters=100, # number of random sampling iterations
    sample_size=500, # sample size for random samples
    n_permutations=1000, # permutations for p-val calculations
)

################################################################################
## try another solution 
https://stackoverflow.com/questions/53014306/error-15-initializing-libiomp5-dylib-but-found-libiomp5-dylib-already-initial
conda install nomkl
scmmd -h 

################################################################################
## run demo https://github.com/calico/scmmd/blob/master/demo/example.sh
$ which scmmd
/Users/xxx/anaconda3/bin/scmmd

wget "https://storage.googleapis.com/calico-website-mca-storage/kang_2017_stim_pbmc.h5ad" -O kang_2017_stim_pbmc.h5ad
scmmd --data kang_2017_stim_pbmc.h5ad --out_path ./ --groupby cell --contrast stim --representation X_pca
# the above command generate two npy files, read them into R for further processing

## inside R
library(zellkonverter)
library(SingleCellExperiment)
inp.sce<-readH5AD('kang_2017_stim_pbmc.h5ad');     
meta=colData(inp.sce)
table(meta$cell) #8 cell types
table(meta$cell,meta$stim)

library(reticulate)
np <- import("numpy")
dist.out=np$load('mmd_distances_cell_stim.npy')
pval.out=np$load('mmd_p_values_cell_stim.npy')
dim(dist.out) #8 100   3 (between, within, within)
dim(pval.out) #8 100   3

dim(dist.out[1,,]) #1st cell type, 100x3



