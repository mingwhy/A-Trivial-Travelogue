
import numpy as np
import pandas as pd
import scanpy as sc
sc.settings.verbosity = 3             # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_header()
sc.settings.set_figure_params(dpi=80, facecolor='white')


results_file = 'write/brain50k.h5ad'  # the file that will store the analysis results

adata = sc.read_10x_mtx(
    '/Users/ming/Documents/single.cell_datasets/fly.brain.atlas/GSE107451_DGRP-551_w1118_WholeBrain_57k_0d_1d_3d_6d_9d_15d_30d_50d_10X_DGEM_MEX.mtx.tsv/',  # the directory with the `.mtx` file
    #var_names='gene_symbols',     # use gene symbols for the variable names (variables-axis index)
    var_names='gene_ids',   
    cache=True)    # write a cache file for faster subsequent reading
adata
# AnnData object with n_obs × n_vars = 56902 × 17473
# var: 'gene_ids'

#adata.var_names_make_unique()  # this is unnecessary if using `var_names='gene_ids'` in `sc.read_10x_mtx`
#adata
adata.var
adata.var_names
tmp=adata.var['gene_symbols'].str.startswith('mt:') 
#tmp=adata.var_names.str.startswith('mt:') 
tmp.unique()
#Index([False, True, nan], dtype='object')
tmp.sum(skipna=True) #37 genes
np.where(tmp.isna()) #locate nan
# (array([13597]),)

adata.var.iloc[[13597]] #slice 1 row pandas dataframe 
#            gene_symbols
#FBgn0036414          NaN
adata.var.iloc[[13597]] = "nan"
adata.var.iloc[[13597]]

type(adata.var_names[13597]) #str or float


tmp=adata[:, adata.var_names.str.match('FBgn0036414')]
#View of AnnData object with n_obs × n_vars = 56902 × 1
#    var: 'gene_ids'
tmp=tmp.to_df()
tmp.sum()
#NaN    7.0
#dtype: float32
tmp1=tmp > 0
tmp1.sum()
#6


# Preprocessing
sc.pl.highest_expr_genes(adata, n_top=20, )


# basic filtering 
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=4)
#filtered out 4917 genes that are detected in less than 4 cells

adata
#AnnData object with n_obs × n_vars = 56902 × 12556

adata.var['mt'] = adata.var['gene_symbols'].str.startswith('mt:')   # annotate the group of mitochondrial genes as 'mt'
adata.var["mt"].sum() #29 mt genes

adata.var["mt"].unique()


sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
max(adata.obs.pct_counts_mt) #13%

#sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],jitter=0.4, multi_panel=True)
sc.pl.scatter(adata, x='total_counts', y='pct_counts_mt')
sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts')

sum(adata.obs.n_genes>200) # 56902, all expr>200 gene
sum(adata.obs.total_counts > 500) #56875
sum(adata.obs.pct_counts_mt < 30) #30%,56902


adata = adata[adata.obs.n_genes >200, :]
adata = adata[adata.obs.total_counts > 500, :]
adata = adata[adata.obs.pct_counts_mt < 30, :]
adata
#View of AnnData object with n_obs × n_vars = 56875 × 12556

#https://github.com/theislab/anndata/issues/274
adata.layers['counts'] = adata.X.copy()

#Total-count normalize (library-size correct) the data matrix 𝐗 to 10,000 reads per cell, so that counts become comparable among cells.
sc.pp.normalize_total(adata, target_sum=1e4)
#Logarithmize the data:
sc.pp.log1p(adata)
#adata.layers['log1p'] = adata.X.copy()

#Identify highly-variable genes.
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
sc.pl.highly_variable_genes(adata)


#create a "fake-sparse" form of the scaled data in .X, so that dimensionality can remain unchanged
#(and lognorm data can live in a layer)
#bdata = adata[:, adata.var['highly_variable']]
#adata.var["highly_variable"].sum() #1029genes

adata.var["highly_variable"].sum() #1029genes
bdata = adata.layers['counts'][:, adata.var['highly_variable']]
bdata
#<56875x1029 sparse matrix of type '<class 'numpy.float32'>'

ndata = adata[:, adata.var['highly_variable']]
ndata.X.A

#construct matrix as lil, which is row based. the HVG index is column based, so two transpositions are in order
adata.shape
#56875, 12556
adata.X
bdata[0:10, 0:10].A #raw umi count

#import scipy.sparse as sps
#X = sps.lil_matrix(adata.shape[::-1])
#X
#12556x56875 sparse matrix
#X[np.arange(adata.shape[1])[adata.var['highly_variable']], :] = bdata.X.T
#adata.X = X.tocsc().T


##############################################
from time import perf_counter as pc
import numpy as np
import scipy.sparse as sps
from sklearn.decomposition import NMF

# NMF 
model = NMF(n_components=20, init='random', max_iter=1000,random_state=0, verbose=True)

start_time = pc()
#W = model.fit_transform(bdata)
ndata.X.A.T.shape #gene by cell
W = model.fit_transform(ndata.X.A.T) 
end_time = pc()

H = model.components_
W.shape
H.shape

print('Used (secs): ', end_time - start_time)
#56875x1029 sparse matrix, Used (secs):  271
print(model.reconstruction_err_)
print(model.n_iter_)

# save mtx format, compatible between R and python
from scipy import sparse
from scipy import io
sW = sparse.csr_matrix(W) 
io.mmwrite('W_matrix.mtx',sW) #basis mat.
sH = sparse.csr_matrix(H) 
io.mmwrite('H_matrix.mtx',sH) #embedding mat
#library(Matrix)
#readMM('W_matrix.mtx')->W

# save  ndata.var and ndata.obs to get cell.barcode and gene names
ndata.var.to_csv('ndata_gene.csv',sep='\t')
ndata.obs.to_csv('ndata_cell.csv',sep='\t')





