##################################################################################
## installation and tutorial
# https://github.com/rajewsky-lab/novosparc/blob/master/reconstruct_drosophila_embryo_tutorial.ipynb
# NovoSpaRc: Nature Protocols: https://www.nature.com/articles/s41596-021-00573-7
# https://github.com/rajewsky-lab/novosparc/blob/master/reconstruct_tissue.py

#https://novosparc.readthedocs.io/tutorials.html#the-drosophila-embryo
#reconstruct_bdtnp_with_markers.py
#https://github.com/rajewsky-lab/novosparc/blob/master/reconstruct_bdtnp_with_markers.py

# install novosparc
pip install novosparc

# install Altair
# https://altair-viz.github.io/getting_started/installation.html
conda install -c conda-forge altair vega_datasets
##################################################################################

#$pwd
#/Users/mingyang/Downloads/novosparc-master 
#!! important for tissue.setup_reconstruction(markers_to_use=markers_to_use, atlas_matrix=atlas_matrix)  
# and you must specify this path in terminal before entering python, have no idea why

import novosparc
import os
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import altair as alt
from scipy.spatial.distance import cdist, squareform, pdist
from scipy.stats import ks_2samp
from scipy.stats import pearsonr

import random
random.seed(0)

pl_genes = ['sna', 'ken', 'eve']

######################################################################################################
## Reading expression data to scanpy AnnData (cells x genes)
data_dir = '/Users/mingyang/Downloads/novosparc-master/novosparc/datasets/drosophila_scRNAseq/'
data_path = os.path.join(data_dir, 'dge_normalized.txt')
dataset = sc.read(data_path).T
gene_names = dataset.var.index.tolist()

num_cells, num_genes = dataset.shape # 1297 cells x 8924 genes

print('number of cells: %d' % num_cells)
print('number of genes: %d' % num_genes)
#number of cells: 1297
#number of genes: 8924

# Preprocess data
# sc.pp.normalize_per_cell(dataset)
# sc.pp.log1p(dataset)

#(Optional): subset cells
#num_cells = 1000
#sc.pp.subsample(dataset, n_obs=num_cells)
#dataset #AnnData object with n_obs × n_vars = 1000 × 8924


## Create a target space
# Alternative 1: using a reference atlas
atlas_dir = '/Users/mingyang/Downloads/novosparc-master/novosparc/datasets/bdtnp/'
target_space_path = os.path.join(atlas_dir, 'geometry.txt')
locations = pd.read_csv(target_space_path, sep=' ')
locations.shape # (6078, 3), x,y,z, coordiantes for 6078 spatial locations
num_locations = 3039
locations_apriori = locations[:num_locations][['xcoord', 'zcoord']].values


## reading reference atlas (Choose a number of markers to use for reconstruction)
locations = locations_apriori
atlas_dir = '/Users/mingyang/Downloads/novosparc-master/novosparc/datasets/bdtnp/'
atlas_path = os.path.join(atlas_dir, 'dge.txt')
atlas = sc.read(atlas_path)
atlas #AnnData object with n_obs × n_vars = 3039 × 84
# this file contains expr.level of 84 marker genes in each of those 3039 locations
atlas_genes = atlas.var.index.tolist()
len(atlas_genes) #84 genes
#atlas.obsm['spatial'] = locations
#novosparc.pl.embedding(atlas, pl_genes)

# params for linear cost 
markers = sorted(list(set(atlas_genes).intersection(gene_names))) #!!!!sorted!!!!
len(markers) #84 markers
markers
sorted(markers)

atlas.to_df()
atlas.to_df()[markers]
atlas_matrix = atlas.to_df()[markers].values
atlas_matrix.shape #(3039, 84)


markers_idx = pd.DataFrame({'markers_idx': np.arange(num_genes)}, index=gene_names)
markers_to_use = np.concatenate(markers_idx.loc[markers].values)
len(markers_to_use) #84
max(markers_to_use) #8898
markers_to_use
len(gene_names) #8924

## Setup and spatial reconstruction 
tissue = novosparc.cm.Tissue(dataset=dataset, locations=locations)
tissue.setup_reconstruction(num_neighbors_s = 5, num_neighbors_t = 5)
tissue.setup_reconstruction(markers_to_use=markers_to_use, atlas_matrix=atlas_matrix) #


# Compute OT of cells to locations with a given alpha parameter
# compute optimal transport of cells to locations
alpha_linear = 0.8
epsilon = 5e-3
tissue.reconstruct(alpha_linear=alpha_linear, epsilon=epsilon)

# calculate spatially informative genes after reconstruction
#tissue.calculate_spatially_informative_genes()  #takes some time


# save the sdge to file
output_folder='/Users/mingyang/Downloads/'
novosparc.io.write_sdge_to_disk(tissue, output_folder)
# a 'sdge_1297_cells_3039_locations.txt' file is saved
# 8924 row by 3039 columns, 8924: gene, 3039: spatial bins. 



# adjust location marginals gradual cell density from the tissue’s center
#rdist = novosparc.gm.prob_dist_from_center(locations)
#atlas.obs['Alternative location marginals'] = rdist
#novosparc.pl.embedding(atlas, ['Alternative location marginals'])
##tissue.reconstruct(alpha_linear=alpha_linear, epsilon=epsilon, p_locations=rdist) 

#Once computed, the transport matrix is available in tissue’s object field tissue.gw (numpy ndarray of dimensions num_cells x num_locations), 
#and the predicted expression in tissue.sdge (numpy ndarray shaped as num_genes x num_locations).
tissue.gw.shape # (1297, 3039)
tissue.sdge.shape # (8924, 3039)

# Validate predicted expression over target space
# reconstructed expression of individual genes
sdge = tissue.sdge
sdge.shape #(8924, 3039)
len(gene_names) #8924
dataset_reconst = sc.AnnData(pd.DataFrame(sdge.T, columns=gene_names))
dataset_reconst.obsm['spatial'] = locations
dataset_reconst
#AnnData object with n_obs × n_vars = 3039 × 8924
#    obsm: 'spatial'

g=pl_genes[0]
dataset_reconst[:,g].X.flatten().size #3039
atlas[:,g].X.flatten().size #3039
pearsonr(dataset_reconst[:,g].X.flatten(),atlas[:,g].X.flatten() )

title = ['%s, corr=%.02f' % (g, pearsonr(dataset_reconst[:,g].X.flatten(),atlas[:,g].X.flatten() )[0] ) for g in pl_genes]
novosparc.pl.embedding(dataset_reconst, pl_genes, title=title)
# atlas[:,g].X.flatten() is a constant, no vairance, so pearsonr is NaN


#Validate localized mapping of individual cells
# probability of individual cells belonging to each location
gw = tissue.gw #(1297 cell, 3039 bins)
gw.sum(0) #col.sum, 3039 bins, each bin.sum equal to 0.00032906
len(gw.sum(1)) #rowSums, 1297

ngw = (gw.T / gw.sum(1)).T #row-wise normalize
ngw.shape #(1297, 3039)
ngw.sum(1) #rowSum, all equal to 1

import numpy as np
a_file=os.path.join(output_folder, 'tissue.gw_1297gene_3039bin.txt')
np.savetxt(a_file, gw) #1297 x 3039 file

a_file=os.path.join(output_folder, 'tissue.ngw_1297gene_3039bin.txt')
np.savetxt(a_file, ngw) #1297 x 3039 file

a_file=os.path.join(output_folder, 'tissue.locaitons_3039xzcoords.txt')
np.savetxt(a_file, locations) #1297 x 3039 file



cell_idx = [1, 12]
cell_prb_cols = ['cell %d' % i for i in cell_idx]
dataset_reconst.obs = pd.DataFrame(ngw.T[:, cell_idx], columns=cell_prb_cols)
dataset_reconst.obs['cell 1'].sum() # equal to 1
dataset_reconst.obs['cell 12'].sum() # equal to 1

title=['Cell %d, entropy=%.02f' % (i, novosparc.an.get_cell_entropy(ngw[i,:])) for i in cell_idx]
novosparc.pl.embedding(dataset_reconst, cell_prb_cols, title=title)

# comparing distributions of entropy for transporting a cell to locations
ent_T, ent_T_unif, ent_T_rproj, ent_T_shuf = novosparc.pl.plot_transport_entropy_dist(tissue)

print(ks_2samp(ent_T, ent_T_rproj))
print(ks_2samp(ent_T, ent_T_shuf))



# Looking for spatially informative genes according to reconstruction in highly variable genes
cyc_genes = [g for g in gene_names if g.startswith('Cyc')]
atlas_genes = list(atlas.var_names)
mI_genes = cyc_genes + atlas_genes

tissue.calculate_spatially_informative_genes(mI_genes)
genes_with_scores = tissue.spatially_informative_genes

genes_with_scores.index = genes_with_scores['genes']

gene_groups = {'Atlas': atlas_genes, 'Cell-cycle': cyc_genes}
novosparc.pl.plot_morans_dists(genes_with_scores, gene_groups)

gene_max_mI = genes_with_scores['genes'].iloc[0]
gene_min_mI = genes_with_scores['genes'].iloc[-1]

title = ['%s, Morans`I=%.02f' % (gene_max_mI, genes_with_scores.loc[gene_max_mI]['mI']), 
         '%s, Morans`I=%.02f' % (gene_min_mI, genes_with_scores.loc[gene_min_mI]['mI'])]

novosparc.pl.embedding(dataset_reconst, [gene_max_mI, gene_min_mI], title=title)

print('Mean Morans I for cell-cycle genes: %.02f' % genes_with_scores.loc[cyc_genes]['mI'].mean())
print('Mean Morans I for atlas genes: %.02f' % genes_with_scores.loc[atlas_genes]['mI'].mean())


# Extract archetypes
# plot spatial expression archtypes
num_clusters = 10
atlas_indices = pd.DataFrame(np.arange(num_genes), index=gene_names)[0].loc[atlas_genes].values
archetypes, clusters, gene_corrs = novosparc.rc.find_spatial_archetypes(num_clusters, sdge[atlas_indices,:])

arch_cols = ['archetype %d'% i for i in np.arange(num_clusters)]
dataset_reconst.obs = pd.DataFrame(index=dataset_reconst.obs.index)
df = pd.DataFrame(archetypes.T, columns=arch_cols)
dataset_reconst.obs = pd.concat((dataset_reconst.obs, df), 1)

novosparc.pl.embedding(dataset_reconst, arch_cols)





# Numpy array dimensions
tissue.gw.shape # (1000, 3039)
tissue.gw[1,].size #3039 potential assigned locations
tissue.gw[1,].sum() #0.01
tissue.gw.sum(axis=1).size #1000, row.sum
tissue.gw.sum(axis=0).size #3039, col.sum
tissue.gw.sum(axis=1) #all equal to 0.001

