# following: https://github.com/czbiohub/tabula-muris-senis/tree/master/2_aging_signature
# /Users/mingyang/Documents/Data_mouse_aging_atlas/github_tms/tabula-muris-senis-master/2_aging_signature/job.downstream/downstream_tissue_cell.dge_summary.ipynb 


import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
sns.set(style="whitegrid")
import numpy as np
import scanpy as sc
from anndata import read_h5ad
from anndata import AnnData
import scipy as sp
import scipy.stats
from gprofiler import GProfiler
import pickle
from adjustText import adjust_text
from matplotlib import gridspec
# Other specific functions 
from itertools import product
from statsmodels.stats.multitest import multipletests
import time
import os

import sys
sys.path.insert(1, '../')
import util # keep util.py in the same working folder 
#/Users/mingyang/Documents/Data_mouse_aging_atlas/github_tms/tabula-muris-senis-master/2_aging_signature/util.py 

# autoreload
%load_ext autoreload
%autoreload 2
# logging
sc.logging.print_versions()


# GLOBAL VARIABLES
DATA_PATH = './'
DGE_RES_PATH = DATA_PATH + '/DGE_result'
DGE_RES_PATH_OLD = DATA_PATH + '/DE_result_old'
ANNO_DATA_PATH = '/Users/mingyang/Documents/Data_mouse_aging_atlas/TMS.gene.data_final/annotation_data'
RESULT_PATH = DATA_PATH + '/result_v1'

METHOD_LIST = ['facs', 'droplet']
DIC_METHOD_NAME = {'facs':'FACS', 'droplet':'droplet'}
CELLCATE_LIST = ['immune', 'stem cell/progenitor', 'stromal', 'endothelial', 'epithelial', 'parenchymal']

## small test
# Load DGE results
version_list = ['1e4', 'n_gene', 'male_only', 'male_only.3_vs_18', 'female_only'] #check out /Users/mingyang/Documents/Data_mouse_aging_atlas/github_tms/tabula-muris-senis-master/2_aging_signature/job.DGE_analysis/DGE_analysis.R 
dic_dge = {}

version='1e4'
df_info_facs,dic_dge_facs = util.load_DGE_res(DATA_PATH, dname='facs.tc', version=version)
len(dic_dge_facs) # 131
df_info_facs # compare with `tabula-muris-senis-facs-official-raw-obj.h5ad`
# there is at least 1 cell in all three age groups
## small test done


# Load DGE results
version_list = ['1e4', 'n_gene', 'male_only', 'male_only.3_vs_18', 'female_only']
dic_dge = {}

for version in version_list:
    df_info_facs,dic_dge_facs = util.load_DGE_res(DATA_PATH, dname='facs.tc', version=version)
    df_info_droplet,dic_dge_droplet = util.load_DGE_res(DATA_PATH, dname='droplet.tc', version=version)

    # Change analyte name
    temp_list = list(dic_dge_facs.keys())
    for analyte in temp_list:
        tissue,cell_type = analyte.split('.')
        cell_type = cell_type.replace('_', ' ')
        dic_dge_facs['%s.%s'%(tissue,cell_type)] = dic_dge_facs[analyte].copy()
        if '%s.%s'%(tissue,cell_type) != analyte: del dic_dge_facs[analyte]

    temp_list = list(dic_dge_droplet.keys())
    for analyte in temp_list:
        tissue,cell_type = analyte.split('.')
        cell_type = cell_type.replace('_', ' ')
        dic_dge_droplet['%s.%s'%(tissue,cell_type)] = dic_dge_droplet[analyte].copy()
        if '%s.%s'%(tissue,cell_type) != analyte:  del dic_dge_droplet[analyte]
        
    # fixit: update bh_p (not sure if this is necessary)
    if version=='1e4':
        dic_dge['facs'] = dic_dge_facs.copy()
        dic_dge['droplet'] = dic_dge_droplet.copy()
    else:
        dic_dge['facs.%s'%version] = dic_dge_facs.copy()
        dic_dge['droplet.%s'%version] = dic_dge_droplet.copy()
        
# Append tissue-level results
df_info_facs_tissue,dic_dge['facs.tissue'] = util.load_DGE_res(DATA_PATH, dname='facs.tissue', version='1e4')
df_info_droplet_tissue,dic_dge['droplet.tissue'] = util.load_DGE_res(DATA_PATH, dname='droplet.tissue',
                                                                     version='1e4')
df_info_bulk_tissue,dic_dge['bulk.tissue'] = util.load_DGE_res(DATA_PATH, dname='bulk.tissue', version='1e4')


df_info_droplet,dic_dge_droplet = util.load_DGE_res(DATA_PATH, dname='droplet.tc', version=version)

# dic_analysis_list and dic_fdr_threshold

# analysis list: facs
min_cell_number = 100
ind_select = (df_info_facs['n_cell_young']>min_cell_number) & (df_info_facs['n_cell_old']>min_cell_number)
analysis_list_facs = list(df_info_facs.index[ind_select])

# analysis list: droplet
min_cell_number = 500
ind_select = (df_info_droplet['n_cell_young']>min_cell_number) & (df_info_droplet['n_cell_old']>min_cell_number)
analysis_list_droplet = list(df_info_droplet.index[ind_select])

dic_analysis_list = {'facs':analysis_list_facs, 'droplet':analysis_list_droplet}
for method in METHOD_LIST:
    print('%s, n_tc=%d'%(method, len(dic_analysis_list[method])))

# thresholds parameters
coef_threshold = 0.005
dic_fdr_threshold = {'facs':0.01, 'droplet':0.01, 'bulk':0.1}

#facs, n_tc=76
#droplet, n_tc=26

## check out `load_DGE_res` funciton in `util.py` script.





