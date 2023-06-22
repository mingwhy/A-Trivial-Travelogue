#https://shendure-web.gs.washington.edu/content/members/DEAP_website/public/

# lineage information
For RNA-based lineage information used to generate the tree in Fig. 3C, [this] is the same RNA annotation table as in the supplementary information, but with an added column for lineage. To help decode the lineage annotations, refer to [this additional file] that translates each lineage into the corresponding germ layer annotation.
table_s1_withLineage.xlsx
celltype_list_rna.txt

Supplemental tables
All supplementary tables are aggregated in the following excel file with one table per spreadsheet: **SupplementaryTables.xlsx.**

Additional processed files
**RNA seurat objects (link). The main seurat object includes all cells that passed QC filters, and includes the time prediction associated with each cell. 
NNv1 is the best performing neural net model with a MSE loss (used for downstream analyses),!!!
NNv2 is the best performing neural net model with the custom loss, and 
lasso refers to the lasso model. 
The [model]_age is the model age precition whereas the [model]_shift refers to the error outside the cells collection window. The time split objects are split into cells divided by collection window (windows), or predicted time windows (pred_windows).**

https://shendure-web.gs.washington.edu/content/members/DEAP_website/public/RNA/update/seurat_objects/
[   ]	14912_29824_subsampled_large_checkpoint_endpoint.Rds	2022-05-19 16:51	706M	 
[   ]	main.Rds	2022-03-24 11:58	1.4G	 
[DIR]	pred_windows/	2022-03-24 11:52	-	 
[DIR]	windows/	2022-03-24 12:02	-	 

I save these files in folder: SupplementaryOnlineMaterial/seurat_objects/


########################################################################
I also download three rds file from:
https://shendure-web.gs.washington.edu/content/members/DEAP_website/public/RNA/tf_correlation_plots/
in RNA_seurat_object/tf_correlation_plots folder



