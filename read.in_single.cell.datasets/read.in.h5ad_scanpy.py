import scanpy 

x=scanpy.read_h5ad('ad_worm_aging.h5ad') #https://scanpy.readthedocs.io/en/stable/generated/scanpy.read_h5ad.html
x
#AnnData object with n_obs × n_vars = 47423 × 20305

x[1:3,1:3].to_df() #https://github.com/scverse/scanpy/issues/262
x[10:30,10:30].to_df()

#https://anndata.readthedocs.io/en/latest/generated/anndata.AnnData.layers.html
x.layers['denoised']

x.layers['denoised'].shape
#Out[9]: (47423, 20305)

#Delete layer: https://anndata.readthedocs.io/en/latest/generated/anndata.AnnData.layers.html
del x.layers

#https://anndata.readthedocs.io/en/latest/generated/anndata.AnnData.write.html

x.write_h5ad("ad_worm_aging_umi.h5ad")

