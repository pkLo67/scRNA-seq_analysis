#!/usr/bin/env python
# coding: utf-8

# In[1]:


import scanpy as sc
import numpy as np
sc.settings.verbosity = 3
sc.settings.set_figure_params(dpi=80)
sc.logging.print_version_and_date()


# In[2]:


sc.logging.print_versions()


# In[3]:


import pandas as pd
path = '/Users/Pang-Kuo/Desktop/10XGenomics/'
features = pd.read_csv(path + 'features.tsv', sep='\t')
features.head(10)


# In[4]:


barcodes = pd.read_csv(path + 'barcodes.tsv', sep='\t')
barcodes.head(10)


# In[5]:


adata = sc.read('matrix.mtx', cache=True).T
adata.var_names = pd.read_csv('features.tsv', header=None, sep='\t')[1]
adata.obs_names = pd.read_csv('barcodes.tsv', header=None)[0]
adata.var_names_make_unique()
adata.to_df()


# In[6]:


adata.shape


# In[7]:


sc.pp.filter_cells(adata, min_genes=1000)
sc.pp.filter_genes(adata, min_cells=6)


# In[8]:


adata.shape


# In[10]:


np.sum(list(filter(lambda x: x > 0, adata.X[0].toarray().tolist()[0])))
adata.X.toarray()


# In[11]:


np.sum(adata.X, axis=1).A1


# In[12]:


mito_genes = [name for name in adata.var_names if name.startswith('MT-')]
# for each cell compute fraction of counts in mito genes vs. all genes
# the '.A1' is only necessary, as X is sparse - it transform to a dense array after summing
adata.obs['percent_mito'] = np.sum(adata[:, mito_genes].X, axis=1).A1 / np.sum(adata.X, axis=1).A1
# add the total counts per cell as observations-annotation to adata
adata.obs['n_counts'] = adata.X.sum(axis=1).A1
sc.pl.violin(adata, ['n_genes', 'n_counts', 'percent_mito'],
                    jitter=0.4, multi_panel=True)


# In[13]:


adata = adata[adata.obs['n_genes'] < 4000, :]
adata = adata[adata.obs['percent_mito'] < 0.15, :]
sc.pl.violin(adata, ['n_genes', 'n_counts', 'percent_mito'],
                    jitter=0.4, multi_panel=True)


# In[14]:


adata.raw = sc.pp.log1p(adata, copy=True)


# In[15]:


sc.pp.normalize_per_cell(adata, counts_per_cell_after=1e4)
filter_result = sc.pp.filter_genes_dispersion(
     adata.X, min_mean=0.0125, max_mean=3, min_disp=0.5)
sc.pl.filter_genes_dispersion(filter_result)


# In[16]:


sum(filter_result.gene_subset) / len(filter_result.gene_subset)*100


# In[17]:


adata = adata[:, filter_result.gene_subset]
adata.shape[1]


# In[18]:


sc.pp.log1p(adata)


# In[22]:


sc.pp.regress_out(adata, ['n_counts', 'percent_mito'])
sc.pp.scale(adata, max_value=10)
adata.write('results_file_1.h5ad')


# In[23]:


sc.tl.pca(adata)
adata.obsm['X_pca'] *= -1     # multiply by -1 to match Seurat
sc.pl.pca_scatter(adata, color='CST3')


# In[24]:


sc.pl.pca_variance_ratio(adata, log=True)


# In[25]:


adata.write('results_file_2.h5ad')
adata


# In[28]:


adata = sc.read_h5ad('results_file.h5ad')
sc.tl.tsne(adata, random_state=2, n_pcs=10)
adata.write('results_file_3.h5ad')


# In[29]:


sc.pl.tsne(adata, color=['CST3', 'NKG7', 'PPBP'])      # plot tSNE of raw data


# In[31]:


sc.pl.tsne(adata, color=['CST3', 'NKG7', 'PPBP'], use_raw=False)     # plot tSNE of filtered data


# In[32]:


adata = sc.read_h5ad('results_file.h5ad')
sc.pp.neighbors(adata, n_neighbors=10)


# In[33]:


sc.tl.umap(adata)
sc.pl.umap(adata, color=['CST3', 'NKG7', 'PPBP'])    #plot UMAP of raw data


# In[34]:


sc.pl.umap(adata, color=['CST3', 'NKG7', 'PPBP'], use_raw=False)     #plot UMAP of filtered data


# In[ ]:




