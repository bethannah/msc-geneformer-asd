
import scanpy as sc
import numpy as np
import scipy.sparse as sp

adata = sc.read_h5ad("preprocessing/vpa_organoids_geneformer_ready.h5ad")

print(adata)
print(adata.obs.columns.tolist())
print(adata.var.columns.tolist())
print(adata.X.data[:20])
print(np.mean(np.isclose(adata.X.data[:100000], np.round(adata.X.data[:100000]))))
print(adata.var[["feature_name", "ensembl_id"]].head(10))
print(adata.obs["cell_subtype"].value_counts())
