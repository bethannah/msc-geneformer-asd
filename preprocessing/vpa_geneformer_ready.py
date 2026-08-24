
from pathlib import Path

import numpy as np
import scanpy as sc
import scipy.sparse as sp


input_file = Path("vpa_organoids.h5ad")
output_file = Path("vpa_organoids_geneformer_ready.h5ad")

# Read original file
adata = sc.read_h5ad(input_file)

print("Original object:")
print(adata)

# Use raw counts stored in adata.raw
if adata.raw is None:
    raise ValueError("adata.raw is missing, so raw counts cannot be found.")

raw_adata = adata.raw.to_adata()

# Keep the cell metadata from the original object
raw_adata.obs = adata.obs.copy()

# Make sure the raw count matrix is stored
if not sp.issparse(raw_adata.X):
    raw_adata.X = sp.csr_matrix(raw_adata.X)
else:
    raw_adata.X = raw_adata.X.tocsr()

# Add total counts per cell for Geneformer
raw_adata.obs["n_counts"] = np.asarray(
    raw_adata.X.sum(axis=1)
).ravel()

# Add Ensembl IDs for Geneformer
raw_adata.var["ensembl_id"] = (
    raw_adata.var_names.astype(str)
    .str.split(".", regex=False)
    .str[0]
)

# Remove cells with zero total counts
raw_adata = raw_adata[
    raw_adata.obs["n_counts"] > 0
].copy()

print("\nPrepared object:")
print(raw_adata)

print("\nCondition values:")
print(raw_adata.obs["condition"].value_counts(dropna=False))

print("\nCell subtype values:")
print(raw_adata.obs["cell_subtype"].value_counts(dropna=False))

print("\nGene metadata:")
print(raw_adata.var.head())

print("\nExample raw count values:")
print(raw_adata.X.data[:30])

print(
    "\nFraction integer-like:",
    np.mean(
        np.isclose(
            raw_adata.X.data[:100000],
            np.round(raw_adata.X.data[:100000])
        )
    )
)

raw_adata.write_h5ad(output_file, compression="gzip")

print(f"\nSaved Geneformer-ready file as: {output_file}")