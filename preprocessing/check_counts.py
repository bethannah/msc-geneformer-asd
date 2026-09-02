

import scanpy as sc
import numpy as np
import scipy.sparse as sp

adata = sc.read_h5ad("preprocessing/vpa_organoids.h5ad")

print("Dataset:")
print(adata)

columns_to_check = [
    "condition",
    "disease",
    "cell_type",
    "cell_subtype",
    "donor_id",
]

for column in columns_to_check:
    print(f"\n{column} values:")
    print(adata.obs[column].value_counts(dropna=False))


def inspect_matrix(matrix, name):
    if sp.issparse(matrix):
        values = matrix.data[:100000]
    else:
        values = np.asarray(matrix).ravel()[:100000]

    print(f"\n{name}")
    print("Example values:", values[:30])
    print(
        "Fraction integer-like:",
        np.mean(np.isclose(values, np.round(values)))
    )
    print("Min:", values.min())
    print("Max:", values.max())


inspect_matrix(adata.X, "adata.X")

if adata.raw is not None:
    inspect_matrix(adata.raw.X, "adata.raw.X")