#!/usr/bin/env python
# coding: utf-8

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

#files
h5ad_file = "tokenizer_input/vpa_organoids_geneformer_ready.h5ad"
candidate_file = "asd_candidate_pool.csv"
output_file = "asd_candidates_detection.csv"

#Load data
adata = ad.read_h5ad(h5ad_file)
candidates = pd.read_csv(candidate_file, dtype=str)

dataset_ensembl_ids = adata.var["ensembl_id"].to_numpy()

#Cell subtypes to screen
cell_subtypes = {
    "RG": "RG",
    "IPC": "IPC/Newborn ExN",
    "Young_dExN": "Young dExN",
}

#screen each subtype
results = candidates.copy()

for label, subtype in cell_subtypes.items():
    mask = (
        (adata.obs["condition"] == "CTRL")
        & (adata.obs["cell_subtype"] == subtype)
    )
    cells = adata[mask]

    detected_counts = []
    detection_percentage = []

    for ensembl_id in candidates["ensembl_id"]:
        gene_position = np.where(dataset_ensembl_ids == ensembl_id)[0]
        if len(gene_position) == 0:
            detected = 0
        else:
            expression = cells[:, gene_position].X

            if sparse.issparse(expression):
                detected = int(
                    (expression.getnnz(axis=1) > 0).sum())
            else:
                detected = int(
                    (np.asarray(expression) > 0)
                    .any(axis=1)
                    .sum())

        detected_counts.append(detected)

        detection_percentage.append(
            100 * detected / cells.n_obs
            if cells.n_obs > 0
            else 0)

    results[f"CTRL_{label}_detected_cells"] = detected_counts
    results[f"CTRL_{label}_detection_percentage"] = detection_percentage

#save results
results.to_csv(output_file, index=False)

print(f"\nSaved to {output_file}")