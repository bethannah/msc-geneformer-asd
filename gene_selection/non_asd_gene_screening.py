#!/usr/bin/env python
# coding: utf-8

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

#files
h5ad_file = "tokenizer_input/vpa_organoids_geneformer_ready.h5ad"
asd_file = "SFARI_Gene_List.csv"
output_file = "non_asd_candidates.csv"
vpa_file = "vpa_reference_genes.csv"

#Load data
adata = ad.read_h5ad(h5ad_file)

#Calculate gene detection for a CTRL cell subtype
def calculate_detection(subtype):
    cells = adata[
        (adata.obs["condition"] == "CTRL") &
        (adata.obs["cell_subtype"] == subtype)]

    if sparse.issparse(cells.X):
        detected_cells = cells.X.getnnz(axis=0)
    else:
        detected_cells = (cells.X > 0).sum(axis=0)

    detection_percentage = (detected_cells / cells.n_obs) * 100

    return(np.asarray(detected_cells).flatten(),
           np.asarray(detection_percentage).flatten())

#Calculate detection in eqach CTRL cell subtype
rg_count, rg_percentage = calculate_detection("RG")
ipc_count, ipc_percentage = calculate_detection("IPC/Newborn ExN")
young_count, young_percentage = calculate_detection("Young dExN")

#Create dataframe containing every gene
genes = pd.DataFrame({
    "gene": adata.var["feature_name"].values,
    "ensembl_id": adata.var["ensembl_id"].values,

    "CTRL_RG_detected_cells": rg_count,
    "CTRL_RG_detection_percentage": rg_percentage,

    "CTRL_IPC_detected_cells": ipc_count,
    "CTRL_IPC_detection_percentage": ipc_percentage,

    "CTRL_Young_detected_cells": young_count,
    "CTRL_Young_detection_percentage": young_percentage
})

#Keep genes detected in at least 50% of CTRL RG cells
genes = genes[
    genes["CTRL_RG_detection_percentage"] >= 50
].copy()

#Load SFARI gene list & remove all SFARI-listed genes
asd = pd.read_csv(asd_file)

genes = genes[
    ~genes["ensembl_id"].isin(asd["ensembl-id"])].copy()

#Load VPA-reference genes and remove
vpa = pd.read_csv(vpa_file)

genes = genes[
    ~genes["ensembl_id"].isin(vpa["ensembl_id"])].copy()

#Remove genes without Ensembl ID
genes = genes.dropna(subset=["ensembl_id"])

#Save
genes.to_csv(output_file, index=False)
print(f"Non-ASD candidate genes: {len(genes)}")
