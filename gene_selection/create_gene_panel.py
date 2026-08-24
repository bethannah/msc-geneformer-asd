import pandas as pd

#files
asd_file = "selected_asd_50.csv"
non_asd_file = "matched_genes.csv"

#load files
asd = pd.read_csv(asd_file)
non_asd = pd.read_csv(non_asd_file)

#add panel labels
asd["panel"] = "ASD"
non_asd["panel"] = "non_ASD"

#add blank SFARI metadata for non-ASD genes
non_asd["sfari_score"] = pd.NA
non_asd["syndromic"] = pd.NA

#Make Young dExN column names consistent with ASD file
non_asd = non_asd.rename(columns={
    "CTRL_Young_detected_cells": "CTRL_Young_dExN_detected_cells",
    "CTRL_Young_detection_percentage": "CTRL_Young_dExN_detection_percentage"
})

#combine both groups
combined = pd.concat([asd, non_asd], ignore_index=True)

##Geneformer input file
gene_panel = combined[
    ["panel", "gene", "ensembl_id"]
]

gene_panel.to_csv("gene_panel.csv", index=False)

##Metadata file
metadata = combined[
    ["gene",
     "ensembl_id",
     "panel",
     "sfari_score",
     "syndromic",
     "CTRL_RG_detected_cells",
     "CTRL_RG_detection_Percentage",
     "CTRL_IPC_detected_cells",
     "CTRL_IPC_detection_percentage",
     "CTRL_Young_detected_Cells",
     "CTRL_Young_detection_percentage"]
]

metadata.to_csv("gene_metadata.csv", index=False)

print(f"Genes in panel: {len(gene_panel)}")