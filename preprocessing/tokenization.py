#!/usr/bin/env python
# coding: utf-8

# ## Tokenizing .loom or .h5ad single cell RNA-seq data to rank value encoding .dataset format

# #### Input data is a directory with .loom or .h5ad files containing raw counts from single cell RNAseq data, including all genes detected in the transcriptome without feature selection. The input file type is specified by the argument file_format in the tokenize_data function.
# 
# #### The discussion below references the .loom file format, but the analagous labels are required for .h5ad files, just that they will be column instead of row attributes and vice versa due to the transposed format of the two file types.
# 
# #### Genes should be labeled with Ensembl IDs (loom row attribute "ensembl_id"), which provide a unique identifer for conversion to tokens. Other forms of gene annotations (e.g. gene names) can be converted to Ensembl IDs via Ensembl Biomart. Cells should be labeled with the total read count in the cell (loom column attribute "n_counts") to be used for normalization.
# 
# #### No cell metadata is required, but custom cell attributes may be passed onto the tokenized dataset by providing a dictionary of custom attributes to be added, which is formatted as loom_col_attr_name : desired_dataset_col_attr_name. For example, if the original .loom dataset has column attributes "cell_type" and "organ_major" and one would like to retain these attributes as labels in the tokenized dataset with the new names "cell_type" and "organ", respectively, the following custom attribute dictionary should be provided: {"cell_type": "cell_type", "organ_major": "organ"}. 
# 
# #### Additionally, if the original .loom file contains a cell column attribute called "filter_pass", this column will be used as a binary indicator of whether to include these cells in the tokenized data. All cells with "1" in this attribute will be tokenized, whereas the others will be excluded. One may use this column to indicate QC filtering or other criteria for selection for inclusion in the final tokenized dataset.
# 
# #### If one's data is in other formats besides .loom or .h5ad, one can use the relevant tools (such as Anndata tools) to convert the file to a .loom or .h5ad format prior to running the transcriptome tokenizer.

# **********************************************************************************************************
# #### OF NOTE: Please ensure the correct token dictionary, gene median file, special token setting, and model input size is used for the correct model version.
# #### Current defaults are for V2 model series. To auto-select the correct settings for V1, set model_version argument to "V1".

# In[ ]:


from geneformer import TranscriptomeTokenizer


# In[ ]:


tk = TranscriptomeTokenizer(
			{"condition": "condition",
			"cell_type": "cell_type",
			"cell_subtype": "cell_subtype",
			"cell_class": "cell_class",
			"donor_id": "donor_id",
			"batch": "batch",
			"sex": "sex",
			"development_stage": "development_stage",
			 "disease": "disease"},
			 nproc=8,
			model_version="V1")  # for V1 model, set model_version="V1"

tk.tokenize_data("tokenizer_input", 
                 "tokenised_data", 
                 "vpa_organoids", 
                 file_format="h5ad")

