#!/usr/bin/env python
# coding: utf-8

# In[ ]:


from geneformer import InSilicoPerturber
from geneformer import InSilicoPerturberStats
#from geneformer import EmbExtractor
from pathlib import Path
import pandas as pd
import pickle

# ### Single gene deletion analysis to test whether ASD- VPA- and control-panel genes shift CTRL radial glia towards the VPA state

# In[ ]:


# first obtain start, goal, and alt embedding positions
# this function was changed to be separate from perturb_data
# to avoid repeating calcuations when parallelizing perturb_data
cell_states_to_model={"state_key": "condition", 
                      "start_state": "CTRL", 
                      "goal_state": "VPA", 
                      "alt_states": []}

filter_data_dict={
		"cell_subtype": ["RG"]}

#embex = EmbExtractor(model_type="Pretrained", # if using previously fine-tuned cell classifier model
                     #num_classes=0,
                     #filter_data=filter_data_dict,
                     #max_ncells=1000,
                     #emb_layer=0,
                     #summary_stat="exact_mean",
                     #forward_batch_size=100,
                     #model_version="V1",  # OF NOTE: SET TO V1 MODEL, PROVIDE V1 MODEL PATH IN SUBSEQUENT CODE
                    # nproc=1)

#state_embs_dict = embex.get_state_embs(cell_states_to_model,
                                      # "/d/projects/u/software/geneformer/content/Geneformer/Geneformer-V1-10M", # example 30M fine-tuned model
                                      # "tokenised_data/vpa_organoids.dataset",
                                      # "results/state_embeddings",
                                      # "RG_CTRL_to_VPA_state_embs")


with open("results/state_embeddings/RG_CTRL_to_VPA_state_embs.pkl","rb") as file:
	state_embs_dict = pickle.load(file)

# In[ ]:

#Paths used for every gene 
model_directory = "/d/projects/u/software/geneformer/content/Geneformer/Geneformer-V1-10M"

input_dataset = "tokenised_data/vpa_organoids.dataset"

#Read the gene panel
gene_panel = pd.read_csv(
	"gene_panel.csv",
	dtype=str,)

#Result directories
perturbation_directory = Path("results/batch_perturbations")

stats_directory = Path("stats/batch_perturbations")

failed_genes = []

#Run each CSV row as a separate single-gene perturbation
for _, row in gene_panel.iterrows():
	panel = row["panel"]
	gene_name = row["gene"]
	ensembl_id = row["ensembl_id"]

	output_prefix = f"{panel}_{gene_name}_RG"

	gene_perturbation_directory = perturbation_directory / output_prefix
	gene_stats_directory = stats_directory / output_prefix

	for directory in (gene_perturbation_directory, gene_stats_directory):
		directory.mkdir(parents=True, exist_ok=True)

	print(f"Running {panel}: {gene_name} ({ensembl_id})")
	try:
		isp =  InSilicoPerturber(
			perturb_type="delete",
			perturb_rank_shift=None,
			genes_to_perturb=[ensembl_id],
			combos=0,
			anchor_gene=None,
			model_type="Pretrained",
			num_classes=0,
			emb_mode="cell",
			cell_emb_style="mean_pool",
			filter_data=filter_data_dict,
			cell_states_to_model=cell_states_to_model,
			state_embs_dict=state_embs_dict,
			max_ncells=1000,
			emb_layer=0,
			forward_batch_size=100,
			model_version="V1",
			nproc=1)

		isp.perturb_data(
			model_directory,
			input_dataset,
			str(gene_perturbation_directory),
			output_prefix)

		ispstats = InSilicoPerturberStats(
			mode="goal_state_shift",
			genes_perturbed=[ensembl_id],
			combos=0,
			anchor_gene=None,
			cell_states_to_model=cell_states_to_model,
			model_version="V1")


		ispstats.get_stats(
			str(gene_perturbation_directory),
			None,
			str(gene_stats_directory),
			output_prefix)
		print(f"Completed {gene_name}, flush=True")

	except Exception as error:
		print(f"Failed {gene_name}: {error}", flush=True,)
		failed_genes.append(
			{"panel": panel,
			"gene": gene_name,
			"ensembl_id": ensembl_id,
			"cell_subtype": cell_subtype,
			"error": str(error),
			}
		)
		continue
	
failed_genes_file = stats_directory / "failed_genes_RG.csv"

pd.DataFrame(failed_genes,
		columns=["panel",
			"gene",
			"ensembl_id",
			"cell_subtype",
			"error"],
		).to_csv(failed_genes_file,index=False)
if failed_genes:
	print(f"{len(failed_genes)} gene(s) failed")
else:
	print("All genes completed")

