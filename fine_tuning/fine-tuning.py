#!/bin/env python3

import datetime
import os

current_date = datetime.datetime.now()
datestamp = f"{str(current_date.year)[-2:]}{current_date.month:02d}{current_date.day:02d}{current_date.hour:02d}{current_date.minute:02d}{current_date.second:02d}"
datestamp_min = f"{str(current_date.year)[-2:]}{current_date.month:02d}{current_date.day:02d}"

output_prefix = "cm_dz_classifier"
output_dir = f"./content/finetuned_models/{datestamp}"
print(datestamp)
os.mkdir(output_dir)

from geneformer import Classifier

training_args = {
    "num_train_epochs": 0.9,
    "learning_rate": 0.000804,
    "lr_scheduler_type": "polynomial",
    "warmup_steps": 1812,
    "weight_decay":0.258828,
    "per_device_train_batch_size": 12,
    "seed": 73,
}
cc = Classifier(classifier="cell",
                cell_state_dict={"state_key": "disease", "states": "all"},
                filter_data={"cell_type": ["Cardiomyocyte"]},
                training_args=training_args,
                max_ncells=None,
                freeze_layers=2,
                num_crossval_splits=1,
                forward_batch_size=32,
                nproc=8,
                model_version="V1" # default is V2, here set to V1 model to fit into Colab 40G GPU resources
)
train_ids = ["1447", "1600", "1462", "1558", "1300", "1508", "1358", "1678", "1561", "1304", "1610", "1430", "1472", "1707", "1726", "1504", "1425", "1617", "1631", "1735", "1582", "1722", "1622", "1630", "1290", "1479", "1371", "1549", "1515"]
eval_ids = ["1422", "1510", "1539", "1606", "1702"]
test_ids = ["1437", "1516", "1602", "1685", "1718"]

train_test_id_split_dict = {"attr_key": "individual",
                            "train": train_ids+eval_ids,
                            "test": test_ids}

train_valid_id_split_dict = {"attr_key": "individual",
                            "train": train_ids,
                            "eval": eval_ids}

cc.prepare_data(input_data_file="./content/input_data/cm_tokenized.dataset",
                output_directory=output_dir,
                output_prefix=output_prefix,
                split_id_dict=train_test_id_split_dict)

os.environ["WANDB_DISABLED"] = "true"

all_metrics = cc.validate(
    model_directory=os.path.abspath("./content/Geneformer/Geneformer-V1-10M/"), # set to V1 model to fit into Colab 40G GPU resources
    prepared_input_data_file=f"{output_dir}/{output_prefix}_labeled_train.dataset",
    id_class_dict_file=f"{output_dir}/{output_prefix}_id_class_dict.pkl",
          output_directory=os.path.abspath(output_dir),
            output_prefix=output_prefix,
            attr_to_split="individual",
    attr_to_balance=["disease","age","sex"],
    n_hyperopt_trials=0 	# Number of trials to run for hyperparameter optimization. Set it to 0 for direct training without hyperparameter optimization.
)

all_metrics_test = cc.evaluate_saved_model(
model_directory=f"{output_dir}/{datestamp_min}_geneformer_cellClassifier_{output_prefix}/ksplit1/",
    	id_class_dict_file=f"{output_dir}/{output_prefix}_id_class_dict.pkl",
     	test_data_file=f"{output_dir}/{output_prefix}_labeled_test.dataset",
     	output_directory=output_dir,
      	output_prefix=output_prefix,
    )
cc.plot_conf_mat(
        conf_mat_dict={"Geneformer": all_metrics_test["conf_matrix"]},
        output_directory=output_dir,
        output_prefix=output_prefix,
        custom_class_order=["nf","hcm","dcm"],	# Disease categories in my data
)

cc.plot_predictions(
    predictions_file=f"{output_dir}/{output_prefix}_pred_dict.pkl",
    id_class_dict_file=f"{output_dir}/{output_prefix}_id_class_dict.pkl",
    title="disease",
    output_directory=output_dir,
    output_prefix=output_prefix,
    custom_class_order=["nf","hcm","dcm"],
)

from geneformer import EmbExtractor

embex = EmbExtractor(
    model_type="CellClassifier", # set to GeneClassifier or Pretrained for those model types
    num_classes=3, # number of classes of fine-tuned model
    filter_data={"cell_type":["Cardiomyocyte"]}, # optionally can extract embeddings from a subset of cells based on the chosen attributes
    max_ncells=10000,
    emb_layer=0, # extracts embeddings from last layer
    emb_label=["disease", "individual"],
    labels_to_plot=["disease", "individual"],
    forward_batch_size=128,
    nproc=8,
    model_version="V1" # default is V2, here set to V1 model to fit into Colab 40G GPU resources
)

embs = embex.extract_embs(f"{output_dir}/{datestamp_min}_geneformer_cellClassifier_{output_prefix}/ksplit1/",
                          "./content/input_data/cm_tokenized.dataset",
                          "./content/embs",
                          "cm_finetuned_embs")
embex.plot_embs(embs=embs,
                plot_style="umap",
                output_directory="./content/embs",
                output_prefix="emb_umap",
                max_ncells_to_plot=10000)
embex.plot_embs(embs=embs,
                plot_style="heatmap",
                output_directory="./content/embs",
                output_prefix="emb_heatmap",
                max_ncells_to_plot=10000)
		
from geneformer import EmbExtractor, InSilicoPerturber

# Perturb genes in cells with "dcm" (dilated cardiomyopathy) state. Quantify the results by the embedding shifts relative to the goal state "nf" (non-failing) and the alternative state "hcm" (hypertrophic cardiomyopathy).
cell_states_to_model = {
    "state_key": "disease",
    "start_state": "dcm",
    "goal_state": "nf",
    "alt_states": ["hcm"] # set to empty list if no alternate states
}
filter_data_dict={"cell_type":["Cardiomyocyte"]} # Optionally extract embeddings and perform perturbations from a subset of the data based on the chosen attribute
embex = EmbExtractor(model_type="CellClassifier",	# type of fine-tuned model
                     num_classes=3,	# number of classes of fine-tuned model
                     filter_data=filter_data_dict,
                     max_ncells=1000,
		     emb_mode="cell",	# set to "cls" for the V2 model
                     emb_layer=0, # use the last layer
                     summary_stat="exact_mean",
                     forward_batch_size=256,
                     nproc=8,
                     model_version="V1" # default is V2, here set to V1 model to fit into Colab 40G GPU resources
                     )
state_embs_dict = embex.get_state_embs(cell_states_to_model,
                    model_directory=f"{output_dir}/{datestamp_min}_geneformer_cellClassifier_{output_prefix}/ksplit1/",
                    input_data_file="./content/input_data/cm_tokenized.dataset",
                    output_directory="./content/isp",
                    output_prefix=output_prefix
                    )
perturbation_type = "delete" # or "overexpress"
genes_to_perturb = "all"

isp = InSilicoPerturber(perturb_type=perturbation_type,
                        genes_to_perturb=genes_to_perturb,
                        model_type="CellClassifier", # set it to "MTLCellClassifier-Quantized" for quantized MTL model, for example.
                        num_classes=3, # number of classes that the classifier was trained to distinguish
                        emb_mode="cell",# set to "cls" for the V2 model
                        filter_data=filter_data_dict,
                        cell_states_to_model=cell_states_to_model,
                        state_embs_dict=state_embs_dict,
                        max_ncells=100, # setting to very low number for the purposes of demonstrating code
                        emb_layer=0,
                        forward_batch_size=256,
                        nproc=8,
			model_version="V1" # default is V2, here set to V1 model to fit into Colab 40G GPU resources
                        )
isp.perturb_data(model_directory = f"{output_dir}/{datestamp_min}_geneformer_cellClassifier_{output_prefix}/ksplit1/",
                 input_data_file="./content/input_data/cm_tokenized.dataset",
                 output_directory="./content/isp",
                 output_prefix=output_prefix)

from geneformer import InSilicoPerturberStats

cell_states_to_model = {
    "state_key": "disease",
    "start_state": "dcm",
    "goal_state": "nf",
    "alt_states": ["hcm"]
}

ispstats = InSilicoPerturberStats(mode="goal_state_shift",
    genes_perturbed="all",
    combos=0,
    anchor_gene=None,
    cell_states_to_model=cell_states_to_model,
    model_version="V1"  # default is V2, here set to V1 model to fit into Colab 40G GPU resources
)
ispstats.get_stats("./content/isp",
                   None,
                   "./content/results",
                   "isd_cm")

import pandas as pd
results = pd.read_csv("./content/results/isd_cm.csv",index_col=0)
print(results)

import gdown
file_id = "1H8x-iBfGKYitayyfmz3x91cztxrJbVuK"
output = "./content/results/isd_cm_all_cells.csv"
gdown.download(f'https://drive.google.com/uc?id={file_id}', output, quiet=False)
import pandas as pd
full_results = pd.read_csv("./content/results/isd_cm_all_cells.csv",index_col=0)
print(full_results)

import matplotlib.pyplot as plt
import seaborn as sns
sns.set_theme(style='white')

# Filter for genes with significant shifts in both directions
df = full_results[(full_results['Goal_end_FDR'] < 0.05) & (full_results['Alt_end_FDR_hcm'] < 0.05)]

# Further filter for genes shifting towards one state and away from the other
df = df[((df['Shift_to_goal_end'] > 0) & (df['Shift_to_alt_end_hcm'] < 0)) |
        ((df['Shift_to_goal_end'] < 0) & (df['Shift_to_alt_end_hcm'] > 0))]

# Filter for number of cells with that gene expressed
df = df[df["N_Detections"] > 100]
n_top_genes = 15
top_genes = df.loc[df[["Shift_to_goal_end", "Shift_to_alt_end_hcm"]].abs().sum(axis=1).sort_values(ascending=False).index]["Gene_name"].tolist()[:n_top_genes]

# List of genes you want to annotate
genes_to_annotate = top_genes + ["GSN", "ESRRG", "HMGB1"]

# Scatter plot of shift to goal state (NF) vs shift to alternate state (HCM)
plt.figure(figsize=(6, 6))
sns.scatterplot(data=df,
                x='Shift_to_goal_end',
                y='Shift_to_alt_end_hcm',
                color='red',  # Color for significant points with desired shifts
                size='N_Detections',  # Size based on the number of detections
                sizes=(20, 200),  # Adjust marker sizes
                legend='brief')
# Add titles and labels
plt.title('Candidate therapeutic targets whose in silico repression \n in DCM cardiomyocytes shifts towards NF, away from HCM', fontsize=16)
plt.xlabel('Shift to goal state (NF)', fontsize=14)
plt.ylabel('Shift to alternate state (HCM)', fontsize=14)

# Add a reference line where shifts to NF and HCM are equal
plt.axline((0, 0), slope=1, color='gray', linestyle='--')

# Add red dotted lines to mark zero positions on both axes
plt.axhline(0, color='red', linestyle=':', linewidth=1.5)
plt.axvline(0, color='red', linestyle=':', linewidth=1.5)
# # Annotate selected genes
for _, row in df.iterrows():
    if row['Gene_name'] in genes_to_annotate:
        plt.text(row['Shift_to_goal_end'], row['Shift_to_alt_end_hcm'], row['Gene_name'],
                 fontsize=10, color='black', weight='bold')

# Show plot
plt.tight_layout()
plt.show()
