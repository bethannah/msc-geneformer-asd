from pathlib import Path
import pandas as pd

project_dir = Path(__file__).resolve().parents[1]

asd_file = project_dir / "gene_selection" / "selected_asd_50.csv"
hpo_file = project_dir / "reference_data" / "phenotype_to_genes.txt"

asd = pd.read_csv(asd_file)
hpo = pd.read_csv(hpo_file, sep="\t", dtype=str)

#clean gene symbols
asd["gene"] = asd["gene"].str.strip()
hpo["gene_symbol"] = hpo["gene_symbol"].str.strip()

#check which ASD genes occur in HPO
hpo_genes = set(hpo["gene_symbol"])

asd["hpo_annotated"] = asd["gene"].isin(hpo_genes)

print(f"ASD Genes: {len(asd)}")
print(f"HPO annotated: {asd['hpo_annotated'].sum()}")
print(f"No HPO annotation: {(~asd['hpo_annotated']).sum()}")

print("\nGenes with no HPO annotation:")
print(asd.loc[~asd["hpo_annotated"], "gene"].to_list())

# HPO phenotype categories to investigate
categories = {
    "seizures": {"HP:0001250"},
    "ID_DD": {"HP:0001249", "HP:0001263"},
    "ADHD": {"HP:0007018"},
    "sleep": {"HP:0002360"},
    "speech_language_delay": {"HP:0000750"},
    "motor_delay": {"HP:0001270"},
    "anxiety": {"HP:0000739"},
    "depression": {"HP:0000716"},
    "aggressive_behaviour": {"HP:0000718"},
}

#Keep HPO records belonging to 50 selected ASD genes
asd_genes = set(asd["gene"])
hpo_asd = hpo[hpo["gene_symbol"].isin(asd_genes)].copy()

#check how many of the ASD genes fall into each category
for category, hpo_ids in categories.items():
    genes = set(hpo_asd.loc[hpo_asd["hpo_id"].isin(hpo_ids), "gene_symbol"])

    print(f"\n{category}: {len(genes)} genes")
    print(sorted(genes))