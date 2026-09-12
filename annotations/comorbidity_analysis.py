from pathlib import Path
import pandas as pd

project_dir = Path(__file__).resolve().parents[1]

comparison_file = project_dir / "results" / "cell_subtype_comparison.csv"
annotation_file = project_dir / "annotations" / "comorbidity_annotations.csv"

comparison = pd.read_csv(comparison_file)
annotations = pd.read_csv(annotation_file)

#keep asd genes only
comparison_asd = comparison[comparison["panel"] == "ASD"].copy()

#merge perturbation rsults awith HPO annotations
combined = comparison_asd.merge(annotations, on="gene", how="left")

output_file = (project_dir / "results" / "asd_integrated_results.csv")

combined.to_csv(output_file, index=False)

#phenotype categories to summarise
phenotypes = [
    "seizures",
    "ID_DD",
    "ADHD",
    "sleep",
    "speech_language_delay",
    "motor_delay",
    "anxiety",
    "aggressive_behaviour"
]

summary = []

for phenotype in phenotypes:
    genes = combined[combined[phenotype] == 1]

    summary.append({
        "phenotype": phenotype,
        "n_genes": len(genes),
        "median_RG_shift": genes["RG_shift"].median(),
        "median_IPC_shift": genes["IPC_shift"].median(),
        "RG_towards_VPA_percent": (genes["RG_shift"] > 0).mean() * 100,
        "IPC_towards_VPA_percent": (genes["IPC_shift"] > 0).mean() * 100,
        "opposite_direction_%": (genes["direction_pattern"] == "opposite_direction").mean() * 100,
        "median_absolute_difference": genes["absolute_difference"].median()
    })

phenotype_summary = pd.DataFrame(summary)

summary_file = project_dir / "results" / "phenotype_summary.csv"
phenotype_summary.to_csv(summary_file, index=False)