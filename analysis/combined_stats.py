from pathlib import Path
import pandas as pd

#Project folders
project_dir = Path(__file__).resolve().parents[1]

#Per-gene perturbation statistics generated on the HPC- this directory isn't included in the Github repository
stats_dir = project_dir / "stats" / "batch_perturbations"
metadata_file = project_dir/ "gene_selection" / "gene_metadata.csv"
output_file = project_dir / "results" / "combined_RG_stats.csv"

results = []

#Go through genes stats folder
for folder in stats_dir.iterdir():
    if not folder.is_dir() or folder.name == "archive":
        continue

    #each csv has same name as folder
    csv_file = folder / f"{folder.name}.csv"

    if not csv_file.exists():
        continue

    stats = pd.read_csv(csv_file)

    #work out panel & gene name from folder name
    if folder.name.startswith("non_ASD_"):
        panel = "non_ASD"
        gene = folder.name.removeprefix("non_ASD_").removesuffix("_RG")

    elif folder.name.startswith("ASD_"):
        panel = "ASD"
        gene = folder.name.removeprefix("ASD_").removesuffix("_RG")

    else:
        continue

    #Extract shift to end goal
    shift = stats["Shift_to_goal_end"].iloc[0]

    results.append({
        "gene": gene,
        "panel": panel,
        "Shift_to_goal_end": shift
    })

#combine all resuts
results = pd.DataFrame(results)

#load metadata
metadata = pd.read_csv(metadata_file)

#merge perturbation results with metadata
combined = metadata.merge(results, on=["gene", "panel"], how="inner")

#add shift direction
combined["direction"] = combined["Shift_to_goal_end"].apply(
    lambda x: "towards_VPA" if x > 0
    else "away_from_VPA" if x < 0
    else "no_shift"
)

#save
combined.to_csv(output_file, index=False)
print(f"combined results: {len(combined)} genes")