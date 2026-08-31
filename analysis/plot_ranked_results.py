import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

cell_subtype = "IPC"

#Project folders
project_dir = Path(__file__).resolve().parents[1]
results_dir = project_dir / "results"
figures_dir = project_dir / "figures"

#files
input_file = results_dir / f"{cell_subtype}_ranked_results.csv"
output_file = figures_dir / f"{cell_subtype}_ranked_gene_plot.png"

#load results
df = pd.read_csv(input_file)

#sort genes from greatest shift away to greatest shift towards
df = df.sort_values("Shift_to_goal_end", ascending=True).reset_index(drop=True)

#add rank position for plotting
df["rank"] = range(1, len(df) + 1)

#separate asd & non-asd genes
asd = df[df["panel"] == "ASD"]
non_asd = df[df["panel"] == "non_ASD"]

fig, ax = plt.subplots(figsize=(12,6))

#plot non-ASD genes
ax.scatter(
    non_asd["rank"],
    non_asd["Shift_to_goal_end"],
    color="hotpink",
    label="Non-ASD genes",
    marker="o",
    s=50,
    alpha=0.8)

#plot ASD genes
ax.scatter(
    asd["rank"],
    asd["Shift_to_goal_end"],
    color="lime",
    label="ASD genes",
    marker="o",
    s=50,
    alpha=0.8)

#add zero reference line
ax.axhline(
    y=0,
    color="black",
    linestyle="--",
    linewidth=1)

#label 5 strongest shifts in each direction
strongest = pd.concat([df.head(3), df.tail(3)])

for _, row in strongest.iterrows():
    ax.annotate(
        row["gene"],
        (row["rank"], row["Shift_to_goal_end"]),
        xytext=(4,5),
        textcoords="offset points",
        fontsize=8)

#labels and titles
ax.set_xlabel("Genes ranked by predicted perturbation shift", fontsize=11)
ax.set_ylabel("Predicted shift towards VPA state", fontsize=11)
ax.set_title("Ranked Geneformer Perturbation Shifts in Control IPC/Newborn ExN cells", fontsize=13)

#display small shift values using scientific notation
ax.ticklabel_format(
    axis="y",
    style="sci",
    scilimits=(0, 0)
)

ax.legend(title="Gene panel", frameon=False)

#tidy up plot
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
plt.tight_layout()

#Save
plt.savefig(output_file, dpi=300, bbox_inches="tight")

plt.show()