import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

cell_subtype = "RG"
comparison = "syndromic"

#Project folders
project_dir = Path(__file__).resolve().parents[1]
results_dir = project_dir / "results"
figures_dir = project_dir / "figures"

input_file = results_dir / f"combined_{cell_subtype}_stats.csv"

#load results
df = pd.read_csv(input_file)

#choose comparison to plot
if comparison == "panel":
    group1 = df[df["panel"] == "ASD"]["Shift_to_goal_end"]
    group2 = df[df["panel"] == "non_ASD"]["Shift_to_goal_end"]

    label1 = "ASD"
    label2 = "Non-ASD"

    output_file = figures_dir / f"{cell_subtype}_ASD_vs_nonASD_plot.png"
    title = f"Geneformer Perturbation Shifts: ASD vs Non-ASD Genes({cell_subtype})"

elif comparison == "syndromic":
    asd = df[df["panel"] == "ASD"]

    group1 = asd[asd["syndromic"] == 1]["Shift_to_goal_end"]
    group2 = asd[asd["syndromic"] == 0]["Shift_to_goal_end"]

    label1 = "Syndromic"
    label2 = "Non-syndromic"

    output_file = figures_dir / f"{cell_subtype}_syndromic_vs_non_syndromic_plot.png"
    title = f"Geneformer Perturbation Shifts: Syndromic vs Non-Syndromic ASD Genes ({cell_subtype})"

#create plot
fig, ax = plt.subplots(figsize=(7.5, 7))

#boxplots
box = ax.boxplot(
    [group1, group2],
    tick_labels=[label1, label2],
    patch_artist=True,
    widths=0.5,
    showfliers=False)

box["boxes"][0].set_facecolor("limegreen")
box["boxes"][1].set_facecolor("hotpink")

for patch in box["boxes"]:
    patch.set_alpha(0.5)

#add individual gene values
rng = np.random.default_rng(42)
group1_x = rng.normal(1, 0.08, len(group1))
group2_x = rng.normal(2, 0.08, len(group2))

ax.scatter(
    group1_x,
    group1,
    color="forestgreen",
    alpha=0.7,
    s=35)

ax.scatter(
    group2_x,
    group2,
    color="deeppink",
    alpha=0.7,
    s=35)

#zero reference line
ax.axhline(
    y=0,
    color="black",
    linestyle="--",
    linewidth=1)

#labels
ax.set_xlabel("ASD gene subgroup", fontsize=11)
ax.set_ylabel("Predicted shift towards VPA state", fontsize=11)
ax.set_title(title, fontsize=13)

#Scientific notation
ax.ticklabel_format(axis="y", style="sci", scilimits=(0,0))

#tidy plot
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
plt.tight_layout()

#Save
plt.savefig(output_file, dpi=300, bbox_inches="tight")

plt.show()
