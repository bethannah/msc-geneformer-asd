import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

cell_subtype = "RG"

#files
input_file = f"combined_{cell_subtype}_stats.csv"
output_file = f"{cell_subtype}_ASD_vs_nonASD_plot.png"

#load results
df = pd.read_csv(input_file)
asd = df[df["panel"] == "ASD"]["Shift_to_goal_end"]
non_asd = df[df["panel"] == "non_ASD"]["Shift_to_goal_end"]

fig, ax = plt.subplots(figsize=(7.5, 7))

#boxplots
box = ax.boxplot(
    [asd, non_asd],
    tick_labels=["ASD", "Non-ASD"],
    patch_artist=True,
    widths=0.5,
    showfliers=False)

box["boxes"][0].set_facecolor("limegreen")
box["boxes"][1].set_facecolor("hotpink")

for patch in box["boxes"]:
    patch.set_alpha(0.5)

#add individual gene values
rng = np.random.default_rng(42)
asd_x = rng.normal(1, 0.08, len(asd))
non_asd_x = rng.normal(2, 0.08, len(non_asd))

ax.scatter(
    asd_x,
    asd,
    color="forestgreen",
    alpha=0.7,
    s=35)

ax.scatter(
    non_asd_x,
    non_asd,
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
ax.set_xlabel("Gene panel", fontsize=11)
ax.set_ylabel("Predicted shift towards VPA state", fontsize=11)
ax.set_title(f"Geneformer Perturbation Shifts: ASD vs non-ASD Genes ({cell_subtype})", fontsize=13)

#Scientific notation
ax.ticklabel_format(axis="y", style="sci", scilimits=(0,0))

#tidy plot
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
plt.tight_layout()

#Save
plt.savefig(output_file, dpi=300, bbox_inches="tight")

plt.show()
