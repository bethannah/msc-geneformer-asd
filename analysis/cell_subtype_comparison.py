from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt


project_dir = Path(__file__).resolve().parents[1]
results_dir = project_dir / "results"
figures_dir = project_dir / "figures"

#load results
cell_subtypes = {
    "RG": "combined_RG_stats.csv",
    "IPC": "combined_IPC_stats.csv"
}

#merge results
comparison = None

for subtype, filename in cell_subtypes.items():
    df = pd.read_csv(results_dir / filename)
    df = df[["gene", "panel", "syndromic", "Shift_to_goal_end"]].rename(columns={"Shift_to_goal_end": f"{subtype}_shift"})

    if comparison is None:
        comparison = df

    else:
        df = df.drop(columns="syndromic")

        comparison = comparison.merge(df, on=["gene", "panel"], how="inner")

#compare direction of effect
def direction_pattern(row):
    if row["RG_shift"] > 0 and row["IPC_shift"] > 0:
        return "consistent_towards"
    elif row["RG_shift"] < 0 and row["IPC_shift"] < 0:
        return "consistent_away"
    else:
        return "opposite_direction"

comparison["direction_pattern"] = comparison.apply(direction_pattern, axis=1)

#difference in magnitutte between subtypes
comparison["absolute_difference"] = (comparison["RG_shift"] - comparison["IPC_shift"]).abs()

comparison.to_csv(results_dir / "cell_subtype_comparison.csv", index=False)

#RG vs IPC scatterplot
fig, ax = plt.subplots(figsize=(7,10))

for panel, colour in [("ASD", "lime"), ("non_ASD", "hotpink")]:
    data = comparison[comparison["panel"] == panel]

    ax.scatter(data["RG_shift"],
               data["IPC_shift"],
               color=colour,
               label=panel,
               alpha=0.8
               )

ax.axhline(0, color="black", linestyle="--", linewidth=1)
ax.axvline(0, color="black", linestyle="--", linewidth=1)

ax.set_xlabel("Predicted shift towards VPA-associated state in RG")
ax.set_ylabel("Predicted shift towards VPA-associated state in IPC/ Newborn ExN")
ax.set_title("Geneformer perturbation shifts: RG vs IPC/Newborn ExN")
ax.legend()

plt.tight_layout()
plt.savefig(figures_dir / "RG_IPC_scatter.png", dpi=300, bbox_inches="tight")

plt.close()


#Gene x cell subtype heatmap
heatmap = comparison.set_index("gene")[["RG_shift", "IPC_shift"]]

#select 25 genes with greatest difference between cell subtypes
difference = (comparison["RG_shift"] - comparison["IPC_shift"]).abs()

top_genes = comparison.loc[difference.nlargest(25).index, "gene"]

heatmap = heatmap.loc[top_genes]

#centre colour scale on zero
limit=abs(heatmap.to_numpy()).max()

fig, ax = plt.subplots(figsize=(5, 9))

image = ax.imshow(heatmap, aspect="auto", cmap="coolwarm",vmin=-limit, vmax=limit)

ax.set_xticks([0, 1])
ax.set_xticklabels(["RG", "IPC/Newborn ExN"])

ax.set_yticks(range(len(heatmap)))
ax.set_yticklabels(heatmap.index, fontsize=8)

ax.set_title("Largest Perturbation Shift Differences Between RG and IPC/Newborn ExN")

fig.colorbar(image, ax=ax, label="Predicted perturbation shift relative to VPA-associated transcriptomic state")

plt.tight_layout()
plt.savefig(figures_dir / "cell_subtype_heatmap.png", dpi=300, bbox_inches="tight")

plt.close()
