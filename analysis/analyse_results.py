import pandas as pd
from pathlib import Path

cell_subtype = "IPC"

#Project folders
project_dir = Path(__file__).resolve().parents[1]
results_dir = project_dir / "results"

#files
input_file = results_dir / f"combined_{cell_subtype}_stats.csv"
summary_file = results_dir / f"{cell_subtype}_summary_stats.csv"
ranked_file = results_dir / f"{cell_subtype}_ranked_results.csv"

#Load combined results
df = pd.read_csv(input_file)

#summarise statistics
def summarise_groups(data, group_column):

    grouped = data.groupby(group_column)["Shift_to_goal_end"]
    
    summary = grouped.agg(
        count="count",
        mean="mean",
        median="median",
        min="min",
        max="max"
    )

    summary["Q1"] = grouped.quantile(0.25)
    summary["Q3"] = grouped.quantile(0.75)
    summary["IQR"] = summary["Q3"] - summary["Q1"]

    #Count shift towards/away from VPA
    direction_counts = pd.crosstab(data[group_column], data["direction"])

    return summary.join(direction_counts)

#save summary statistics
summary = summarise_groups(df, "panel")
summary.to_csv(summary_file)

#rank genes by predicted shifts
ranked = df.sort_values("Shift_to_goal_end", ascending=False)

ranked.to_csv(ranked_file, index=False)


#identify ASD genes with greatest predicted shifts
asd_ranked = ranked[ranked["panel"] == "ASD"].copy()

towards_vpa = asd_ranked.head(10)

away_vpa = asd_ranked.sort_values("Shift_to_goal_end", ascending=True).head(10)

top_asd_shifts = pd.concat([towards_vpa, away_vpa])

top_asd_shifts.to_csv(results_dir / f"{cell_subtype}_top_ASD_shifts.csv", index=False)


#Syndromic vs non-syndromic comparison
asd = df[df["panel"] == "ASD"].copy()

asd["syndromic_group"] = asd["syndromic"].map({
    1: "Syndromic",
    0: "Non-syndromic"
})

syndromic_summary = summarise_groups(asd, "syndromic_group")

syndromic_summary.to_csv(results_dir / f"{cell_subtype}_syndromic_summary.csv")