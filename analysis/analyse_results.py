import pandas as pd
from pathlib import Path

cell_subtype = "RG"

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
summary = df.groupby("panel")["Shift_to_goal_end"].agg(
    count="count",
    mean="mean",
    median="median",
    min="min",
    max="max"
)

summary["Q1"] = df.groupby("panel")["Shift_to_goal_end"].quantile(0.25)
summary["Q3"] = df.groupby("panel")["Shift_to_goal_end"].quantile(0.75)
summary["IQR"] = summary["Q3"] - summary["Q1"]

#Count shift towards/away from VPA
direction_counts = pd.crosstab(df["panel"], df["direction"])

#add direction counts to summary table
summary = summary.join(direction_counts)

#save summary statistics
summary.to_csv(summary_file)

#rank genes by shifts
ranked = df.sort_values("Shift_to_goal_end", ascending=False)

ranked.to_csv(ranked_file, index=False)

print("\nSummary statistics:")
print(summary)

print("\nStrongest shifts towards VPA:")
print(ranked[["gene", "panel", "Shift_to_goal_end"]].head(10))

print("\nStrongest shifts away from VPA:")
print(ranked[["gene", "panel", "Shift_to_goal_end"]].tail(10))