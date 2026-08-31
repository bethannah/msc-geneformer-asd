import pandas as pd
from scipy.stats import mannwhitneyu
from pathlib import Path

#cell subtype
cell_subtype = "IPC"

#project folders
project_dir = Path(__file__).resolve().parents[1]
results_dir = project_dir / "results"

#files
input_file = results_dir / f"combined_{cell_subtype}_stats.csv"
output_file = results_dir / f"{cell_subtype}_statistical_comparison.csv"

#load results
df = pd.read_csv(input_file)

#function for statistical comparison
def compare_stats(group1, group2, comparison, group1_name, group2_name):

    u_stat, p_value = mannwhitneyu(group1, group2, alternative="two-sided")

    rank_biserial = ((2 * u_stat) / (len(group1) * len(group2))) - 1


    return {
        "comparison": comparison,
        "group1": group1_name,
        "group2": group2_name,
        "group1_n": len(group1),
        "group2_n": len(group2),
        "group1_median": group1.median(),
        "group2_median": group2.median(),
        "median_difference": group1.median() - group2.median(),
        "Mann_Whitney_U": u_stat,
        "p_value": p_value,
        "rank_biserial_correlation": rank_biserial
    }

#ASD vs non-ASD
asd = df[df["panel"] == "ASD"]["Shift_to_goal_end"]
non_asd = df[df["panel"] == "non_ASD"]["Shift_to_goal_end"]

panel_result = compare_stats(asd, non_asd, "ASD vs non_ASD", "ASD", "Non-ASD")

#Syndromic vs non-syndromic ASD
asd_df = df[df["panel"] == "ASD"]

syndromic = asd_df[asd_df["syndromic"] == 1]["Shift_to_goal_end"]
non_syndromic = asd_df[asd_df["syndromic"] == 0]["Shift_to_goal_end"]

syndromic_result = compare_stats(syndromic, non_syndromic, "Syndromic vs non-syndromic", "Syndromic", "Non-Syndromic")

#Combine and save results
results = pd.DataFrame([panel_result, syndromic_result])

results.to_csv(output_file, index=False)

