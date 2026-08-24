import pandas as pd

#Files
asd_file = "selected_asd_50.csv"
non_asd_file = "non_asd_candidates.csv"
output_file =  "matched_genes.csv"

#Load data
asd = pd.read_csv(asd_file)
non_asd = pd.read_csv(non_asd_file)

#Sort ASD genes by detection %
asd = asd.sort_values("CTRL_RG_detection_percentage", ascending=False)

matched_rows = []

#Match each ASD gene to the closest unused non-ASD gene
for _, asd_gene in asd.iterrows():
    non_asd["detection_difference"] = abs(
        non_asd["CTRL_RG_detection_percentage"] - asd_gene["CTRL_RG_detection_percentage"])

    best_match = non_asd.loc[
        non_asd["detection_difference"].idxmin()].copy()

    best_match["matched_ASD_gene"] = asd_gene["gene"]
    best_match["matched_ASD_detection_percentage"] = (asd_gene["CTRL_RG_detection_percentage"])

    matched_rows.append(best_match)

    #remove selected gene so can't be used again
    non_asd = non_asd.drop(best_match.name)

#combine matches
matched = pd.DataFrame(matched_rows)

#Save
matched.to_csv(output_file, index=False)

print(f"Mean detection difference: {matched['detection_difference'].mean():.2f}%")
print(f"Maximum detection difference: {matched['detection_difference'].max():.2f}%")