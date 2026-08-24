import pandas as pd

df = pd.read_csv("asd_candidates_detection.csv")

#keep ASD genes detected in at least 50% of CTRL RG cells
eligible = df[df["CTRL_RG_detection_percentage"] >= 50].copy()

#randomly select 50 for reproducibility
selected = eligible.sample(n=50, random_state=42)

#sort for viewing
selected = selected.sort_values(
    "CTRL_RG_detection_percentage",
    ascending=False)

selected.to_csv("selected_asd_50.csv", index=False)
