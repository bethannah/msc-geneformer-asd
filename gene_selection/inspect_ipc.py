import pandas as pd

df = pd.read_csv("gene_selection/gene_metadata.csv")

# Summary of IPC detection by panel
print("\nIPC detection summary:")
print(
    df.groupby("panel")["CTRL_IPC_detection_percentage"].describe()
)

# Genes detected in fewer than 50% of CTRL IPC cells
low_detection = df[
    df["CTRL_IPC_detection_percentage"] < 50
][["gene", "panel", "CTRL_IPC_detection_percentage"]]

print("\nGenes with <50% IPC detection:")
print(low_detection.to_string(index=False))

# Number below 50% in each panel
print("\nNumber of genes with <50% IPC detection:")
print(
    df.assign(
        below_50=df["CTRL_IPC_detection_percentage"] < 50
    ).groupby("panel")["below_50"].sum()
)

# Lowest-detected genes
print("\n10 lowest-detected genes in IPC:")
print(
    df[
        ["gene", "panel", "CTRL_IPC_detection_percentage"]
    ]
    .sort_values("CTRL_IPC_detection_percentage")
    .head(10)
    .to_string(index=False)
)