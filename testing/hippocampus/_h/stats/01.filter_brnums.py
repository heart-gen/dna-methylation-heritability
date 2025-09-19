import pandas as pd

# Load main dataset (tab-delimited)
main_df = pd.read_csv("/projects/p32505/users/elisa/inputs/phenotypes/_m/phenotypes-AA.tsv", sep="\t")

# Load brnums from the filter file (assuming it's tab-delimited or space-delimited)
filter_df = pd.read_csv("/projects/p32505/users/elisa/testing/hippocampus/_m/samples.txt", sep="\t", header=None)

# First column of filter file is brnum
brnums_to_keep = filter_df[0].tolist()

# Filter main dataset
filtered = main_df[main_df["brnum"].isin(brnums_to_keep)]

# Remove duplicate brnums (keep the first occurrence)
filtered = filtered.drop_duplicates(subset=["brnum"], keep="first")

# Save the result
filtered.to_csv("filtered_file.txt", sep="\t", index=False)

print(f"Filtered file saved with {len(filtered)} unique brnums.")
