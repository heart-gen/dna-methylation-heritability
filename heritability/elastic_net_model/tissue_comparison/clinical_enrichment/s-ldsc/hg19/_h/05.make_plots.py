import pandas as pd
import os

#################################################
# 1. Combine all files into one long-format table
#################################################

# Read file paths from text file
with open("ldsc_file_list.txt", "r") as f:
    file_paths = f.read().strip().split()

dfs = []

for f in file_paths:
    basename = os.path.basename(f)  # e.g., ad_caudate_non_heritable_hg19.results
    # Remove file extension
    name_no_ext = os.path.splitext(basename)[0]  # ad_caudate_non_heritable_hg19
    
    parts = name_no_ext.split("_")
    
    disorder = parts[0]         # "ad"
    brain_region = parts[1]     # "caudate"
    # Join the remaining parts except the last if it's hg19
    if parts[-1] == "hg19":
        h2_metric = "_".join(parts[2:-1])  # "non_heritable" or "heritable"
    else:
        h2_metric = "_".join(parts[2:])    # fallback
    
    # Read LDSC file (tab-delimited)
    df = pd.read_csv(f, sep="\t")
    
    # Add metadata
    df["Disorder"] = disorder
    df["BrainRegion"] = brain_region
    df["h2Metric"] = h2_metric
    
    dfs.append(df)

# Combine all dataframes
combined_df = pd.concat(dfs, ignore_index=True)

# Reorder columns
cols = ["Disorder", "BrainRegion", "h2Metric"] + [c for c in combined_df.columns if c not in ["Disorder", "BrainRegion", "h2Metric"]]
combined_df = combined_df[cols]

# Save tab-delimited
combined_df.to_csv("combined_ldsc_results.tsv", sep="\t", index=False)

print("Combined LDSC results saved! Shape:", combined_df.shape)

######################################
# 2. Filter for significant enrichment
######################################

df = pd.read_csv("combined_ldsc_results.tsv", sep="\t")

print("Total rows:", df.shape[0])
print("Rows with valid Enrichment_p:", df['Enrichment_p'].notna().sum())

from statsmodels.stats.multitest import multipletests

# Load combined LDSC results
df = pd.read_csv("combined_ldsc_results.tsv", sep="\t")

# Keep only rows with valid enrichment p-values and make a copy
df_valid = df[df["Enrichment_p"].notna()].copy()

# Compute FDR per group
df_valid["FDR"] = df_valid.groupby(["Disorder", "BrainRegion", "h2Metric"])["Enrichment_p"].transform(
    lambda x: multipletests(x, method="fdr_bh")[1]
)

# Filter significant enrichment
sig_df = df_valid[(df_valid["FDR"] < 0.05) & (df_valid["Enrichment"] > 1)]

# Save
sig_df.to_csv("ldsc_significant_enrichment.tsv", sep="\t", index=False)

print("Rows passing filter:", sig_df.shape[0])

# Remove the outlier category
filtered_df = sig_df[sig_df['Category'] != "MAF_Adj_ASMCL2_0"]

##############################################
# 3. Aggregate/enumerate top/shared categories
##############################################

# Example: filtered_df contains filtered significant enrichment
# Ensure your functional category column is named consistently
category_col = "Category"  # or "Functional_Category" depending on your file

# Count how many disorders each category is enriched in
disorders_per_category = filtered_df.groupby(category_col)["Disorder"].nunique().reset_index()
disorders_per_category.rename(columns={"Disorder": "Disorders_Enriched"}, inplace=True)

# Count how many brain regions each category is enriched in
regions_per_category = filtered_df.groupby(category_col)["BrainRegion"].nunique().reset_index()
regions_per_category.rename(columns={"BrainRegion": "BrainRegions_Enriched"}, inplace=True)

# Compute average enrichment per category
avg_enrichment = filtered_df.groupby(category_col)["Enrichment"].mean().reset_index()
avg_enrichment.rename(columns={"Enrichment": "Average_Enrichment"}, inplace=True)

# Merge all summaries into one table
summary_df = disorders_per_category.merge(regions_per_category, on=category_col)
summary_df = summary_df.merge(avg_enrichment, on=category_col)

# Sort by number of disorders enriched and average enrichment
summary_df = summary_df.sort_values(["Disorders_Enriched", "Average_Enrichment"], ascending=[False, False])

# Save summary table
summary_df.to_csv("ldsc_summary_table.tsv", sep="\t", index=False)

print("Summary table saved! Preview:")
print(summary_df.head(10))

#################
# 4. Make heatmap
#################

import seaborn as sns
import matplotlib.pyplot as plt

# Assume filtered_df is your filtered significant enrichment
# Ensure columns: 'Category', 'Disorder', 'BrainRegion', 'h2Metric', 'Enrichment', 'FDR'

# Select top N categories by average enrichment for clarity
top_categories = filtered_df.groupby('Category')['Enrichment'].mean().sort_values(ascending=False).head(20).index
heatmap_df = filtered_df[filtered_df['Category'].isin(top_categories)]

# Pivot to MultiIndex: rows = Category, columns = Disorder, BrainRegion, h2Metric
heatmap_matrix = heatmap_df.pivot_table(
    index='Category',
    columns=['Disorder','BrainRegion','h2Metric'],
    values='Enrichment',
    fill_value=0
)

# Flatten MultiIndex: join with underscores for column labels
heatmap_matrix.columns = ['_'.join(col) for col in heatmap_matrix.columns]

# --- Create annotation matrix with asterisks based on FDR ---
annot_matrix = heatmap_df.pivot_table(
    index='Category',
    columns=['Disorder','BrainRegion','h2Metric'],
    values='FDR',
    fill_value=1  # non-significant default
)
annot_matrix.columns = ['_'.join(col) for col in annot_matrix.columns]

# Function to convert FDR to asterisks
def fdr_to_stars(fdr):
    if fdr < 0.0001:
        return '****'
    elif fdr < 0.001:
        return '***'
    elif fdr < 0.01:
        return '**'
    elif fdr < 0.05:
        return '*'
    else:
        return ''

# Apply function to entire annotation matrix
annot_matrix = annot_matrix.applymap(fdr_to_stars)

# --- Plot heatmap with asterisks ---
plt.figure(figsize=(16, 10))
sns.heatmap(
    heatmap_matrix,
    cmap='Purples',
    linewidths=0.5,
    annot=annot_matrix,  # asterisks instead of numbers
    fmt='',
    cbar_kws={'label': 'Enrichment'}
)
plt.title("LDSC Partitioned Heritability Enrichment Heatmap")
plt.xlabel("Disorder_BrainRegion_h2Metric")
plt.ylabel("Functional Category")

plt.tight_layout()

# Save as PNG
plt.savefig("ldsc_enrichment_heatmap_multiindex_asterisks.png", dpi=300)
plt.show()



