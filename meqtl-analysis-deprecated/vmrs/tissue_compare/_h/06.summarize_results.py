import pandas as pd
import glob
import os
import seaborn as sns
import matplotlib.pyplot as plt

# Paths to your gzipped input files
lfsr_folder = "./local/temp/"
posterior_folder = "./local/temp/"
intermediate_folder = "filtered_rows/"
os.makedirs(intermediate_folder, exist_ok=True)

# Define brain regions to filter by
brain_regions = ["Caudate", "Dlpfc", "Hippocampus"]  # change this to the regions you want

# Get all gzipped files
lfsr_files = sorted(glob.glob(os.path.join(lfsr_folder, "lfsr_*_ancestry.txt.gz")))
posterior_files = sorted(glob.glob(os.path.join(posterior_folder, "posterior_mean_*_ancestry.txt.gz")))

# Ensure files correspond
assert len(lfsr_files) == len(posterior_files), "Mismatch in number of files"

all_results = []

for lfsr_file, posterior_file in zip(lfsr_files, posterior_files):
    print(f"\nProcessing {lfsr_file} and {posterior_file}...")

    # Load gzipped input files
    lfsr_df = pd.read_csv(lfsr_file, sep="\t", compression='gzip')
    posterior_df = pd.read_csv(posterior_file, sep="\t", compression='gzip')
    
    # Filter rows where all selected brain regions ≤ 0.05
    if not all(region in lfsr_df.columns for region in brain_regions):
        raise ValueError(f"One or more specified brain regions not in LFSR file columns: {brain_regions}")
    
    lfsr_filtered = lfsr_df[(lfsr_df[brain_regions] <= 0.05).all(axis=1)]
    
    # Keep only corresponding posterior rows
    posterior_filtered = posterior_df[posterior_df['effect'].isin(lfsr_filtered['effect'])]
    
    # Print filtered rows preview
    print("Filtered rows preview:")
    print(posterior_filtered.head(10))
    print(f"Total filtered rows: {len(posterior_filtered)}")
    
    # Save intermediate filtered posterior rows as plain TSV
    base_name = os.path.basename(posterior_file).replace(".tsv.gz", "_filtered.tsv")
    intermediate_path = os.path.join(intermediate_folder, base_name)
    posterior_filtered.to_csv(intermediate_path, sep="\t", index=False)
    
    # Average by phenotype_id
    avg_posterior = posterior_filtered.groupby('phenotype_id').mean(numeric_only=True).reset_index()
    
    # Include 'effect' column as a list
    avg_posterior['effects'] = posterior_filtered.groupby('phenotype_id')['effect'].agg(list).values
    
    # Store source filename
    avg_posterior['source_file'] = os.path.basename(posterior_file)
    
    all_results.append(avg_posterior)

# Concatenate all averaged results
final_df = pd.concat(all_results, ignore_index=True)

# Save final averaged results as plain TSV
final_df.to_csv("local_all_averaged_posterior_filtered.tsv", sep="\t", index=False)

print("\nAll files processed, intermediate filtered rows saved as TSV, and final averaged results saved as TSV.")

# Load your TSV
df = pd.read_csv("local_all_averaged_posterior_filtered.tsv", sep="\t")

# Set phenotype_id as the index
df_heat = df.set_index("phenotype_id")[["Caudate", "Dlpfc", "Hippocampus"]]

# Optional: sort by mean across regions for nicer visualization
df_heat["mean"] = df_heat.mean(axis=1)
df_heat = df_heat.sort_values("mean", ascending=False).drop(columns="mean")

# Set plot style
sns.set(style="white")

plt.figure(figsize=(10, max(6, 0.5*len(df_heat))))  # adjust height based on number of VMRs

# Create heatmap
sns.heatmap(
    df_heat,
    cmap="coolwarm",          # blue-red color scale
    center=0,                 # white at 0
    annot=True,               # optional: shows the numbers
    fmt=".3f",
    linewidths=.5,
    cbar_kws={'label': 'Posterior Mean'}
)

plt.title("Posterior Means by VMR and Brain Region")
plt.xlabel("Brain Region")
plt.ylabel("VMR (Phenotype ID)")
plt.tight_layout()

# Save as PNG
plt.savefig("local_posterior_means_heatmap.png", dpi=300)

plt.show()

print("Heatmap saved as local_posterior_means_heatmap.png")

