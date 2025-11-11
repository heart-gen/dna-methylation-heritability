import pandas as pd
import os
import argparse
from pyliftover import LiftOver

parser = argparse.ArgumentParser(description="Separate TSV into groups with genomic sorting and liftover from hg38 to hg19")
parser.add_argument("--input_file", required=True, help="Path to input TSV file")
parser.add_argument("--output_dir", required=True, help="Directory to save output files")
parser.add_argument("--chain_file", required=True, help="Path to hg38ToHg19.over.chain.gz")
args = parser.parse_args()

input_file = args.input_file
output_dir = args.output_dir
chain_file = args.chain_file

os.makedirs(output_dir, exist_ok=True)

# Load the data
df = pd.read_csv(input_file, sep='\t')

# Subset groups
heritable = df[(df["h2_unscaled"] >= 0.1) & (df["r_squared_cv"] > 0.75)]
non_heritable = df[(df["h2_unscaled"] < 0.1) & (df["r_squared_cv"] > 0.75)]

cols_to_keep = ["chrom", "start", "end"]
heritable = heritable[cols_to_keep].copy()
non_heritable = non_heritable[cols_to_keep].copy()

# Optional: define low_prediction if you plan to use it later
# low_prediction = df[(df["r_squared_cv"] <= 0.75)][cols_to_keep].copy()

# Convert to "chr" format
for df_subset in [heritable, non_heritable]:
    df_subset["chrom"] = df_subset["chrom"].astype(str).apply(lambda x: f"chr{x}" if not str(x).startswith("chr") else x)

# Perform hg38 → hg19 liftover
lo = LiftOver(chain_file)

def liftover_coordinate(chrom, pos):
    try:
        new = lo.convert_coordinate(chrom, int(pos))
        if new:
            return int(new[0][1])
        return None
    except Exception:
        return None

for df_subset in [heritable, non_heritable]:
    df_subset["start_hg19"] = df_subset.apply(lambda x: liftover_coordinate(x["chrom"], x["start"]), axis=1)
    df_subset["end_hg19"] = df_subset.apply(lambda x: liftover_coordinate(x["chrom"], x["end"]), axis=1)
    # Drop unmapped entries
    df_subset.dropna(subset=["start_hg19", "end_hg19"], inplace=True)
    # Replace start/end with hg19 positions
    df_subset["start"] = df_subset["start_hg19"].astype(int)
    df_subset["end"] = df_subset["end_hg19"].astype(int)
    df_subset.drop(columns=["start_hg19", "end_hg19"], inplace=True)

# Sort by chromosome and position
chrom_order = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
for df_subset in [heritable, non_heritable]:
    df_subset["chrom"] = pd.Categorical(df_subset["chrom"], categories=chrom_order, ordered=True)
    df_subset.sort_values(by=["chrom", "start"], ascending=[True, True], inplace=True)

# Save BED files
heritable.to_csv(os.path.join(output_dir, "heritable_hg19.bed"), sep='\t', index=False, header=False)
non_heritable.to_csv(os.path.join(output_dir, "non_heritable_hg19.bed"), sep='\t', index=False, header=False)

print("Separation and liftover complete! Files saved without headers.")
print(f"Heritable group: {len(heritable)} rows → {os.path.join(output_dir, 'heritable_hg19.bed')}")
print(f"Non-heritable group: {len(non_heritable)} rows → {os.path.join(output_dir, 'non_heritable_hg19.bed')}")