import pandas as pd
import os
import argparse

parser = argparse.ArgumentParser(description="Separate TSV into groups with genomic sorting")
parser.add_argument("--input_file", required=True, help="Path to input TSV file")
parser.add_argument("--output_dir", required=True, help="Directory to save output files")
args = parser.parse_args()

input_file = args.input_file
output_dir = args.output_dir

os.makedirs(output_dir, exist_ok=True)

df = pd.read_csv(input_file, sep='\t')

heritable = df[(df["h2_unscaled"] > 0.1) & (df["r_squared_cv"] > 0.75)]
non_heritable = df[(df["h2_unscaled"] < 0.1) & (df["r_squared_cv"] > 0.75)]
low_prediction = df[df["r_squared_cv"] < 0.75]

cols_to_keep = ["chrom", "start", "end"]
heritable = heritable[cols_to_keep]
non_heritable = non_heritable[cols_to_keep]
low_prediction = low_prediction[cols_to_keep]

for df_subset in [heritable, non_heritable, low_prediction]:
    df_subset["chrom"] = df_subset["chrom"].astype(str).apply(lambda x: f"chr{x}")

chrom_order = [f"chr{i}" for i in range(1, 23)]

for df_subset in [heritable, non_heritable, low_prediction]:
    df_subset["chrom"] = pd.Categorical(df_subset["chrom"], categories=chrom_order, ordered=True)
    df_subset.sort_values(by=["chrom", "start"], ascending=[True, True], inplace=True)

heritable.to_csv(os.path.join(output_dir, "heritable.bed"), sep='\t', index=False, header=False)
non_heritable.to_csv(os.path.join(output_dir, "non_heritable.bed"), sep='\t', index=False, header=False)
low_prediction.to_csv(os.path.join(output_dir, "low_prediction.bed"), sep='\t', index=False, header=False)

print("Separation complete! Files saved without headers.")
print(f"Heritable group: {len(heritable)} rows → {os.path.join(output_dir, 'heritable.bed')}")
print(f"Non-heritable group: {len(non_heritable)} rows → {os.path.join(output_dir, 'non_heritable.bed')}")
print(f"Low-prediction group: {len(low_prediction)} rows → {os.path.join(output_dir, 'low_prediction.bed')}")