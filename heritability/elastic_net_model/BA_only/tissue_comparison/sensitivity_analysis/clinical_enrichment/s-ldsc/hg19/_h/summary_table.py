import pandas as pd
import glob
import os

# Base directory containing all results
base_dir = "./results"

# Diseases to keep
allowed_diseases = {"height", "smoking", "substance_abuse", "ad", "pd", "scz"}

# Find all result files recursively
file_pattern = "**/*.results"

rows = []

for filepath in glob.glob(os.path.join(base_dir, file_pattern), recursive=True):
    try:
        parts = filepath.split(os.sep)

        # Extract metadata
        disease = parts[-4]
        tissue = parts[-3]
        h2 = parts[-2]

        # Skip unwanted diseases
        if disease not in allowed_diseases:
            continue

        # Read file
        df = pd.read_csv(filepath, sep="\t")

        # Filter for L2_1
        l2_row = df[df.iloc[:, 0] == "L2_1"]

        if not l2_row.empty:
            l2_row = l2_row.copy()

            # Keep only Coefficient_z-score + metadata
            keep_cols = ["Coefficient_z-score"]
            if "Coefficient_z-score" not in l2_row.columns:
                print(f"Warning: 'Coefficient_z-score' not found in {filepath}")
                continue

            l2_row = l2_row[keep_cols]
            l2_row["disease"] = disease
            l2_row["tissue"] = tissue
            l2_row["heritability"] = h2
            l2_row["source_file"] = os.path.basename(filepath)

            rows.append(l2_row)
        else:
            print(f"No L2_1 row in {filepath}")

    except Exception as e:
        print(f"Error processing {filepath}: {e}")

# Combine all rows
if rows:
    combined_df = pd.concat(rows, ignore_index=True)

    # Reorder columns (metadata first)
    cols = ["disease", "tissue", "heritability", "source_file", "Coefficient_z-score"]
    combined_df = combined_df[cols]

    # Save
    combined_df.to_csv("combined_L2_1_Coeff_zscore.tsv", sep="\t", index=False)
    print("Done! Saved combined_L2_1_Coeff_zscore.tsv")
else:
    print("No matching L2_1 rows found.")