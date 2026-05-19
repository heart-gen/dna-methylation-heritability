import pandas as pd
import glob
import os
from scipy.stats import norm
from statsmodels.stats.multitest import multipletests

# Base directory containing all results
base_dir = "./results"

# Diseases to keep
allowed_diseases = {"smoking", "ad", "pd", "scz"}

# Find all result files recursively
file_pattern = "**/*.results"

rows = []

for filepath in glob.glob(os.path.join(base_dir, file_pattern), recursive=True):
    try:
        parts = filepath.split(os.sep)

        # Extract metadata
        disease = parts[-3]
        tissue = parts[-2]

        # Skip unwanted diseases
        if disease not in allowed_diseases:
            continue

        # Read file
        df = pd.read_csv(filepath, sep="\t")

        # Filter for L2_1
        l2_row = df[df.iloc[:, 0].str.contains(r"^ANNOT_Q\d+L2_1$", regex=True, na=False)]

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
            l2_row["annotation"] = l2_row.iloc[:, 0].str.extract(r"(ANNOT_Q\d+L2_1)")

            rows.append(l2_row)
        else:
            print(f"No ANNOT_Q#L2_1 rows in {filepath}")

    except Exception as e:
        print(f"Error processing {filepath}: {e}")

# Combine all rows
if rows:
    combined_df = pd.concat(rows, ignore_index=True)

    # Reorder columns (metadata first)
    cols = ["disease", "tissue", "annotation", "Coefficient_z-score"]
    combined_df = combined_df[cols]

    # Save
    combined_df.to_csv("combined_ANNOT_Q_L2_1_Coeff_zscore.tsv", sep="\t", index=False)
    print("Done! Saved combined_ANNOT_Q_L2_1_Coeff_zscore.tsv")
else:
    print("No matching ANNOT_Q#L2_1 rows found.")

# --- Step 1: Load your summary table ---
input_file = "combined_ANNOT_Q_L2_1_Coeff_zscore.tsv"  # change this to your file path
df = pd.read_csv(input_file, sep="\t")

# --- Step 2: Convert Coefficient_z_score to two-tailed p-values ---
df['p_value'] = 2 * (1 - norm.cdf(df['Coefficient_z-score'].abs()))

# --- Step 3: Apply FDR correction (Benjamini-Hochberg) ---
df['FDR_q'] = multipletests(df['p_value'], method='fdr_bh')[1]

# --- Step 4: Save the new table ---
output_file = "summary_table_fdr.tsv"
df.to_csv(output_file, sep="\t", index=False)

print(f"FDR-corrected table saved to {output_file}")