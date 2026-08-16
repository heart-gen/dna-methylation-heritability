import pandas as pd
import glob
import os
from scipy.stats import norm
from statsmodels.stats.multitest import multipletests

# Base directory containing all results
base_dir = "./EA"

# Diseases to keep
allowed_diseases = {"smoking", "ad", "pd", "scz"}

# Find all result files recursively
file_pattern = "**/*.results"

rows = []

for filepath in glob.glob(os.path.join(base_dir, file_pattern), recursive=True):
    try:
        parts = filepath.split(os.sep)

        # Expected structure:
        # ./$ancestry/results/$disease/$tissue/file.results

        # Extract metadata
        disease = parts[-3]
        tissue = parts[-2]
        ancestry = parts[-5]

        # Skip unwanted diseases
        if disease not in allowed_diseases:
            continue

        # Read file
        df = pd.read_csv(filepath, sep="\t")

        # Filter rows like ANNOT_Q#L2_1
        row_mask = df.iloc[:, 0].str.contains(
            r"^ANNOT_Q\d+L2_1$", regex=True, na=False
        )

        l2_rows = df[row_mask].copy()

        if not l2_rows.empty:

            # Check required column
            if "Coefficient_z-score" not in l2_rows.columns:
                print(f"Warning: 'Coefficient_z-score' not found in {filepath}")
                continue

            # Extract annotation names from matching rows only
            l2_rows["annotation"] = (
                l2_rows.iloc[:, 0]
                .str.extract(r"(ANNOT_Q\d+L2_1)", expand=False)
            )

            # Keep desired columns
            l2_rows = l2_rows[["annotation", "Coefficient_z-score"]]

            # Add metadata
            l2_rows["disease"] = disease
            l2_rows["tissue"] = tissue
            l2_rows["ancestry"] = ancestry

            rows.append(l2_rows)

        else:
            print(f"No ANNOT_Q#L2_1 rows in {filepath}")

    except Exception as e:
        print(f"Error processing {filepath}: {e}")

# Combine all rows
if rows:
    combined_df = pd.concat(rows, ignore_index=True)

    # Reorder columns
    cols = [
        "disease",
        "tissue",
        "ancestry",
        "annotation",
        "Coefficient_z-score"
    ]
    combined_df = combined_df[cols]

    # Save combined table
    combined_output = "combined_Coeff_zscore.tsv"
    combined_df.to_csv(combined_output, sep="\t", index=False)

else:
    print("No matching ANNOT_Q#L2_1 rows found.")
    exit()

# ---------------------------------------------------------
# Calculate p-values and FDR separately for each ancestry
# ---------------------------------------------------------

results = []

for ancestry, subdf in combined_df.groupby("ancestry"):

    subdf = subdf.copy()

    # Convert z-scores to two-tailed p-values
    subdf["p_value"] = 2 * (
        1 - norm.cdf(subdf["Coefficient_z-score"].abs())
    )

    # FDR correction within ancestry
    subdf["FDR_q"] = multipletests(
        subdf["p_value"],
        method="fdr_bh"
    )[1]

    results.append(subdf)

    # Save ancestry-specific output
    ancestry_output = f"summary_table_fdr_{ancestry}.tsv"
    subdf.to_csv(ancestry_output, sep="\t", index=False)

    print(f"Saved ancestry-specific table: {ancestry_output}")