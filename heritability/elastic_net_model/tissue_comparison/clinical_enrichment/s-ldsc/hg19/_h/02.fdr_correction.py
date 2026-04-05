import pandas as pd
from scipy.stats import norm
from statsmodels.stats.multitest import multipletests

# --- Step 1: Load your summary table ---
input_file = "combined_L2_1_Coeff_zscore.tsv"  # change this to your file path
df = pd.read_csv(input_file, sep="\t")

# --- Step 2: Convert Coefficient_z_score to two-tailed p-values ---
df['p_value'] = 2 * (1 - norm.cdf(df['Coefficient_z-score'].abs()))

# --- Step 3: Apply FDR correction (Benjamini-Hochberg) ---
df['FDR_q'] = multipletests(df['p_value'], method='fdr_bh')[1]

# --- Step 4: Save the new table ---
output_file = "summary_table_fdr.tsv"
df.to_csv(output_file, sep="\t", index=False)

print(f"FDR-corrected table saved to {output_file}")