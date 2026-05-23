import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

# =========================
# Load data
# =========================
df = pd.read_csv("summary_table_fdr_BA_only.tsv", sep="\t")

# =========================
# Simplify annotation labels
# ANNOT_Q1L2_1 -> Q1
# =========================
df["annotation_short"] = (
    df["annotation"]
    .str.extract(r'ANNOT_(Q\d)')[0]
)

# =========================
# Create row labels
# =========================
df["row_label"] = df["disease"] + "_" + df["tissue"]

# =========================
# Create pivot tables
# =========================
heatmap_data = df.pivot(
    index="row_label",
    columns="annotation_short",
    values="Coefficient_z-score"
)

# Ensure Q1-Q5 ordering
quintile_order = ["Q1", "Q2", "Q3", "Q4", "Q5"]
heatmap_data = heatmap_data[quintile_order]

pvals = df.pivot(
    index="row_label",
    columns="annotation_short",
    values="p_value"
)

pvals = pvals[quintile_order]

# =========================
# Convert p-values to stars
# =========================
def pval_to_star(p):
    if p < 0.001:
        return "***"
    elif p < 0.01:
        return "**"
    elif p < 0.05:
        return "*"
    else:
        return ""

stars = pvals.applymap(pval_to_star)

# =========================
# Custom red-blue colormap
# =========================
custom_cmap = LinearSegmentedColormap.from_list(
    "custom_red_blue",
    ["#2166AC", "white", "#B2182B"]
)

# =========================
# Plot heatmap
# =========================
plt.figure(figsize=(8, 8))

ax = sns.heatmap(
    heatmap_data,
    cmap=custom_cmap,
    center=0,
    linewidths=0.5,
    linecolor="black",
    annot=stars,
    fmt="",
    cbar_kws={"label": "Coefficient z-score"}
)

# Labels
plt.title("BA-only Coefficient Z-scores")
plt.xlabel("Quintile")
plt.ylabel("Disease_Tissue")

plt.xticks(rotation=0)
plt.yticks(rotation=0)

plt.tight_layout()

# =========================
# Save outputs
# =========================
plt.savefig("BA_only.pdf", bbox_inches="tight")
plt.savefig("BA_only.png", dpi=300, bbox_inches="tight")

plt.show()

print("Saved:")
print(" - BA_only.pdf")
print(" - BA_only.png")