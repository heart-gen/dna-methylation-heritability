"""
LDSC Partitioned Heritability Analysis and Visualization

This module combines LDSC results, filters for significant enrichment,
aggregates categories, and generates heatmap visualizations.
"""
import os
import session_info
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from statsmodels.stats.multitest import multipletests
from scipy.stats import norm


def combine_ldsc_files(file_list_path, output_path="combined_ldsc_results.tsv"):
    """
    Combine only the first data row (after header) from each LDSC result file
    into one long-format table.
    """
    with open(file_list_path, "r") as f:
        file_paths = f.read().strip().split()

    dfs = []
    for fpath in file_paths:
        basename = os.path.basename(fpath)
        name_no_ext = os.path.splitext(basename)[0]
        parts = name_no_ext.split("_")

        disorder = parts[0]
        brain_region = parts[1]
        if parts[-1] == "hg19":
            h2_metric = "_".join(parts[2:-1])
        else:
            h2_metric = "_".join(parts[2:])

        # ✅ Read only first row after header
        df = pd.read_csv(fpath, sep="\t", nrows=1)

        df["Disorder"] = disorder
        df["BrainRegion"] = brain_region
        df["h2Metric"] = h2_metric
        dfs.append(df)

    combined_df = pd.concat(dfs, ignore_index=True)

    cols = ["Disorder", "BrainRegion", "h2Metric"] + [
        c for c in combined_df.columns
        if c not in ["Disorder", "BrainRegion", "h2Metric"]
    ]
    combined_df = combined_df[cols]
    combined_df.to_csv(output_path, sep="\t", index=False)

    print("Combined LDSC results saved! Shape:", combined_df.shape)
    return combined_df


def filter_significant_enrichment(
    input_path, output_path="ldsc_significant_enrichment.tsv",
    fdr_threshold=0.05, exclude_categories=None
):
    """
    Filter LDSC results for significant enrichment.

    Parameters
    ----------
    input_path : str
        Path to combined LDSC results.
    output_path : str
        Path to save filtered results.
    fdr_threshold : float
        FDR threshold for significance.
    exclude_categories : list, optional
        Categories to exclude from results.

    Returns
    -------
    pd.DataFrame
        Filtered significant enrichment results.
    """

    df = pd.read_csv(input_path, sep="\t")
    print("Total rows:", df.shape[0])
    print("Rows with valid Enrichment_p:", df['Enrichment_p'].notna().sum())

    df_valid = df[df["Enrichment_p"].notna()].copy()
    df_valid["FDR"] = df_valid.groupby(
        ["Disorder", "BrainRegion", "h2Metric"]
    )["Enrichment_p"].transform(
        lambda x: multipletests(x, method="fdr_bh")[1]
    )

    sig_df = df_valid[(df_valid["FDR"] < fdr_threshold) & (df_valid["Enrichment"] > 1)]
    sig_df.to_csv(output_path, sep="\t", index=False)
    print("Rows passing filter:", sig_df.shape[0])

    filtered_df = sig_df[~sig_df['Category'].isin(exclude_categories)]
    return filtered_df


def filter_significant_tau(
    input_path, output_path="ldsc_significant_tau.tsv",
    fdr_threshold=0.05, exclude_categories=None
):
    """
    Filter LDSC results for significant tau (Coefficient) effects using
    two-sided p-values derived from Coefficient_z-score.
    """
    if exclude_categories is None:
        exclude_categories = ["MAF_Adj_ASMCL2_0"]

    df = pd.read_csv(input_path, sep="\t")
    print("Total rows:", df.shape[0])
    print("Rows with valid Coefficient_z-score:", df['Coefficient_z-score'].notna().sum())

    df_valid = df[df["Coefficient_z-score"].notna()].copy()
    df_valid["Tau_p"] = 2 * norm.sf(df_valid["Coefficient_z-score"].abs())
    df_valid["Tau_FDR"] = df_valid.groupby(
        ["Disorder", "BrainRegion", "h2Metric"]
    )["Tau_p"].transform(
        lambda x: multipletests(x, method="fdr_bh")[1]
    )

    sig_df = df_valid[df_valid["Tau_FDR"] < fdr_threshold]
    sig_df.to_csv(output_path, sep="\t", index=False)
    print("Rows passing tau filter:", sig_df.shape[0])

    filtered_df = sig_df[~sig_df['Category'].isin(exclude_categories)]
    return filtered_df


def create_summary_table(filtered_df, output_path="ldsc_summary_table.tsv",
                         category_col="Category"):
    """
    Aggregate and enumerate top/shared categories.

    Parameters
    ----------
    filtered_df : pd.DataFrame
        Filtered significant enrichment data.
    output_path : str
        Path to save summary table.
    category_col : str
        Name of the category column.

    Returns
    -------
    pd.DataFrame
        Summary table with category statistics.
    """
    disorders_per_category = filtered_df.groupby(category_col)["Disorder"].nunique().reset_index()
    disorders_per_category.rename(columns={"Disorder": "Disorders_Enriched"}, inplace=True)

    regions_per_category = filtered_df.groupby(category_col)["BrainRegion"].nunique().reset_index()
    regions_per_category.rename(columns={"BrainRegion": "BrainRegions_Enriched"}, inplace=True)

    avg_enrichment = filtered_df.groupby(category_col)["Enrichment"].mean().reset_index()
    avg_enrichment.rename(columns={"Enrichment": "Average_Enrichment"}, inplace=True)

    summary_df = disorders_per_category.merge(regions_per_category, on=category_col)
    summary_df = summary_df.merge(avg_enrichment, on=category_col)
    summary_df = summary_df.sort_values(
        ["Disorders_Enriched", "Average_Enrichment"],
        ascending=[False, False]
    )

    summary_df.to_csv(output_path, sep="\t", index=False)
    print("Summary table saved! Preview:")
    print(summary_df.head(10))

    return summary_df


def fdr_to_stars(fdr):
    """Convert FDR value to significance asterisks."""
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


def make_heatmap(
    filtered_df, output_path="ldsc_enrichment_heatmap_multiindex_asterisks.png",
    top_n=20, figsize=(16, 10)
):
    """
    Create enrichment heatmap with significance asterisks.

    Parameters
    ----------
    filtered_df : pd.DataFrame
        Filtered significant enrichment data.
    output_path : str
        Path to save heatmap image.
    top_n : int
        Number of top categories to include.
    figsize : tuple
        Figure size (width, height).
    """
    top_categories = (
        filtered_df.groupby('Category')['Enrichment']
        .mean()
        .sort_values(ascending=False)
        .head(top_n)
        .index
    )
    heatmap_df = filtered_df[filtered_df['Category'].isin(top_categories)]
    heatmap_matrix = heatmap_df.pivot_table(
        index='Category',
        columns=['Disorder', 'BrainRegion', 'h2Metric'],
        values='Enrichment',
        fill_value=0
    )
    heatmap_matrix.columns = ['_'.join(col) for col in heatmap_matrix.columns]

    annot_matrix = heatmap_df.pivot_table(
        index='Category',
        columns=['Disorder', 'BrainRegion', 'h2Metric'],
        values='FDR',
        fill_value=1
    )
    annot_matrix.columns = ['_'.join(col) for col in annot_matrix.columns]
    annot_matrix = annot_matrix.map(fdr_to_stars)

    plt.figure(figsize=figsize)
    sns.heatmap(
        heatmap_matrix, cmap='Purples', linewidths=0.5,
        annot=annot_matrix, fmt='', cbar_kws={'label': 'Enrichment'}
    )
    plt.title("LDSC Partitioned Heritability Enrichment Heatmap")
    plt.xlabel("Disorder_BrainRegion_h2Metric")
    plt.ylabel("Functional Category")
    plt.tight_layout()

    plt.savefig(output_path, dpi=300)
    plt.show()
    print(f"Heatmap saved to {output_path}")


def main():
    """Run the full LDSC analysis and visualization pipeline."""
    # Combine all LDSC files
    combined_df = combine_ldsc_files("ldsc_file_list.txt")

    # Filter for significant enrichment
    filtered_df = filter_significant_enrichment("combined_ldsc_results.tsv")

    # Filter for significant tau (Coefficient) effects
    filtered_tau_df = filter_significant_tau("combined_ldsc_results.tsv")

    # Create summary table
    summary_df = create_summary_table(filtered_df)

    # Generate heatmap
    make_heatmap(filtered_df)

    # Session information
    session_info.show()


if __name__ == "__main__":
    main()
