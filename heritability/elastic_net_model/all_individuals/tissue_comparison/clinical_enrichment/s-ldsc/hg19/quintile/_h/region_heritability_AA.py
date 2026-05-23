import argparse
import session_info
import pandas as pd
from pathlib import Path
from pyliftover import LiftOver


def liftover_coordinate(lo, chrom, pos):
    try:
        new = lo.convert_coordinate(chrom, int(pos))
        if new:
            return int(new[0][1])
        return None
    except Exception:
        return None


def process_subset(df_subset, lo):
    df_subset["chrom"] = df_subset["chrom"].astype(str).apply(
        lambda x: f"chr{x}" if not str(x).startswith("chr") else x
    )
    df_subset["start_hg19"] = df_subset.apply(
        lambda x: liftover_coordinate(lo, x["chrom"], x["start"]), axis=1
    )
    df_subset["end_hg19"] = df_subset.apply(
        lambda x: liftover_coordinate(lo, x["chrom"], x["end"]), axis=1
    )
    df_subset.dropna(subset=["start_hg19", "end_hg19"], inplace=True)

    df_subset["start"] = df_subset["start_hg19"].astype(int)
    df_subset["end"] = df_subset["end_hg19"].astype(int)

    df_subset.drop(columns=["start_hg19", "end_hg19"], inplace=True)
    df_subset = df_subset[df_subset["start"] <= df_subset["end"]]

    return df_subset


def sort_by_genomic_position(df_subset):
    chrom_order = [f"chr{i}" for i in range(1, 23)]
    df_subset["chrom"] = pd.Categorical(
        df_subset["chrom"], categories=chrom_order, ordered=True
    )
    df_subset.sort_values(by=["chrom", "start"], inplace=True)
    return df_subset

import numpy as np

def create_quintiles(df, n_quintiles=5):
    # Compute quantile breaks
    breaks = np.quantile(
        df["h2_unscaled_AA"].dropna(),
        q=np.linspace(0, 1, n_quintiles + 1)
    )

    # Remove duplicate breaks (same behavior as R `unique()`)
    breaks = np.unique(breaks)
    n_bins = len(breaks) - 1

    if n_bins < 1:
        print("Warning: Insufficient unique h2 values to compute quantile bins")
        return pd.DataFrame()

    if n_bins < n_quintiles:
        print(f"Warning: Reduced h2 bins from {n_quintiles} to {n_bins}")

    # Assign bins
    df["h2_quintile_AA"] = pd.cut(
        df["h2_unscaled_AA"],
        bins=breaks,
        labels=[f"Q{i}" for i in range(1, n_bins + 1)],
        include_lowest=True
    )

    # Drop NA bins
    df = df.dropna(subset=["h2_quintile_AA"])

    return df

def main():
    parser = argparse.ArgumentParser(
        description="Create continuous annotation BED with liftover from hg38 to hg19"
    )
    parser.add_argument("--input_file", required=True, help="Path to input TSV file")
    parser.add_argument("--output_dir", required=True, help="Directory to save output files")
    parser.add_argument("--chain_file", required=True, help="Path to hg38ToHg19.over.chain.gz")
    args = parser.parse_args()

    input_file = Path(args.input_file)
    output_dir = Path(args.output_dir)
    chain_file = Path(args.chain_file)

    output_dir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(input_file, sep='\t')

    # Remove low prediction VMRs
    df = df[df["r_squared_cv_AA"] >= 0.3]

    # Keep only necessary columns (continuous annotation retained)
    cols_to_keep = ["chrom", "start", "end", "h2_unscaled_AA"]
    df = df[cols_to_keep].copy()

    lo = LiftOver(str(chain_file))

    df = process_subset(df, lo)

    df = create_quintiles(df)

    # Create separate BED files per quintile
    for quintile, subset in df.groupby("h2_quintile_AA"):
        # Convert category to string just in case
        quintile_str = str(quintile)

        output_file = output_dir / f"{quintile_str}.bed"

        # BED format: chrom, start, end
        subset[["chrom", "start", "end"]].to_csv(
            output_file,
            sep="\t",
            index=False,
            header=False
        )

        print(f"Saved {quintile_str}: {len(subset)} rows → {output_file}")

if __name__ == "__main__":
    main()