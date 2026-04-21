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

    # Keep only necessary columns (continuous annotation retained)
    cols_to_keep = ["chrom", "start", "end", "h2_unscaled"]
    df = df[cols_to_keep].copy()

    lo = LiftOver(str(chain_file))

    # Process full dataset (no partitioning)
    df = process_subset(df, lo)
    df = sort_by_genomic_position(df)

    output_bed = output_dir / "continuous_annotation_hg19.bed"

    df.to_csv(output_bed, sep='\t', index=False, header=False)

    print("Liftover complete! Continuous annotation file saved.")
    print(f"Total rows: {len(df)} → {output_bed}")

    session_info.show()


if __name__ == "__main__":
    main()