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
    df_subset.drop(df_subset[df_subset["start"] > df_subset["end"]].index, inplace=True)
    return df_subset


def sort_by_genomic_position(df_subset):
    chrom_order = [f"chr{i}" for i in range(1, 23)]
    df_subset["chrom"] = pd.Categorical(
        df_subset["chrom"], categories=chrom_order, ordered=True
    )
    df_subset.sort_values(by=["chrom", "start"], ascending=[True, True], inplace=True)
    return df_subset


def main():
    parser = argparse.ArgumentParser(
        description="Separate TSV into groups with genomic sorting and liftover from hg38 to hg19"
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

    heritable = df[(df["h2_unscaled"] >= 0.3) & (df["r_squared_cv"] >= 0.75)]
    non_heritable = df[(df["h2_unscaled"] < 0.3) & (df["r_squared_cv"] >= 0.75)]
    low_prediction = df[df["r_squared_cv"] < 0.75]

    cols_to_keep = ["chrom", "start", "end", "h2_unscaled"]
    heritable = heritable[cols_to_keep].copy()
    non_heritable = non_heritable[cols_to_keep].copy()
    low_prediction = low_prediction[cols_to_keep].copy()

    lo = LiftOver(str(chain_file))

    heritable = process_subset(heritable, lo)
    non_heritable = process_subset(non_heritable, lo)
    low_prediction = process_subset(low_prediction, lo)

    heritable = sort_by_genomic_position(heritable)
    non_heritable = sort_by_genomic_position(non_heritable)
    low_prediction = sort_by_genomic_position(low_prediction)

    heritable_bed = output_dir / "heritable_hg19.bed"
    non_heritable_bed = output_dir / "non_heritable_hg19.bed"
    low_prediction_bed = output_dir / "low_prediction_hg19.bed"

    heritable.to_csv(heritable_bed, sep='\t', index=False, header=False)
    non_heritable.to_csv(non_heritable_bed, sep='\t', index=False, header=False)
    low_prediction.to_csv(low_prediction_bed, sep='\t', index=False, header=False)

    print("Separation and liftover complete! Files saved without headers.")
    print(f"Heritable group: {len(heritable)} rows → {heritable_bed}")
    print(f"Non-heritable group: {len(non_heritable)} rows → {non_heritable_bed}")
    print(f"Low prediction group: {len(low_prediction)} rows → {low_prediction_bed}")

    session_info.show()


if __name__ == "__main__":
    main()
