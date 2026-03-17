"""
This script calculates the percent overlap of VMRs
across brain regions
"""
import pandas as pd
from pathlib import Path
import session_info

def get_bed(bed_file):
    cols = ["chromA", "startA", "endA", "chromB", "startB", "endB", "overlap"]
    return pd.read_csv(bed_file, sep="\t", header=None, names=cols)

def annot_overlap(bed_df):
    # Percent overlap for each VMR
    bed_df["lenA"] = bed_df.endA - bed_df.startA
    bed_df["lenB"] = bed_df.endB - bed_df.startB
    bed_df["pctA"] = (bed_df.overlap / bed_df.lenA) * 100
    bed_df["pctB"] = (bed_df.overlap / bed_df.lenB) * 100

    # Reciprocal overlap
    bed_df["reciprocal"] = bed_df[["pctA", "pctB"]].min(axis=1)
    return bed_df

def main():
    for threshold_dir in Path(".").glob("*_0.*"):
        input_dir = threshold_dir / "sets"
        out_dir = threshold_dir / "percent_overlap"

        out_dir.mkdir(exist_ok=True)

        for bed_file in input_dir.glob("*overlap*.bed"):
            if "3tissues" in bed_file.name:
                print(f"{bed_file} is not pairwise. Skipping.")
                continue

            output_fn = out_dir / bed_file.with_suffix(".tsv").name
            
            # Read in bed 
            df = get_bed(bed_file)

            # Calculate percent overlap
            df = annot_overlap(df)

            # Write dataframe
            df.to_csv(output_fn, sep="\t", index=False)

    # Session information
    session_info.show()

if __name__ == "__main__":
    main()