import argparse
from pathlib import Path
import session_info
import pandas as pd
import numpy as np

def filter_vmr(enet):
    vmr = enet.dropna().copy()

    vmr.loc[(vmr["h2_unscaled"] >= 0.1) & (vmr["r_squared_cv"] > 0.3), "h2_category"] = "heritable"
    vmr.loc[(vmr["h2_unscaled"] < 0.1) & (vmr["r_squared_cv"] > 0.3), "h2_category"] = "non-heritable"
    vmr.loc[(vmr["r_squared_cv"] <= 0.3), "h2_category"] = "low_prediction"
    
    return vmr
    
def sort_by_genomic_position(df):
    return df.sort_values(by=["chrom", "start"], ascending=[True, True])

def main():
    parser = argparse.ArgumentParser(
        description="Separate TSV into groups with genomic sorting and liftover from hg38 to hg19"
    )
    parser.add_argument("--input_file", required=True, help="Path to input TSV file")
    parser.add_argument("--output_dir", required=True, help="Directory to save output files")
    parser.add_argument("--region", required=True, help="Brain region")
    args = parser.parse_args()

    input_file = Path(args.input_file)
    output_dir = Path(args.output_dir)
    region = args.region

    output_dir.mkdir(parents=True, exist_ok=True)

    enet = pd.read_csv(input_file, sep='\t')
    vmr_all = filter_vmr(enet)
    h2_categories = ["heritable", "non-heritable", "low_prediction"]

    for h2_cat in h2_categories:
        cols_to_keep = ["chrom", "start", "end"]
        vmr = vmr_all[vmr_all['h2_category'] == h2_cat]
        vmr = vmr[cols_to_keep].copy()
        vmr = sort_by_genomic_position(vmr)
        out_fn = output_dir / f"{region}_{h2_cat}.bed"
        vmr.to_csv(out_fn, sep='\t', index=False, header=False)

    session_info.show()

if __name__ == "__main__":
    main()
