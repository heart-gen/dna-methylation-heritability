"""
This script prepares meQTL results (localQTL) for mash modeling
by combining per-chromosome parquet files for selected brain regions.
"""
import session_info
import argparse, os
import pandas as pd
from glob import glob
from pathlib import Path

# Brain regions to include
BRAIN_REGIONS = ["caudate", "dlpfc", "hippocampus"]

def extract_shared_vmrs(shared_key_fn):
    """
    Extract shared VMRs across all brain regions.
    """
    return pd.read_csv(shared_key_fn, sep="\t", 
                       usecols=["shared_feature_id", "feature_id"])

def load_meqtl_parquet(region, shared, localqtl=False):
    """
    Load all chromosome-specific parquet files for a given brain region and feature.
    """
    base_path = f"../../{region}/local_analysis/_m/"
    if localqtl:
        prefix = "TOPMed_LIBD.haps"
    else:
        prefix = "TOPMed_LIBD"
    pattern = os.path.join(base_path, f"{prefix}.chr*.parquet")
    files = sorted(glob(pattern))

    if not files:
        raise FileNotFoundError(f"No parquet files found for {region} at {pattern}")

    dfs = []
    for f in files:
        df = pd.read_parquet(f, columns=["phenotype_id", "variant_id",
                                         "slope", "slope_se"])
        df - df.merge(shared, left_on = "phenotype_id", 
                      right_on = "phenotype_id", how = "inner")

        dfs.append(df)

    return pd.concat(dfs, ignore_index=True)

def extract_meqtls(shared, localqtl):
    """
    Extract meQTLs from each brain region.
    """
    data = {}
    for region in BRAIN_REGIONS:
        data[region] = load_meqtl_parquet(region, shared, localqtl)
    return data


def extract_dataframe(region_data, variable, label, outdir):
    """
    Combine the selected variable across regions into a single dataframe 
    for only shared VMRs.
    """
    dfs = []
    for region, df in region_data.items():
        renamed = df.loc[:, ["phenotype_id", "variant_id", variable]] \
                    .rename(columns={variable: region.capitalize()})
        dfs.append(renamed)

    # Merge on phenotype_id and variant_id
    result = dfs[0]
    for other in dfs[1:]:
        result = result.merge(other, on=["phenotype_id", "variant_id"])

    output_file = outdir / f"{label}_nominal_3regions_AA.txt.gz"
    result.to_csv(output_file, sep='\t', index=False)
    print(f"Saved: {output_file}")


def main():
    # Parser
    parser = argparse.ArgumentParser()
    parser.add_argument('--localqtl', action='store_true')
    parser.add_argument('--outdir', type=Path, default=Path("./"))
    args = parser.parse_args()

    # Extract data
    shared_data_fn = f"../../shared_vmrs/_m/TOPMed_LIBD_shared_vmr_key.tsv"
    shared = extract_shared_vmrs(shared_data_fn)
    region_data = extract_meqtls(shared, args.localqtl)

    extract_dataframe(region_data, "slope", "bhat", args.outdir)
    extract_dataframe(region_data, "slope_se", "shat", args.outdir)

    # Session information
    session_info.show()


if __name__ == '__main__':
    main()
