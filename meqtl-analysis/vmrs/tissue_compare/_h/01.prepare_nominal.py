"""
This script prepares meQTL results (localQTL) for mash modeling
by combining per-chromosome parquet files for selected brain regions.
"""
import session_info
import argparse, os
import pandas as pd
from glob import glob

# Brain regions to include
BRAIN_REGIONS = ["caudate", "dlpfc", "hippocampus"]

def load_meqtl_parquet(region, localqtl=False):
    """
    Load all chromosome-specific parquet files for a given brain region and feature.
    """
    base_path = f"../../{region}/local_analysis/_m/"
    if localqtl:
        prefix = "TOPMed_LIBD.haps"
    else:
        prefix = "TOPMed_LIBD"
    pattern = os.path.join(base_path, f"{prefix}.cis_qtl_pairs.chr*.parquet")
    files = sorted(glob(pattern))

    if not files:
        raise FileNotFoundError(f"No parquet files found for {region} at {pattern}")

    dfs = []
    for f in files:
        df = pd.read_parquet(f, columns=["phenotype_id", "variant_id",
                                         "slope", "slope_se"])
        dfs.append(df)

    return pd.concat(dfs, ignore_index=True)


def extract_meqtls(localqtl):
    """
    Extract meQTLs from each brain region.
    """
    data = {}
    for region in BRAIN_REGIONS:
        data[region] = load_meqtl_parquet(region, localqtl)
    return data


def extract_dataframe(region_data, variable, label):
    """
    Combine the selected variable across regions into a single dataframe.
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

    output_file = f"{label}_nominal_3regions_AA.txt.gz"
    result.to_csv(output_file, sep='\t', index=False)
    print(f"Saved: {output_file}")


def main():
    # Parser
    parser = argparse.ArgumentParser()
    parser.add_argument('--localqtl', action='store_true')
    args = parser.parse_args()

    # Extract data
    region_data = extract_meqtls(args.localqtl)
    extract_dataframe(region_data, "slope", "bhat")
    extract_dataframe(region_data, "slope_se", "shat")

    # Session information
    session_info.show()

if __name__ == '__main__':
    main()
