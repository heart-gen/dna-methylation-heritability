import os
import session_info
import pandas as pd
from pyhere import here
from pathlib import Path

from genboostgpu.vmr_runner import run_single_window
from genboostgpu.hyperparams import enet_from_targets
from genboostgpu.orchestration import run_windows_with_dask
from genboostgpu.snp_processing import count_snps_in_window
from genboostgpu.data_io import load_genotypes, load_phenotypes
from genboostgpu.tuning import select_tuning_windows, global_tune_params

def construct_data_path(num_samples, dtype="phen"):
    base_dir = Path(here("inputs/simulated-data/_m")) / f"sim_{num_samples}_indiv"
    if dtype.lower() == "plink":
        return base_dir / "plink_sim" / "simulated.bed"
    elif dtype.lower() == "phen":
        return base_dir / "simulated.phen"
    else:
        raise ValueError(f"Unknown dtype: {dtype}")


def get_pheno_loc(num_samples):
    base_dir = Path(here("inputs/simulated-data/_m")) / f"sim_{num_samples}_indiv"
    mapping_file = os.path.join(base_dir, "snp_phenotype_mapping.tsv")

    if not os.path.exists(mapping_file):
        raise FileNotFoundError(f"Mapping file not found: {mapping_file}")

    mapped_df = pd.read_csv(mapping_file, sep="\t", header=0)

    return mapped_df

def build_windows(num_samples):
    # load genotype + bim/fam
    # geno_path = construct_data_path(num_samples, "plink")

    # Load phenotypes
    pheno_path = construct_data_path(num_samples, "phen")

    # Build windows config list
    pheno_loc_df = get_pheno_loc(num_samples)
    windows = []
    for _, row in pheno_loc_df[0:10].iterrows(): ## Test the first 10
        windows.append({
            "chrom": row["chrom"],
            "start": row["start"],
            "end": row["end"],
            "pheno_id": row["phenotype_id"],
            "pheno_path": pheno_path,
        })
    return windows


def tune_windows(num_samples):
    geno_path = construct_data_path(num_samples, "plink")
    geno_arr, bim, fam = load_genotypes(str(geno_path))
    N = len(fam)
    windows = build_windows(num_samples)
    tw = select_tuning_windows(
        windows, bim, frac=0.05, n_min=60, n_max=300, window_size=500_000,
        use_window=False, n_bins=3, per_chrom_min=1, seed=13
    )
    # Tune alpha scale only
    best = global_tune_params(
        tuning_windows=tw,
        geno_arr=geno_arr, bim=bim, fam=fam,
        window_size=500_000, by_hand=False, use_window=False,
        
    # Limit grid to alpha-scale:
    grid={
        "c_lambda":       [0.5, 0.7, 1.0, 1.4, 2.0],
        "c_ridge":        [1.0],      # keep EN balance fixed here
        "subsample_frac": [0.7],      # keep fixed for now
        "batch_size":     [4096],     # keep fixed for now
    },
    early_stop={"patience": 5, "min_delta": 1e-4, "warmup": 5},
    batch_size=4096
    )
    print("Best alpha scale:", best)

    return best


def main():
    num_samples = os.environ.get("NUM_SAMPLES")

    if not num_samples:
        raise ValueError("NUM_SAMPLES environment variable must be set")

    # Run with dask orchestration
    df = run_windows_with_dask(
        windows, outdir="results", window_size=500_000,
        n_iter=100, n_trials=10, use_window=False,
        save=True, prefix="simu_100"
    )
    print(f"Completed {len(df)} VMRs")
    print(df.head())

    # Session information
    session_info.show()

if __name__ == "__main__":
    main()
