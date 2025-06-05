import argparse
import numpy as np
import pandas as pd
import session_info
from tqdm import tqdm
from pathlib import Path
from scipy.stats import norm
from random import seed, shuffle
from scipy.linalg import toeplitz


def get_chromosome_size(chrom):
    """
    Given chromosomes. Return chromosome size.
    """
    fn = "/projects/b1213/resources/genomes/human/gencode-v47/fasta/" +\
        "chromosome_sizes.txt"
    chrom_sizes = pd.read_csv(fn, sep="\t", header=None, names=["chr", "pos"])
    chrom_sizes = chrom_sizes[(chrom_sizes["chr"].str.startswith("chr"))].copy()
    chrom_size_dict = chrom_sizes.set_index("chr").to_dict()['pos']
    return chrom_size_dict[f"chr{chrom}"]


def generate_phenotypes(num_pheno, num_chroms, snp_window):
    """
    Generates phenotype definitions with causal SNPs and heritability.

    Parameters:
        num_pheno (int): Total phenotypes to generate
        num_chroms (int): Number of chromosomes (1-22)
        snp_window (int): Window size around causal SNPs (±bp)

    Returns:
        tuple: (phenotype definitions, simulation regions)
    """
    pheno_h2_flags = [True]*int(num_pheno * 0.25) + \
        [False]*int(num_pheno * 0.75)
    shuffle(pheno_h2_flags)
    # Store info about each phenotype's causal SNPs and H2
    pheno_dict = []; all_regions = {}
    # Assign causal SNPs to phenotypes
    print(" Assigning causal SNPs to phenotypes...")
    for pheno_idx in range(num_pheno):
        is_heritable = pheno_h2_flags[pheno_idx]
        h2 = np.random.uniform(0.1, 0.8) if is_heritable else np.random.uniform(0.0, 0.099)
        num_causal = np.random.randint(1, 6)
        causal_snps = []
        for _ in range(num_causal):
            chrom = np.random.randint(1, num_chroms + 1)
            max_bp_pos = get_chromosome_size(chrom)
            # Ensure causal_pos allows for a full window
            causal_pos = np.random.randint(snp_window + 1, max_bp_pos - snp_window)
            # Define the region for this specific causal SNP
            region_key = (chrom, causal_pos)
            if region_key not in all_regions:
                all_regions[region_key] = (
                    max(1, causal_pos - snp_window), # ensure start >=1
                    min(max_bp_pos, causal_pos + snp_window)
                )
                causal_snps.append({'chrom': chrom, 'pos': causal_pos})
        pheno_dict.append({
            'pheno_id': f"pheno_{pheno_idx+1}",
            'causal_snps_info': causal_snps,
            'heritability': h2
        })
    return pheno_dict, all_regions


def simulate_ld_genotypes(num_snps, num_samples, ld_decay):
    """
    Simulates genotypes for each SNP using the inverse CDF of the normal
    distribution, ensuring Hardy-Weinberg equilibrium genotype frequencies.

    Parameters:
        num_snps (int): Number of SNPs to simulate.
        num_samples (int): Number of individuals/samples.
        ld_decay (float): Exponential decay factor for LD.

    Returns:
        np.ndarray: Genotype matrix (num_snps x num_samples) with values 0, 1, or 2.
    """
    # Create LD correlation matrix
    ld_vector = ld_decay ** np.arange(num_snps)
    corr_matrix = toeplitz(ld_vector).astype(np.float32)
    try:
        L = np.linalg.cholesky(corr_matrix + 1e-5 * np.eye(num_snps, dtype=np.float32))
    except np.linalg.LinAlgError:
        L = np.eye(num_snps, dtype=np.float32)
    raw = np.random.normal(0, 1, size=(num_snps, num_samples)).astype(np.float32)
    correlated = L @ raw
    # Simulate minor allele frequencies (MAFs)
    mafs = np.random.uniform(0.05, 0.5, size=num_snps).astype(np.float32)
    genotypes = np.zeros((num_snps, num_samples), dtype=np.int8)
    for i, maf in enumerate(mafs):
        # Hardy-Weinberg genotype frequencies: (1-p)^2, 2p(1-p), p^2
        p = maf; q = 1 - p
        # Cumulative probabilities for AA, Aa, aa
        cum_prob_AA = q**2
        cum_prob_Aa = q**2 + 2*p*q
        # Set normal distribution thresholds for AA, Aa, aa
        threshold_AA = norm.ppf(cum_prob_AA)
        threshold_Aa = norm.ppf(cum_prob_Aa)
        # Assign genotypes
        geno = correlated[i, :]
        genotypes[i, :] = np.where(geno <= threshold_AA, 0,
                                   np.where(geno <= threshold_Aa, 1, 2))
    return genotypes


def simulate_genotypes(simu_regions, min_snps, max_snps, num_samples, ld_decay):
    """
    Simulates genotypes with LD and creates PLINK BIM records.

    Parameters:
        simu_regions (dict): Genomic regions to simulate
        min_snps (int): Minimum SNPs per region
        max_snps (int): Maximum SNPs per region
        num_samples (int): Number of individuals
        ld_decay (float): LD decay rate

    Returns:
        tuple: (genotype matrix, SNP location map, BIM records)
    """
    bim_records = []; simu_snp_loc = {}; geno_lt = []; id_counter = 1
    sorted_keys = sorted(simu_regions.keys())
    print(f" Simulating genotypes for {len(sorted_keys)} unique causal SNP regions...")
    for region_idx, (chrom, causal_pos) in enumerate(tqdm(sorted_keys, desc="Simulating regions", unit="region")):
        region_start, region_end = simu_regions[(chrom, causal_pos)]
        num_snps = np.random.randint(min_snps, max_snps + 1)
        region_span = region_end - region_start + 1
        if region_span < num_snps:
            snp_locs = np.arange(region_start, region_end + 1)
        else:
            causal = np.array([causal_pos])
            candidate_locs = np.arange(region_start, region_end + 1)
            candidate_locs = candidate_locs[candidate_locs != causal_pos]
            if len(candidate_locs) < num_snps - 1:
                chosen_locs = candidate_locs
            else:
                chosen_locs = np.random.choice(candidate_locs, num_snps - 1, replace=False)
            snp_locs = np.sort(np.append(chosen_locs, causal))
        snps_in_region = len(snp_locs)
        if snps_in_region == 0:
            continue
        region_genotypes = simulate_ld_genotypes(snps_in_region, num_samples, ld_decay)
        if region_genotypes.shape[0] == 0:
            continue
        # Fast batch appending
        for i, pos in enumerate(snp_locs):
            snp_id = f"rs{id_counter}"
            bim_records.append((chrom, snp_id, 0, pos, 'A', 'G'))
            simu_snp_loc[(chrom, pos)] = (snp_id, len(geno_lt))
            geno_lt.append(region_genotypes[i])
            id_counter += 1
    geno_mat = np.vstack(geno_lt) if geno_lt else np.empty((0, num_samples), dtype=np.int8)
    return geno_mat, simu_snp_loc, bim_records


def cal_pheno(pheno_dict, simu_snp_loc, geno_mat, num_samples, num_pheno):
    """
    Standardize genetic_value_for_pheno to have variance = h2
    This makes E(Var(G)) = h2 if E(Var(P)) = 1.
    Let P = G + E. Var(P) = Var(G) + Var(E). h2 = Var(G) / Var(P).
    If we set Var(G) = h2, then Var(E) = (1-h2). Then Var(P) = 1.
    """
    output_mat = np.zeros((num_samples, num_pheno))
    map_records = [] # For the TSV file
    print(" Calculating phenotype values...")
    for pheno_idx, pheno_def in enumerate(pheno_dict):
        h2 = pheno_def['heritability']
        pheno_id = pheno_def['pheno_id']
        genetic_val = np.zeros(num_samples)
        causal_snp = []
        for causal_snp_info in pheno_def['causal_snps_info']:
            chrom = causal_snp_info['chrom']
            pos = causal_snp_info['pos']
            if (chrom, pos) in simu_snp_loc:
                rs_id, geno_idx = simu_snp_loc[(chrom, pos)]
                geno_casual = geno_mat[geno_idx, :]
                effect_size = np.random.normal(0, 0.5) # Scaled
                genetic_val += geno_casual * effect_size
                causal_snp.append({'chrom': chrom, 'pos': pos, 'rs_id': rs_id})
        # Scale genetic component and add environmental noise
        if len(pheno_def['causal_snps_info']) > 0 and np.var(genetic_val) > 1e-6:
            genetic_var = np.var(genetic_val)
            if genetic_var > 0:
                genetic_val_scaled = (genetic_val - np.mean(genetic_val)) / np.sqrt(genetic_var) * np.sqrt(h2)
            else:
                genetic_val_scaled = np.zeros(num_samples)
            environ_var_scale = np.sqrt(1 - h2) if (1-h2) > 0 else 0
            environ_comp = np.random.normal(0, environ_var_scale, num_samples)
            output_mat[:, pheno_idx] = genetic_val_scaled + environ_comp
        else:
            environ_comp = np.random.normal(0, 1, num_samples)
            output_mat[:, pheno_idx] = environ_comp
            if h2 > 0 and not causal_snp :
                pheno_def['heritability'] = 0.0
        # Store mapping info
        if causal_snp:
            chrom_lt = [str(cs['chrom']) for cs in causal_snp]
            pos_lt = [str(cs['pos']) for cs in causal_snp]
            rsid_lt = [cs['rs_id'] for cs in causal_snp]
            map_records.append({
                "phenotype_id": pheno_id, "num_causal_snps": len(causal_snp),
                "chroms": ",".join(chrom_lt), "positions": ",".join(pos_lt),
                "rsids": ",".join(rsid_lt), "target_heritability": h2,
                "simulated_heritability": np.var(genetic_val_scaled) if len(pheno_def['causal_snps_info']) > 0 and np.var(genetic_val) > 1e-6 and 'genetic_val_scaled' in locals() else 0
            })
        else:
            map_records.append({
                "phenotype_id": pheno_id, "num_causal_snps": 0,
                "chroms": "", "positions": "", "rsids": "",
                "target_heritability": h2, "simulated_heritability": 0.0
            })
    return map_records, output_mat


def save_phen(output_mat, fam_df, num_pheno, outdir):
    pheno_file_path = outdir / "simulated.phen"
    with open(pheno_file_path, 'w') as f:
        # Write header
        f.write("FID\tIID\t" + "\t".join(f"pheno_{i+1}" for i in range(num_pheno)) + "\n")
        # Write data row by row
        for i in range(output_mat.shape[0]):
            fid = fam_df["FID"].iloc[i]
            iid = fam_df["IID"].iloc[i]
            row_data = "\t".join(map(str, output_mat[i, :]))
            f.write(f"{fid}\t{iid}\t{row_data}\n")
    print(f"  Saved Phenotype file: {pheno_file_path}")


def _save_fam(fam_data, plink_dir):
    fam_file_path = plink_dir / "simulated.fam"
    with open(fam_file_path, 'w') as f:
        for row in fam_data:
            f.write("\t".join(map(str, row)) + "\n")
    print(f"  Saved FAM file: {fam_file_path}")
    # Return DataFrame if needed downstream
    return pd.DataFrame(fam_data, columns=["FID", "IID", "FatherID", "MotherID", "Sex", "Phenotype"])


def _save_bim(bim_records, plink_dir):
    bim_file_path = plink_dir / "simulated.bim"
    if bim_records:
        with open(bim_file_path, 'w') as f:
            for row in bim_records:
                f.write("\t".join(map(str, row)) + "\n")
        print(f"  Saved BIM file: {bim_file_path}")
    else:
        print("  Warning: No SNP records generated for BIM file.")


def _save_bed(geno_mat, plink_dir):
    bed_file_path = plink_dir / "simulated.bed"
    if geno_mat.size > 0: # Check if there's data
        num_total_snps = geno_mat.shape[0]
        num_total_samples = geno_mat.shape[1]
        with open(bed_file_path, 'wb') as f_bed:
            # Write PLINK magic numbers: 0x6c (l), 0x1b (k), 0x01 (SNP-major mode)
            f_bed.write(bytes([0x6c, 0x1b, 0x01]))
            # Iterate over each SNP (row in geno_mat)
            for snp_idx in range(num_total_snps):
                byte_val = 0
                bit_pos_in_byte = 0 # from 0 to 7 (or 0 to 3 pairs)
                for sample_idx in range(num_total_samples):
                    genotype = geno_mat[snp_idx, sample_idx] # This is 0, 1, or 2
                    plink_geno_code = 0b00 # Default to AA (00)
                    if genotype == 1: # AG
                        plink_geno_code = 0b10 # Heterozygous
                    elif genotype == 2: # GG
                        plink_geno_code = 0b11 # Homozygous for A2
                    byte_val |= (plink_geno_code << (bit_pos_in_byte))
                    bit_pos_in_byte += 2
                    if bit_pos_in_byte == 8: # Byte is full
                        f_bed.write(bytes([byte_val]))
                        byte_val = 0
                        bit_pos_in_byte = 0
                # Write any remaining part of the last byte for this SNP
                if bit_pos_in_byte > 0:
                    f_bed.write(bytes([byte_val]))
        print(f"  Saved BED file: {bed_file_path}")
    else:
        print("  Warning: No genotype data generated for BED file. BED file not written.")


def save_plink(fam_data, bim_records, geno_mat, plink_dir):
    # FAM file (.fam)
    fam_df = _save_fam(fam_data, plink_dir)
    # BIM File (.bim)
    _save_bim(bim_records, plink_dir)
    # BED File (.bed) - Binary Genotype File
    _save_bed(geno_mat, plink_dir)
    print(f"PLINK files are in: {PLINK_OUTPUT_DIR}")
    return fam_df


def main(NUM_PHENOTYPES, NUM_SAMPLES, NUM_CHROMOSOMES, SNP_WINDOW_SIZE,
         MIN_SNPS, MAX_SNPS, LD_DECAY, OUTPUT_DIR, PLINK_OUTPUT_DIR):
    # Global Setup & Covariates
    print("Step I: Global Setup & Covariates")
    sex_lt = np.random.choice([1, 2], NUM_SAMPLES) # PLINK typically 1=male, 2=female
    age_lt = np.random.normal(50, 10, NUM_SAMPLES)
    batch_lt = np.random.randint(1, 4, NUM_SAMPLES) # Simulating 3 batches

    fam_data = []
    for i in range(NUM_SAMPLES):
        fam_data.append([f"ID{i+1}", f"ID{i+1}", "0", "0", sex_lt[i], -9])

    # Define Causal Loci and Simulation Regions
    print("Step II: Define Causal Loci and Simulation Regions")
    pheno_dict, all_regions = generate_phenotypes(
        NUM_PHENOTYPES, NUM_CHROMOSOMES, SNP_WINDOW_SIZE)

    # Genotype Simulation (PLINK .bim, .bed content)
    print("Step III: Genotype Simulation")
    geno_mat, simu_snp_loc, bim_records = simulate_genotypes(
        all_regions, MIN_SNPS, MAX_SNPS, NUM_SAMPLES, LD_DECAY)

    # Phenotype Value Calculation
    print("Step IV: Phenotype Value Calculation")
    map_records, output_mat = cal_pheno(pheno_dict, simu_snp_loc, geno_mat,
                                        NUM_SAMPLES, NUM_PHENOTYPES)

    # Save Output Files
    print("Step V: Saving Output Files")
    fam_df = save_plink(fam_data, bim_records, geno_mat, PLINK_OUTPUT_DIR)
    save_phen(output_mat, fam_df, NUM_PHENOTYPES, OUTPUT_DIR)

    # SNP to Phenotype Mapping File (.tsv)
    if map_records:
        mapping_df = pd.DataFrame(map_records)
        mapping_file_path = OUTPUT_DIR / "snp_phenotype_mapping.tsv"
        mapping_df.to_csv(mapping_file_path, sep="\t", index=False, na_rep="NA")
        print(f"  Saved SNP-Phenotype mapping file: {mapping_file_path}")
    else:
        print("  Warning: No SNP-Phenotype mapping records generated.")

    # Categorical Covariates (.covar)
    covar_df = pd.DataFrame({"FID": fam_df["FID"], "IID": fam_df["IID"],
                             "sex": sex_lt, "batch": batch_lt})
    covar_file_path = OUTPUT_DIR / "simulated.covar"
    covar_df.to_csv(covar_file_path, sep="\t", index=False, header=True)
    print(f"  Saved Categorical Covariate file: {covar_file_path}")
    # Quantitative Covariates (.qcovar)
    qcovar_df = pd.DataFrame({"FID": fam_df["FID"], "IID": fam_df["IID"],
                              "age": age_lt})
    qcovar_file_path = OUTPUT_DIR / "simulated.qcovar"
    qcovar_df.to_csv(qcovar_file_path, sep="\t", index=False, header=True)

    print(f"  Saved Quantitative Covariate file: {qcovar_file_path}")
    print("\nSimulation complete.")
    print(f"Output files are in: {OUTPUT_DIR}")

    # Example of how many SNPs were simulated in total
    if 'geno_mat' in locals() and geno_mat.size > 0:
        print(f"Total SNPs simulated: {geno_mat.shape[0]}")
    elif bim_records:
        print(f"Total SNPs simulated (from bim_records): {len(bim_records)}")
    else:
        print("Total SNPs simulated: 0")


def parse_args():
    parser = argparse.ArgumentParser(
        description="Simulate phenotypes and genotypes with heritability and LD.")
    parser.add_argument("--num_phenotypes", type=int, default=1000,
                        help="Number of phenotypes to simulate")
    parser.add_argument("--num_samples", type=int, default=500,
                        help="Number of individuals")
    parser.add_argument("--num_chromosomes", type=int, default=22,
                        help="Number of chromosomes (default 22)")
    parser.add_argument("--snp_window", type=int, default=500_000,
                        help="Window size around causal SNP (±500kb)")
    parser.add_argument("--min_snps", type=int, default=150,
                        help="Minimum SNPs per window")
    parser.add_argument("--max_snps", type=int, default=5000,
                        help="Maximum SNPs per window")
    parser.add_argument("--ld_decay", type=float, default=0.8,
                        help="LD decay rate (e.g. 0.8)")
    parser.add_argument("--output_dir", type=str, default="./simulation_output",
                        help="Output directory path")
    parser.add_argument("--seed", type=int, default=13,
                        help="Random seed")
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()

    # Constants from arguments
    NUM_PHENOTYPES = args.num_phenotypes
    NUM_SAMPLES = args.num_samples
    NUM_CHROMOSOMES = args.num_chromosomes
    SNP_WINDOW_SIZE = args.snp_window
    MIN_SNPS = args.min_snps
    MAX_SNPS = args.max_snps
    LD_DECAY = args.ld_decay

    # Set output directories
    OUTPUT_DIR = Path(args.output_dir)
    OUTPUT_DIR.mkdir(exist_ok=True)
    PLINK_OUTPUT_DIR = OUTPUT_DIR / "plink_sim"
    PLINK_OUTPUT_DIR.mkdir(exist_ok=True)

    # Set seeds
    np.random.seed(args.seed)
    seed(args.seed)

    # Call main
    main(NUM_PHENOTYPES, NUM_SAMPLES, NUM_CHROMOSOMES, SNP_WINDOW_SIZE,
         MIN_SNPS, MAX_SNPS, LD_DECAY, OUTPUT_DIR, PLINK_OUTPUT_DIR)

    # Session information
    session_info.show()
