"""This script runs localQTL -- local ancestry aware."""
import session_info, pandas as pd
from localqtl.cis import CisMapper
from localqtl.haplotypeio import RFMixReader
from localqtl import PgenReader, read_phenotype_bed
from localqtl.cis.postproc import get_significant_pairs

def get_genotypes():
    plink_prefix_path = "/projects/b1213/users/alexis/projects/dna-methylation-heritability/meqtl-analysis/vmrs/dlpfc/_m/protected_data/TOPMed_LIBD"
    pgr = PgenReader(plink_prefix_path)
    variant_df = pgr.variant_df
    variant_df.loc[:, "chrom"] = "chr" + variant_df.chrom
    mask = ~variant_df.duplicated(subset=["chrom", "pos"], keep="first")
    return pgr.load_genotypes()[mask], variant_df[mask]


def get_covars(feature = "vmrs"):
    covar_file = f"/projects/b1213/users/alexis/projects/dna-methylation-heritability/meqtl-analysis/vmrs/dlpfc/_m/{feature}.combined_covariates.txt"
    return pd.read_csv(covar_file, sep='\t', index_col=0).T


def get_phenotype(feature = "vmrs"):
    meth_bed = f"/projects/b1213/users/alexis/projects/dna-methylation-heritability/meqtl-analysis/vmrs/dlpfc/_m/{feature}.methylation.bed.gz"
    return read_phenotype_bed(meth_bed)


def get_haplotypes(select_samples):
    """Make sure that samples match between haplotypes and genotypes."""
    rfmix_data_path = "/projects/b1213/resources/processed-data/" +\
        "local-ancestry/rfmix-version/_m/"
    binary_path = f"{rfmix_data_path}/binary_files"
    rfr = RFMixReader(prefix_path=rfmix_data_path,
                      binary_path=binary_path,
                      select_samples=select_samples)
    return rfr.load_haplotypes(), rfr.loci_df


def main():
    # Load data
    feature = "vmrs"
    covars_df = get_covars(feature)
    phenotype_df, phenotype_pos_df = get_phenotype(feature)
    genotype_df, variant_df = get_genotypes()
    cols = ["phenotype_id","variant_id", "slope", "slope_se",
            "pval_nominal", "start_distance","end_distance"]

    # Permutation without haplotypes (tensorQTL-style)
    prefix = "TOPMed_LIBD"
    mapper = CisMapper(
        genotype_df=genotype_df, variant_df=variant_df,
        phenotype_df=phenotype_df, phenotype_pos_df=phenotype_pos_df,
        covariates_df=covars_df, window=500_000, maf_threshold=0.01,
        device="auto", out_dir="./", out_prefix=prefix,
        tensorqtl_flavor=True
    )
    perm_df = mapper.map_permutations(nperm=1_000, beta_approx=True,
                                      seed=13131313)
    perm_df = mapper.calculate_qvalues(perm_df, fdr=0.05)
    perm_df.to_csv(f"{prefix}.permutation.txt.gz", sep='\t', index=False)

    signif_pairs = get_significant_pairs(
        res_df=perm_df, nominal_files=f"./{prefix}*.parquet",
        fdr=0.05, columns=cols # pre-select columns for speed
    )
    signif_pairs.to_csv(f"{prefix}.signif_variants.txt.gz", sep="\t",
                        index=False)

    # Permutation with haplotypes
    prefix = "TOPMed_LIBD.haps"
    select_samples = genotype_df.columns.values
    H, loci_df = get_haplotypes(select_samples)
    mapper = CisMapper(
        genotype_df=genotype_df, variant_df=variant_df,
        phenotype_df=phenotype_df, phenotype_pos_df=phenotype_pos_df,
        covariates_df=covars_df, haplotypes=H, loci_df=loci_df,
        window=500_000, maf_threshold=0.01, device="auto",
        out_dir="./", out_prefix=prefix,
    )
    perm_df = mapper.map_permutations(nperm=1_000, beta_approx=True,
                                      seed=13131313)
    perm_df = mapper.calculate_qvalues(perm_df, fdr=0.05)
    perm_df.to_csv(f"{prefix}.permutation.txt.gz", sep='\t',
                   index=False)

    signif_pairs = get_significant_pairs(
        res_df=perm_df, nominal_files=f"./{prefix}*.parquet",
        fdr=0.05, columns=cols # pre-select columns for speed
    )
    signif_pairs.to_csv(f"{prefix}.signif_variants.txt.gz", sep="\t",
                        index=False)

    # Session information
    session_info.show()


if __name__ == "__main__":
    main()
