"""This script runs localQTL."""
import session_info, pandas as pd
from localqtl.cis import CisMapper
from localqtl.haplotypeio import RFMixReader
from localqtl import PgenReader, read_phenotype_bed

def get_genotypes():
    plink_prefix_path = "../../_m/protected_data/TOPMed_LIBD"
    pgr = PgenReader(plink_prefix_path)
    variant_df = pgr.variant_df
    variant_df.loc[:, "chrom"] = "chr" + variant_df.chrom
    mask = ~variant_df.duplicated(subset=["chrom", "pos"], keep="first")
    return pgr.load_genotypes()[mask], variant_df[mask]


def get_covars(feature = "genes"):
    covar_file = f"../../_m/{feature}.combined_covariates.txt"
    return pd.read_csv(covar_file, sep='\t', index_col=0).T


def get_phenotype(feature = "genes"):
    expr_bed = f"../../_m/{feature}.expression.bed.gz"
    return read_phenotype_bed(expr_bed)


def get_haplotypes(select_samples):
    """Make sure that samples match between haplotypes and genotypes."""
    rfmix_data_path = "/projects/b1213/resources/processed-data/" +\
        "local-ancestry/rfmix-version/_m/"
    # rfmix_data_path = "/ocean/projects/bio250020p/shared/resources/" +\
    #     "processed-data/local-ancestry/rfmix-version/_m/"
    binary_path = f"{rfmix_data_path}/binary_files"
    rfr = RFMixReader(prefix_path=rfmix_data_path,
                      binary_path=binary_path,
                      select_samples=select_samples)
    return rfr.load_haplotypes(), rfr.loci_df


def main():
    # Load data
    feature = "genes"
    covars_df = get_covars(feature)
    phenotype_df, phenotype_pos_df = get_phenotype(feature)
    genotype_df, variant_df = get_genotypes()

    # Nominal without haplotypes (tensorQTL-style)
    prefix = "TOPMed_LIBD"
    mapper = CisMapper(
        genotype_df=genotype_df, variant_df=variant_df,
        phenotype_df=phenotype_df, phenotype_pos_df=phenotype_pos_df,
        covariates_df=covars_df,  window=500_000, maf_threshold=0.01,
        device="auto", out_dir="./", out_prefix=prefix,
        tensorqtl_flavor=True,
    )
    mapper.map_nominal()

    # Nominal with haplotypes
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
    mapper.map_nominal()

    # Session information
    session_info.show()


if __name__ == "__main__":
    main()
