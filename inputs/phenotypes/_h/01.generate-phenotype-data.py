import session_info
import polars as pl
from pyhere import here


def load_methyl_pheno():
    # Replace to load general information
    fn = "../_m/phenotypes-DNAm.tsv"
    return pl.read_csv(fn, separator="\t")


def load_snp_pcs():
    # Data with genotypes
    fn = here("inputs/genotypes/TOPMed_LIBD.AA.eigenvec")
    df = pl.read_csv(fn, separator="\t")
    rename_dict = {
        col: "snpPC" + col[2:]
        for col in df.columns if col.startswith("PC")
    }
    return df.rename(rename_dict)\
             .select(pl.exclude("IID"))


def load_protected_data():
    # Load detailed medical history (select)
    cols = ["brnumerical","smoking","codeine", "morphine", "cocaine", "ethanol",
            "amphetamines", "education", "marital_status", "antipsychotics",
            "past_suicide_attempts", "psychosis", "hallucinations", "delusions",
            "fentanyl", "ph", "hx_military_service", "past_self_mutilation",
            "lifetime_antipsych", "hx_other_trauma", "hx_physical_abuse",
            "hx_sexual_abuse", "nicotine", "cerad", "braak", "age_onset_schizo",
            "fsiq", "manner_of_death"]
    fn = "../_h/protected-data/phenotype.csv"
    return pl.read_csv(fn, infer_schema_length=5000)\
             .select(cols)\
             .with_columns([
                 pl.format("Br{}", pl.col("brnumerical")).alias("BrNum")
             ])


def main():
    # Main
    df = load_methyl_pheno()\
        .join(load_snp_pcs(), left_on="brnum", right_on="#FID", how="left")\
        .join(load_protected_data(), left_on="brnum", right_on="BrNum")
    print(df.group_by(["region", "primarydx", "race"])\
          .agg([pl.len().alias("n")]).sort("region"))

    # Save data
    df.write_csv("phenotypes-AA.tsv", separator="\t")

    # Session information
    session_info.show()


if __name__ == "__main__":
    main()
