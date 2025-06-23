import session_info
import pandas as pd
from pyhere import here
import os
import re
from pathlib import Path

def process_hsq_file(filename):
    greml_df = pd.read_csv(filename, sep="\t")\
                 .melt(id_vars="Source", value_vars=["Variance", "SE"],
                       var_name="Stat", value_name="Value")\
                 .pivot_table(index=None, columns=["Source", "Stat"],
                              values="Value").dropna(axis=1)
    greml_df.columns = [f"{src.replace('/', '_')}_{stat}" for src, stat in greml_df.columns]
    return greml_df

def main():
    ## Main
    basepath = here("simulation-analysis", "gcta", "_m", "h2")
    output_path = here("simulation-analysis", "gcta", "_m", "summary")
    output_path.mkdir(exist_ok=True)

    dt = pd.DataFrame()
    pattern_dir = re.compile(r"sim_(\d+)_indiv")
    pattern_file = re.compile(r"greml_pheno(\d+)\.hsq")
    for fn in Path(basepath).glob("*/*hsq"):
        greml_df = process_hsq_file(fn)

        # extract sample size from filename
        match_dir = pattern_dir.search(str(fn))
        if match_dir:
            n_indiv = int(match_dir.group(1))
            greml_df["N"] = n_indiv
        else:
            print(f"Number of individuals not found in filename: {fn}")
            continue

        # extract pheno number from filename
        match_file = pattern_file.search(fn.name)
        if match_file:
            pheno_str = f"pheno_{match_file.group(1)}"
            greml_df["pheno"] = pheno_str
        else:
            print(f"Pheno number not found in filename: {fn}")
            continue



        # Move sample size to the beginning of df
        cols = ["N", "pheno"] + [col for col in greml_df.columns if col not in ["N", "pheno"]]
        greml_df = greml_df[cols]

        # Add to data frame and sort
        dt = pd.concat([dt, greml_df], axis=0)

    dt["pheno_num"] = dt["pheno"].str.extract(r"pheno_(\d+)", expand=False).astype(int)
    dt = dt.sort_values(by="pheno_num").drop(columns="pheno_num")
    
    ## Save
    dt.to_csv(os.path.join(output_path, "greml_summary.tsv"), sep='\t', index=False)

    ## Session information
    session_info.show()

if __name__ == '__main__':
    main()
