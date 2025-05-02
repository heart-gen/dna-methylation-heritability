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
    basepath = here("heritability", "gcta", "caudate", "_m", "h2")
    output_path = here("heritability", "gcta", "caudate", "_m", "summary")
    output_path.mkdir(exist_ok=True)

    dt = pd.DataFrame()
    pattern = re.compile(r"chr_(\d+)/TOPMed_LIBD\.AA\.(\d+)_(\d+)\.hsq")
    for fn in Path(basepath).glob("*/*hsq"):
        greml_df = process_hsq_file(fn)

        # extract chromosome and position from filename
        match = pattern.search(str(fn))
        if match:
            chrom, start, end = match.groups()
            greml_df["chr"] = chrom
            greml_df["start"] = start
            greml_df["end"] = end
        else:
            print(f"Filename pattern not matched: {fn}")
            continue

        # Move chr, start, end to the beginning of df
        cols = ["chr", "start", "end"] + [col for col in greml_df.columns if col not in ["chr", "start", "end"]]
        greml_df = greml_df[cols]

        # Add to data frame and sort
        dt = pd.concat([dt, greml_df], axis=0)
        dt = dt.astype({"chr": int}).sort_values(by=["chr", "start", "end"])
    
    ## Save
    dt.to_csv(os.path.join(output_path, "greml_summary.tsv"), sep='\t', index=False)

    ## Session information
    session_info.show()

if __name__ == '__main__':
    main()
