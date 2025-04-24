import glob
import pandas as pd
import os
from pyhere import here
import session_info

def merge_csv(input_dir, output_dir):

                                        # recursively find all .csv files
    csv_files = glob.glob(os.path.join(input_dir, "**", "*csv"), recursive=True)

                                        # read each csv into list
    df = []
    for file in csv_files:
        df_tmp = pd.read_csv(file, sep='\t')

                                        # merge all files into df
    combined_csv = pd.concat(df, ignore_index=True)

                                        # write to file
    combined_csv.to_csv(os.path.join(output_dir, 'greml_summary.csv'), 
                        sep=',', index=False)

input_dir = here('testing', 'caudate', '_m', 'summary')
output_dir = here('testing', 'caudate', '_m', 'summary')
merge_csv(input_dir, output_dir)

session_info.show()