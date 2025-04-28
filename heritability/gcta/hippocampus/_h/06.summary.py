import pandas as pd
from pyhere import here
import os
import argparse 
import session_info

def process_hsq_file(input, chr, start_pos, end_pos, output):

                                        # read in .hsq file
    df = pd.read_csv(input, sep='\t', header=None, 
                     names=['source', 'variance', 'standard_error'])

                                        # exclude header names from output
    df = df[~df['source'].str.contains('^Source$', case=False, na=False)]

                                        # extract results with only one col
    one_col = {'logL', 'logL0', 'LRT', 'df', 'Pval', 'n'}
    df1 = df[~df['source'].isin(one_col)]
    df2 = df[df['source'].isin(one_col)]

    column_names = []
    extracted_data = []

                                        # extract results with both cols
    for _, values in df1.iterrows():
        variable_name = values['source']
        variance = f'Variance_{variable_name.replace("/", "_")}'
        se = f'SE_{variable_name.replace("/", "_")}'
        column_names.extend([variance, se])
        extracted_data.extend([values['variance'], values['standard_error']])

                                        # add one col data
    one_col_val = {}
    for var in one_col:
        match = df2[df2['source'] == var]
        one_col_val[var] = match.iloc[0, 1] if not match.empty else None
        column_names.append(var)
        extracted_data.append(one_col_val[var])

                                        # add chr, start pos, end pos 
                                        # as the first columns
    greml = [chr, start_pos, end_pos] + extracted_data
    greml_df = pd.DataFrame(
            [greml],
            columns=['Chr', 'Start_Pos', 'End_Pos'] + column_names
        )
                                        # write to csv
    out_hsq = f'greml_{start_pos}_{end_pos}.csv'
    greml_df.to_csv(os.path.join(output, out_hsq), sep='\t', index=False)

    return greml_df

def main():

                                        # define parser for chr, start_pos
                                        # and end_pos inputs
    parser = argparse.ArgumentParser(
        description='Write GREML results to csv for each region')
    parser.add_argument('--chr', type=int, required=True)
    parser.add_argument('--start_pos', type=int, required=True)
    parser.add_argument('--end_pos', type=int, required=True)
    args = parser.parse_args()
    
    input = here('heritability', 'gcta', 'hippocampus', '_m', 'h2', f'chr_{args.chr}', f'TOPMed_LIBD.AA.{args.start_pos}_{args.end_pos}.hsq')

                                        # define output directory
    output = here('heritability', 'gcta', 'hippocampus', '_m', 'summary', f'chr_{args.chr}')

                                        # create output directory if it 
                                        # doesn't exist
    os.makedirs(output, exist_ok=True)

    greml_df = process_hsq_file(input, args.chr, args.start_pos, 
                                 args.end_pos, output)

                                        # reproducibility information
    session_info.show()

if __name__ == '__main__':
    main()
