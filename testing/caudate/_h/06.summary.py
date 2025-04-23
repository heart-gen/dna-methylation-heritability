import pandas as pd
from pyhere import here
import os
import argparse 
import session_info

def process_hsq_file(input, chr, start_pos, end_pos, output):
    # Read in .hsq file
    df = pd.read_csv(input, sep='\t', header=None, 
                     names=['source', 'variance', 'standard_error'])
    extracted_data = pd.concat([
        df.iloc[1:12, -2:],
        df.iloc[16:17, [1]]
    ])
    transposed_data = extracted_data.values.flatten()

    # Add Chromosome, start position, and end position as the first columns
    greml = [chr, start_pos, end_pos] + list(transposed_data)
    greml_df = pd.DataFrame(
            [greml],
            columns=[
                'Chr', 'Start_Pos', 'End_Pos',
                'Variance_Vg1', 'SE_Vg1',
                'Variance_Vg2', 'SE_Vg2',
                'Variance_Vg3', 'SE_Vg3',
                'Variance_Vg4', 'SE_Vg4',
                'Variance_Ve', 'SE_Ve',
                'Variance_Vp', 'SE_Vp',
                'Variance_Vg1_Vp', 'SE_Vg1_Vp',
                'Variance_Vg2_Vp', 'SE_Vg2_Vp',
                'Variance_Vg3_Vp', 'SE_Vg3_Vp',
                'Variance_Vg4_Vp', 'SE_Vg4_Vp',
                'Variance_Sum_Vg_Vp', 'SE_Sum_Vg_Vp',
                'Pval', 'n'
            ]
        )
    out_hsq = f'greml_{start_pos}_{end_pos}.csv'
    greml_df.to_csv(os.path.join(output, out_hsq), sep='\t', index=False)

    return greml_df

def main():
    parser = argparse.ArgumentParser(
        description='Write GREML results to csv for each region')
    parser.add_argument('--chr', required=True)
    parser.add_argument('--start_pos', required=True)
    parser.add_argument('--end_pos', required=True)
    args = parser.parse_args()
    
    input = here('testing', 'caudate', '_m', 'h2', f'chr_{args.chr}', f'TOPMed_LIBD.AA.{args.start_pos}_{args.end_pos}.hsq')
    output = here('testing', 'caudate', '_m', 'summary', f'chr_{args.chr}')
    os.makedirs(output, exist_ok=True)

    greml_df = process_hsq_file(input, args.chr, args.start_pos, 
                                 args.end_pos, output)
    session_info.show()

if __name__ == '__main__':
    main()