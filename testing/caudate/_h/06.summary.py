import pandas as pd
from pyhere import here
import os
import argparse 
import session_info

def process_hsq_file(input, chr, start_pos, end_pos, output):
    # Read in .hsq file
    df = pd.read_csv(input, sep='\t', header=None, names=['source', 'variance', 'standard_error'])
    extracted_data = df.iloc[1:5, -2:]
    transposed_data = extracted_data.values.flatten()

    # Add Chromosome, start position, and end position as the first columns
    greml = [chr, start_pos, end_pos] + list(transposed_data)
    greml_df = pd.DataFrame([greml], columns=['Chr', 'Start_Pos', 'End_Pos', 'Variance_Vg', 'Standard_Error_Vg', 
                                                'Variance_Ve', 'Standard_Error_Ve',
                                                'Variance_Vp', 'Standard_Error_Vp',
                                                'Variance_Vg_Vp', 'Standard_Error_Vg_Vp'])
    greml_df.to_csv(output, sep='\t', index=False)

    return greml_df

def main():
    parser = argparse.ArgumentParser(description='Write GREML results to csv for all regions')
    parser.add_argument('--chr', required=True)
    parser.add_argument('--start_pos', required=True)
    parser.add_argument('--end_pos', required=True)
    args = parser.parse_args()
    
    input = here('testing', 'caudate', '_m', 'h2')
    output = here('testing', 'caudate', '_m', 'summary')

    if not os.path.exists('./summary'):
        os.makedirs('./summary')

    result_df = process_hsq_file(input, args.chr, args.start_pos, args.end_pos, output)
    session_info.show()

if __name__ == '__main__':
    main()