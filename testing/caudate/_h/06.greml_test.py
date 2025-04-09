import pandas as pd

def process_hsq_file(input_file_path, chr, start_pos, end_pos, output_file_path):
    # Read in .hsq file
    df = pd.read_csv(input_file_path, sep='\t', header=None, names=['source', 'variance', 'standard_error'])
    extracted_data = df.iloc[1:5, -2:]
    transposed_data = extracted_data.values.flatten()

    # Add Chromosome, start position, and end position as the first columns
    greml = [chr, start_pos, end_pos] + list(transposed_data)
    greml_df = pd.DataFrame([greml], columns=['Chr', 'Start_Pos', 'End_Pos', 'Variance_Vg', 'Standard_Error_Vg', 
                                                'Variance_Ve', 'Standard_Error_Ve',
                                                'Variance_Vp', 'Standard_Error_Vp',
                                                'Variance_Vg_Vp', 'Standard_Error_Vg_Vp'])
    greml_df.to_csv(output_file_path, sep='\t', index=False)

    return greml_df

# Test run
input_file_path = '/projects/p32505/projects/dna-methylation-heritability/testing/_m/TOPMed_LIBD.AA.VMR1.hsq'
output_file_path = '/projects/p32505/projects/dna-methylation-heritability/testing/_m/TOPMed_LIBD.AA.VMR1.txt'
chr = 'chr1'
start_pos = 403969
end_pos = 1404084

result_df = process_hsq_file(input_file_path, chr, start_pos, end_pos, output_file_path)
print(result_df)
