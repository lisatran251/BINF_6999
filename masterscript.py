import pandas as pd
import numpy as np 
import time
import argparse
import os
import re
import sys
import subprocess 
import traceback
from Bio import SeqIO, Entrez 
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# How to run:      
# python3 -m venv ~/my_venv
# source ~/my_venv/bin/activate
# pip install pandas numpy biopython argparse
# python3 masterscript.py fasta_file primer_file 

# Set the display.max_rows option to print all rows 
pd.set_option('display.max_rows', None)

def extract_products(input_file, primer_file):
    df = pd.read_csv(primer_file)

    # Create lists of unique forward and reverse primers 
    F_primers = list(set(df['F_truseq']))
    R_primers = list(set(df['R_truseq']))

    # Create lists of reverse complements of each primer 
    F_reverse_complement = [str(Seq(primer).reverse_complement()) for primer in F_primers]
    R_reverse_complement = [str(Seq(primer).reverse_complement()) for primer in R_primers]

    # Start time tracking 
    start_time = time.time()

    # Create empty dictionary to store results 
    record_sequences = {}

    # Open the input file and parse it as a fasta file 
    with open(input_file, 'r') as handle:
        records = SeqIO.parse(handle, 'fasta')

        # Loop through ea record 
        for record in records:
            sequence = str(record.seq)
            # Create a list to store the position of ea primer in the sequence
            primers_positions = [] 

            # For ea primer, find all positions in the sequence where the primer can be found 
            for primer in F_primers + R_primers + F_reverse_complement + R_reverse_complement:
                for match in re.finditer(primer, sequence):
                    pos = match.start() # Start position of the match
                    # Append a tuple with the primer and its position to list 
                    primers_positions.append((primer, pos))

            # If there are at least 2 primer positions found 
            if len(primers_positions) >= 2:
                record_id = record.id   # Find the ID of the record
                # Sort the list nof primers and positions by postions  
                primers_positions.sort(key=lambda x: x[1])

                # Create a list to store info for ea seq that was found 
                sequence_info = []

                # For ea pair of primers 
                for i in range(len(primers_positions) - 1):
                    for j in range(i + 1, len(primers_positions)):
                        # Extract the start and end primer and their pos
                        start_primer, start_pos = primers_positions[i]
                        end_primer, end_pos = primers_positions[j]
                        
                        # Extract the product seq btw the start and end primers
                        product = sequence[start_pos:end_pos + len(end_primer)]
                        
                        # Caculate legnth of product 
                        length = len(product)
                        combination = f"{start_primer}-{end_primer}"
                        
                        # Find source of ea primer 
                        start_primer_source = None
                        end_primer_source = None

                        if start_primer in F_primers:
                            start_primer_source = "F_primers"
                        elif start_primer in R_primers:
                            start_primer_source = "R_primers"
                        elif start_primer in F_reverse_complement:
                            start_primer_source = "F_reverse_complement"
                        elif start_primer in R_reverse_complement:
                            start_primer_source = "R_reverse_complement"

                        if end_primer in F_primers:
                            end_primer_source = "F_primers"
                        elif end_primer in R_primers:
                            end_primer_source = "R_primers"
                        elif end_primer in F_reverse_complement:
                            end_primer_source = "F_reverse_complement"
                        elif end_primer in R_reverse_complement:
                            end_primer_source = "R_reverse_complement"

                        info = {
                            'Start Primer': start_primer,
                            'End Primer': end_primer,
                            'Start Position': start_pos,
                            'End Position': end_pos + len(end_primer),
                            'Length': length,
                            'Combination': f"{start_primer_source}-{end_primer_source}",
                            'Product': product,
                        }
                        sequence_info.append(info)

                # Add the info for ea record to the record_sequences dictionary
                record_sequences[record_id] = sequence_info

    # End the time tracking 
    end_time = time.time()

    # Print out the total execution time 
    print("Total Execution Time:", end_time - start_time)

    return record_sequences

def main(input_file, primer_file):
    # Extract products using the input file and primer file 
    record_sequences = extract_products(input_file, primer_file)
    
    # Create a list to store the result dictionaries 
    results = []

    # For ea record in the record sequences dictionary 
    for record_id, sequence_info in record_sequences.items():
        for info in sequence_info:
            # Create result dictionary and add it to the results list 
            result = {
                "Record ID": record_id,
                "Start Primer": info["Start Primer"],
                "End Primer": info["End Primer"],
                "Start Position": info["Start Position"],
                "End Position": info["End Position"],
                "Length": info["Length"],
                "Combination": info["Combination"],
                "Product": info["Product"],
            }
            results.append(result)

    # Convert the results list to df 
    df = pd.DataFrame(results)

    # If result file already exists, append results to it, otherwise create a new file with header 
    if os.path.isfile('raw_results.csv'):
        df.to_csv('raw_results.csv', mode='a', header=False, index=False)
    else:
        df.to_csv('raw_results.csv', index=False)

# Call main functon 
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('input_file', type=str)
    parser.add_argument('primer_file', type=str)
    args = parser.parse_args()
    input_file = args.input_file
    primer_file = args.primer_file
    main(input_file, primer_file)

# Function to get the reverse complement of a sequence
def reverse_complement(seq):
    complement_dict = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return "".join(complement_dict[base] for base in reversed(seq))

primers = pd.read_csv(primer_file)
results = pd.read_csv('raw_results.csv')

# Assuming df is your DataFrame
primers['target_gene'] = primers['target_gene'].str.replace('^cfrB', 'cfr_B', regex=True)

# Fill 'F_primer' and 'R_primer' columns based on the value of 'Combination' column
for i, row in results.iterrows():
    if row['Combination'] == 'F_primers-R_primers':
        results.at[i, 'F_primer'] = row['Start Primer']
        results.at[i, 'R_primer'] = row['End Primer']
        results.at[i, 'F_Product'] = row['Product']

    elif row['Combination'] == 'R_primers-F_primers':
        results.at[i, 'F_primer'] = row['End Primer']
        results.at[i, 'R_primer'] = row['Start Primer']
        results.at[i, 'F_Product'] = row['Product'][::-1]
    
    elif row['Combination'] == 'R_reverse_complement-F_reverse_complement':
        results.at[i, 'F_primer'] = reverse_complement(row['End Primer'])
        results.at[i, 'R_primer'] = reverse_complement(row['Start Primer'])
        results.at[i, 'F_Product'] = reverse_complement(row['Product'][::-1])

    elif row['Combination'] == 'F_reverse_complement-R_reverse_complement':
        results.at[i, 'F_primer'] = reverse_complement(row['Start Primer'])
        results.at[i, 'R_primer'] = reverse_complement(row['End Primer'])
        results.at[i, 'F_Product'] = reverse_complement(row['Product'])
        
# Create 'Target' column with default value 'N/A'
results['Target'] = 'N/A'

# Create 'Qualified?' column with default value 'N/A'
results['Qualified?'] = 'No'

# Specify the valid combinations
qualified_combinations = ['F_primers-R_primers', 'R_primers-F_primers', 'R_reverse_complement-F_reverse_complement', 
'F_reverse_complement-R_reverse_complement']

# Filter results
results = results[results['Combination'].isin(qualified_combinations)]

# Set 'Target' column if 'F_primer' and 'R_primer' match with 'F_truseq' and 'R_truseq' in the primers df, respectively
for i, row in results.iterrows():
    matched_rows = primers[(primers['F_truseq'] == row['F_primer']) & (primers['R_truseq'] == row['R_primer'])]
    if len(matched_rows) > 0:
        results.at[i, 'Target'] = matched_rows.iloc[0]['target_gene']

# Set 'Qualified?' column based on product length and 'Target' gene
results.loc[(results['Length'] >= 150) & (results['Length'] <= 400) & (results['Target'] != 'N/A'), 'Qualified?'] = 'Yes'
results.loc[results['Qualified?'] != 'Yes', 'Qualified?'] = 'No'

# Print out final result
results.to_csv('final_result.csv', index=False)

# Split the df into two based on 'Qualified?' column
match = results[results['Qualified?'] == 'Yes']
nonmatch = results[results['Qualified?'] == 'No']

# Save to separate csv files without the 'Qualified?' column, are these neccessary?
match.drop(columns=['Qualified?']).to_csv('match.csv', index=False)
nonmatch.drop(columns=['Qualified?']).to_csv('nonmatch.csv', index=False)

# Set the display.max_rows option to print all rows 
pd.set_option('display.max_rows', None)

# Provide email to NCBI, maybe turn this into input value later
Entrez.email = "thuyduye@uoguelph.ca"

# Define the function that used to run multiple Shell script 
def run_command(command):
    # Run the command 
    result = subprocess.run(command, capture_output=True, text=True)
    # Check the return code to determine if the command was successful or not 
    if result.returncode != 0:
        # If the return code is non-zero, command failed
        print(f"Command '{' '.join(command)}' failed with error code {result.returncode}, stderr: {result.stderr}")
    else:
        # If the return code is 0, command succeeded 
        print(result.stdout)

# Run abricate: abricate --db resfinder --quiet contigs_ex.fasta > abricate_results.csv
run_command(['./run_abricate.sh'])

# Load the data
results = pd.read_csv('final_result.csv')
ab_results = pd.read_csv('abricate_results.csv', delimiter='\t')

# Create a df for the final_result
filtered_result = results[['Record ID', 'Length', 'Combination', 'F_Product', 'F_primer', 'R_primer', 'Target']]
filtered_result = filtered_result.dropna(subset=['F_Product'])
filtered_result = filtered_result.dropna(subset=['Target'])

# Create a length column for abricate data
ab_results['LENGTH'] = ab_results['END'] - ab_results['START']

# Create a df for abricate_result
ab_results = ab_results[['SEQUENCE', 'LENGTH', 'STRAND', 'GENE', '%COVERAGE', '%IDENTITY', 'ACCESSION', 'PRODUCT']]

# Define a function to fetch sequences
def fetch_sequence(accession):
    try:
        handle = Entrez.efetch(db='nucleotide', id=accession, rettype='fasta', retmode='text')
        record = SeqIO.read(handle, 'fasta')
        return str(record.seq)
    except Exception as e:
        print(f"Could not fetch the sequence for accession {accession}: {e}")
        # traceback.print_exc()
        return None
    
# Apply the function to fetch sequences and add them to the df 
ab_results['originalSEQUENCE'] = ab_results['ACCESSION'].apply(fetch_sequence)

# Define a function that remove the last part from the primers 
def shorten_target(s):
    if ':' in s:
        return s
    s = s.rsplit('_', 1)[0]
    if not s[-1].isdigit():
        s = s.rsplit('_', 1)[0]
    return s

# Use the function to create the new column with shorter name 
filtered_result['shortTarget'] = filtered_result['Target'].apply(shorten_target)

# Function to create 'shortProduct'
def create_short_product(x):
    if ':' in x:
        return x.split(':')[0].split('_')[0]
    elif re.search(r'_\d+$', x):  # If ends with underscore and digits
        return re.sub(r'_\d+$', '', x)
    elif re.search(r'_-\w+_\d+$', x):  # If has pattern '_-word_digit'
        return re.sub(r'_\d+$', '', x)
    elif re.search(r'\d+$', x):  # If ends with digits
        return re.sub(r'\d+$', '', x)
    else:
        return x

# Create 'shortProduct' column
filtered_result['shortProduct'] = filtered_result['shortTarget'].apply(create_short_product)

# Create a copy of primer df  
filtered_primer = primers.copy()

# Create a new column that convert the name of the primer to shorter name 
filtered_primer['shortTarget'] = filtered_primer['target_gene'].apply(shorten_target)

# Define a function to replace every "'" with "_"
def replace_chars(s):
    return re.sub(r"[\(''\)]", "_", s)

# Apply the function to each string in the 'GENE' column
ab_results['prgGENE'] = ab_results['GENE'].apply(replace_chars)

# Apply the function to each string in the 'PRODUCT' column 
ab_results['prgPRODUCT'] = ab_results['PRODUCT'].apply(replace_chars)

# Create a df that merges the abricate results and original results
merged_program_ab = filtered_result.merge(ab_results, how='outer', left_on=['Record ID','shortProduct'], right_on=['SEQUENCE','prgPRODUCT'], indicator=True)

# Create 'Gene Match?' column
merged_program_ab['Gene Match?'] = merged_program_ab.apply(lambda row: 'Yes' if row['shortProduct'] == row['prgPRODUCT'] else 'No', axis=1)

# Create 'Variant Match?' column
def check_variant_match(row):
    if row['Gene Match?'] == 'Yes':
        # Extract the trailing number from 'shortTarget' and 'prgGENE'
        short_target_number = re.search(r'\d+$', row['shortTarget'])
        prg_gene_number = re.search(r'\d+$', row['prgGENE'])
        
        if short_target_number is not None and prg_gene_number is not None and short_target_number.group() == prg_gene_number.group():
            return 'Yes'
        else:
            return 'No'
    else:
        return 'No'

merged_program_ab['Variant Match?'] = merged_program_ab.apply(check_variant_match, axis=1)
# Final report after merging
merged_program_ab.to_csv('merged_prog_ab.csv', index=False)

# Rows in program not in abricate
program_only = merged_program_ab[merged_program_ab['_merge'] == 'left_only']
program_only = program_only.iloc[:, :9]
program_only.to_csv('program_only.csv', index=False)

# Rows in abricate not in program
abricate_only = merged_program_ab[merged_program_ab['_merge'] == 'right_only']
abricate_only = abricate_only.iloc[:, 9:20]
abricate_only.to_csv('abricate_only.csv', index=False)

# Rows found in both progsram ands abricate
both_program_ab_result = merged_program_ab[merged_program_ab['_merge'] == 'both']
both_program_ab_result.drop('_merge', axis=1).to_csv('both_program_ab_result.csv', index=False)

## Scenario 1: target gene found - explore variants, artifacts etc. 
# Select rows from merged_program_ab where 'Gene Match?' = 'Yes' and 'Variant Match?' = 'No'
filtered_merged_program_ab = merged_program_ab[(merged_program_ab['Gene Match?'] == 'Yes') & (merged_program_ab['Variant Match?'] == 'No')]
filtered_merged_program_ab.to_csv('filtered_merged_program_ab.csv', index=False)

# Define a function to get the sequence from the original contigs 
def get_sequence(record_id):    
    with open(input_file, 'r') as f:
        for record in SeqIO.parse(f, 'fasta'):
            if record.id == record_id:
                return str(record.seq)
    return None 

# Create a new column in df to contain the contigs, can create a copy to avoid warning (opt)
filtered_merged_program_ab['Contigs'] = filtered_merged_program_ab['Record ID'].apply(get_sequence)

# Write a fasta include Record ID and Contigs
with open('variants_originalCONTIG.fasta', 'w') as f:
    for index, row in filtered_merged_program_ab.iterrows():
        f.write(f">seq_{index}_originalCONTIGS\n{row['Contigs']}\n")

# Write a fasta include Record ID and original sequence from Abricate program
with  open('variants_originalSEQUENCE.fasta', 'w') as f:
    for index, row in filtered_merged_program_ab.iterrows():
        f.write(f">seq_{index}_originalSEQUENCE\n{row['originalSEQUENCE']}\n")

## Run cdhit in Unix: cd-hit-est -i variants.fasta -o cdhit_rs -c 0.9 -n 10 -d 0
# The main reason is the product is very short while the sequence achieved using accession number is much longer. Therefore, when run cd-hit which based on the similarity between two sequences and the identiy % much be higher than 70%, therefore we may want to run BLAST instead.  

# ## Scenario 2: found in program not in abricate- re-run with lower %identity on abricate
# Create a set of all record IDs from the 'Record ID' column 
record_ids = program_only['Record ID'].tolist()

fasta_seqs = SeqIO.parse(open(input_file), 'fasta')

# Create new fasta file that contains seqs only found by the program
with open("program_only.fasta", "w") as out_file:
    for fasta in fasta_seqs:
        name, sequence = fasta.id, str(fasta.seq)
        if name in record_ids:
            SeqIO.write(fasta, out_file, "fasta")

# Re-run Abricate with lower coverage: abricate --db resfinder --quiet --mincov 10 program_only.fasta > program_only_ab10_results.csv
run_command(['./run_abricate_10.sh'])

# Load files
program_only_ab10_results = pd.read_csv('program_only_ab10_results.csv', delimiter='\t')

# Apply the function to each string in the 'GENE' column
program_only_ab10_results['prgGENE'] = program_only_ab10_results['GENE'].apply(replace_chars)

# Apply the function to each string in the 'PRODUCT' column 
program_only_ab10_results['prgPRODUCT'] = program_only_ab10_results['PRODUCT'].apply(replace_chars)

# Merge to find similar seqs between the ones program found and abricate (10% coverage) found
merged_program_only_ab10_results = program_only.merge(program_only_ab10_results, how='outer', left_on=['Record ID','shortProduct'], right_on=['SEQUENCE','prgPRODUCT'], indicator=True)

# Create 'Gene Match?' column
merged_program_only_ab10_results['Gene Match?'] = merged_program_only_ab10_results.apply(lambda row: 'Yes' if row['shortProduct'] == row['prgPRODUCT'] else 'No', axis=1)

# Create 'Variant Match?' column 
merged_program_only_ab10_results['Variant Match?'] = merged_program_only_ab10_results.apply(check_variant_match, axis=1)

# Export the results to CSV file 
merged_program_only_ab10_results.to_csv('merged_program_only_ab10_both_id_and_gene.csv', index=False)

# Print out df includes matches 
merged_program_only_ab10_match = merged_program_only_ab10_results[merged_program_only_ab10_results['_merge'] == 'both']
merged_program_only_ab10_match.to_csv('merged_program_only_ab10_match.csv', index=False)

# Print out df includes non-matches 
merged_program_only_ab10_nonmatch = merged_program_only_ab10_results[merged_program_only_ab10_results['_merge'].isin(['left_only', 'right_only'])]
merged_program_only_ab10_nonmatch.to_csv('merged_program_only_ab10_nonmatch.csv', index=False)

# Eventually should create a table to tell the genes, variants found

# Scenario 3: found in abricate but not in the program- extracting contigs using Unix, identifying associated primers, and then using NCBI for alignment of the primers with the extracted contigs
# Put the list of seqs that needed to be compare into a .txt file: abricateOnly.txt
# seqtk subseq contigs_ex.fasta abricateOnly.txt > abricateOnly.fasta # This will print out all the seqs for the comparison 

# Create 'Reason' column in abricate_only df
# If 'prgGENE' is found in 'filtered_primer', fill with 'Mismatches', else 'Primers not found'
abricate_only['Reason'] = abricate_only['prgGENE'].apply(lambda x: 'Mismatches' if x in filtered_primer['shortTarget'].values else 'Primers not found')
abricate_only.to_csv('abricate_only_and_primers.csv', index=False)

# Filter out rows with 'Reason' == 'Mismatches' in abricate_only dataframe
mismatches_primers_df = abricate_only[abricate_only['Reason'] == 'Mismatches']

# # Get corresponding 'F_truseq' and 'R_truseq' from filtered_primer dataframe for each 'prgGENE' in mismatches_primers_df
mismatches_primers_df = mismatches_primers_df.merge(filtered_primer[['target_gene', 'F_truseq', 'R_truseq', 'shortTarget']], left_on='prgGENE', right_on='shortTarget', how='inner')

# Write all the primers in a single line in the text file
with open('mismatches_primers.txt', 'w') as file:
    mismatches_primers = ', '.join([f"{row['F_truseq']}, {row['R_truseq']}" for _, row in mismatches_primers_df.iterrows()])
    file.write(mismatches_primers)

# Create a set of all record IDs from the 'Record ID' column 
record_ids = mismatches_primers_df['SEQUENCE'].tolist()

# Read the contigs
fasta_seqs = SeqIO.parse(open(input_file), 'fasta')

# Create new fasta file that contains only seqs found in program
with open("abricate_only.fasta", "w") as out_file:
    for fasta in fasta_seqs:
        name, sequence = fasta.id, str(fasta.seq)
        if name in record_ids:
            SeqIO.write(fasta, out_file, "fasta")

# # Run the primers supposed to be found by the program against the contigs to check for mismatches 
run_command(['./run_blast.sh'])
