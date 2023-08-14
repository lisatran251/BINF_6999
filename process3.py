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
# pip install pandas numpy biopython
# python3 process3.py email_address primer_file result_file original_fasta_file 
# e.g: python3 process3.py lisaa.tran2501@gmail.com DARTE-QM_primer_design.csv final_result.csv contigs_ex.fasta 

# Set the display.max_rows option to print all rows 
pd.set_option('display.max_rows', None)

# Provide email to NCBI, maybe turn this into input value later
Entrez.email = sys.argv[1]

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
# run_command(['./run_abricate.sh'])

# Load the data
primers_file = sys.argv[2]
results_file = sys.argv[3]
primers = pd.read_csv(primers_file) 
results = pd.read_csv(results_file)
ab_results = pd.read_csv('abricate_results.csv', delimiter='\t')
# Read the contigs
contigs = sys.argv[4]

# Fill in the original contigs in result df
sequence_mapping = {record.id: str(record.seq) for record in SeqIO.parse(contigs, "fasta")}

# Map the sequences to the df to create the new column 
results['originalCONTIG'] = results['Record ID'].map(sequence_mapping)

# Create a df for the final_result
filtered_result = results[['Record ID', 'Length', 'Combination', 'F_Product', 'F_primer', 'R_primer', 'Target', 'originalCONTIG']]
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
# Replace values in 'shortProduct' column
filtered_result['shortProduct'] = filtered_result['shortProduct'].replace(r'^aac_3_-IV.*$', 'aac_3_-IV', regex=True)

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

# # Sort before merging 
# merged_program_ab = merged_program_ab.sort_values(by='Record ID')

# Final report after merging
merged_program_ab.to_csv('merged_prog_ab.csv', index=False)

# Rows in program not in abricate
program_only = merged_program_ab[merged_program_ab['_merge'] == 'left_only']
program_only = program_only.iloc[:, :10]
# program_only.to_csv('program_only.csv', index=False)

# Rows in abricate not in program
abricate_only = merged_program_ab[merged_program_ab['_merge'] == 'right_only']
abricate_only = abricate_only.iloc[:, 10:20]
# abricate_only.to_csv('abricate_only.csv', index=False)

# Rows found in both progsram ands abricate
both_program_ab_result = merged_program_ab[merged_program_ab['_merge'] == 'both']
# both_program_ab_result.drop('_merge', axis=1).to_csv('both_program_ab_result.csv', index=False)

## Scenario 1: target gene found - explore variants, artifacts etc. 
# Select rows from merged_program_ab where 'Gene Match?' = 'Yes' and 'Variant Match?' = 'No'
# filtered_merged_program_ab = merged_program_ab[(merged_program_ab['Gene Match?'] == 'Yes') & (merged_program_ab['Variant Match?'] == 'No')]
# # filtered_merged_program_ab.to_csv('filtered_merged_program_ab.csv', index=False)

# # Define a function to get the sequence from the original contigs 
# def get_sequence(record_id):    
#     with open(contigs, 'r') as f:
#         for record in SeqIO.parse(f, 'fasta'):
#             if record.id == record_id:
#                 return str(record.seq)
#     return None 

# # Create a new column in df to contain the contigs, can create a copy to avoid warning (opt)
# filtered_merged_program_ab['Contigs'] = filtered_merged_program_ab['Record ID'].apply(get_sequence)

# # Write a fasta include Record ID and Contigs
# with open('variants_originalCONTIG.fasta', 'w') as f:
#     for index, row in filtered_merged_program_ab.iterrows():
#         f.write(f">seq_{index}_originalCONTIGS\n{row['Contigs']}\n")

# # Write a fasta include Record ID and original sequence from Abricate program
# with  open('variants_originalSEQUENCE.fasta', 'w') as f:
#     for index, row in filtered_merged_program_ab.iterrows():
#         f.write(f">seq_{index}_originalSEQUENCE\n{row['originalSEQUENCE']}\n")

## Run cdhit in Unix: cd-hit-est -i varsiants.fasta -o cdhit_rs -c 0.9 -n 10 -d 0
# The main reason is tsshe product is very short while the sequence achieved using accession number is much longer. Therefore, when run cd-hit which based on the similarity between two sequences and the identiy % much be higher than 70%, therefore we may want to run BLAST instead.  

## Scenario 2: found in program not in abricate- re-run with lower %identity on abricate
# Create a set of all record IDs from the 'Record ID' column 
record_ids = program_only['Record ID'].tolist()

fasta_seqs = SeqIO.parse(open(contigs), 'fasta')

# Create new fasta file that contains seqs only found by the program
with open("program_only.fasta", "w") as out_file:
    for fasta in fasta_seqs:
        name, sequence = fasta.id, str(fasta.seq)
        if name in record_ids:
            SeqIO.write(fasta, out_file, "fasta")

# Re-run Abricate with lower coverage: abricate --db resfinder --quiet --mincov 10 program_only.fasta > program_only_ab10_results.csv
# run_command(['./run_abricate_10.sh'])

# Load files
program_only_ab10_results = pd.read_csv('program_only_ab10_results.csv', delimiter='\t')

# Create a length column for abricate data
program_only_ab10_results['LENGTH'] = program_only_ab10_results['END'] - program_only_ab10_results['START']

# Create a df for abricate_result
program_only_ab10_results = program_only_ab10_results[['SEQUENCE', 'LENGTH', 'STRAND', 'GENE', '%COVERAGE', '%IDENTITY', 'ACCESSION', 'PRODUCT']]

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
# merged_program_only_ab10_match.to_csv('merged_program_only_ab10_match.csv', index=False)

# Print out df includes non-matches 
merged_program_only_ab10_nonmatch = merged_program_only_ab10_results[merged_program_only_ab10_results['_merge'].isin(['left_only', 'right_only'])]
# merged_program_only_ab10_nonmatch.to_csv('merged_program_only_ab10_nonmatch.csv', index=False)

# Scenario 3: found in abricate but not in the program- extracting contigs using Unix, identifying associated primers, and then using NCBI for alignment of the primers with the extracted contigs
# Put the list of seqs that needed to be compare into a .txt file: abricateOnly.txt
# seqtk subseq contigs_ex.fasta abricateOnly.txt > abricateOnly.fasta # This will print out all the seqs for the comparison 

# Create 'Reason' column in abricate_only df
# If 'prgGENE' is found in 'filtered_primer', fill with 'Mismatches', else 'Primers not found'
abricate_only['Reason'] = abricate_only['prgGENE'].apply(lambda x: 'Mismatches' if x in filtered_primer['shortTarget'].values else 'Primers not found')

# Filter out rows with 'Reason' == 'Mismatches' in abricate_only dataframe
mismatches_primers_df = abricate_only[abricate_only['Reason'] == 'Mismatches']

# Merge to get corresponding 'F_truseq' and 'R_truseq' from filtered_primer dataframe for each 'prgGENE' in mismatches_primers_df
mismatches_primers_df = mismatches_primers_df.merge(filtered_primer[['target_gene', 'F_truseq', 'R_truseq', 'shortTarget']], left_on='prgGENE', right_on='shortTarget', how='inner')

# Calculate reverse complements for each row in mismatches_primers_df
mismatches_primers_df['F_truseq_reverse_complement'] = mismatches_primers_df['F_truseq'].apply(lambda x: str(Seq(x).reverse_complement()))
mismatches_primers_df['R_truseq_reverse_complement'] = mismatches_primers_df['R_truseq'].apply(lambda x: str(Seq(x).reverse_complement()))

# Write all the primers (original and reverse complemented) in a single line in the text file
with open('mismatches_primers.txt', 'w') as file:
    mismatches_primers = ', '.join([f"{row['F_truseq']}, {row['R_truseq']}, {row['F_truseq_reverse_complement']}, {row['R_truseq_reverse_complement']}" for _, row in mismatches_primers_df.iterrows()])
    file.write(mismatches_primers)

# Create a set of all record IDs from the 'Record ID' column 
record_ids = mismatches_primers_df['SEQUENCE'].tolist()

# Read the contigs
fasta_seqs = SeqIO.parse(open(contigs), 'fasta')

# Create new fasta file that contains only seqs found in program
with open("abricate_only.fasta", "w") as out_file:
    for fasta in fasta_seqs:
        name, sequence = fasta.id, str(fasta.seq)
        if name in record_ids:
            SeqIO.write(fasta, out_file, "fasta")

# # Run the primers supposed to be found by the program against the contigs to check for mismatches 
run_command(['./run_blast.sh'])

## Printing final report
# Create copies of the dataframes
df1_copy = merged_program_ab.copy()
df2_copy = merged_program_only_ab10_results.copy()
df3_copy = abricate_only.copy()

# Add the 'extra note' column to each copy without specifying a name
df1_copy[""] = 'high coverage'
df2_copy[""] = 'low coverage'
df3_copy[""] = 'only found in Abricate'

# Append the dataframes one after the other
final_df = df1_copy.append(df2_copy).append(df3_copy).reset_index(drop=True)

# # Export to HTML
# html = final_df.to_html()

# # Save to a file
# with open('output.html', 'w') as f:
#     f.write(html)

final_df = final_df.drop(['originalCONTIG', 'originalSEQUENCE'], axis=1)
final_df.to_csv('final_report.csv', index=False)

#############################################################################################################