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
# python3 process3.py email_address original_fasta_file primer_file result_file 
# e.g: python3 process3.py lisaa.tran2501@gmail.com *trimmed.contigs.fa DARTE-QM_primer_design.csv final_result.csv 

# Set the display.max_rows option to print all rows 
pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)
pd.set_option('display.width', None)
pd.set_option('display.max_colwidth', None)

# Provide email to NCBI, maybe turn this into input value later
Entrez.email = sys.argv[1]
contigs = sys.argv[2]

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

# Run abricate at 10% cov
run_command(["./run_abricate.sh", contigs])

# Load the data
primers_file = sys.argv[3]
results_file = sys.argv[4]
primers = pd.read_csv(primers_file) 
results = pd.read_csv(results_file)
ab_results = pd.read_csv('abricate_results.csv', delimiter='\t')

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
ab_results = ab_results[['SEQUENCE', 'LENGTH', 'STRAND', 'GENE', '%COVERAGE', '%IDENTITY', 'DATABASE', 'ACCESSION', 'PRODUCT']]

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

# Create a new copy of primers 
primer_info = primers.copy()

# Write a FASTA file 
def write_to_fasta(filename, names, seqs):
    with open(filename, 'w') as f:
        for name, seq in zip(names, seqs):
            f.write(f">{name}\n{seq}\n")

# Rows where 'Record ID' contains ':DARTE'
darte_rows = primer_info[primer_info['target_gene'].str.contains(':DARTE')]

# Rows where 'Record ID' doesn't contain ':DARTE'
non_darte_rows = primer_info[~primer_info['target_gene'].str.contains(':DARTE')]

# Write to 'naming_card.fa'
write_to_fasta('naming_card.fa', darte_rows['target_gene'], darte_rows['product'])

# Write to 'naming_res.fa'
write_to_fasta('naming_res.fa', non_darte_rows['target_gene'], non_darte_rows['product'])

# Trimming name
def trim_name(name):
    # If name contains ':DARTE' or doesn't contain '_', return the original name
    if ':DARTE' in name or '_' not in name:
        return name
    # Split by underscore and join all parts except the last one
    return '_'.join(name.split('_')[:-1])

# Apply the function to the 'name' column to create 'shortName'
primer_info['short_target_gene'] = primer_info['target_gene'].apply(trim_name)

# Create a dictionary of exceptions
name_exceptions = {
    'tet_44__1_NZ_ABDU01000081': 'tet44',
    'tetB_P__1_NC_010937': 'tetbp',
    'erm_X__4_NC_005206': 'ermx',
    'aac_6__-Iaa_1_NC_003197': 'aac6iaa',
    'erm_44_v_1_LN623525': 'erm44',
    'aac_6__-IIc_1_NC_012555': 'aac6iic'
}

def process_card_name(name):
    if pd.isna(name) or not isinstance(name, str):
        return name 
    # Bypass rule for specific names using the dictionary
    if name in name_exceptions:
        return name_exceptions[name]
    if 'DARTE' in name:
        return name.split('_')[0].replace('-', '').replace('(','').replace(')','').lower()
    else:
        if name.startswith('bla'):
            # Remove 'bla' from the start of the name
            name_without_bla = name[3:]
            # Keep only the uppercase characters from the remaining string
            uppercase_only = ''.join([char for char in name_without_bla if char.isupper()])
            # Convert the uppercase characters to lowercase
            lowercase = uppercase_only.lower()
            # Split at the last underscore, take the first part, remove all underscores and hyphens, and convert to lowercase
            remaining = name_without_bla.rsplit('_', 1)[0].replace('_', '').replace('-', '').replace('(','').replace(')','').lower()
            return lowercase + remaining[len(lowercase):]  # Combine the lowercase letters with the rest of the string
        else:
            # Split at the last underscore and take the first part
            first_split = name.rsplit('_', 1)[0]
            # Split the first part again at the second-to-last underscore
            second_split = first_split.rsplit('_', 1)[0]
            # Remove all remaining underscores and hyphens, and convert to lowercase
            cleaned_name = second_split.replace('_', '').replace('-', '').replace('(','').replace(')','').lower()
            return cleaned_name

def process_res_name(name):
    if pd.isna(name) or not isinstance(name, str):
        return name.lower() 

    # Bypass rule for specific names using the dictionary
    if name in name_exceptions:
        return name_exceptions[name]
    
    if 'DARTE' in name:
        return name.split('_')[0].replace('-', '').replace('(','').replace(')','').lower()
    else:
        # Split at the last underscore and take the first part
        first_split = name.rsplit('_', 1)[0]
        # Split the first part again at the second-to-last underscore
        second_split = first_split.rsplit('_', 1)[0]
        # Remove all remaining underscores and hyphens, and convert to lowercase
        cleaned_name = second_split.replace('_', '').replace('-', '').replace('(','').replace(')','').lower()
        return cleaned_name

# Apply the function to the 'short_target_gene' column to create 'card_name'
primer_info['card_gene'] = primer_info['short_target_gene'].apply(process_card_name)
primer_info['res_gene'] = primer_info['short_target_gene'].apply(process_res_name)

filtered_result['Card Name'] = filtered_result['Target'].apply(process_card_name)
filtered_result['Res Name'] = filtered_result['Target'].apply(process_res_name)

def process_res_name_ab(name):
    # Remove everything after the last '_'
    short_name = name.rsplit('_', 1)[0].lower()
    # Remove specified characters
    for char in ['_', '-', '(', ')', "'"]:
        short_name = short_name.replace(char, '').lower()
    return short_name

def process_card_name_ab(name):
    # If there are two underscores, take only the part after the last '_'
    if name.count('_') == 2:
        name = name.rsplit('_', 1)[-1]
    # Remove specified characters
    for char in ['_', '-', '(', ')', "'"]:
        name = name.replace(char, '').lower()
    return name

# General function to check if one column value is a substring of another
def is_substring(row, col1, col2):
    # Check if either value is NaN
    if pd.isna(row[col1]) or pd.isna(row[col2]):
        return False
    # Ensure both values are strings
    return str(row[col1]) in str(row[col2])

# Function to apply the custom merge logic
def apply_custom_merge(df, left_col, right_col):
    df['_merge'] = df.apply(lambda row: 'both' if is_substring(row, right_col, left_col) else row['_merge'], axis=1)
    return df

# Apply the function to each string in the 'GENE' column
ab_results['Gene Check'] = ab_results['GENE'].apply(process_res_name_ab)

# Create a df that merges the abricate results and original results
merged_program_ab = filtered_result.merge(ab_results, how='outer', left_on=['Record ID','Res Name'], right_on=['SEQUENCE','Gene Check'], indicator=True)
# merged_program_ab.to_csv('merged_program_ab.csv', index=False)
merged_program_ab = apply_custom_merge(merged_program_ab, 'Res Name', 'Gene Check')

program_only = merged_program_ab[merged_program_ab['_merge'] == 'left_only']
program_only = program_only.iloc[:, :10]

# Create a set of all record IDs from the 'Record ID' column 
record_ids = program_only['Record ID'].tolist()

fasta_seqs = SeqIO.parse(open(contigs), 'fasta')

# Create new fasta file that contains seqs only found by the program
with open("program_only.fasta", "w") as out_file:
    for fasta in fasta_seqs:
        name, sequence = fasta.id, str(fasta.seq)
        if name in record_ids:
            SeqIO.write(fasta, out_file, "fasta")

# Run abricate at 80% cov and 80% identity using CARD db
run_command(['./run_abricate_card.sh'])

# Read ABR results of CARD db 
ab_results_card = pd.read_csv('abricate_results_card.csv', delimiter='\t')

# Create a length column for abricate data
ab_results_card['LENGTH'] = ab_results_card['END'] - ab_results_card['START']

# Create a df for abricate_result
ab_results_card = ab_results_card[['SEQUENCE', 'LENGTH', 'STRAND', 'GENE', '%COVERAGE', '%IDENTITY', 'DATABASE', 'ACCESSION', 'PRODUCT']]

# Apply the function to each string in the 'GENE' column
ab_results_card['Gene Check'] = ab_results_card['GENE'].apply(process_card_name_ab)

# Matching 'No' from Resfinder but may be found in CARD
matched_df = program_only.merge(ab_results_card, how='outer', left_on=['Record ID','Card Name'], right_on=['SEQUENCE','Gene Check'], indicator=True)
matched_df = apply_custom_merge(matched_df, 'Card Name', 'Gene Check')

# Dropping rows where 'Gene Match?' is 'No'
merged_program_ab = merged_program_ab[merged_program_ab['_merge'] != 'left_only']

# Concat to the final df 
merged_program_ab = pd.concat([merged_program_ab, matched_df], ignore_index=True)
# merged_program_ab.to_csv('merged_program_ab.csv', index=False)

# Rows in abricate not in program
abricate_only = merged_program_ab[merged_program_ab['_merge'] == 'right_only']
abricate_only = abricate_only.iloc[:, 10:20]
# abricate_only.to_csv('abricate_only.csv', index=False)

# Create 'Reason' column in abricate_only df
# If 'prgGENE' is found in 'primer_info', fill with 'Mismatches', else 'Primers not found'
abricate_only['Reason'] = abricate_only['Gene Check'].apply(lambda x: 'Check Report' if x in primer_info['res_gene'].values or x in primer_info['card_gene'].values else 'Primer not available')

# Filter out rows with 'Reason' == 'Mismatches' in abricate_only dataframe
mismatches_primers_df = abricate_only[abricate_only['Reason'] == 'Check Report']

# Split the dataframe based on database values
resfinder_df = mismatches_primers_df[mismatches_primers_df['DATABASE'] == 'resfinder']
card_df = mismatches_primers_df[mismatches_primers_df['DATABASE'] == 'card']

# Merge based on database
merged_resfinder = resfinder_df.merge(primer_info[['target_gene', 'F_truseq', 'R_truseq', 'res_gene', 'card_gene']], 
                                      left_on='Gene Check', right_on='res_gene', how='inner')

merged_card = card_df.merge(primer_info[['target_gene', 'F_truseq', 'R_truseq', 'res_gene', 'card_gene']], 
                            left_on='Gene Check', right_on='card_gene', how='inner')

# Concatenate the merged results
mismatches_primers_df = pd.concat([merged_resfinder, merged_card], ignore_index=True)

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

# Run the primers supposed to be found by the program against the contigs to check for mismatches 
run_command(['./run_blast.sh'])

# Analyze the results of a BLAST search
def parse_blast_report(blast_report):
    blast_results = {}
    with open(blast_report, 'r') as file:
        current_seq = None
        current_query = None  # Variable to store the current Query sequence
        for line in file:
            if line.startswith('>'):
                current_seq = line.split()[1]  
                blast_results[current_seq] = {'seqs': [], 'queries': [], 'identities': []}
            elif "Identities" in line and current_seq:
                identity_match = re.search(r'Identities = (\d+/\d+) \((\d+)%\)', line)
                if identity_match:
                    identity = identity_match.group(2) + "%"  # Keeping only the percentage
                    blast_results[current_seq]['identities'].append(identity)
            elif 'Query' in line and current_seq:
                query_match = re.search(r'Query\s+\d+\s+([A-Za-z-]+)\s+\d+', line)
                if query_match:
                    current_query = query_match.group(1)  # Extract the Query sequence
            elif 'Sbjct' in line and current_seq and current_query:
                sbjct_match = re.search(r'Sbjct\s+\d+\s+([A-Za-z-]+)\s+\d+', line)
                if sbjct_match:
                    sbjct_seq = sbjct_match.group(1)
                    blast_results[current_seq]['seqs'].append(sbjct_seq)
                    blast_results[current_seq]['queries'].append(current_query)  # Add the corresponding Query sequence
    return blast_results

# Update more info of mismatches
def update_primers_df(blast_results, mismatches_df):
    # Initialize new columns
    mismatches_df['primer 1'] = ''
    mismatches_df['primer identities 1'] = ''
    mismatches_df['primer 2'] = ''
    mismatches_df['primer identities 2'] = ''

    for index, row in mismatches_df.iterrows():
        sequence_id = row['SEQUENCE']
        if sequence_id in blast_results:
            for sbjct_seq, identity, query_seq in zip(blast_results[sequence_id]['seqs'], 
                                                      blast_results[sequence_id]['identities'], 
                                                      blast_results[sequence_id]['queries']):
                # Check for the first primer
                if sbjct_seq == row['F_truseq']:
                    mismatches_df.at[index, 'primer 1'] = query_seq
                    mismatches_df.at[index, 'primer identities 1'] = identity
                    if row['R_truseq'] in blast_results[sequence_id]['seqs']:
                        pair_index = blast_results[sequence_id]['seqs'].index(row['R_truseq'])
                        mismatches_df.at[index, 'primer 2'] = blast_results[sequence_id]['queries'][pair_index]
                        mismatches_df.at[index, 'primer identities 2'] = blast_results[sequence_id]['identities'][pair_index]
                elif sbjct_seq == row['R_truseq']:
                    mismatches_df.at[index, 'primer 1'] = query_seq
                    mismatches_df.at[index, 'primer identities 1'] = identity
                    if row['F_truseq'] in blast_results[sequence_id]['seqs']:
                        pair_index = blast_results[sequence_id]['seqs'].index(row['F_truseq'])
                        mismatches_df.at[index, 'primer 2'] = blast_results[sequence_id]['queries'][pair_index]
                        mismatches_df.at[index, 'primer identities 2'] = blast_results[sequence_id]['identities'][pair_index]
                # Similar checks for reverse complements
                elif sbjct_seq == row['F_truseq_reverse_complement']:
                    mismatches_df.at[index, 'primer 1'] = query_seq
                    mismatches_df.at[index, 'primer identities 1'] = identity
                    if row['R_truseq_reverse_complement'] in blast_results[sequence_id]['seqs']:
                        pair_index = blast_results[sequence_id]['seqs'].index(row['R_truseq_reverse_complement'])
                        mismatches_df.at[index, 'primer 2'] = blast_results[sequence_id]['queries'][pair_index]
                        mismatches_df.at[index, 'primer identities 2'] = blast_results[sequence_id]['identities'][pair_index]
                elif sbjct_seq == row['R_truseq_reverse_complement']:
                    mismatches_df.at[index, 'primer 1'] = query_seq
                    mismatches_df.at[index, 'primer identities 1'] = identity
                    if row['F_truseq_reverse_complement'] in blast_results[sequence_id]['seqs']:
                        pair_index = blast_results[sequence_id]['seqs'].index(row['F_truseq_reverse_complement'])
                        mismatches_df.at[index, 'primer 2'] = blast_results[sequence_id]['queries'][pair_index]
                        mismatches_df.at[index, 'primer identities 2'] = blast_results[sequence_id]['identities'][pair_index]
    return mismatches_df

# Parsing the blast report
blast_results = parse_blast_report("abricate_only_blast_report.txt")
mismatches_primers_df = update_primers_df(blast_results, mismatches_primers_df)
# print(mismatches_primers_df)

# Removing the percentage sign from the 'primer identities' columns
mismatches_primers_df['primer identities 1'] = mismatches_primers_df['primer identities 1'].str.replace('%', '')
mismatches_primers_df['primer identities 2'] = mismatches_primers_df['primer identities 2'].str.replace('%', '')

# Converting the 'primer identities' columns to numeric.
mismatches_primers_df['primer identities 1'] = pd.to_numeric(mismatches_primers_df['primer identities 1'], errors='coerce')
mismatches_primers_df['primer identities 2'] = pd.to_numeric(mismatches_primers_df['primer identities 2'], errors='coerce')

# Add 'fix_name' and 'reason' columns to a df based on primer identity values
def add_fix_name_and_reason(df):
    def determine_fix_name_and_reason(row):
        # Extract primer identity and primer values from the row
        primer1_identity = row['primer identities 1']
        primer2_identity = row['primer identities 2']
        primer1 = row['primer 1']
        primer2 = row['primer 2']
        fix_name = 'N/A'
        reason = 'N/A'

        # Convert identity values to integers, except NaN 
        if not pd.isna(primer1_identity):
            primer1_identity = int(float(primer1_identity))
        if not pd.isna(primer2_identity):
            primer2_identity = int(float(primer2_identity))

        # Determine 'fix_name' and 'reason' based on primer identity values and database 
        if primer1_identity == 100 and primer2_identity == 100:
            # If both primer identities are 100% -> partial primers
            fix_name = row['res_gene'] if row['DATABASE'] == 'resfinder' else (row['card_gene'] if row['DATABASE'] == 'card' else 'N/A')
            reason = 'Partial primers'

        elif primer1_identity >= 80 and primer2_identity >= 80:
            # If both primer identities are >= 80% -> mismatches
            fix_name = row['res_gene'] if row['DATABASE'] == 'resfinder' else (row['card_gene'] if row['DATABASE'] == 'card' else 'N/A')
            reason = 'Mismatches'

        elif (primer1_identity >= 80 and pd.isna(primer2_identity)) or (primer2_identity >= 80 and pd.isna(primer1_identity)):
            # If one primer identity is >= 80% and the other is NaN -> missing primer
            fix_name = row['res_gene'] if row['DATABASE'] == 'resfinder' else (row['card_gene'] if row['DATABASE'] == 'card' else 'N/A')
            reason = 'Missing primer'

        return fix_name, reason

    # Apply the function to each row of the df 
    df[['fix_name', 'reason']] = df.apply(determine_fix_name_and_reason, axis=1, result_type='expand')
    return df

# Add 'fix_name' and 'reason' columns to the df
df = add_fix_name_and_reason(mismatches_primers_df)
df.to_csv('abricate_only_report.csv', index=False)

# Create copies of the dataframes
df1_copy = merged_program_ab[merged_program_ab['_merge'] == 'both'].copy()
df2_copy = merged_program_ab[merged_program_ab['_merge'] == 'left_only'].copy()
df3_copy = abricate_only.copy()

# Add the 'extra note' column to each copy without specifying a name
df1_copy["Location"] = 'found in both program and ABR'
df2_copy["Location"] = 'only found in program'
df3_copy["Location"] = 'only found in ABR'

# Assign 'right_only' to the '_merge' column of df3_copy
df3_copy['_merge'] = 'right_only'

# Append the dfs 
final_df = df1_copy.append(df2_copy).append(df3_copy).reset_index(drop=True)

# Convert to final reports
final_df.to_csv('final_report.csv', index=False)
#############################################################################################################


# find . -type d -name 'SRR*' -exec sh -c 'find "$0" -name "final_report.csv" -print0 | xargs -0 awk "{print}"' {} \; > combined_final_report.csv
# find . -type d -name 'SRR*' -exec sh -c 'find "$0" -name "abricate_only_report.csv" -print0 | xargs -0 awk "{print}"' {} \; > combined_abricate_only_report.csv

