import pandas as pd
import numpy as np 
import time
import argparse
import os
import re
import sys
import subprocess 
import traceback
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
from Bio import SeqIO, Entrez 
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# How to run:      
# python3 -m venv ~/my_venv
# source ~/my_venv/bin/activate
# pip install pandas numpy biopython argparse
# python3 masterscript.py fasta_file primer_file 
# e.g: python3 masterscript.py testing2.fasta DARTE-QM_primer_design.csv

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

# Provide email to NCBI, maybe turn this into input value later
Entrez.email = "lisaa.tran2501@gmail.com"
contigs = input_file

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
results = pd.read_csv('final_result.csv')
primers = pd.read_csv(primer_file) 
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

# Open the FASTA file
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

# Apply for individual assemblies
# find . -type d -name 'SRR*' -exec sh -c 'find "$0" -name "final_report.csv" -print0 | xargs -0 awk "{print}"' {} \; > combined_final_report.csv
# find . -type d -name 'SRR*' -exec sh -c 'find "$0" -name "abricate_only_report.csv" -print0 | xargs -0 awk "{print}"' {} \; > combined_abricate_only_report.csv

# Depend on the type of result, could be individual file or requires to be combine before using
report = pd.read_csv('final_report.csv')

# Read the primer file 
ori_primer = pd.read_csv(primer_file)

output_file = open('analysis_result.txt', 'w')
## Apply this for multiple individual assemblies
## Prepare an empty DataFrame to hold all combined data
# report = pd.DataFrame()

## Path to the parent directory containing all your 'SRR' directories
# parent_directory = '.'

## Loop through each directory in the parent directory
# for directory in os.listdir(parent_directory):
#     if directory.startswith('SRR'):
#         # Construct the full path to the 'final_report.csv' file
#         file_path = os.path.join(parent_directory, directory, 'final_report.csv')
        
#         # Check if the file exists to avoid errors
#         if os.path.exists(file_path):
#             # Read the file into a DataFrame
#             df = pd.read_csv(file_path)
#             # Add a new column with the name of the directory (sample name)
#             df['Sample'] = directory
#             # Append this DataFrame to the combined DataFrame
#             report = report.append(df, ignore_index=True)

## Identify rows where the first column equals 'Record ID'
# rows_to_remove = report.iloc[:, 0] == 'Record ID'

## Remove these rows
# report = report[~rows_to_remove]

# Now proceed with the rest of your operations
report = report.drop('originalCONTIG', axis=1)
# report.rename(columns={report.columns[20]: 'Location'}, inplace=True)
# report.to_csv('short_combined_final_report.csv', index=False)
report_for_plot = report.copy()

# Write a .txt file results
output_file = open('analysis_result.txt', 'w')

# Count total products found 
total_products = len(report.index)

# Count the number of products found by both program and ABR
num_gene = report[report['_merge'] == 'both'].shape[0]

# Count the number of products only found by the program
num_left_only = report[report['_merge'] == 'left_only'].shape[0]

# Count the number of products only found by ABR (w or w/o primers)
num_unavai_primer = report[(report['_merge'] == 'right_only') & (report['Reason'] == 'Primer not available')].shape[0]
num_miss = report[(report['_merge'] == 'right_only') & (report['Reason'] != 'Primer not available')].shape[0]

# Categories
categories = ['Both', 'Program Only', 'ABR Only (primers available)', 'ABR Only (primers not available)']

# Counts
counts = [num_gene, num_left_only, num_unavai_primer, num_miss]

# Creating the bar chart
plt.figure(figsize=(10, 6))
plt.bar(categories, counts, color=['blue', 'green', 'red', 'orange'])

# print the result
print(f"Analysis Summarization", file=output_file)
print(f"Number of products found: {total_products}", file=output_file)
print(f"Number of target genes found (both program and ABR): {num_gene}", file=output_file)
print(f"Number of products found only by the program: {num_left_only}", file=output_file)
print(f"Number of products found only by the ABR (primers available): {num_miss}", file=output_file)
print(f"Number of products found only by the ABR (primers not available): {num_unavai_primer}\n", file=output_file)

# Concatenate the values from F_primer and R_primer columns, excluding the header row
info = pd.concat([report['F_primer'][1:], report['R_primer'][1:]])

# Count the frequency of each primer
primer_freq = info.value_counts()
primer_freq_df = primer_freq.reset_index() 

# Print primers found and their counts
primer_freq_df.columns = ['Primer', 'Frequency']
print(primer_freq_df, file=output_file)
# primer_freq_df.to_csv('primer_freq_df.csv', index=False)

# Combine 'F_truseq' and 'R_truseq' columns
combined_primers = pd.concat([ori_primer['F_truseq'], ori_primer['R_truseq']]).unique()

# Convert the primer_freq_df 'Primer' column to a set for faster lookup
primer_freq_set = set(primer_freq_df['Primer'])

# Find primers not found in primer_freq_df
not_found_primers = [primer for primer in combined_primers if primer not in primer_freq_set]

# Print the number of primers not found and the primers themselves
print(f"\nNumber of primers not found: {len(not_found_primers)}", file=output_file)

# Join the primers into a single string separated by commas
primers_string = ', '.join(not_found_primers)
print("Unfound primers:", file=output_file)
print(primers_string, file=output_file)

# Count the frequency of each primer in F_primer and R_primer columns
f_primer_freq = report['F_primer'].value_counts().reset_index()
f_primer_freq.columns = ['Primer', 'F_Frequency']

r_primer_freq = report['R_primer'].value_counts().reset_index()
r_primer_freq.columns = ['Primer', 'R_Frequency']

# For comparison
primer_freq_combined = pd.merge(f_primer_freq, r_primer_freq, on='Primer', how='outer')
print(f"\n", file=output_file)
print(primer_freq_combined, file=output_file)

# Count the frequency of each target gene
target_freq = report['Gene Check'].value_counts()
target_freq_df = target_freq.reset_index()
print(f"\n", file=output_file)
target_freq_df.columns = ['Target Gene', 'Frequency']

# Print target genes and their frequencies
print(target_freq_df, file=output_file)

# Visualization 
# Store counts in a dictionary
product_counts = {
    'Category': ['Both', 'Program Only', 'ABR Only (primers available)', 'ABR Only (primers not available)'],
    'Count': [num_gene, num_left_only, num_miss, num_unavai_primer]
}

# Convert dictionary to df
df_product_counts = pd.DataFrame(product_counts) 
 
# Plotting the product counts
plt.figure(figsize=(15, 10))  # Increase figure width as needed
plt.bar(categories, counts, color=['blue', 'green', 'red', 'orange'])
plt.title('Comparison of Product Findings')
plt.xlabel('Category')
plt.ylabel('Counts')
plt.savefig('product_count.png')
plt.close()  

# Assign unique identifiers to each primer
f_primer_ids = {primer: f'F{i}' for i, primer in enumerate(report['F_primer'].unique())}
r_primer_ids = {primer: f'R{i}' for i, primer in enumerate(report['R_primer'].unique())}

report['F_primer_id'] = report['F_primer'].map(f_primer_ids)
report['R_primer_id'] = report['R_primer'].map(r_primer_ids)

# Filter the data for 'found in both program and ABR'
filtered_data = report[report['Location'] == 'found in both program and ABR']

# Create a pivot table
pivot_table = pd.pivot_table(filtered_data, values='Record ID', 
                             index='F_primer_id', columns='R_primer_id', 
                             aggfunc='count', fill_value=0)  

primer_freq_combined.fillna(0, inplace=True) 

# Plotting a bar plot with counts as the x-axis
plt.figure(figsize=(10, 6))
plt.bar(target_freq_df['Target Gene'], target_freq_df['Frequency'], alpha=0.7)

# Remove the x-tick labels to avoid overlaps
plt.xticks([])
# Labels and title
plt.xlabel('Target Genes')
plt.ylabel('Counts')
plt.title('Bar Plot of Target Gene Counts')
plt.savefig('target_frq.png')
plt.close() 

# Assuming combined_df is your DataFrame
# gene_counts = report.groupby(['Gene Check', 'Sample']).size().unstack(fill_value=0)

# Sort the DataFrame to get the top 30 primers by frequency
top_genes = target_freq_df.sort_values(by='Frequency', ascending=False).head(31)

# Generate a list of colors based on the 'Frequency' value
num_unique_freqs = top_genes['Frequency'].nunique()
palette = sns.color_palette("Spectral", num_unique_freqs)

# Map each frequency to a color
freq_to_color = {freq: color for freq, color in zip(sorted(top_genes['Frequency'].unique()), palette)}

# Apply the mapping to assign a color to each row in the df
top_genes['Color'] = top_genes['Frequency'].map(freq_to_color)

# Create the scatter plot with manual color assignment
plt.figure(figsize=(15, 8))
for _, row in top_genes.iterrows():
    plt.scatter(row['Target Gene'], row['Frequency'], color=row['Color'], label=row['Target Gene'])

# Set the y-axis to start from 20 with increments of 5
y_ticks = np.arange(0, top_genes['Frequency'].max() + 5, 5)
plt.yticks(y_ticks)

# Remove the x-axis labels
plt.xticks([])

# Create custom legend handles
legend_handles = [mpatches.Patch(color=freq_to_color[freq], label=primer) 
                  for freq, primer in zip(top_genes['Frequency'], top_genes['Target Gene'])]

# Create the legend, placing it outside of the plot to the right
plt.legend(handles=legend_handles, title='Target Genes', bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)

plt.title('Top 30 Most Found Target Genes')
plt.ylabel('Frequency')
plt.tight_layout()
plt.savefig('target_genes.png', bbox_inches='tight')
plt.close()

# Plotting the heatmap
plt.figure(figsize=(12, 10))
sns.heatmap(pivot_table, annot=True, cmap='Paired')
plt.title('Heatmap of Primer Combinations Leading to Correct Target Genes')
plt.xticks([])  # Remove x-axis labels
plt.yticks([])  # Remove y-axis labels
plt.ylabel('Forward Primer Identifier')
plt.xlabel('Reverse Primer Identifier')
plt.savefig('heatmap.png')
plt.close()

## Apply for individual assembly
## Path to the parent directory containing all your 'SRR' directories
# parent_directory = '.'

# all_df_product_counts = []

# # Loop through each directory in the parent directory
# for directory in os.listdir(parent_directory):
#     if directory.startswith('SRR'):
#         # Construct the full path to the result file in the directory
#         file_path = os.path.join(parent_directory, directory, 'final_report.csv')
        
#         # Check if the file exists to avoid errors
#         if os.path.exists(file_path):
#             # Read the file into a df
#             report = pd.read_csv(file_path)

#             # Count the number of products found by both program and ABR
#             num_gene = report[report['_merge'] == 'both'].shape[0]

#             # Count the number of products only found by the program
#             num_left_only = report[report['_merge'] == 'left_only'].shape[0]

#             # Count the number of products only found by ABR (w or w/o primers)
#             num_unavai_primer = report[(report['_merge'] == 'right_only') & (report['Reason'] == 'Primer not available')].shape[0]
#             num_miss = report[(report['_merge'] == 'right_only') & (report['Reason'] != 'Primer not available')].shape[0]

#             # Store counts in a dictionary
#             product_counts = {
#                 'Category': ['Both', 'Program Only', 'ABR Only (primers available)', 'ABR Only (primers not available)'],
#                 'Count': [num_gene, num_left_only, num_miss, num_unavai_primer]
#             }

#             # Convert dictionary to df
#             df_product_counts = pd.DataFrame(product_counts) 

#             # Add the sample name (directory name) as a column
#             df_product_counts['Sample'] = directory

#             # Add the df to the list
#             all_df_product_counts.append(df_product_counts)

# # Concatenate all dfs in the list into one df
# final_df_product_counts = pd.concat(all_df_product_counts)

# # Reset index of the final df
# final_df_product_counts.reset_index(drop=True, inplace=True)

# # print(final_df_product_counts)
# # Create a pivot table for the stacked bar chart
# pivot_df = final_df_product_counts.pivot(index='Sample', columns='Category', values='Count')

# # Increase the figure size for better visibility
# plt.figure(figsize=(14, 8))  # You can adjust this size as needed

# # Plot stacked bar chart
# pivot_df.plot(kind='bar', stacked=True, figsize=(14, 8))  

# # Title and labels
# # plt.title('Stacked Bar Chart of Categories per Sample')
# plt.xlabel('Samples')
# plt.ylabel('Counts')

# # Rotate x-axis labels to prevent overlap and improve readability
# plt.xticks(rotation=90)  # Rotate by 90 degrees

# # Move the legend outside of the plot to the right side
# plt.legend(title='Category', bbox_to_anchor=(1, 1), loc='upper left')

# # Adjust layout to make room for the x-axis labels and the legend
# plt.tight_layout(rect=[0, 0, 0.75, 1])  # You might need to tweak these values

# plt.savefig('barchart.png')
# plt.close()
print(f"\n", file=output_file)
print("Figure saved as 'product_count.png' showing the comparison of product found in each category.", file=output_file)
print("Figure saved as 'target_frq.png' showing the bar plot of target genes countss.", file=output_file)
print("Figure saved as 'target_genes.png' showing top 30 target genes found.", file=output_file)
print("Figure saved as 'heatmap.png' showing the primer combinations leading to correct target genes.", file=output_file)
output_file.close()

