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
# python3 process2.py DARTE-QM_primer_design.csv raw_results.csv 

# Set the display.max_rows option to print all rows 
pd.set_option('display.max_rows', None)

# Function to get the reverse complement of a sequence
def reverse_complement(seq):
    complement_dict = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return "".join(complement_dict[base] for base in reversed(seq))

# Load the data
primers_file = sys.argv[1] 
results_file = sys.argv[2] 
primers = pd.read_csv(primers_file)
results = pd.read_csv(results_file)

# Assuming df is your DataFrame
primers['target_gene'] = primers['target_gene'].str.replace('^cfrB', 'cfr_B_', regex=True)

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
# match.drop(columns=['Qualified?']).to_csv('match.csv', index=False)
# nonmatch.drop(columns=['Qualified?']).to_csv('nonmatch.csv', index=False)


