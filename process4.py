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
# e.g: python3 process4.py 

# Set the display.max_rows option to print all rows 
pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)
pd.set_option('display.width', None)
pd.set_option('display.max_colwidth', None)

# Load the DataFrame
report = pd.read_csv('combined_final_report.csv')

# Identify rows where the first column equals 'Record ID'
rows_to_remove = report.iloc[:, 0] == 'Record ID'

# Remove these rows
report = report[~rows_to_remove]

# Now proceed with the rest of your operations
report = report.drop('originalCONTIG', axis=1)
report.rename(columns={report.columns[21]: 'Location'}, inplace=True)
report.to_csv('short_combined_final_report.csv', index=False)

# Concatenate the values from F_primer and R_primer columns, excluding the header row
info = pd.concat([report['F_primer'][1:], report['R_primer'][1:]])

# Count the frequency of each primer
primer_freq = info.value_counts()

# Convert to DataFrame, if necessary
primer_freq_df = primer_freq.reset_index()
primer_freq_df.columns = ['Primer', 'Frequency']

print(primer_freq_df)

# Count the frequency of each primer in F_primer and R_primer columns
f_primer_freq = report['F_primer'].value_counts().reset_index()
f_primer_freq.columns = ['Primer', 'F_Frequency']

r_primer_freq = report['R_primer'].value_counts().reset_index()
r_primer_freq.columns = ['Primer', 'R_Frequency']

# Optionally, merge the two frequency DataFrames for comparison
primer_freq_combined = pd.merge(f_primer_freq, r_primer_freq, on='Primer', how='outer')

print(primer_freq_combined)


