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
from Bio import SeqIO, Entrez 
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# How to run:
# python3 -m venv ~/my_venv
# source ~/my_venv/bin/activate
# pip install pandas numpy biopython matplotlib
# python3 process3.py email_address original_fasta_file primer_file result_file 
# e.g: python3 process4.py 

# Set the display.max_rows option to print all rows 
pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)
pd.set_option('display.width', None)
pd.set_option('display.max_colwidth', None)

# Depend on the type of result, could be individual file or requires to be combine before using
result = sys.argv[1]

# Load the DataFrame
report = pd.read_csv(result)

# Identify rows where the first column equals 'Record ID'
rows_to_remove = report.iloc[:, 0] == 'Record ID'

# Remove these rows
report = report[~rows_to_remove]

# Now proceed with the rest of your operations
report = report.drop('originalCONTIG', axis=1)
report.rename(columns={report.columns[20]: 'Location'}, inplace=True)
report.to_csv('short_combined_final_report.csv', index=False)
# print(report.head())

# count total products found 
total_products = len(report.index)

# count the number of products found by both program and ABR
num_gene = report[report['_merge'] == 'both'].shape[0]

# count the number of products only found by the program
num_left_only = report[report['_merge'] == 'left_only'].shape[0]

# count the number of products only found by ABR (w or w/o primers)
num_unavai_primer = report[(report['_merge'] == 'right_only') & (report['Reason'] == 'Primer not available')].shape[0]
num_miss = report[(report['_merge'] == 'right_only') & (report['Reason'] != 'Primer not available')].shape[0]

# print the result
print(f"\nSummarization")
print(f"Number of products found: {total_products}")
print(f"Number of target genes found (both program and ABR): {num_gene}")
print(f"Number of products found only by the program: {num_left_only}")
print(f"Number of products found only by the ABR (primers available): {num_miss}")
print(f"Number of products found only by the ABR (primers not available: {num_unavai_primer}\n")

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

# for comparison
primer_freq_combined = pd.merge(f_primer_freq, r_primer_freq, on='Primer', how='outer')

# replace NaN with an empty string
# primer_freq_combined.fillna('', inplace=True)

# # convert frequency columns to integers where not empty
# primer_freq_combined['F_Frequency'] = primer_freq_combined['F_Frequency'].apply(lambda x: int(x) if x != '' else x)
# primer_freq_combined['R_Frequency'] = primer_freq_combined['R_Frequency'].apply(lambda x: int(x) if x != '' else x)

print(primer_freq_combined)

# Count the frequency of each target gene
target_freq = report['Target'].value_counts()

# Convert to DataFrame, if necessary
target_freq_df = target_freq.reset_index()
target_freq_df.columns = ['Target Gene', 'Frequency']

print(target_freq_df)

# Visualization 
# store counts in a dictionary
product_counts = {
    'Category': ['Both', 'Program Only', 'ABR Only (primers available)', 'ABR Only (primers not available)'],
    'Count': [num_gene, num_left_only, num_miss, num_unavai_primer]
}

# convert dictionary to df
df_product_counts = pd.DataFrame(product_counts)

# Plotting the product counts
plt.figure(figsize=(15, 6))  
plt.bar(df_product_counts['Category'], df_product_counts['Count'], color='blue')
plt.xlabel('Category')
plt.ylabel('Count')
plt.title('Product Counts by Category')
plt.xticks(rotation=0, fontsize=10) 
# plt.show()
plt.savefig('product_count.png')
plt.close()  

primer_freq_combined.fillna(0, inplace=True)

# Scatter Plot
plt.scatter(primer_freq_combined['F_Frequency'], primer_freq_combined['R_Frequency'])
plt.xlabel('F_Primer Frequency')
plt.ylabel('R_Primer Frequency')
plt.title('Scatter Plot of Primer Frequencies')
plt.savefig('fr_primer_frq.png')
plt.close()  

# plotting the line chart of target gene frequency
plt.xticks([])  # remove the x-tick labels to avoid overlaps
plt.plot(target_freq_df['Target Gene'], target_freq_df['Frequency'])
plt.xlabel('Target Gene')
plt.ylabel('Frequency')
plt.title('Line Chart of Target Gene Frequencies')
plt.savefig('target_frq.png')
plt.close()  

