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
# pip install pandas numpy biopython matplotlib
# python3 process4.py final_report_file primer_file
# e.g: python3 process4.py final_report.csv DARTE-QM_primer_design.csv

# Set the display.max_rows option to print all rows 
pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)
pd.set_option('display.width', None)
pd.set_option('display.max_colwidth', None)

# Depend on the type of result, could be individual file or requires to be combine before using
result = sys.argv[1]
report = pd.read_csv(result)

# Read the primer file 
primer = sys.argv[2]
ori_primer = pd.read_csv(primer)

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


# Write first 500 results from sample dataset into a FASTA file for testing purposes 
# awk -F, 'NR>1 {print ">"$1"\n"$8}' combined_final_report.csv | head -n1000 > testing.fasta 

