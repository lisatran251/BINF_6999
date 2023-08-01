import pandas as pd
import numpy as np 
import sys
from Bio import SeqIO 
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import re
import subprocess 

# How to run:
# python3 -m venv ~/my_venv
# source ~/my_venv/bin/activate
# pip install pandas numpy biopython

# Set the display.max_rows option to print all rows 
pd.set_option('display.max_rows', None)

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
primers = pd.read_csv('DARTE-QM_primer_design.csv')
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

# Create a copy of primer df  
filtered_primer = primers.copy()

# Create a new column that convert the name of the primer to shorter name 
filtered_primer['shortTarget'] = filtered_primer['target_gene'].apply(shorten_target)

# Define a function to replace every "'" with "_"
def replace_chars(s):
    return re.sub(r"[\(''\)]", "_", s)

# Apply the function to each string in the 'GENE' column
ab_results['prgGENE'] = ab_results['GENE'].apply(replace_chars)

# Create a df that merges the abricate results and original results
merged_program_ab = filtered_result.merge(ab_results, how='outer', left_on=['Record ID','shortTarget'], right_on=['SEQUENCE','prgGENE'], indicator=True)
merged_program_ab.to_csv('merged_prog_ab.csv', index=False)

# Rows in program not in abricate
program_only = merged_program_ab[merged_program_ab['_merge'] == 'left_only']
# program_only.drop('_merge', axis=1).to_csv('program_only.csv', index=False)

# Rows in abricate not in program
abricate_only = merged_program_ab[merged_program_ab['_merge'] == 'right_only']
# abricate_only.drop('_merge', axis=1).to_csv('abricate_only.csv', index=False)

# Rows found in both progsram not in abricate
both_program_ab_result = merged_program_ab[merged_program_ab['_merge'] == 'both']
# both_program_ab_result.drop('_merge', axis=1).to_csv('both_program_ab_result.csv', index=False)

# ## Scenario 1: target gene found - explore variants, artifacts etc. 
# # Save the F_product to a fasta file
# with open('F_product.fasta', 'w') as f:
#     for index, row in results.iterrows():
#         # Include target in identifier if it exists and is not NaN
#         target = f"_{row['Target']}" if pd.notna(row['Target']) else ""
#         f.write(f">{row['Record ID']}_{target}\n{row['F_Product']}\n")

# # # Run cdhit in Unix: cd-hit-est -i F_product.fasta -o cdhit_rs -c 0.9 -n 10 -d 0

# ## Scenario 2: found in program not in abricate- re-run with lower %identity on abricate
# # Create a set of all record IDs from the 'Record ID' column 
# record_ids = program_only['Record ID'].tolist()

# # Read the contigs
# fasta_seqs = SeqIO.parse(open('contigs_ex.fasta'), 'fasta')

# # Create new fasta file that contains seqs only found by the program
# with open("program_only.fasta", "w") as out_file:
#     for fasta in fasta_seqs:
#         name, sequence = fasta.id, str(fasta.seq)
#         if name in record_ids:
#             SeqIO.write(fasta, out_file, "fasta")

# # # Re-run Abricate with lower coverage: abricate --db resfinder --quiet --mincov 10 program_only.fasta > program_only_ab10_results.csv
# # # run_command(['./run_abricate_10.sh'])

# # Load files
# program_only_ab10_results = pd.read_csv('program_only_ab10_results.csv', delimiter='\t')

# # Apply the function to each string in the 'GENE' column
# program_only_ab10_results['prgGENE'] = program_only_ab10_results['GENE'].apply(replace_chars)

# # Merge to find similar seqs between the ones program found and abricate (10% coverage) found
# merged_program_only_ab10 = program_only[['Record ID','Target','shortTarget']].merge(program_only_ab10_results[['SEQUENCE','GENE','prgGENE']], how='outer', left_on=['Record ID','shortTarget'], right_on=['SEQUENCE','prgGENE'], indicator=True)
# # merged_program_only_ab10.to_csv('merged_program_only_ab10_both_id_and_gene.csv', index=False)

# # Print out df includes matches 
# merged_program_only_ab10_match = merged_program_only_ab10[merged_program_only_ab10['_merge'] == 'both']
# # merged_program_only_ab10_match.to_csv('merged_program_only_ab10_match.csv', index=False)

# # Print out df includes non-matches 
# merged_program_only_ab10_nonmatch = merged_program_only_ab10[merged_program_only_ab10['_merge'].isin(['left_only', 'right_only'])]
# # merged_program_only_ab10_nonmatch.to_csv('merged_program_only_ab10_nonmatch.csv', index=False)

# # ## Scenario 3: found in abricate but not in the program- extracting contigs using Unix, identifying associated primers, and then using NCBI for alignment of the primers with the extracted contigs
# # # Put the list of seqs that needed to be compare into a .txt file: abricateOnly.txt
# # # seqtk subseq contigs_ex.fasta abricateOnly.txt > abricateOnly.fasta # This will print out all the seqs for the comparison 

# # Create 'Reason' column in abricate_only df
# # If 'prgGENE' is found in 'filtered_primer', fill with 'Mismatches', else 'Primers not found'
# abricate_only = abricate_only[['SEQUENCE','GENE','prgGENE']]
# abricate_only['Reason'] = abricate_only['prgGENE'].apply(lambda x: 'Mismatches' if x in filtered_primer['shortTarget'].values else 'Primers not found')
# # abricate_only.to_csv('abricate_only_and_primers.csv', index=False)

# # Filter out rows with 'Reason' == 'Mismatches' in abricate_only dataframe
# mismatches_primers_df = abricate_only[abricate_only['Reason'] == 'Mismatches']

# # # Get corresponding 'F_truseq' and 'R_truseq' from filtered_primer dataframe for each 'prgGENE' in mismatches_primers_df
# mismatches_primers_df = mismatches_primers_df.merge(filtered_primer[['shortTarget', 'F_truseq', 'R_truseq']], left_on='prgGENE', right_on='shortTarget', how='inner')

# # Write all the primers in a single line in the text file
# with open('mismatches_primers.txt', 'w') as file:
#     mismatches_primers = ', '.join([f"{row['F_truseq']}, {row['R_truseq']}" for _, row in mismatches_primers_df.iterrows()])
#     file.write(mismatches_primers)

# # Create a set of all record IDs from the 'Record ID' column 
# record_ids = mismatches_primers_df['SEQUENCE'].tolist()

# # Read the contigs
# fasta_seqs = SeqIO.parse(open('contigs_ex.fasta'), 'fasta')

# # Create new fasta file that contains only seqs found in program
# with open("abricate_only.fasta", "w") as out_file:
#     for fasta in fasta_seqs:
#         name, sequence = fasta.id, str(fasta.seq)
#         if name in record_ids:
#             SeqIO.write(fasta, out_file, "fasta")

# # # Run the primers supposed to be found by the program against the contigs to check for mismatches 
# run_command(['./blast.sh'])

# #############################################################################################################