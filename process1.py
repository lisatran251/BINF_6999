import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
import time
import argparse
import os
import re

# How to run:
# python3 -m venv ~/my_venv
# source ~/my_venv/bin/activate
# pip install pandas numpy biopython argparse
# python3 fasta_file primer_file 

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