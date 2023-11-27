#!/bin/bash

module load StdEnv/2020  gcc/9.3.0 abricate/1.0.0

# Take the fasta file name from the command line argument
fasta_file="$1"

# Run abricate with the provided fasta file
abricate --db resfinder --quiet --mincov 80 --minid 80 "$fasta_file" > abricate_results.csv


