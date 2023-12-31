#!/bin/bash

module load nixpkgs/16.09  gcc/5.4.0 blast+/2.6.0

# Run blastn with similarities of at least 10 nt, gap allowed 5 nt 
blastn -query mismatches_primers.txt -subject abricate_only.fasta -word_size 10 -evalue 0.01 -gapopen 5 -gapextend 2 -out abricate_only_blast_report.txt 