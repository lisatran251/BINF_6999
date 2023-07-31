#!/bin/bash

module load StdEnv/2020  gcc/9.3.0 abricate/1.0.0

abricate --db resfinder --quiet contigs_ex.fasta > abricate_results.csv 

