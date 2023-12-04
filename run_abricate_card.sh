#!/bin/bash

module load StdEnv/2020  gcc/9.3.0 abricate/1.0.0

abricate --db card --quiet --mincov 80 --minid 80 program_only.fasta > abricate_results_card.csv 

