#!/bin/bash

module load StdEnv/2020  gcc/9.3.0 abricate/1.0.0

abricate --db resfinder --quiet --mincov 10 program_only.fasta > program_only_ab10_results.csv
