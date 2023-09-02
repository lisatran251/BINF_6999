#!/bin/bash

module load StdEnv/2020  gcc/9.3.0 abricate/1.0.0

abricate --db resfinder --quiet --mincov 10 originalProduct.fasta | awk 'BEGIN {OFS = ","} {$1=$1; print}' > res_10cov.csv


