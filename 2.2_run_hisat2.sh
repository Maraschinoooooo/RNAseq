#!/bin/bash

# Define an array of sample names
samples=("F1" "F3" "L_backup" "L1" "L2")

# Loop through each sample and run the hisat2 command
for sample in "${samples[@]}"; do
  hisat2 -x /grch38_snp_tran/genome_snp_tran \
    -1 "231_${sample}_1.fq.gz" \
    -2 "231_${sample}_2.fq.gz" \
    -S "${sample}_alignment.sam"
done