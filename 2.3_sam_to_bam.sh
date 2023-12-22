#!/bin/bash

# Define an array of sample names
samples=("SRR20810435" "SRR20810436" "F_backup" "F1" "F3" "L_backup" "L1" "L2")

# Loop through each sample and run the samtools commands
for sample in "${samples[@]}"; do
  samtools view -S "${sample}_alignment.sam" -b | \
    samtools sort -@ 8 -o "${sample}_sort.bam" | \
    samtools index - "${sample}_sort.bam.index"
done




# $ samtools view -S SRR20810434_alignment.sam -b | samtools sort -@ 8 -o SRR20810434_sort.bam | samtools index - SRR20810434_sort.bam.index

# $ samtools view -S SRR20810435_alignment.sam -b | samtools sort -@ 8 -o SRR20810435_sort.bam | samtools index - SRR20810435_sort.bam.index

# $ samtools view -S SRR20810436_alignment.sam -b | samtools sort -@ 8 -o SRR20810436_sort.bam | samtools index - SRR20810436_sort.bam.index

# $ samtools view -S F_backup_alignment.sam -b | samtools sort -@ 8 -o F_backup_sort.bam | samtools index - F_backup.bam.index

# $ samtools view -S F1_alignment.sam -b | samtools sort -@ 8 -o F1_sort.bam | samtools index - F1.bam.index

# $ samtools view -S F3_alignment.sam -b | samtools sort -@ 8 -o F3_sort.bam | samtools index - F3.bam.index

# $ samtools view -S L_backup_alignment.sam -b | samtools sort -@ 8 -o L_backup_sort.bam | samtools index - L_backup.bam.index

# $ samtools view -S L1_alignment.sam -b | samtools sort -@ 8 -o L1_sort.bam | samtools index - L1.bam.index

# $ samtools view -S L2_alignment.sam -b | samtools sort -@ 8 -o L2_sort.bam | samtools index - L2.bam.index