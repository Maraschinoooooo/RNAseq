#!/bin/sh
#
# RNA-Seq Data Analysis Workflow

# ===============================
# Step 1: Download SRA data
# ===============================
# Example** Download SRA data for experiment SRX16830050
# $prefetch SRX16830050
# $fastq-dump SRX16830050
#
# Example**download fastq.gz data. (By adding the --gzip option, the fastq-dump command will generate compressed FASTQ files (fq.gz) instead of the uncompressed ones.)
$fastq-dump --gzip SRR20810434

# Example**take a peek at first or last few lines of the contents 
$zcat SRR20810434.fastq.gz | head
$tail SRR20810434.fastq

# Example**use the vdb-dump tool from the SRA Toolkit to get information about the run
$vdb-dump --info SRR20810434


# ===============================
# Step 2: Split and convert SRA to fastq
# ===============================
$fastq-dump --gzip --split-3 $sample

# ===============================
# Step 3: Quality Control with fastQC
# ===============================
$fastqc "${sample}_1.fastq.gz" "${sample}_2.fastq.gz"

# ===============================
# Step 4: Reads mapping with hisat2
# ===============================
$hisat2 -x $genome -1 "${sample}_1.fastq.gz" -2 "${sample}_2.fastq.gz" -S "${sample}_alignment.sam"

# ===============================
# Step 5: Transfer SAM to BAM, sort, and index
# ===============================
$samtools view -S "${sample}_alignment.sam" -b | samtools sort -@ 8 -o "${sample}_sort.bam" && samtools index "${sample}_sort.bam"

# ===============================
# Step 6: FeatureCount - Count the reads number
# ===============================
$featureCounts -p -t -g gene_id -M -T 8 -a $genomegt -o all.featurecounts.txt ./*_sort.bam

# ===============================
# Step 7: MultiQC
# ===============================
$multiqc all.featurecounts.txt.summary -o all.counts.summary













