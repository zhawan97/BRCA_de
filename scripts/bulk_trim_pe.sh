#!/bin/bash

INPUT="data"
OUT_PAIRED="data/trimmed/paired"
OUT_UNPAIRED="data/trimmed/unpaired"
ADAPTERS="/opt/anaconda3/envs/bioinfo/share/trimmomatic-0.39-2/adapters/TruSeq3-PE.fa"

mkdir -p data/trimmed
mkdir -p data/trimmed/paired
mkdir -p data/trimmed/unpaired

# Loop over forward reads
for R1 in "$INPUT"/*_1.fastq; do

    # Get Matching reverse read
    R2="${R1/_1.fastq/_2.fastq}"

    # Extract basename
    BASENAME=$(basename "$R1" _1.fastq)

    OUT1="$OUT_PAIRED/${BASENAME}_1_paired.fastq"
    OUT1_U="$OUT_UNPAIRED/${BASENAME}_1_unpaired.fastq"
    OUT2="$OUT_PAIRED/${BASENAME}_2_paired.fastq"
    OUT2_U="$OUT_UNPAIRED/${BASENAME}_2_unpaired.fastq"


    # Run Trimmomatic
    trimmomatic PE -threads 4 -phred33 \
        "$R1" "$R2" \
        "$OUT1" "$OUT1_U" \
        "$OUT2" "$OUT2_U" \
        ILLUMINACLIP:$ADAPTERS:2:30:10 \
        LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

    # Run fastqc
    fastqc -o results/FASTQC/trimmed "$OUT1"
    fastqc -o results/FASTQC/trimmed "$OUT2"

done

