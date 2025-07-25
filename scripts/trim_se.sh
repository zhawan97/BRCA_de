#!/bin/bash


# File: trim.sh
# Author: Zain Awan
# Purpose: run trimmomatic SE on input file and output to specified directory

usage() {
    echo "Usage: $0 <reads.fastq> <output_dir>"
    exit 1
}

# Get args

READS=$1
OUTPUT_DIR=$2

# Check args

if [[ -z "$READS" || -z "$OUTPUT_DIR" ]]; then
    usage
fi

if [[ ! -f "$READS" ]]; then
    echo "ERROR: Reads file not found: $READS"
    exit 1
fi

if [[ ! -d "$OUTPUT_DIR" ]]; then
    echo "Creating output directory: $OUTPUT_DIR"
    mkdir -p "$OUTPUT_DIR"
fi


BASE=$(basename "$READS" .fastq)
OUTPUT_FILE="$OUTPUT_DIR/${BASE}_trimmed.fastq"

# Run Trimmomatic

ADAPTERS=/opt/anaconda3/envs/bioinfo/share/trimmomatic-0.39-2/adapters/TruSeq3-PE.fa
echo "Running Trimmomatic SE..."
trimmomatic SE -phred33 \
    "$READS" \
    "$OUTPUT_FILE" \
    ILLUMINACLIP:$ADAPTERS:2:30:10 \
    LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36

if [[ -f "$OUTPUT_FILE" ]]; then
    echo "Trimmed reads succesfully!"
fi
