#!/bin/bash


# File: trim.sh
# Author: Zain Awan
# Purpose: run trimmomatic on input file and output to specified directory

usage() {
    echo "Usage: $0 <forward_reads.fastq> <reverse_reads.fastq> <output_dir>"
    exit 1
}

# Get args

FORWARD_READS=$1
REVERSE_READS=$2
OUTPUT_DIR=$3

# Check args

if [[ -z "$FORWARD_READS" || -z "$REVERSE_READS" || -z "$OUTPUT_DIR" ]]; then
    usage
fi

if [[ ! -f "$FORWARD_READS" ]]; then
    echo "ERROR: Forward reads file not found: $FORWARD_READS"
    exit 1
fi

if [[ ! -f "$REVERSE_READS" ]]; then
    echo "ERROR: Reverse reads file not found: $REVERSE_READS"
    exit 1
fi

if [[ ! -d "$OUTPUT_DIR" ]]; then
    echo "Creating output directory: $OUTPUT_DIR"
    mkdir -p "$OUTPUT_DIR"
fi

# Get basenames

BASE1=$(basename "$FORWARD_READS" .fastq)
BASE2=$(basename "$REVERSE_READS" .fastq)

# Output files

PAIRED1="$OUTPUT_DIR/${BASE1}_paired.fastq"
UNPAIRED1="$OUTPUT_DIR/${BASE1}_unpaired.fastq"
PAIRED2="$OUTPUT_DIR/${BASE2}_paired.fastq"
UNPAIRED2="$OUTPUT_DIR/${BASE2}_unpaired.fastq"

# Run Trimmomatic
echo "Running Trimmomatic PE..."
trimmomatic PE -phred33 \
    "$FORWARD_READS" "$REVERSE_READS" \
    "$PAIRED1" "$UNPAIRED1" \
    "$PAIRED2" "$UNPAIRED2" \
    LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36

if [[ -f "$PAIRED1" && "$PAIRED2" && "$UNPAIRED1" && "$UNPAIRED2" ]]; then
    echo "Trimmed reads succesfully"
fi
