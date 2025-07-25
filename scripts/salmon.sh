#!/bin/bash

usage() {
    echo "Usage: $0 <reference.fa> <genetransfer.gtf>"
    echo "Please provide the basename of the reference and gtf"
    exit 1
}


# Get Args

REFERENCE=$1
GTF=$2

if [[ -z "$REFERENCE" || -z "$GTF" ]]; then
    usage
fi

# Create tx2gene.tsv from GTF
gawk '$3 == "transcript" { 
    match($0, /transcript_id "([^"]+)"/, tx); 
    match($0, /gene_id "([^"]+)"/, gene); 
    if (tx[1] && gene[1]) print tx[1] "\t" gene[1]; 
}' "$GTF" > tx2gene.tsv


# Run Salmon in data directory
cd data
echo "Running Salmon..."
salmon index -t reference/"$REFERENCE" -i salmon_index

salmon quant -i salmon_index -l A \
  -1 trimmed/SRR12192507_1_paired.fastq \
  -2 trimmed/SRR12192507_2_paired.fastq \
  -p 4 --validateMappings -o ../results/salmon


# Go back to parent directory
cd ..
