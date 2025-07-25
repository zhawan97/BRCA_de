#!/bin/bash


usage() {
    echo "Error: Missing required args."
    echo "Usage: "$0" <reference.fasta> <path/to/indexdir> <sample.fastq> <aligned.sam>"
    echo "Ensure to include relative paths for all args."
    echo "Additional Note: The index arg should be a directory, not a file"
    exit 1
}

REFERENCE=$1
INDEX=$2
SAMPLE=$3
SAM=$4


if [[ -z "$REFERENCE" || -z "$INDEX" || -z "$SAMPLE" || -z "$SAM" ]]; then
    usage
fi

if [[ ! -d "$INDEX" ]]; then
    mkdir -p "$INDEX"
fi

bowtie2-build "$REFERENCE" "$INDEX/reference_index"

bowtie2 -x "$INDEX/reference_index" -U "$SAMPLE" -S "$SAM"


