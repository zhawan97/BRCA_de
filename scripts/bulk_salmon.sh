#!/bin/bash


INPUT="data/trimmed/paired"
INDEX="results/salmon/index"
REF="data/reference/Homo_sapiens.GRCh38.cdna.all.fa"

mkdir -p results/salmon

# Build index if not there
if [ ! -d "$INDEX" ] || [ ! -f "$INDEX"/versioninfo.json ]; then
    salmon index -t "$REF" -i "$INDEX"
else
    echo "Index already built, skipping step..."
fi

# Loop over all trimmed reads
for R1 in "$INPUT"/*_1_paired.fastq; do

    # Get basename for output fodler
    OUT_DIR=$(basename "$R1" _1_paired.fastq)

    # Get matching reverse read
    R2="${R1/_1_paired.fastq/_2_paired.fastq}"

    # Run salmon
    salmon quant -i "$INDEX" \
        -l A \
        -1 "$R1" \
        -2 "$R2" \
        -p 4 \
        -o results/salmon/"$OUT_DIR"
done

