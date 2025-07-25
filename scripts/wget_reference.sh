#!/bin/usr bash


# Download reference and GTF into data/reference

cd data/reference

wget -O Homo_sapiens.GRCh38.transcriptome.fa.gz \
ftp://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz
gunzip Homo_sapiens.GRCh38.cdna.all.fa.gz

wget -O Homo_sapiens.GRCh38.110.gtf.gz \
ftp://ftp.ensembl.org/pub/release-110/gtf/homo_sapiens/Homo_sapiens.GRCh38.110.gtf.gz
gunzip Homo_sapiens.GRCh38.110.gtf.gz

# Create tx2gene.tsv
gawk '$3 == "transcript" {
    match($0, /transcript_id "([^"]+)"/, tx);
    match($0, /gene_id "([^"]+)"/, gene);
    if (tx[1] && gene[1]) print tx[1] "\t" gene[1];
}' Homo_sapiens.GRCh38.110.gtf > tx2gene.tsv


# Return to parent directory
cd ../..
