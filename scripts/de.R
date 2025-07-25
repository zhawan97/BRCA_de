#!/usr/bin/env Rscript

# Load required libraries
suppressMessages(library(tximport))
suppressMessages(library(readr))
suppressMessages(library(DESeq2))

# Load tx2gene and sample files
salmon_out <- "results/salmon"
tx2gene_file <- "data/reference/tx2gene.tsv"
tx2gene <- read_tsv(tx2gene_file, col_names = c("TXNAME", "GENEID"), show_col_types=FALSE)

samples <- list.files(path=salmon_out, 
                      pattern="^SRR", 
                      recursive=FALSE, 
                      include.dirs=FALSE)

print("----- Loading in the following samples for differential expression analysis -----")
print(samples)
sample_paths <- file.path(salmon_out, samples, "quant.sf")
print("----- Located the following quant.sf files -----")
print(sample_paths)
names(sample_paths) <- samples

sample_table <- data.frame(
                           sample = samples,
                           condition = c("case1", "case2", "case3", 
                                         "case4", "case5", "case6", 
                                         "case7", "case8", "case9"), 
                           row.names = samples)

txi <- tximport(sample_paths, type = "salmon", tx2gene = tx2gene, ignoreTxVersion=TRUE)

dds <- DESeqDataSetTximport(txi, colData = sample_table, design = ~ condition)

dds <- DESeq(dds)

res <- results(dds)

write.csv(as.data.frame(res), file = "results/deseq2_results.csv")

