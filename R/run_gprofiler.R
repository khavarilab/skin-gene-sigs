#!/usr/bin/env Rscript

# run gprofiler
# requires gene set and background gene set
library(gProfileR)

# folders files etc
args <- commandArgs(trailingOnly=TRUE)
gene_list_file <- args[1]
background_list_file <- args[2]
out_dir <- args[3]

# read in gene list and background gene list
gene_list <- read.table(gzfile(gene_list_file), sep="\t", header=TRUE, stringsAsFactors=FALSE)
gene_list <- gene_list$ensembl_gene_id

background_list <- read.table(gzfile(background_list_file), sep="\t", header=TRUE, stringsAsFactors=FALSE)
background_list <- background_list$gene_id

prefix <- sub('\\.txt.gz$', "", basename(gene_list_file))
padj_cutoff <- 0.1

# run gProfileR
results <- gprofiler(
    gene_list,
    ordered_query=TRUE,
    organism="hsapiens",
    custom_bg=background_list)
out_file <- paste(out_dir, "/", prefix, ".go_gprofiler.txt", sep="")
write.table(results, out_file, quote=FALSE, sep="\t")
