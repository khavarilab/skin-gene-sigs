#!/usr/bin/env Rscript

# plot distribution of gene expr values
library(reshape2)
library(ggplot2)

args <- commandArgs(trailingOnly=TRUE)
tpm_file <- args[1]
plot_file <- args[2]

# read data
tpms <- read.table(gzfile(tpm_file), header=TRUE, row.names=1)

# remove zeros and convert to log2
tpms <- tpms[rowSums(tpms) != 0,]
print(dim(tpms))
tpms_log2 <- log2(tpms)
is.na(tpms_log2) <- sapply(tpms_log2, is.infinite)
tpms_log2[is.na(tpms_log2)] <- 0
tpms_flat <- melt(tpms_log2)

# plot
ggplot(tpms_flat, aes(x=value, colour=variable)) + geom_density() + geom_vline(xintercept=1)
ggsave(plot_file)
