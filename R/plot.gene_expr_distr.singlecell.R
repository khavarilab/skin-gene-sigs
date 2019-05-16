#!/usr/bin/env Rscript

# plot distribution of gene expr values
library(reshape2)
library(ggplot2)

args <- commandArgs(trailingOnly=TRUE)
tpm_file <- args[1]
plot_file <- args[2]

# read data
tpms <- read.table(tpm_file, header=TRUE, row.names=1)

# remove zeros and convert to log2
tpms <- tpms[rowSums(tpms) != 0,]
print(dim(tpms))
#tpms_log2 <- log2(tpms)
#is.na(tpms_log2) <- sapply(tpms_log2, is.infinite)
#tpms_log2[is.na(tpms_log2)] <- 0
#tpms <- tpms_log2
tpms_flat <- melt(tpms)

# plot
x_line <- 0.0001
ggplot(tpms_flat, aes(x=value, colour=variable)) + geom_density() + geom_vline(xintercept=x_line) + xlim(0,0.005)
ggsave(plot_file)


# if cut here, what's left
tpms <- tpms[apply(tpms, 1, max) > 0.0001,]
print(dim(tpms))
