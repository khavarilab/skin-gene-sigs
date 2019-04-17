#!/usr/bin/env Rscript

# Description: calculates deseq2 pairwise sequentially
# through a count matrix

library(DESeq2)

print('running DESeq2...')

# folders files etc
args <- commandArgs(trailingOnly=TRUE)
in_file <- args[1]
fdr_cutoff <- as.numeric(args[2])
prefix <- args[3]
out_file <- args[4]

# make count matrix
count_data <- as.matrix(read.table(
    gzfile(in_file),
    sep='\t',
    header=TRUE,
    row.names=1))

# run pairs sequentially (assumes count matrix is ordered)
# also assumes single replicates only
for (col_idx in seq(2, ncol(count_data))) {

    # make smaller count matrix
    desired_columns <- c(1, col_idx)
    print(desired_columns)
    pairwise_count_data <- count_data[,desired_columns]
    
    # make condition table
    conditions <- colnames(pairwise_count_data)
    unique_conditions <- unique(conditions)
    cond_baseline <- unique_conditions[1]
    cond_compare <- unique_conditions[2]
    cond_table <- data.frame(condition=conditions)
    rownames(cond_table) <- colnames(pairwise_count_data)
    
    # make DESeq2 dataset
    dds <- DESeqDataSetFromMatrix(
        countData=pairwise_count_data,
        colData=cond_table,
        design = ~ condition)
    
    # build DESeq model
    dds <- DESeq(dds)

    # get results and save out to tables
    res <- results(
        dds,
        contrast=c('condition', cond_compare, cond_baseline),
        alpha=fdr_cutoff)

    # remove NA and filter for cutoff
    res_noNA <- res[!is.na(res$padj),]
    res_filt <- res_noNA[res_noNA$padj < fdr_cutoff,]

    # separate up and down sets
    res_filt_up <- res_filt[res_filt$log2FoldChange > 0,]
    res_filt_down <- res_filt[res_filt$log2FoldChange < 0,]

    # adjust for file names
    cond_compare <- gsub(".", "-", cond_compare, fixed=TRUE)
    cond_baseline <- gsub(".", "-", cond_baseline, fixed=TRUE)

    # write out all results (incl non significant)
    write.table(
        res_noNA,
        file=gzfile(
            paste(prefix, '.', cond_compare, '_over_', cond_baseline, '_resultsAll.txt.gz', sep='')),
        quote=FALSE, sep='\t')

    # write out all sig
    write.table(
        res_filt,
        file=gzfile(
            paste(prefix, '.', cond_compare, '_over_', cond_baseline, '_sigResultsAll.txt.gz', sep='')),
        quote=FALSE, sep='\t')

    # write out sig up
    write.table(
        rownames(res_filt_up),
        file=gzfile(
            paste(prefix, '.', cond_compare, '_over_', cond_baseline, '_sigResultsUp.txt.gz', sep='')),
        quote=FALSE, row.names=FALSE, col.names=FALSE,
        sep='\t')

    # write out sig down
    write.table(
        rownames(res_filt_down),
        file=gzfile(
            paste(prefix, '.', cond_compare, '_over_', cond_baseline, '_sigResultsDown.txt.gz', sep='')),
        quote=FALSE, row.names=FALSE, col.names=FALSE,
        sep='\t')

    # require at least log2 FC > 1 and write out again
    res_filt <- res_filt[abs(res_filt$log2FoldChange) >= 1,]
    res_filt_up <- res_filt[res_filt$log2FoldChange >= 1,]
    res_filt_down <- res_filt[res_filt$log2FoldChange <= -1,]

    write.table(
        res_filt,
        file=gzfile(
            paste(prefix, '.', cond_compare, '_over_', cond_baseline, '_sigResultsAll.log2_thresh.txt.gz', sep='')),
        quote=FALSE, sep='\t')
    
    write.table(
        rownames(res_filt_up),
        file=gzfile(
            paste(prefix, '.', cond_compare, '_over_', cond_baseline, '_sigResultsUp.log2_thresh.txt.gz', sep='')),
        quote=FALSE, row.names=FALSE, col.names=FALSE,
        sep='\t')

    write.table(
        rownames(res_filt_down),
        file=gzfile(
            paste(prefix, '.', cond_compare, '_over_', cond_baseline, '_sigResultsDown.log2_thresh.txt.gz', sep='')),
        quote=FALSE, row.names=FALSE, col.names=FALSE,
        sep='\t')

}

if (FALSE) {
# now merge sigResults to get list of full regions that are significant
    merge_sigresults <- paste("zcat ",
                              dirname(prefix),
                              "/*sigResultsAll.txt.gz | ",
                              "awk -F '\t' '{ print $1 }' | ",
                              "grep -v baseMean | ",
                              "sort | ",
                              "uniq | ",
                              "gzip -c > ",
                              out_file,
                              sep="")
    system(merge_sigresults)
}
