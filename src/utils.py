
import os

import numpy as np
import pandas as pd


def make_matrix(
        quant_dirs,
        quant_main_dir,
        out_dir,
        quant_type="rsem"):
    """merge quant files
    """
    if quant_type == "rsem":
        quant_file = "Quant.genes.results"

    for file_idx in range(len(quant_dirs)):
        quant_dir = quant_dirs[file_idx]
        full_path_file = "{}/{}/{}".format(quant_main_dir, quant_dir, quant_file)
        col_header = "{}.{}".format(
            quant_dir.split(".")[0].split("-")[-1],
            quant_dir.split(".")[1].split("_")[0])

        # read in with pandas
        file_data = pd.read_csv(full_path_file, sep="\t")

        # make a tpm file
        tpms = file_data[["gene_id", "TPM"]]
        tpms.columns = ["gene_id", col_header]
        if file_idx == 0:
            all_tpms = tpms.copy()
        else:
            all_tpms = all_tpms.merge(tpms, on="gene_id")
        
        # make a counts file
        counts = file_data[["gene_id", "expected_count"]]
        counts.columns = ["gene_id", col_header]
        if file_idx == 0:
            all_counts = counts.copy()
        else:
            all_counts = all_counts.merge(counts, on="gene_id")

    # move gene ids to index and save out
    all_tpms["gene_id"] = all_tpms["gene_id"].str.split(".").str[0]
    all_tpms = all_tpms.set_index("gene_id")
    all_counts["gene_id"] = all_counts["gene_id"].str.split(".").str[0]
    all_counts = all_counts.set_index("gene_id").astype(int)

    # save out
    tpm_file = "{}/tpms.mat.txt.gz".format(out_dir)
    all_tpms.to_csv(tpm_file, sep="\t", compression="gzip")
    counts_file = "{}/counts.mat.txt.gz".format(out_dir)
    all_counts.to_csv(counts_file, sep="\t", compression="gzip")
    
    return tpm_file, counts_file


def filter_by_ids(matrix_file, filter_list_file):
    """filter file if rowname is filter list file
    """
    # load files
    keep_ids = pd.read_csv(filter_list_file, header=None).iloc[:,0].values
    data = pd.read_csv(matrix_file, index_col=0, sep="\t")

    # filter
    data = data[data.index.isin(keep_ids)]

    # save out
    pc_file = "{}.pc.mat.txt.gz".format(matrix_file.split(".mat")[0])
    data.to_csv(pc_file, sep="\t")

    return pc_file


def filter_for_expressed(mat_pc_files, threshold=1):
    """vis in R (to confirm empirical cutoff) and cut here
    threshold is 1 in the log2(TPM) space.
    """
    # plot to confirm
    plot_file = "{}.log2.expr_distr.pdf".format(mat_pc_files[0].split(".mat")[0])
    plot_cmd = "plot.gene_expr_distr.R {} {}".format(mat_pc_files[0], plot_file)
    print plot_cmd
    os.system(plot_cmd)

    # load data
    tpm_data = pd.read_csv(mat_pc_files[0], sep="\t", index_col=0)
    count_data = pd.read_csv(mat_pc_files[1], sep="\t", index_col=0)

    # filter
    tpm_data_log2 = np.log2(tpm_data)
    tpm_data_filt = tpm_data[np.max(tpm_data_log2.values, axis=1) >= threshold]
    count_data_filt = count_data[np.max(tpm_data_log2.values, axis=1) >= threshold]

    # save out
    tpm_filt_file = "{}.expr_filt.mat.txt.gz".format(mat_pc_files[0].split(".mat")[0])
    tpm_data_filt.to_csv(tpm_filt_file, sep="\t")
    count_filt_file = "{}.expr_filt.mat.txt.gz".format(mat_pc_files[1].split(".mat")[0])
    count_data_filt.to_csv(count_filt_file, sep="\t")
    
    return tpm_filt_file, count_filt_file


def run_sequential_deseq2(counts_file, out_prefix, fdr_cutoff=0.20):
    """go to R
    """
    diff_cmd = "deseq2.sequential.R {} {} {}".format(
        counts_file,
        fdr_cutoff,
        out_prefix)
    print diff_cmd
    os.system(diff_cmd)
    
    return


def make_gene_signature_file(
        deseq_results_file,
        gene_mappings_file,
        out_file,
        background_file=None,
        filt_file=None,
        sort_ascending=True):
    """take deseq2 results file and condense to relevant info
    """
    # pull in results, filter if needed
    deseq2_results = pd.read_csv(deseq_results_file, sep="\t")
    if filt_file is not None:
        filt_list = pd.read_csv(filt_file, header=None).iloc[:,0].values
        deseq2_results = deseq2_results[deseq2_results.index.isin(filt_list)]
    deseq2_results["ensembl_gene_id"] = deseq2_results.index.values
        
    # pull in mapping ids
    mappings = pd.read_csv(gene_mappings_file, sep="\t")
    results = deseq2_results.merge(mappings, on="ensembl_gene_id", how="left")
    
    # keep only needed columns
    results = results[[
        "ensembl_gene_id",
        "hgnc_symbol",
        "log2FoldChange",
        "lfcSE",
        "padj"]]

    # reduce duplicates
    results = results.drop_duplicates()

    # sort by logFC
    results = results.sort_values("log2FoldChange", ascending=sort_ascending)

    # save out
    results.to_csv(out_file, sep="\t", compression="gzip", index=False)

    # check each signature with gprofiler
    if background_file is not None:
        gprofiler_cmd = "run_gprofiler.R {} {} {}".format(
            out_file, background_file, os.path.dirname(out_file))
        print gprofiler_cmd
        os.system(gprofiler_cmd)
    
    return
