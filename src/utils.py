
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
    all_tpms = all_tpms.set_index("gene_id")
    all_counts = all_counts.set_index("gene_id").astype(int)

    # save out
    tpm_file = "{}/tpms.mat.txt.gz".format(out_dir)
    all_tpms.to_csv(tpm_file, sep="\t", compression="gzip")
    counts_file = "{}/counts.mat.txt.gz".format(out_dir)
    all_counts.to_csv(counts_file, sep="\t", compression="gzip")
    
    return
