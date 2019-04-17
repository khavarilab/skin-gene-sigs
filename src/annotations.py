

import os


def get_proteincoding_gene_list_from_gtf(gtf_file, out_file):
    """Using a gtf file of protein coding genes, extract ids
    Note that it leaves the decimals on there

    Args:
      gtf_file: gtf annotation file with protein coding annotations
      out_file: output file of gene ids
    """
    assert gtf_file.endswith("gtf")
    assert out_file.endswith(".gz")
    extract_gene_ids = (
        "cat {0} | "
        "grep 'protein_coding' | "
        "awk -F '\"' '{{ print $2 }}' | "
        "awk -F '.' '{{ print $1 }}' | "
        "sort | "
        "uniq | "
        "gzip -c > {1}").format(
            gtf_file,
            out_file)
    os.system(extract_gene_ids)
    
    return None



def get_ensembl_to_geneid_mapping(ensembl_id_file, out_mapping_file):
    """uses biomart to produce a mapping table
    """
    get_mappings = ("annot.ensembl_to_mappings.R {} {}").format(
        ensembl_id_file, out_mapping_file)
    os.system(get_mappings)
    
    return None

