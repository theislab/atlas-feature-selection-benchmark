#!/usr/bin/env python

"""
Evaluate integration using cell cycle conservation

Usage:
    metric-cellCycle.py --dataset=<str> --method=<str> --integration=<str> --exprs=<file> --cc-genes=<file> --out-file=<path> <file>

Options:
    -h --help            Show this screen.
    --dataset=<str>      Name of the dataset to calculate the metric for.
    --method=<str>       Name of the method to calculate the metric for.
    --integration=<str>  Name of the integration to calculate the metric for.
    --exprs=<file>       Path to H5AD file containing the expression matrix.
    --cc-genes=<file>    Path to TSV file containing cell cycle genes.
    --out-file=<path>    Path to output file.
"""


def calculate_cellCycleConservation(adata, exprs, s_genes, g2m_genes):
    """
    Calculate the cell cycle conservation score for an integrated dataset.

    Parameters
    ----------
    adata
        AnnData object containing the integrated dataset
    exprs
        AnnData object containing the unintegrated dataset with the expression matrix
    s_genes
        List containing S phase genes
    g2m_genes
        List containing G2M phase genes

    Returns
    -------
    Cell cycle conservation score
    """
    from scib.metrics import cell_cycle
    from scanpy.tools import score_genes_cell_cycle

    if adata.uns["Species"].lower() not in ["human", "mouse"]:
        from warnings import warn

        warn(
            f"'{adata.uns['Species']}' is not a valid species ('human', 'mouse'). A score of 'NA' will be returned."
        )
        return "NA"

    print("Calculating batch cell cycle scores...")
    batches = exprs.obs["Batch"].unique()
    for batch in batches:
        print(f"Batch '{batch}'")
        exprs_batch = exprs[exprs.obs["Batch"] == batch].copy()
        score_genes_cell_cycle(exprs_batch, s_genes, g2m_genes)
        exprs.obs.loc[exprs_batch.obs_names, "S_score"] = exprs_batch.obs["S_score"]
        exprs.obs.loc[exprs_batch.obs_names, "G2M_score"] = exprs_batch.obs["G2M_score"]
        adata.obs.loc[exprs_batch.obs_names, "S_score"] = exprs_batch.obs["S_score"]
        adata.obs.loc[exprs_batch.obs_names, "G2M_score"] = exprs_batch.obs["G2M_score"]

    print("Calculating cell cycle conservation score...")
    score = cell_cycle(
        adata_pre=exprs,
        adata_post=adata,
        batch_key="Batch",
        embed="X_emb",
        organism=adata.uns["Species"].lower(),
        recompute_cc=False,
        verbose=True,
    )
    print(f"Final score: {score}")

    return score


def read_cellcycle_genes(file, var_names, species):
    """
    Calculate the cell cycle conservation score for an integrated dataset.

    Parameters
    ----------
    file
        Path to TSV file containing cell cycle genes
    var_names
        List of the variable names for the dataset that will be scored
    species
        The species of the data to be scored, either 'human' or 'mouse'

    Returns
    -------
    A tuple of lists of S phase genes and G2M phase genes
    """

    from pandas import read_csv

    print(f"Reading cell cycle genes from '{file}'...")
    genes_df = read_csv(file, sep="\t")

    print(f"Selecting genes for species '{species}'...")
    if species.lower() == "human":
        genes_df = genes_df[genes_df["Species"] == "Human"]
        ensembl_prefix = "ENSG"
    elif species.lower() == "mouse":
        genes_df = genes_df[genes_df["Species"] == "Mouse"]
        ensembl_prefix = "ENSMUSG"
    else:
        from warnings import warn

        warn(
            f"'{species}' is not a valid species ('human', 'mouse'). Returning empty gene sets."
        )
        return ([], [])

    is_ensembl = all(var_name.startswith(ensembl_prefix) for var_name in var_names)
    if is_ensembl:
        print(f"ENSEMBL var names detected, using ENSEMBL IDs")
        id_col = "ENSEMBL"
    else:
        print("ENSEMBL var names not found, using gene names")
        id_col = "Gene"

    print("Selecting gene sets...")
    genes_df = genes_df.dropna(subset=[id_col])
    s_genes = list(genes_df[id_col][genes_df["Phase"] == "S"])
    s_genes = [gene for gene in s_genes if gene in var_names]
    print(f"Found {len(s_genes)} 'S' phase genes")
    g2m_genes = list(genes_df[id_col][genes_df["Phase"] == "G2M"])
    g2m_genes = [gene for gene in g2m_genes if gene in var_names]
    print(f"Found {len(g2m_genes)} 'G2M' phase genes")

    return (s_genes, g2m_genes)


def main():
    """The main script function"""
    from docopt import docopt
    from scanpy import read_h5ad
    from functions.metrics import format_metric_results

    args = docopt(__doc__)

    file = args["<file>"]
    dataset = args["--dataset"]
    method = args["--method"]
    integration = args["--integration"]
    exprs_file = args["--exprs"]
    cc_file = args["--cc-genes"]
    out_file = args["--out-file"]

    print(f"Reading data from '{file}'...")
    input = read_h5ad(file)
    print("Read data:")
    print(input)
    print(f"Reading expression data from '{exprs_file}'...")
    exprs = read_h5ad(exprs_file)
    print("Read expression data:")
    print(exprs)
    s_genes, g2m_genes = read_cellcycle_genes(
        cc_file, exprs.var_names, exprs.uns["Species"]
    )
    score = calculate_cellCycleConservation(input, exprs, s_genes, g2m_genes)
    output = format_metric_results(
        dataset, method, integration, "IntegrationBio", "CellCycle", score
    )
    print(output)
    print(f"Writing output to '{out_file}'...")
    output.to_csv(out_file, sep="\t", index=False)
    print("Done!")


if __name__ == "__main__":
    main()
