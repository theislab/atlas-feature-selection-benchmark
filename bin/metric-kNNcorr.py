#!/usr/bin/env python

"""
Evaluate mapping using the query kNN correlation metric

Usage:
    metric-kNNcorr.py --dataset=<str> --method=<str> --integration=<str> --exprs=<path> --out-file=<path> <file>

Options:
    -h --help            Show this screen.
    --dataset=<str>      Name of the dataset to calculate the metric for.
    --method=<str>       Name of the method to calculate the metric for.
    --integration=<str>  Name of the integration to calculate the metric for.
    --exprs=<path>       Path to the query expression H5AD file.
    --out-file=<path>    Path to output file.
"""


def calculate_kNNcorr(input, exprs):
    """
    Calculate the query kNN correlation metric

    Parameters
    ----------
    input
        AnnData object containing the integrated dataset
    exprs
        AnnData object containing the unintegrated dataset with the expression matrix

    Returns
    -------
    The kNN correlation score
    """

    from scanpy.preprocessing import normalize_total, log1p, neighbors
    from scanpy.tools import pca
    from sklearn.metrics.pairwise import euclidean_distances
    from scipy.stats import spearmanr
    import numpy as np

    # Make sure cell order is the same
    if not (input.obs_names == exprs.obs_names).all():
        exprs = exprs[input.obs_names].copy()

    exprs.obsm["X_emb"] = input.obsm["X_emb"].copy()

    n_neighbors = 100
    batches = input.obs["Batch"].unique()
    batch_scores = []
    for batch in batches:
        print(f"Calculating kNN correlation for batch '{batch}'...")
        adata_batch = exprs[exprs.obs["Batch"] == batch].copy()
        normalize_total(adata_batch, target_sum=10000)
        log1p(adata_batch)
        pca(adata_batch, n_comps=20)
        neighbors(adata_batch, n_neighbors=n_neighbors, use_rep="X_pca")

        cell_scores = []
        for cell_idx in range(adata_batch.n_obs):
            batch_distances = adata_batch.obsp["distances"][cell_idx, :]
            neighbours = np.nonzero(batch_distances)[1]
            cell_embedding = adata_batch.obsm["X_emb"][cell_idx, :]
            neighbours_embedding = adata_batch.obsm["X_emb"][neighbours, :]
            embedding_distances = euclidean_distances(
                neighbours_embedding,
                cell_embedding.reshape(1, -1)
            )
            batch_distances = batch_distances[:, neighbours].toarray()
            cor, _ = spearmanr(batch_distances.flatten(), embedding_distances.flatten())
            cell_scores.append(cor)

        # Change range from [-1, 1] to [0, 1]
        cell_scores = [(score + 1) / 2 for score in cell_scores]
        batch_scores.append(np.mean(cell_scores))

    print("Calculating final score...")
    score = np.mean(batch_scores)

    return score


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
    out_file = args["--out-file"]

    print(f"Reading data from '{file}'...")
    input = read_h5ad(file)
    print("Read data:")
    print(input)
    print(f"Reading expression data from '{exprs_file}'...")
    exprs = read_h5ad(exprs_file)
    print("Read expression data:")
    print(exprs)
    score = calculate_kNNcorr(input, exprs)
    output = format_metric_results(
        dataset, method, integration, "Mapping", "kNNcorr", score
    )
    print(output)
    print(f"Writing output to '{out_file}'...")
    output.to_csv(out_file, sep="\t", index=False)
    print("Done!")


if __name__ == "__main__":
    main()
