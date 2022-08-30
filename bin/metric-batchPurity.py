#!/usr/bin/env python

"""
Evaluate integration using a cluster batch purity metric (NOTE: This is an example only and should not be used for a real analysis)

Usage:
    metric-batchPurity.py --dataset=<str> --method=<str> --integration=<str> --out-file=<path> [options] <file>

Options:
    -h --help            Show this screen.
    --dataset=<str>      Name of the dataset to calculate the metric for.
    --method=<str>       Name of the method to calculate the metric for.
    --integration=<str>  Name of the integration to calculate the metric for.
    --out-file=<path>    Path to output file.
"""


def calculate_batch_purity(adata):
    """
    Calculate the batch purity metric for an integrated dataset

    Parameters
    ----------
    adata
        AnnData object containing the integrated dataset

    Returns
    -------
    The batch purity score
    """
    from scanpy.preprocessing import neighbors
    from scanpy.tools import leiden

    print("Calculating nearest neighbours...")
    neighbors(adata, use_rep="X_emb")

    print("Calculating clusters...")
    leiden(adata, key_added="Cluster")

    print("Calculating cluster batch purities...")
    cluster_purities = []
    for cluster in adata.obs["Cluster"].cat.categories:
        adata_cluster = adata[adata.obs["Cluster"] == cluster]
        counts = adata_cluster.obs["Batch"].value_counts().tolist()
        props = [count / adata_cluster.n_obs for count in counts]
        purity = max(props)
        print(f"Cluster {cluster}: {purity}")
        cluster_purities.append(purity)

    print("Calculating final score...")
    # Subtract from 1 so that 0 is the worst score and 1 the best
    score = 1 - (sum(props) / len(props))
    print(f"Final score: {score}")

    return score


def main():
    """The main script function"""
    from docopt import docopt
    from scanpy import read_h5ad
    from _functions import format_metric_results

    args = docopt(__doc__)

    file = args["<file>"]
    dataset = args["--dataset"]
    method = args["--method"]
    integration = args["--integration"]
    out_file = args["--out-file"]

    print(f"Reading data from '{file}'...")
    input = read_h5ad(file)
    print("Read data:")
    print(input)
    score = calculate_batch_purity(input)
    output = format_metric_results(
        dataset, method, integration, "Integration", "BatchPurity", score
    )
    print(output)
    print(f"Writing output to '{out_file}'...")
    output.to_csv(out_file, sep="\t", index=False)
    print("Done!")


if __name__ == "__main__":
    main()
