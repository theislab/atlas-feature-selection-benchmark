#!/usr/bin/env python

"""
Calculate the cell mapping distance metric

Usage:
    metric-cellDist.py --dataset=<str> --method=<str> --integration=<str> --reference=<path> --out-file=<path> [options] <file>

Options:
    -h --help            Show this screen.
    --dataset=<str>      Name of the dataset to calculate the metric for.
    --method=<str>       Name of the method to calculate the metric for.
    --integration=<str>  Name of the integration to calculate the metric for.
    --reference=<path>   Path to the reference H5AD file.
    --out-file=<path>    Path to output file.
"""


def calculate_unseen_cell_distance(reference, query):
    """
    Calculate the cell Mahalanobis unseen population distance metric

    Parameters
    ----------
    reference
        H5AD object containing the integrated reference dataset
    query
        H5AD object containing the mapped query dataset

    Returns
    -------
    The cell mapping distance score
    """

    from scipy.spatial.distance import mahalanobis
    from functions.distances import get_inverse_covariances, get_centroids
    from numpy import mean, quantile

    reference_coords = reference.obsm["X_emb"]
    reference_labels = reference.obs["Label"].tolist()

    print("Calculating inverse covariance matrices...")
    inverse_covariances = get_inverse_covariances(reference_coords, reference_labels)
    print("Calculating centroid positions...")
    centroids = get_centroids(reference_coords, reference_labels)

    print("Calculating reference cell 90th quantile distances...")
    ref_quantiles = {}
    for label in set(reference_labels):
        label_distances = []
        label_reference = reference[reference.obs["Label"] == label].copy()

        for idx in range(label_reference.n_obs):
            label_coord = label_reference.obsm["X_emb"][idx, :]

            distance = mahalanobis(
                label_coord, centroids[label], inverse_covariances[label]
            )
            label_distances.append(distance)

        ref_quantiles[label] = quantile(label_distances, 0.9)

    print("Selecting unseen query cells...")
    query = query[query.obs["Unseen"], :].copy()

    print("Calculating cells outside 90th quantile...")
    is_outside = []
    for idx in range(query.n_obs):
        query_label = query.obs["Label"][idx]
        query_coord = query.obsm["X_emb"][idx, :]

        min_distance = 1e6
        min_label = None
        for label in set(reference_labels):
            distance = mahalanobis(
                query_coord, centroids[label], inverse_covariances[label]
            )
            if distance < min_distance:
                min_distance = distance
                min_label = label

        is_outside.append(min_distance > ref_quantiles[min_label])

    print("Calculating final score...")
    # The score is the proportion of unseen cells outside the 90th quantile
    score = mean(is_outside)

    return score


def main():
    """The main script function"""
    from docopt import docopt
    from scanpy import read_h5ad
    from functions.metrics import format_metric_results

    args = docopt(__doc__)

    query_file = args["<file>"]
    dataset = args["--dataset"]
    method = args["--method"]
    integration = args["--integration"]
    reference_file = args["--reference"]
    out_file = args["--out-file"]

    print(f"Reading query data from '{query_file}'...")
    query = read_h5ad(query_file)
    query.obs["Dataset"] = "Query"
    print("Read data:")
    print(query)
    print(f"Reading reference data from '{reference_file}'...")
    reference = read_h5ad(reference_file)
    reference.obs["Dataset"] = "Reference"
    print("Read data:")
    print(reference)
    score = calculate_unseen_cell_distance(reference, query)
    output = format_metric_results(
        dataset, method, integration, "Unseen", "unseenCellDist", score
    )
    print(output)
    print(f"Writing output to '{out_file}'...")
    output.to_csv(out_file, sep="\t", index=False)
    print("Done!")


if __name__ == "__main__":
    main()
