#!/usr/bin/env python

"""
Calculate the label unseen population distance metric

Usage:
    metric-unseenLabelDist.py --dataset=<str> --method=<str> --integration=<str> --reference=<path> --out-file=<path> [options] <file>

Options:
    -h --help            Show this screen.
    --dataset=<str>      Name of the dataset to calculate the metric for.
    --method=<str>       Name of the method to calculate the metric for.
    --integration=<str>  Name of the integration to calculate the metric for.
    --reference=<path>   Path to the reference H5AD file.
    --out-file=<path>    Path to output file.
"""


def calculate_unseen_label_distance(reference, query):
    """
    Calculate the label Mahalanobis unseen population distance metric

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
    from numpy import mean

    reference_coords = reference.obsm["X_emb"]
    reference_labels = reference.obs["Label"].tolist()

    print("Selecting unseen query cells...")
    query = query[query.obs["Unseen"], :].copy()

    query_coords = query.obsm["X_emb"]
    query_labels = query.obs["Label"].tolist()

    print("Calculating inverse covariance matrices...")
    inverse_covariances = get_inverse_covariances(query_coords, query_labels)
    print("Calculating query centroid positions...")
    query_centroids = get_centroids(query_coords, query_labels)
    print("Calculating reference centroid positions...")
    reference_centroids = get_centroids(reference_coords, reference_labels)

    print("Calculating maximum query cell distances...")
    max_distances = {}
    for label in set(query_labels):
        label_distances = []
        label_query = query[query.obs["Label"] == label].copy()

        for idx in range(label_query.n_obs):
            label_coord = label_query.obsm["X_emb"][idx, :]

            distance = mahalanobis(
                label_coord, query_centroids[label], inverse_covariances[label]
            )
            label_distances.append(distance)

        max_distances[label] = max(label_distances)

    print("Calculating label Mahalonobis distances...")
    distances = []
    for label in set(query_labels):
        print(f"Calculating distances for label '{label}'...")

        # Skip labels with less than 2 * dim samples
        if query_labels.count(label) < 2 * query_coords.shape[1]:
            print(f"Label '{label}' has less than 2 * dim samples, skipping...")
            continue

        label_distances = []
        for ref_label in set(reference_labels):
            distance = mahalanobis(
                reference_centroids[ref_label],
                query_centroids[label],
                inverse_covariances[label],
            )
            label_distances.append(distance)

        scaled_distance = min(label_distances) / max_distances[label]
        distances.append(scaled_distance)

    print("Calculating final score...")
    # Set distances greated then the maximum to 1
    distances = [1.0 if d > 1.0 else d for d in distances]

    score = mean(distances)

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
    score = calculate_unseen_label_distance(reference, query)
    output = format_metric_results(
        dataset, method, integration, "Unseen", "unseenLabelDist", score
    )
    print(output)
    print(f"Writing output to '{out_file}'...")
    output.to_csv(out_file, sep="\t", index=False)
    print("Done!")


if __name__ == "__main__":
    main()
