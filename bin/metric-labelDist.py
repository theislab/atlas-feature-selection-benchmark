#!/usr/bin/env python

"""
Calculate the label mapping distance metric

Usage:
    metric-labelDist.py --dataset=<str> --method=<str> --integration=<str> --reference=<path> --out-file=<path> [options] <file>

Options:
    -h --help            Show this screen.
    --dataset=<str>      Name of the dataset to calculate the metric for.
    --method=<str>       Name of the method to calculate the metric for.
    --integration=<str>  Name of the integration to calculate the metric for.
    --reference=<path>   Path to the reference H5AD file.
    --out-file=<path>    Path to output file.
"""


def calculate_label_distance(reference, query):
    """
    Calculate the label Mahalanobis mapping distance metric

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
    from scipy.stats import chi2
    from functions.distances import get_inverse_covariances, get_centroids
    from numpy import mean

    reference_coords = reference.obsm["X_emb"]
    reference_labels = reference.obs["Label"].tolist()

    query_coords = query.obsm["X_emb"]
    query_labels = query.obs["Label"].tolist()

    print("Calculating inverse covariance matrices...")
    inverse_covariances = get_inverse_covariances(query_coords, query_labels)
    print("Calculating query centroid positions...")
    query_centroids = get_centroids(query_coords, query_labels)
    print("Calculating reference centroid positions...")
    reference_centroids = get_centroids(reference_coords, reference_labels)

    print("Calculating label Mahalonobis distances...")
    distances = []
    for label in set(query_labels):
        # Skip labels with less than 2 * dim samples
        if query_labels.count(label) < 2 * query_coords.shape[1]:
            print(f"Label '{label}' has less than 2 * dim samples, skipping...")
            continue

        # Skip labels that are not in the reference
        if label not in reference_centroids.keys():
            print(f"Label '{label}' is not in the reference, skipping...")
            continue

        distance = mahalanobis(
            reference_centroids[label],
            query_centroids[label],
            inverse_covariances[label],
        )
        distances.append(distance)

    print("Calculating p-values...")
    df = query_coords.shape[1]
    p_vals = [1 - chi2.cdf(dist, df=df) for dist in distances]

    print("Calculating final score...")
    score = mean(p_vals)

    return score


def main():
    """The main script function"""
    from docopt import docopt
    from scanpy import read_h5ad
    from functions.functions import format_metric_results

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
    score = calculate_label_distance(reference, query)
    output = format_metric_results(
        dataset, method, integration, "Mapping", "labelDist", score
    )
    print(output)
    print(f"Writing output to '{out_file}'...")
    output.to_csv(out_file, sep="\t", index=False)
    print("Done!")


if __name__ == "__main__":
    main()
