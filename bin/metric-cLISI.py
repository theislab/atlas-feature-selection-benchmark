#!/usr/bin/env python

"""
Evaluate integration using Cell-type LISI (cLISI)

Usage:
    metric-cLISI.py --dataset=<str> --method=<str> --integration=<str> --out-file=<path> <file>

Options:
    -h --help            Show this screen.
    --dataset=<str>      Name of the dataset to calculate the metric for.
    --method=<str>       Name of the method to calculate the metric for.
    --integration=<str>  Name of the integration to calculate the metric for.
    --out-file=<path>    Path to output file.
"""


def calculate_cLISI(adata):
    """
    Calculate the Cell-type LISI (cLISI) score for an integrated dataset.

    Parameters
    ----------
    adata
        AnnData object containing the integrated dataset

    Returns
    -------
    cLISI score
    """

    from scib.metrics import clisi_graph
    from scanpy.preprocessing import neighbors
    from igraph import Graph

    # Reduce the k0 parameter for small datasets to avoid theislab/scib/issues/374
    # This should only happen for tiny test datasets
    if adata.n_obs < 500:
        import warnings

        warnings.warn("Less than 500 cells, setting k0=50")
        k0 = 50
    else:
        k0 = 90

    adata.obs["Label"] = adata.obs["Label"].cat.remove_unused_categories()

    # Remove cells with fewer than k0 neighbours to avoid theislab/scib/issues/374

    # Calculate the neighbourhood graph using same settings as the metric
    neighbors(adata, n_neighbors=15, use_rep="X_emb")

    # Create a graph object from the neighbourhood graph
    graph = Graph.Weighted_Adjacency(adata.obsp["connectivities"])

    # Get the connected components
    components = graph.connected_components()
    component_sizes = [len(component) for component in components]

    # Check which components have fewer than k0 cells
    connected = []
    n_unconnected = 0
    for idx in range(len(components)):
        if component_sizes[idx] >= k0:
            connected += components[idx]
        else:
            n_unconnected += component_sizes[idx]

    if n_unconnected > 0:
        from warnings import warn
        warn(f"Found {n_unconnected} cells with fewer than {k0} neighbours. These cells will be skipped.")

    # Delete the neighbourhood graph so it's not used by the metric
    del adata.uns["neighbors"]
    del adata.obsp["distances"]
    del adata.obsp["connectivities"]

    print("Calculating cLISI score...")
    score = clisi_graph(
        adata[connected, :],
        label_key="Label",
        type_="embed",
        use_rep="X_emb",
        k0=k0,
        subsample=None,
        scale=True,
        n_cores=2,
        verbose=True,
    )
    print(f"Final score: {score}")

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
    out_file = args["--out-file"]

    print(f"Reading data from '{file}'...")
    input = read_h5ad(file)
    print("Read data:")
    print(input)
    score = calculate_cLISI(input)
    output = format_metric_results(
        dataset, method, integration, "IntegrationBio", "cLISI", score
    )
    print(output)
    print(f"Writing output to '{out_file}'...")
    output.to_csv(out_file, sep="\t", index=False)
    print("Done!")


if __name__ == "__main__":
    main()
