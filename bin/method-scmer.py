#!/usr/bin/env python

"""
Select features using scmer

Usage:
    method-scmer.py --out-file=<path> [options] <file>

Options:
    -h --help            Show this screen.
    -o --out-file=<path>    Path to output file.
    -l --lasso              Lasso strength [default: 0]
    -n --n-features=<int>   Number of features to select [default: 1000].
    -b --batch              Apply to each batch. Requires column 'Batch' in .obs.
"""

def select_features_scanpy(adata, lasso, n, batch):
    """
    Select features using scanpy.pp.highly_variable_genes

    Parameters
    ----------
    adata
        AnnData object
    lasso
        Strength of L1 regularization in elastic net
    n
    	Number of features to select
    batch
    	Flag whether to use 'Batch' as batch_key

    Returns
    ----------
    DataFrame containing the selected features
    """

    import scanpy as sc
    from scmer import UmapL1

    if adata.n_vars < n:
        import warnings

        warnings.warn(
            "Number of features to select is greater than the number present, setting n to adata.n_vars"
        )
        n = adata.n_vars
    
    if batch:
        batches=adata.obs["Batch"].values
    else:
        batches=None
    
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    model = UmapL1(lasso=lasso).fit(X=adata.X, batches=batches)
    adata = model.transform(adata)
    
    model = UmapL1().tune(X=adata.X, target_n_features=n)
    adata = model.transform(adata)
    adata.var.shape

    adata.var["Feature"] = adata.var.index
    selected_features = adata.var
    
    return selected_features


def main():
    """The main script function"""
    from docopt import docopt
    from scanpy import read_h5ad

    args = docopt(__doc__)

    file = args["<file>"]
    lasso = int(args["--lasso"])
    n_features = int(args["--n-features"])
    batch = args["--batch"]
    out_file = args["--out-file"]

    print(f"Reading data from '{file}'...")
    input = read_h5ad(file)
    print("Read data:")
    print(input)
    output = select_features_scanpy(input, lasso, n_features, batch)
    print(f"Writing output to '{out_file}'...")
    output.to_csv(out_file, sep="\t", index=False)
    print("Done!")


if __name__ == "__main__":
    main()
