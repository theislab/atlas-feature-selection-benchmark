#!/usr/bin/env python

"""
Select features using scmer

Usage:
    method-scmer.py --out-file=<path> [options] <file>

Options:
    -h --help            Show this screen.
    --out-file=<path>    Path to output file.
    --n-features=<int>   Number of features to select [default: 1000].
    --flavor=<str>       Flavor of feature selection method.
    One of: seurat, cell_ranger, seurat_v3 [default: seurat].    
    --batch=<bool>        Whether to apply to each batch.
    Requires Batch in .obs [default: False].
"""

def select_features_scanpy(adata, n, flavor, batch):
    """
    Select features using scanpy.pp.highly_variable_genes

    Parameters
    ----------
    adata
        AnnData object
    n
    	Number of features to select
    flavor
    	Method flavor ('seurat', 'cell_ranger', 'seurat_v3')
    batch
    	Boolean whether to use 'Batch' as batch_key

    Returns
    ----------
    DataFrame containing the selected features
    """

    if adata.n_vars < n:
        import warnings

        warnings.warn(
            "Number of features to select is greater than the number present, setting n to adata.n_vars"
        )
        n = adata.n_vars

    import scanpy as sc
    
    if batch:
        batch_key='Batch'
    else:
        batch_key=None
    
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=n, flavor=flavor, batch_key=batch_key)
    
    adata.var["Feature"] = adata.var.index
    selected_features = adata.var[adata.var["highly_variable"] == True]
    
    return selected_features


def main():
    """The main script function"""
    from docopt import docopt
    from scanpy import read_h5ad

    args = docopt(__doc__)

    file = args["<file>"]
    n_features = int(args["--n-features"])
    flavor = args["--flavor"]
    batch = args["--batch"]
    out_file = args["--out-file"]

    print(f"Reading data from '{file}'...")
    input = read_h5ad(file)
    print("Read data:")
    print(input)
    output = select_features_scanpy(input, n_features, flavor, batch)
    print(f"Writing output to '{out_file}'...")
    output.to_csv(out_file, sep="\t", index=False)
    print("Done!")


if __name__ == "__main__":
    main()
