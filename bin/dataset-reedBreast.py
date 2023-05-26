#!/usr/bin/env python

"""
Download the Reed Human Breast Cell Atlas dataset (DOI: 10.1101/2023.04.21.537845)

Usage:
    dataset-reedBreast.py --out-file=<path> [options]

Options:
    -h --help            Show this screen.
    --out-file=<path>    Path to output file.
"""


def get_reedBreast(temp_dir):
    """
    Get the Reed Human Breast Cell Atlas dataset

    Parameters
    ----------
    temp_dir
        Temporary directory object

    Returns
    -------
    AnnData containing the Reed HBCA dataset
    """

    from cellxgene_census import download_source_h5ad
    from anndata import read_h5ad
    from os.path import join

    print("Downloading Reed breast dataset from cellxgene...")
    dataset_id = "0ba636a1-4754-4786-a8be-7ab3cf760fd6"
    print(f"Using dataset ID '{dataset_id}'")
    temp_h5ad = join(temp_dir.name, "temp.h5ad")
    download_source_h5ad(dataset_id, to_path = temp_h5ad)

    print("Reading downloaded H5AD...")
    adata = read_h5ad(temp_h5ad, backed=True)

    print("Setting counts matrix...")
    adata.X = adata.raw.X.to_memory()

    print("Subsetting to WT and BRCA1 cells and removing doublets...")
    is_wt_brca = adata.obs["brca_status"].isin(["WT", "assume_WT", "BRCA1"])
    is_doublet = adata.obs["level2"] == "Doublet"
    adata = adata[is_wt_brca & ~is_doublet]

    return adata


def main():
    """The main script function"""
    from docopt import docopt
    from tempfile import TemporaryDirectory

    args = docopt(__doc__)

    out_file = args["--out-file"]

    temp_dir = TemporaryDirectory()

    output = get_reedBreast(temp_dir)
    print("Read dataset:")
    print(output)
    print(f"Writing output to '{out_file}'...")
    output.write_h5ad(out_file)
    print("Cleaning up temporary directory...")
    temp_dir.cleanup()
    print("Done!")


if __name__ == "__main__":
    main()
