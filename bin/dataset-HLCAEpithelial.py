#!/usr/bin/env python

"""
Download the HLCA (Epithelial) dataset (DOI: 10.1101/2022.03.10.483747)

Usage:
    dataset-HLCAEpithelial.py --out-file=<path>

Options:
    -h --help            Show this screen.
    --out-file=<path>    Path to output file.
"""


def get_HLCAEpithelial(temp_dir):
    """
    Get the HLCA (Epithelial) dataset

    Parameters
    ----------
    temp_dir
        Temporary directory object

    Returns
    -------
    AnnData containing the HLCA (Epithelial) dataset
    """

    from cellxgene_census import download_source_h5ad
    from anndata import read_h5ad
    from os.path import join

    print("Downloading HLCA dataset from cellxgene...")
    dataset_id = "066943a2-fdac-4b29-b348-40cede398e4e"
    print(f"Using dataset ID '{dataset_id}'")
    temp_h5ad = join(temp_dir.name, "temp.h5ad")
    download_source_h5ad(dataset_id, to_path=temp_h5ad)

    print("Reading downloaded H5AD...")
    adata = read_h5ad(temp_h5ad, backed=True)

    print("Setting counts matrix...")
    adata.X = adata.raw.X.to_memory()

    print("Setting var names...")
    adata.var_names = adata.var["feature_name"].astype(str)
    adata.var_names_make_unique()

    print("Subsetting to epithelial cells...")
    is_epithelial = adata.obs["ann_level_1"] == "Epithelial"
    is_non_epithelial_dataset = adata.obs["dataset"].isin(
        [
            "Lafyatis_Rojas_2019_10Xv1",
        ]
    )
    adata = adata[is_epithelial & ~is_non_epithelial_dataset]

    return adata


def main():
    """The main script function"""
    from docopt import docopt
    from tempfile import TemporaryDirectory

    args = docopt(__doc__)

    out_file = args["--out-file"]

    temp_dir = TemporaryDirectory()

    output = get_HLCAImmune(temp_dir)
    print("Read dataset:")
    print(output)
    print(f"Writing output to '{out_file}'...")
    output.write_h5ad(out_file)
    print("Cleaning up temporary directory...")
    temp_dir.cleanup()
    print("Done!")


if __name__ == "__main__":
    main()
