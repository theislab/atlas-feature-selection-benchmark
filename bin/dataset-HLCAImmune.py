#!/usr/bin/env python

"""
Download the HLCA (Immune) dataset (DOI: 10.1101/2022.03.10.483747)

Usage:
    dataset-HLCAImmune.py --out-file=<path>

Options:
    -h --help            Show this screen.
    --out-file=<path>    Path to output file.
"""


def get_HLCAImmune(temp_dir):
    """
    Get the HLCA (Immune) dataset

    Parameters
    ----------
    temp_dir
        Temporary directory object

    Returns
    -------
    AnnData containing the HLCA dataset
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

    print("Subsetting to immune cells...")
    is_immune = adata.obs["ann_level_1"] == "Immune"
    is_non_immune_dataset = adata.obs["dataset"].isin(
        [
            "Barbry_Leroy_2020",
            "Jain_Misharin_2021_10Xv2",
            "Jain_Misharin_2021_10Xv1",
            "Seibold_2020_10Xv2",
            "Seibold_2020_10Xv3",
        ]
    )
    adata = adata[is_immune & ~is_non_immune_dataset]

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
