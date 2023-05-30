#!/usr/bin/env python

"""
Download the human endoderm atlas (DOI: 10.1016/j.cell.2021.04.028)

Usage:
    dataset-humanEndoderm.py --out-file=<path> [options]

Options:
    -h --help               Show this screen.
    -o --out-file=<path>    Path to output file.
"""


def get_humanEndoderm():
    """
    Get the human endoderm atlas

    Returns
    -------
    AnnData containing the dataset
    """

    import scanpy as sc
    import pandas as pd
    import os
    from os.path import join
    from os.path import isfile
    from os.path import dirname
    from tempfile import TemporaryDirectory

    # Download archive
    url = "https://md-datasets-cache-zipfiles-prod.s3.eu-west-1.amazonaws.com/x53tts3zfr-1.zip"
    temp_dir = TemporaryDirectory()
    dest = join(temp_dir.name, "human-endoderm-atlas.zip")

    if isfile(dest):
        print(f"File '{dest}' already exists...")
    else:
        print(f"Reading dataset from '{url}'...")
        call = "curl" + " " + url + " -o " + dest
        os.system(call)

    print(f"Unpacking '{dest}'...")
    call = "unzip" + " " + dest + " -d " + dirname(dest)
    os.system(call)

    dest = dest.replace(".zip", "")
    datapath = join(dest, "2. Count matrix and meta data/Fetal atlas")

    print("Create adata from counts")
    adata = sc.read_mtx(join(datapath, "Table_fetal_atlas_count.mtx.gz"))
    adata = adata.T

    print("Adding .obs")
    coldata = pd.read_csv(join(datapath, "Table_fetal_atlas_meta_info.csv.gz"))
    adata.obs = coldata

    print("Adding .var")
    genes = pd.read_csv(join(datapath, "Table_fetal_atlas_gene_symbol.csv.gz"), header=None)
    adata.var["name"] = genes.values
    adata.var.index = adata.var["name"].values
    adata.var_names_make_unique()

    print("Cleaning up temporary directory...")
    temp_dir.cleanup()

    return adata


def main():
    """The main script function"""
    from docopt import docopt

    args = docopt(__doc__)

    out_file = args["--out-file"]

    output = get_humanEndoderm()
    print("Read dataset:")
    print(output)
    print(f"Writing output to '{out_file}'...")
    output.write_h5ad(out_file)
    print("Done!")


if __name__ == "__main__":
    main()
