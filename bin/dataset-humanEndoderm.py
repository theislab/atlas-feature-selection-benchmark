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

    from scanpy import read_mtx
    from pandas import read_csv
    from os.path import join
    from urllib.request import urlopen
    from zipfile import ZipFile
    from tempfile import TemporaryDirectory
    import certifi
    import ssl

    # Download archive
    url = "https://prod-dcd-datasets-cache-zipfiles.s3.eu-west-1.amazonaws.com/x53tts3zfr-5.zip"
    temp_dir = TemporaryDirectory()
    zip_file = join(temp_dir.name, "human-endoderm-atlas.zip")

    print(f"Downloading dataset from '{url}'...")
    ssl_context = ssl.create_default_context(cafile=certifi.where())
    with urlopen(url, context=ssl_context) as url_file:
        with open(zip_file, "wb") as out_file:
            out_file.write(url_file.read())

    unpack_dir = join(temp_dir.name, "human-endoderm-atlas")
    print(f"Unpacking '{zip_file}' to '{unpack_dir}'...")
    with ZipFile(zip_file, "r") as zip:
        zip.extractall(unpack_dir)

    data_dir = join(unpack_dir, "2. Count matrix and meta data/Fetal atlas")
    counts_file = join(data_dir, "Table_fetal_atlas_count.mtx.gz")
    print(f"Reading counts from '{counts_file}'...")
    adata = read_mtx(counts_file)
    adata = adata.T

    obs_file = join(data_dir, "Table_fetal_atlas_meta_info.csv.gz")
    print(f"Reading obs from '{obs_file}'...")
    adata.obs = read_csv(obs_file)

    var_file = join(data_dir, "Table_fetal_atlas_gene_symbol.csv.gz")
    print(f"Reading var from '{var_file}'...")
    genes = read_csv(var_file, header=None)
    adata.var["name"] = genes.values
    adata.var.index = adata.var["name"].values
    adata.var_names_make_unique()

    print("Removing undefined cells...")
    is_undefined = adata.obs["Cell_type"].isin(["Undefined"])
    adata = adata[~is_undefined, :]

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
