#!/usr/bin/env python

"""
Download the Ren et al. COVID-19 dataset (DOI: 10.1016/j.cell.2021.01.053)

Usage:
    dataset-COVID19.py --out-file=<path>

Options:
    -h --help            Show this screen.
    --out-file=<path>    Path to output file.
"""


def get_COVID19():
    """
    Get the scIB Pancreas dataset

    Returns
    -------
    AnnData containing the COVID-19 dataset
    """

    import scanpy as sc
    import numpy as np
    import os
    from os.path import join
    from tempfile import TemporaryDirectory

    FGGID="1IwWcn4W-YKgNbz4DpNweM2cKxlx1hbM0"
    FNAME="COVID19_ALL.h5ad.tar.gz"
    temp_dir = TemporaryDirectory()
    url = "'https://drive.google.com/uc?export=download&id=" + FGGID +"&confirm=yes' -O " + FNAME

    print("Reading dataset from " + url + "...")
    os.system("cd " + temp_dir.name)
    os.system("wget " + url)
    os.system("tar -xvf " + FNAME)
    adata = read(join(temp_dir.name, "COVID19_ALL.h5ad"))

    print("Cleaning up temporary directory...")
    temp_dir.cleanup()
    
    adataPBMC = adata[adata.obs['Sample type']=="fresh PBMC"]
    adataPBMCSmall = adataPBMC[[city in ["Shenzhen","Beijing","Huanggang","Harbin"] for city in adataPBMC.obs["City"]]].copy()
    adataPBMCSmall.obs = adataPBMCSmall.obs[['celltype','sampleID']]

    return adataPBMCSmall


def main():
    """The main script function"""
    from docopt import docopt

    args = docopt(__doc__)

    out_file = args["--out-file"]

    output = get_COVID19()
    print("Read dataset:")
    print(output)
    print("Writing output to '{out_file}'...")
    output.write_h5ad(out_file)
    print("Done!")


if __name__ == "__main__":
    main()
