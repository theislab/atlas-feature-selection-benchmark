#!/usr/bin/env python

"""
Download the HLCA dataset (DOI: 10.1101/2022.03.10.483747)

Usage:
    dataset-HLCA.py --out-file=<path>

Options:
    -h --help            Show this screen.
    --out-file=<path>    Path to output file.
"""


def get_HLCA():
    """
    Get the HLCA dataset

    Returns
    -------
    AnnData containing the HLCA dataset
    """

    import scanpy as sc
    import numpy as np
    import os
    from os.path import join
    from tempfile import TemporaryDirectory

    originalFolder = os.getcwd()
    temp_dir = TemporaryDirectory()
    url = 'curl -o local.h5ad "https://corpora-data-prod.s3.amazonaws.com/13825e35-ea32-4104-a0b7-1986673cd3fc/local.h5ad?X-Amz-Algorithm=AWS4-HMAC-SHA256&X-Amz-Credential=ASIATLYQ5N5XY4TTPYHT%2F20220922%2Fus-west-2%2Fs3%2Faws4_request&X-Amz-Date=20220922T075845Z&X-Amz-Expires=604800&X-Amz-SignedHeaders=host&X-Amz-Security-Token=IQoJb3JpZ2luX2VjEN3%2F%2F%2F%2F%2F%2F%2F%2F%2F%2FwEaCXVzLXdlc3QtMiJGMEQCIGqdCwi5z3MfPgw%2FPGFzc26pTh4zaVSZ1%2FjfikjPiRnnAiACDfwMzCp0PNQMqCvC5ouWEM4VirhZzoRd2zwbrZ9TzSr0AwiG%2F%2F%2F%2F%2F%2F%2F%2F%2F%2F8BEAEaDDIzMTQyNjg0NjU3NSIM53jPy5QnbUzV0VJ%2FKsgDWTwdKbB8sjZh4tbclURG%2BAOM7ABYDklCwI%2Fg1iS7MYuKqxeR9CIk6cKtU0P0%2FQ57MuqRO2CeqME4c%2FpUgJIT5Aw1Xy43SqvU09bV2dll9BOs3R6YFCWKPkwyxTCfnJv0iN%2FvUFnq2FeCIl9%2FJi7LKSOE8y7oyGIR29zd8ij9Tumw5enbkZmeooA2w3bxJEbqduqacugc32m6YCVknMdH0prnx%2BXYBP%2FOLyqSqAPt0AK7Tw74E1rUSmEDr8FCEQ%2FS657zXKWOtBXBKbF1nTIcp%2B5egKxP6jBl%2BYPchCBMXQ%2ByFsgvkDBNYQHajrmCWhG267iQyROZVBk5gItkjQhmbeNIdyGnqlhJAxF1dNiRGhtaiE%2BbilYoBJokIR%2BYkIF0J%2BCf5NF6j5PjCL8XVdZBHbcqLCTgECaKNK%2BHz%2F7kYowGeMJ3s7C4Zy46X0D4bFyjRMJZCoc8LPmqmWeYRwbp%2BulwFHV18md1nw8pq56aiqFaSWPg%2BPOQNeV6xz1Hx%2FU%2BFAUqO7xzwiaLU0N%2BvO%2BTEHSU7eacnKboF2sZ3dtBfXs%2FwGPBnRbC3rqkB0w9KPagpKlOATsyK8yJbb2IU5SltCoqISCUFL78MPPYr5kGOqYB33Ne3CMNPbVbLkz%2Bt5npWYhBb482b1gVL35CS4ryaQbt4aiQbhn5QS%2FX70eY9XIrpIVYxUgDQ0ebqYBGUCq9lVu0RQDyC6Uwqnuw8pDWG%2BDPkCR65hv%2FyZVBnZ%2Fei4IVp7E%2FtAPk%2F%2Bqpx3qNbtnxZ10RIudsMzds0u%2FzQAXBwotI5FaoP%2Bd6VML8cbb6PIoJl%2BS%2BgsW7JORsfS%2B4sfLAynUHzObanA%3D%3D&X-Amz-Signature=3f5c9f6174dda96fac6f7df7af1199919a5bcda3eadf357d3bf0a2bee688191e"'

    print("Reading dataset from " + url + "...")
    os.chdir(temp_dir.name)
    os.system(url)

    adata = sc.read_h5ad(join(temp_dir.name, "local.h5ad"))
    
    print("Cleaning up temporary directory...")
    temp_dir.cleanup()
    
    adata = adata[adata.obs['ann_level_1']=="Immune"]
    adataSmall = adata[[study in ["Banovich_Kropski_2020","Misharin_2021","Krasnow_2020","Meyer_2019"] for study in adata.obs["study"]]].copy()
    adataSmall.obs = adataSmall.obs[['ann_finest_level','sample']]
    
    os.chdir(originalFolder)

    return adataSmall


def main():
    """The main script function"""
    from docopt import docopt

    args = docopt(__doc__)

    out_file = args["--out-file"]

    output = get_HLCA()
    print("Read dataset:")
    print(output)
    print(f"Writing output to '{out_file}'...")
    output.write_h5ad(out_file)
    print("Done!")


if __name__ == "__main__":
    main()
