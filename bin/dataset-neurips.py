def get_neurips():
    """
    Get the neurips2021 open problem cite-seq dataset (RNA)
    Returns
    -------
    AnnData containing the BMMC dataset
    """
    import os
    import gzip
    from scanpy import read
    from tempfile import TemporaryDirectory
    from os.path import join
    
    FNAME="GSE194122_openproblems_neurips2021_cite_BMMC_processed.h5ad.gz"
    temp_dir = TemporaryDirectory()
    url = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE194nnn/GSE194122/suppl/GSE194122_openproblems_neurips2021_cite_BMMC_processed.h5ad.gz"


    print("Reading dataset from " + url + "...")
    os.system("cd " + temp_dir.name)
    os.system("wget " + url)
    os.system("tar -xvf " + FNAME)
    adata = read("GSE194122_openproblems_neurips2021_cite_BMMC_processed.h5ad")
    
    print("Cleaning up temporary directory...")
    temp_dir.cleanup()

    adata = adata[:,adata.var["feature_types"] == "GEX"]
    return adata

def main():
    """The main script function"""
    from docopt import docopt

    args = docopt(__doc__)

    out_file = args["--out-file"]

    output = get_neurips()
    print("Read dataset:")
    print(output)
    print("Writing output to '{out_file}'...")
    output.write_h5ad(out_file)
    print("Done!")


if __name__ == "__main__":
    main()