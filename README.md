# Atlas feature selection benchmarking

[![PAPER](https://img.shields.io/badge/Nature%20Methods-10.1038%2Fs41592--025--02624--3-eb5b25)](https://doi.org/10.1038/s41592-025-02624-3)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.13995812.svg)](https://doi.org/10.5281/zenodo.13995812)
[![DATA](https://img.shields.io/badge/figshare-10.6084%2Fm9.figshare.c.7521966-a60845?logo=figshare)](https://doi.org/10.6084/m9.figshare.c.7521966)

This repository contains code for benchmarking the effect of feature selection on scRNA-seq atlas construction and use.

For more information please refer to the documentation on the [wiki](https://github.com/theislab/atlas-feature-selection-benchmark/wiki).

## Directory structure

* `analysis/` - Notebooks used to perform analysis of the results
  * `R/` - R functions used in the analysis notebooks
* `bin/` - Scripts used in **Nextflow** workflows
  * `functions/` - Functions used across multiple scripts
* `conf/` - **Nextflow** configuration files
* `data/` - Data files used in the analysis. Metric score files can downloaded from [figshare](https://figshare.com/articles/dataset/Metric_scores/26390911).
  * `selected-features/` - Selected features files for each dataset. Can be downloaded from [figshare](https://figshare.com/articles/dataset/Selected_features/26390944).
* `envs/` - **conda** environment YAML files
* `output/` - Output from Nextflow workflows (not included in **git**)
* `reports/` - **RMarkdown** files and functions for output reports generated by the pipeline
* `work/` - The **Nextflow** working directory (not included in **git**)
* `workflows/` - **Nextflow** workflow files
* `LICENSE` - The project license
* `main.nf` - Main **Nextflow** workflow file
* `nextflow.config` - Main **Nextflow** config file
* `README.md` - This README
* `style_bin.sh` - A script for styling the files in `bin/`

## Data files

Data files for this project are available from [figshare](https://doi.org/10.6084/m9.figshare.c.7521966).
The analysis reports in the `analysis/` directory assume:

* [Metric files](https://figshare.com/articles/dataset/Metric_scores/26390911) are downloaded to the `data/` directory
* [Selected features files](https://figshare.com/articles/dataset/Selected_features/26390944) are downloaded to the `data/selected-features/` directory

## Citation

To refer to this repository, please cite:

> Zappia, L., Ramírez-Suástegui, C., Kfuri-Rubens, R., Vornholz, L., Wang, W., Dietrich, O., Frishberg, A., D Luecken, M. & J Theis, F. [_Feature selection methods affect the performance of scRNA-seq data integration and querying._](https://doi.org/10.1038/s41592-025-02624-3) Nature methods 1–11 (2025). doi:10.1038/s41592-025-02624-3

```bibtex
@ARTICLE{Zappia2025,
  title     = "{Feature selection methods affect the performance of scRNA-seq
               data integration and querying}",
  author    = "Zappia, Luke and Richter, Sabrina and Ramírez-Suástegui, Ciro and
               Kfuri-Rubens, Raphael and Vornholz, Larsen and Wang, Weixu and
               Dietrich, Oliver and Frishberg, Amit and Luecken, Malte D and
               Theis, Fabian J",
  journal   = "Nature methods",
  publisher = "Nature Publishing Group",
  pages     = "1--11",
  month     =  mar,
  year      =  2025,
  doi       = "10.1038/s41592-025-02624-3",
  pmid      =  40082610,
  issn      = "1548-7091,1548-7105",
  language  = "en"
}
```
