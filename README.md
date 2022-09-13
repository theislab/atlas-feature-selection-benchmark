# Atlas feature selection benchmarking

This repository contains code for benchmarking the effect of feature selection on scRNA-seq atlas construction and use.

For more information please refer to the documentation on the [wiki](https://github.com/theislab/atlas-feature-selection-benchmark/wiki).

## Directory structure

* `bin/` - Scripts used in **Nextflow** workflows
* `conf/` - **Nextflow** configuration files
* `envs/` - **conda** environment YAML files
* `output/` - Output from Nexflow workflows (not included in **git**)
* `reports/` - **RMarkdown** files and functions for output reports
* `work/` - The **Nextflow** working directory (not included in **git**)
* `workflows/` - **Nextflow** workflow files
* `LICENSE` - The project license
* `main.nf` - Main **Nextflow** workflow file
* `nextflow.config` - Main **Nextflow** config file
* `README.md` - This README
* `style_bin.sh` - A script for styling the files in `bin/`
