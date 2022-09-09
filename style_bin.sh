#!/usr/bin/env bash

echo "Styling Python files using Black..."
black bin

echo "Styling R files using {styler} and {biocthis}..."
R --slave -e 'styler::style_dir("bin", transformers = biocthis::bioc_style())'

echo "Done!"
