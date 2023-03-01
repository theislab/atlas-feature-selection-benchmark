#!/usr/bin/env bash

echo "Styling Python files using Black..."
echo "Styling scripts..."
black bin
echo "Styling functions..."
black bin/functions

echo "Styling R files using {styler} and {biocthis}..."
echo "Styling scripts..."
R --slave -e 'styler::style_dir("bin", transformers = biocthis::bioc_style())'
echo "Styling functions..."
R --slave -e 'styler::style_dir("bin/functions", transformers = biocthis::bioc_style())'

echo "Done!"
