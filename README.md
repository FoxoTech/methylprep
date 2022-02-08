`methylprep` is a python package for processing Illumina methylation array data.
View on [ReadTheDocs.](https://life-epigenetics-methylprep.readthedocs-hosted.com/en/latest/)

[![tests](https://github.com/FoxoTech/methylprep/workflows/tests/badge.svg)](https://github.com/FoxoTech/methylprep/actions/workflows/ci.yml) [![Readthedocs](https://readthedocs.com/projects/life-epigenetics-methylprep/badge/?version=latest)](https://life-epigenetics-methylprep.readthedocs-hosted.com/en/latest/) [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT) [![CircleCI](https://circleci.com/gh/FoxoTech/methylprep.svg?style=shield)](https://circleci.com/gh/FoxoTech/methylprep) [![Codacy Badge](https://app.codacy.com/project/badge/Grade/e7228cfdfd714411bda7d6f8d6656630)](https://www.codacy.com/gh/FoxoTech/methylprep/dashboard?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=FoxoTech/methylprep&amp;utm_campaign=Badge_Grade) [![Coverage Status](https://coveralls.io/repos/github/FoxoTech/methylprep/badge.svg?t=mwigt8)](https://coveralls.io/github/FoxoTech/methylprep) [![PyPI-Downloads](https://img.shields.io/pypi/dm/methylprep.svg?label=pypi%20downloads&logo=PyPI&logoColor=white)](https://pypi.org/project/methylprep/)

## Methylprep is part of the methylsuite

![](https://raw.githubusercontent.com/FoxoTech/methylprep/master/docs/methyl-suite.png)

`methylprep` is part of the [methylsuite](https://pypi.org/project/methylsuite/) of python packages that provide functions to process and analyze DNA methylation data from Illumina's Infinium arrays (27k, 450k, and EPIC, as well as mouse arrays). The `methylprep` package contains functions for processing raw data files from arrays and downloading/processing public data sets from GEO (the NIH Gene Expression Omnibus database repository), or from ArrayExpress. It contains both a command line interface (CLI) for processing data from local files, and a set of functions for building a custom pipeline in a jupyter notebook or python scripting environment. The aim is to offer a standard process, with flexibility for those who want it.

`methylprep` data processing has also been tested and benchmarked to match the outputs of two popular R packages: [sesame](https://bioconductor.org/packages/release/bioc/html/sesame.html) (v1.10.4) and [minfi](https://bioconductor.org/packages/release/bioc/html/minfi.html) (v1.38).

## Methylsuite package components

You should install all three components, as they work together. The parts include:

- `methylprep`: (this package) for processing `idat` files or downloading GEO datasets from NIH. Processing steps include
   - infer type-I channel switch
   - NOOB (normal-exponential convolution on out-of-band probe data)
   - poobah (p-value with out-of-band array hybridization, for filtering lose signal-to-noise probes)
   - qualityMask (to exclude historically less reliable probes)
   - nonlinear dye bias correction (AKA signal quantile normalization between red/green channels across a sample)
   - calculate beta-value, m-value, or copy-number matrix
   - large batch memory management, by splitting it up into smaller batches during processing

- `methylcheck`: for quality control (QC) and analysis, including
   - functions for filtering out unreliable probes, based on the published literature
      - Note that `methylprep process` will exclude a set of unreliable probes by default. You can disable that using the --no_quality_mask option from CLI.
   - sample outlier detection
   - array level QC plots, based on Genome Studio functions
   - a python clone of Illumina's Bead Array Controls Reporter software (QC)
   - data visualization functions based on `seaborn` and `matplotlib` graphic libraries.
   - predict sex of human samples from probes
   - interactive method for assigning samples to groups, based on array data, in a Jupyter notebook

- `methylize` provides more analysis and interpretation functions
   - differentially methylated probe statistics (between treatment and control samples)
   - volcano plots (which probes are the most different?)
   - manhattan plots (where in genome are the differences?)

## Installation

`methylprep` maintains configuration files for your Python package manager of choice: [pipenv](https://pipenv.readthedocs.io/en/latest/) or [pip](https://pip.pypa.io/en/stable/). Conda install is coming soon.

```shell
>>> pip install methylprep
```

or if you want to install all three packages at once:
```shell
>>> pip install methylsuite
```

## Tutorials and Guides
If you're new to DNA methylation analysis, we recommend reading through [this introduction](docs/introduction/introduction.md) in order get the background knowledge needed to best utilize `methylprep` effectively. Otherwise, you're ready to use `methylprep` for:
<br>

- processing [your own methylation data](docs/general_walkthrough.md#processing-your-own-data)
- downloading [unprocessed data](docs/general_walkthrough.md#downloading-from-geo) (like IDAT files) from GEO.
- downloading [preprocessed data](docs/special_cases.md#using-beta-bake-for-preprocessed-data) (like beta values) from GEO.
- building a composite dataset [using control samples](docs/special_cases.md#building-a-composite-dataset-using-meta-data) from GEO.
- building a composite dataset from GEO data [with any keyword you choose](docs/special_cases.md#building-a-composite-dataset-with-alert-and-composite) (e.g. combining all GEO datasets that have methylation data from patients with brain cancer).

<!-- Add link to methods paper when available -->
