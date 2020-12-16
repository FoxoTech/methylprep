`methylprep` is a python package for processing Illumina methylation array data.
View on [ReadTheDocs.](https://life-epigenetics-methylprep.readthedocs-hosted.com/en/latest/)

[![Readthedocs](https://readthedocs.com/projects/life-epigenetics-methylprep/badge/?version=latest)](https://life-epigenetics-methylprep.readthedocs-hosted.com/en/latest/) [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT) [![CircleCI](https://circleci.com/gh/FOXOBioScience/methylprep.svg?style=shield)](https://circleci.com/gh/FOXOBioScience/methylprep) [![Build status](https://ci.appveyor.com/api/projects/status/jqhqss0ks58kt4mh?svg=true)](https://ci.appveyor.com/project/life_epigenetics/methpype-ck8v2)
[![Codacy Badge](https://api.codacy.com/project/badge/Grade/5112cd82685548ffb8c64961e286180b)](https://www.codacy.com/app/marcmaxmeister/methylprep?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=FOXOBioScience/methylprep&amp;utm_campaign=Badge_Grade)
[![Coverage Status](https://coveralls.io/repos/github/FOXOBioScience/methylprep/badge.svg?t=mwigt8)](https://coveralls.io/github/FOXOBioScience/methylprep) ![PyPI-Downloads](https://img.shields.io/pypi/dm/methylprep.svg?label=pypi%20downloads&logo=PyPI&logoColor=white)

## Methylprep is part of the methyl-suite

![](https://github.com/FOXOBioScience/methylprep/blob/dev/docs/methyl-suite.png?raw=true)

`methylprep` is part of a methyl-suite of python packages that provide functions to process and analyze DNA methylation data from Illumina arrays (27, 450k, and EPIC/850k supported). The `methylprep` package contains functions for processing raw data files from arrays, or downloading (and processing) public data sets from GEO (the NIH Gene Expression Omnibus is a database repository), or from ArrayExpress. It contains both a command line interface (CLI) for processing data from local files, and a set of functions for building a custom pipeline in a jupyter notebook or python scripting environment. The aim is to offer a standard process, with flexibility for those who want it.

## Related packages

You should install all three components, as they work together.

- `methylcheck` includes
   - quality control (QC) functions for filtering out unreliable probes, based on the published literature and outlier detection.
   - sample outlier detection
   - array level QC plots, based on Genome Studio functions
   - data visualization functions based on seaborn and matplotlib graphic libraries.
   - predict sex of human samples from probes
   - interactive method for assigning samples to groups, based on array data, in a Jupyter notebook
- `methylize` provides analysis functions
   - differentially methylated probe statistics (between treatment and control samples)
   - volcano plots (which probes are the most different)
   - manhattan plot (where in genome are the differences)

## Installation

methylprep maintains configuration files for your Python package manager of choice: [pipenv](https://pipenv.readthedocs.io/en/latest/) or [pip](https://pip.pypa.io/en/stable/). Conda install is coming soon.

```python
pip install methylprep
```

---

## Command line data processing

The most common use case is processing `.idat` files on a computer within a command line interface. This can also be done in a Jupyter notebook, but large data sets take hours to run and Jupyter will take longer to run these than command line.

![processing pipeline](https://github.com/FOXOBioScience/methylprep/blob/master/docs/methylprep-processing-pipeline.png?raw=true)

### `process`

```shell
python -m methylprep -v process -d <filepath> --all
```
The `--all` option applies the most common settings. Here are some specific options:

Argument | Type | Default | Description
--- | --- | --- | ---
`data_dir` | `str`, `Path` | **REQUIRED** | Base directory of the sample sheet and associated IDAT files
`array_type` | `str` | `None` | Code of the array type being processed. Possible values are `custom`, `27k`, `450k`, `epic`, and `epic+`. If not provided, the pacakage will attempt to determine the array type based on the number of probes in the raw data. If the batch contains samples from different array types, this may not work. Our data `download` function attempts to split different arrays into separate batches for processing to accommodate this.
`manifest_filepath` | `str`, `Path` | `None` | File path for the array's manifest file. If not provided, this file will be downloaded from a Life Epigenetics archive.
`no_sample_sheet` | `bool` | `None` | pass in "--no_sample_sheet" from command line to trigger sample sheet auto-generation. Sample names will be based on idat filenames. Useful for public GEO data sets that lack sample sheets.
`sample_sheet_filepath` | `str`, `Path` | `None` | File path of the project's sample sheet. If not provided, the package will try to find one based on the supplied data directory path.
`sample_name` | `str` to list | `None` | List of sample names to process, in the CLI format of `-n sample1 sample2 sample3 etc`. If provided, only those samples specified will be processed. Otherwise all samples found in the sample sheet will be processed.
`export` | `bool` | `False` | Add flag to export the processed data to CSV.
`betas` | `bool` | `False` | Add flag to output a pickled dataframe of beta values of sample probe values.
`m_value` | `bool` | `False` | Add flag to output a pickled dataframe of m_values of samples probe values.
`batch_size` | `int` | `None` | Optional: splits the batch into smaller sized sets for processing. Useful when processing hundreds of samples that can't fit into memory. Produces multiple output files. This is also used by the package to process batches that come from different array types.

`data_dir` is the one required field. If you do not provide the file path for the project's sample_sheet, it will find one based on the supplied data directory path. It will also auto detect the array type and download the corresponding manifest file for you.


### `run_pipeline` (within a python interpreter, such as IDLE or Jupyter)

Run the complete methylation processing pipeline for the given project directory, optionally exporting the results to file.

Returns: A collection of DataContainer objects for each processed sample

```python
from methylprep import run_pipeline

data_containers = run_pipeline(data_dir, array_type=None, export=False, manifest_filepath=None, sample_sheet_filepath=None, sample_names=None)
```
Note: All the same input parameters from command line apply to `run_pipeline`, except `--all`. Type `dir(methylprep.run_pipeline)` in an interactive python session to see details.

Note: By default, if `run_pipeline` is called as a function in a script, a list of SampleDataContainer objects is returned. However, if you specify `betas=True` or `m_value=True`, a dataframe of beta values or m-values is returned instead. All `methylcheck` functions are designed to work on a dataframe or a folder to the processed data generated by `run_pipeline`.

### Getting help from command line

methylprep provides a command line interface (CLI) so the package can be used directly in bash/batchfile or windows/cmd scripts as part of building your custom processing pipeline.

All invocations of the methylprep CLI will provide contextual help, supplying the possible arguments and/or options available based on the invoked command. If you specify verbose logging the package will emit log output of DEBUG levels and above.

```Shell
python -m methylprep

usage: methylprep [-h] [-v] {process,sample_sheet} ...

Utility to process methylation data from Illumina IDAT files

positional arguments:
  {process,sample_sheet}
    process             process help
    sample_sheet        sample sheet help

optional arguments:
  -h, --help            show this help message and exit
  -v, --verbose         Enable verbose logging
```

---

### Other commands

The methylprep cli provides these top-level commands, which make it easier to use GEO datasets:

- `process` is the main function for processing methylation data from `idat` files. Covered already.
- `download` download and process public data sets in NIH GEO or ArrayExpress collections. Provide the public Accession ID and it will handle the rest.
- `beta_bake` combines `download`, `meta_data`, and file format conversion functions to produce a package
    that can be processed (with `process`) or loaded with `methylcheck.load` for analysis.
- `sample_sheet` will find/read/validate/create a sample sheet for a data set, or display its contents
    This is part of `process` and be applied using the `--no_sample_sheet` flag.
- `alert` scan GEO database and construct a CSV / dataframe of sample meta data and phenotypes for all studies matching a keyword
- `composite` download a bunch of datasets from a list of GEO ids, process them all, and combine into a large dataset
- `meta_data` will download just the meta data for a GEO dataset and convert it to a samplesheet CSV

### `download`

There are thousands of publically accessible DNA methylation data sets available via the GEO (US NCBI NIH) https://www.ncbi.nlm.nih.gov/geo/ and ArrayExpress (UK) https://www.ebi.ac.uk/arrayexpress/ websites. This function makes it easy to import them and build a reference library of methylation data.

The CLI now includes a `download` option. Supply the GEO ID or ArrayExpress ID and it will locate the files, download the idats, process them, and build a dataframe of the associated meta data. This dataframe format should be compatible with methylcheck and methylize.

Argument | Type | Default | Description
--- | --- | --- | ---  
-h, --help ||        show this help message and exit
  -d , --data_dir | `str` | [required path] | path to where the data series will be saved. Folder must exist already.
  -i ID, --id ID | `str` | [required ID] |The dataset's reference ID (Starts with `GSE` for GEO or `E-MTAB-` for ArrayExpress)
  -l LIST, --list LIST | `multiple strings` | optional | List of series IDs (can be either GEO or ArrayExpress), for partial downloading
  -o, --dict_only | `True` | pass flag only | If passed, will only create dictionaries and not process any samples
  -b BATCH_SIZE, --batch_size BATCH_SIZE | `int` | optional | Number of samples to process at a time, 100 by default.

When processing large batches of raw `.idat` files, specify `--batch_size` to break the processing up into smaller batches so the computer's memory won't overload. This is off by default when using `process` but is ON when using `download` and set to batch_size of 100. Set to 0 to force processing everything as one batch. The output files will be split into multiple files afterwards, and you can recomine them using `methylcheck.load`.

### `beta_bake`

Covered under [Public GEO Datasets](https://life-epigenetics-methylprep.readthedocs-hosted.com/en/latest/docs/example_download.html).

### `sample_sheet`

Find and parse the sample sheet in a given directory and emit the details of each sample. This is not required for actually processing data.

```Shell
>>> python -m methylprep sample_sheet
```

usage: methylprep sample_sheet -d DATA_DIR

optional arguments:

  Argument | Type | Description
  --- | --- | ---   
  -h, --help || show this help message and exit
  -d, --data_dir | string | Base directory of the sample sheet and associated IDAT files
  -c, --create | bool | If specified, this creates a sample sheet from idats instead of parsing an existing sample sheet. The output file will be called "samplesheet.csv".
  -o OUTPUT_FILE, --output_file OUTPUT_FILE | string | If creating a sample sheet, you can provide an optional output filename (CSV).                        


#### example of creating a sample sheet

```shell
~/methylprep$ python -m methylprep -v sample_sheet -d ~/GSE133062/GSE133062 --create
INFO:methylprep.files.sample_sheets:[!] Created sample sheet: ~/GSE133062/GSE133062/samplesheet.csv with 70 GSM_IDs
INFO:methylprep.files.sample_sheets:Searching for sample_sheet in ~/GSE133062/GSE133062
INFO:methylprep.files.sample_sheets:Found sample sheet file: ~/GSE133062/GSE133062/samplesheet.csv
INFO:methylprep.files.sample_sheets:Parsing sample_sheet
200861170112_R01C01
200882160083_R03C01
200861170067_R02C01
200498360027_R04C01
200498360027_R08C01
200861170067_R01C01
200861170072_R05C01
200498360027_R06C01
200861170072_R01C01
200861170067_R03C01
200882160070_R02C01
...
```

### `composite`

A tool to build a data set from a list of public datasets.

optional arguments:

Argument | Type | Description
--- | --- | ---  
  -h, --help || show this help message and exit
  -l LIST, --list LIST | filepath | A text file containg several GEO/ArrayExpress series ids. One ID per line in file. Note: The GEO Accession Viewer lets you export search results in this format.
  -d DATA_DIR, --data_dir DATA_DIR | filepath | Folder where to save data (and read the ID list file).
  -c, --control | bool | If flagged, this will only save samples that have the word "control" in their meta data.
  -k KEYWORD --keyword KEYWORD | string | Only retain samples that include this keyword (e.g. blood) somewhere in their meta data.
  -e, --export | bool | If passed, saves raw processing file data for each sample. (unlike meth-process, this is off by default)
  -b, --betas  | bool | If passed, output returns a dataframe of beta values for samples x probes. Local file beta_values.npy is also created.
  -m, --m_value  | bool | If passed, output returns a dataframe of M-values for samples x probes. Local file m_values.npy is also created.

### `alert`

Function to check for new datasets on GEO and update a csv each time it is run. Usable as a weekly cron command line function. Saves data to a local csv to compare with old datasets in <pattern>_meta.csv. Saves the dates of each dataset from GEO; calculates any new ones as new rows. updates csv.

optional arguments:

   Argument | Type | Description
   --- | --- | ---  
   keyword | string | Specify a word or phrase to narrow the search, such as "spleen blood".

### `meta_data`

Provides a more feature-rich meta data parser for public MINiML (formatted) GEO datasets. Run this
after downloading the dataset using `download` command. This reads all the
meta data from MINiML into a samplesheet.csv and meta data dataframe.

#### Sample exclusion filtering

You can use `meta_data` to identify 'control' or samples containing a specific keyword (e.g. blood,
tumor, etc) and remove any samples from sheet that lack these criteria, and
delete the associated idats that don't have these keywords. After, run
`process` on the rest, saving time. You can effectively ignore the parts of
datasets that you don't need based on the associated meta data.

optional arguments:

Argument | Type | Description
--- | --- | ---
  -h, --help  ||   show this help message and exit
  -i ID, --id ID |str|   Unique ID of the series (the GEO GSExxxx ID)
  -d DATA_DIR, --data_dir DATA_DIR |str or path|
                        Directory to search for MINiML file.
  -c, --control  |str| [experimental]: If flagged, this will look at the sample sheet and only save samples that appear to be"controls".
  -k KEYWORD, --keyword KEYWORD |str| [experimental]: Retain samples that include this keyword (e.g. blood, case insensitive) somewhere in samplesheet values.
  -s, --sync_idats | bool |      [experimental]: If flagged, this will scan the `data_dir` and remove all idat files that are not in the filtered samplesheet, so they won't be processed.
  -o, --dont_download  | bool | By default, this will first look at the local filepath (--data-dir) for `GSE..._family.xml` files. IF this is specified, it wont later look online to download the file. Sometimes a series has multiple files and it is easier to download, extract, and point this parser to each file instead.
