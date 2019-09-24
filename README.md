`methylprep` is a python package for processing Illumina methylation array data.
View on [ReadTheDocs.](https://life-epigenetics-methylprep.readthedocs-hosted.com/en/latest/)

[![Readthedocs](https://readthedocs.com/projects/life-epigenetics-methylprep/badge/?version=latest)](https://life-epigenetics-methylprep.readthedocs-hosted.com/en/latest/) [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT) [![CircleCI](https://circleci.com/gh/LifeEGX/methylprep.svg?style=shield)](https://circleci.com/gh/LifeEGX/methylprep) [![Build status](https://ci.appveyor.com/api/projects/status/jqhqss0ks58kt4mh?svg=true)](https://ci.appveyor.com/project/life_epigenetics/methpype-ck8v2)
[![Codacy Badge](https://api.codacy.com/project/badge/Grade/5112cd82685548ffb8c64961e286180b)](https://www.codacy.com/app/marcmaxmeister/methylprep?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=LifeEGX/methylprep&amp;utm_campaign=Badge_Grade)
[![Coverage Status](https://coveralls.io/repos/github/LifeEGX/methylprep/badge.svg?t=mwigt8)](https://coveralls.io/github/LifeEGX/methylprep)

## methylprep Package

The methylprep package contains both high-level APIs for processing data from local files and low-level functionality allowing you to customize the flow of data and how it is processed.

## Installation

methylprep maintains configuration files for your Python package manager of choice: [conda](https://conda.io), [pipenv](https://pipenv.readthedocs.io/en/latest/), and [pip](https://pip.pypa.io/en/stable/).

```python
pip install methylprep
```

---

## High-Level Processing

The primary methylprep API provides methods for the most common data processing and file retrieval functionality.

### `run_pipeline`

Run the complete methylation processing pipeline for the given project directory, optionally exporting the results to file.

Returns: A collection of DataContainer objects for each processed sample

```python
from methylprep import run_pipeline

data_containers = run_pipeline(data_dir, array_type=None, export=False, manifest_filepath=None, sample_sheet_filepath=None, sample_names=None)
```

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

Note: By default, if `run_pipeline` is called as a function in a script, a list of SampleDataContainer objects is returned.

### methylprep Command Line Interface (CLI)

methylprep provides a command line interface (CLI) so the package can be used directly in bash/batchfile scripts as part of building your custom processing pipeline.

All invocations of the methylprep CLI will provide contextual help, supplying the possible arguments and/or options available based on the invoked command. If you specify verbose logging the package will emit log output of DEBUG levels and above.

```Shell
>>> python -m methylprep

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

### Commands

The methylprep cli provides two top-level commands:

- `process` to process methylation data
- `download` script to download and process public data sets in NIH GEO or ArrayExpress collections. Provide the public Accession ID and it will handle the rest.
- `sample_sheet` to find/read/validate a sample sheet and output its contents

### `process`

Process the methylation data for a group of samples listed in a single sample sheet.

If you do not provide the file path for the project's sample_sheet the module will try to find one based on the supplied data directory path.
You must supply either the name of the array being processed or the file path for the array's manifest file. If you only specify the array type, the array's manifest file will be downloaded from a Life Epigenetics archive.

```Shell
>>> python -m methylprep process

usage: methylprep idat [-h] -d DATA_DIR [-a {custom,27k,450k,epic,epic+}]
                       [-m MANIFEST] [-s SAMPLE_SHEET] [--no_sample_sheet]
                       [-n [SAMPLE_NAME [SAMPLE_NAME ...]]] [-e] [-b]
                       [--m_value] [--batch_size BATCH_SIZE]

Process Illumina IDAT files

optional arguments:
  -h, --help            show this help message and exit
  -d DATA_DIR, --data_dir DATA_DIR
                        Base directory of the sample sheet and associated IDAT
                        files. If IDAT files are in nested directories, this
                        will discover them.
  -a {custom,27k,450k,epic,epic+}, --array_type {custom,27k,450k,epic,epic+}
                        Type of array being processed. If omitted, this will
                        autodetect it.
  -m MANIFEST, --manifest MANIFEST
                        File path of the array manifest file. If omitted, this
                        will download the appropriate file from `s3`.
  -s SAMPLE_SHEET, --sample_sheet SAMPLE_SHEET
                        File path of the sample sheet. If omitted, this will
                        discover it. There must be only one CSV file in the
                        data_dir for discovery to work.
  --no_sample_sheet     If your dataset lacks a sample sheet csv file, specify
                        --no_sample_sheet to have it create one on the fly.
                        This will read .idat file names and ensure processing
                        works. If there is a matrix file, it will add in
                        sample names too.
  -n [SAMPLE_NAME [SAMPLE_NAME ...]], --sample_name [SAMPLE_NAME [SAMPLE_NAME ...]]
                        Sample(s) to process. You can pass multiple sample
                        names with multiple -n params.
  -e, --no_export       Default is to export data to csv in same folder where
                        IDAT file resides. Pass in --no_export to suppress
                        this.
  -b, --betas           If passed, output returns a dataframe of beta values
                        for samples x probes. Local file beta_values.npy is
                        also created.
  --m_value             If passed, output returns a dataframe of M-values for
                        samples x probes. Local file m_values.npy is also
                        created.
  --batch_size BATCH_SIZE
                        If specified, samples will be processed and saved in
                        batches no greater than the specified batch size
```

### `download`

There are thousands of publically accessible DNA methylation data sets available via the GEO (US NCBI NIH) https://www.ncbi.nlm.nih.gov/geo/ and ArrayExpress (UK) https://www.ebi.ac.uk/arrayexpress/ websites. This function makes it easy to import them and build a reference library of methylation data.

Argument | Type | Default | Description
--- | --- | --- | ---  
  -d , --data_dir | `str` | [required path] | path to where the data series will be saved. Folder must exist already.
  -i ID, --id ID | `str` | [required ID] |The dataset's reference ID (Starts with `GSM` for GEO or `E-MTAB-` for ArrayExpress)
  -l LIST, --list LIST | `multiple strings` | optional | List of series IDs (can be either GEO or ArrayExpress), for partial downloading
  -o, --dict_only | `True` | pass flag only | If passed, will only create dictionaries and not process any samples
  -b BATCH_SIZE, --batch_size BATCH_SIZE | `int` | optional | Number of samples to process at a time, 100 by default. Set to 0 for processing everything as one batch. Regardless of this number, the resulting file structure will be the same. But most machines cannot process more than 200 samples in memory at once, so this helps the user set the memory limits for their machine.

### `sample_sheet`

Find and parse the sample sheet in a given directory and emit the details of each sample. This is not required for actually processing data.

```Shell
>>> python -m methylprep sample_sheet

usage: methylprep sample_sheet [-h] -d DATA_DIR

Process Illumina sample sheet file

optional arguments:
  -h, --help            show this help message and exit
  -d, --data_dir        Base directory of the sample sheet and associated IDAT
                        files
  -c, --create          If specified, this creates a sample sheet from idats
                        instead of parsing an existing sample sheet. The
                        output file will be called "samplesheet.csv".
  -o OUTPUT_FILE, --output_file OUTPUT_FILE
                        If creating a sample sheet, you can provide an
                        optional output filename (CSV).                        
```

#### example of creating a sample sheet
```bash
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

#### `download`
The CLI now includes a `download` option. Supply the GEO ID or ArrayExpress ID and it will locate the files, download the idats, process them, and build a dataframe of the associated meta data. This dataframe format should be compatible with methylcheck and methylize. 

##### optional arguments:

Argument | Type | Description
--- | --- | ---
  -h, --help ||        show this help message and exit
  -d DATA_DIR, --data_dir DATA_DIR | path (required) | Directory to download series to
  -i ID, --id ID | string | Unique ID of the series (either GEO or ArrayExpress ID)
  -l LIST, --list LIST | multiple string arguments | List of series IDs (can be either GEO or ArrayExpress)
  -o, --dict_only | no args | If passed, will only create dictionaries and not process any samples
  -b BATCH_SIZE, --batch_size BATCH_SIZE | number | Number of samples to process at a time, 100 by default

- When processing large batches of raw `.idat` files, specify `--batch_size` to break the processing up into smaller batches so the computer's memory won't overload. This is off by default when using `process` but is ON when using `download` and set to batch_size of 100.


---

## Low-Level Processing

These are some functions that you can use within methylprep. `run_pipeline` calls them for you as needed.

### `get_sample_sheet`

Find and parse the sample sheet for the provided project directory path.

Returns: A SampleSheet object containing the parsed sample information from the project's sample sheet file

```python
from methylprep import get_sample_sheet

sample_sheet = get_sample_sheet(dir_path, filepath=None)
```

Argument | Type | Default | Description
--- | --- | --- | ---
`data_dir` | `str`, `Path` | - | Base directory of the sample sheet and associated IDAT files
`sample_sheet_filepath` | `str`, `Path` | `None` | File path of the project's sample sheet. If not provided, the package will try to find one based on the supplied data directory path.

### `get_manifest`

Find and parse the manifest file for the processed array type.

Returns: A Manifest object containing the parsed probe information for the processed array type

```python
from methylprep import get_manifest

manifest = get_manifest(raw_datasets, array_type=None, manifest_filepath=None)
```

Argument | Type | Default | Description
--- | --- | --- | ---
`raw_datasets` | `RawDataset` collection | - | Collection of RawDataset objects containing probe information from the raw IDAT files.
`array_type` | `str` | `None` | Code of the array type being processed. Possible values are `custom`, `450k`, `epic`, and `epic+`. If not provided, the pacakage will attempt to determine the array type based on the provided RawDataset objects.
`manifest_filepath` | `str`, `Path` | `None` | File path for the array's manifest file. If not provided, this file will be downloaded from a Life Epigenetics archive.

### `get_raw_datasets`

Find and parse the IDAT files for samples within a project's sample sheet.

Returns: A collection of RawDataset objects for each sample's IDAT file pair.

```python
from methylprep import get_raw_datasets

raw_datasets = get_raw_datasets(sample_sheet, sample_names=None)
```

Argument | Type | Default | Description
--- | --- | --- | ---
`sample_sheet` | `SampleSheet` | - | A SampleSheet instance from a valid project sample sheet file.
`sample_names` | `str` collection | `None` | List of sample names to process. If provided, only those samples specified will be processed. Otherwise all samples found in the sample sheet will be processed.


