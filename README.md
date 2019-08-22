`methylprep` is a python package for processing Illumina methylation array data.
View on [ReadTheDocs.](https://life-epigenetics-methylprep.readthedocs-hosted.com/en/latest/)

[![Readthedocs](https://readthedocs.com/projects/life-epigenetics-methylprep/badge/?version=latest)](https://life-epigenetics-methylprep.readthedocs-hosted.com/en/latest/) [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT) [![CircleCI](https://circleci.com/gh/LifeEGX/methylprep.svg?style=shield)](https://circleci.com/gh/LifeEGX/methylprep) [![Build status](https://ci.appveyor.com/api/projects/status/jqhqss0ks58kt4mh?svg=true)](https://ci.appveyor.com/project/life_epigenetics/methpype-ck8v2) [![Codacy Badge](https://api.codacy.com/project/badge/Grade/5112cd82685548ffb8c64961e286180b)](https://www.codacy.com/app/marcmaxmeister/methylprep?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=LifeEGX/methylprep&amp;utm_campaign=Badge_Grade) [![Coverage Status](https://coveralls.io/repos/github/LifeEGX/methylprep/badge.svg?t=mwigt8)](https://coveralls.io/github/LifeEGX/methylprep)

## Methpype Package

The MethPype package contains both high-level APIs for processing data from local files and low-level functionality allowing you to customize the flow of data and how it is processed.

## Installation

MethPype maintains configuration files for your Python package manager of choice: [conda](https://conda.io), [pipenv](https://pipenv.readthedocs.io/en/latest/), and [pip](https://pip.pypa.io/en/stable/).

```python
pip install methylprep
```

---

## High-Level Processing

The primary Methpype API provides methods for the most common data processing and file retrieval functionality.

### `run_pipeline`

Run the complete methylation processing pipeline for the given project directory, optionally exporting the results to file.

Returns: A collection of DataContainer objects for each processed sample

```python
from methylprep import run_pipeline

data_containers = run_pipeline(data_dir, array_type=None, export=False, manifest_filepath=None, sample_sheet_filepath=None, sample_names=None)
```

Argument | Type | Default | Description
--- | --- | --- | ---
`data_dir` | `str`, `Path` | - | Base directory of the sample sheet and associated IDAT files
`array_type` | `str` | `None` | Code of the array type being processed. Possible values are `custom`, `450k`, `epic`, and `epic+`. If not provided, the pacakage will attempt to determine the array type based on the number of probes in the raw data.
`export` | `bool` | `False` | Whether to export the processed data to CSV
`manifest_filepath` | `str`, `Path` | `None` | File path for the array's manifest file. If not provided, this file will be downloaded from a Life Epigenetics archive.
`sample_sheet_filepath` | `str`, `Path` | `None` | File path of the project's sample sheet. If not provided, the package will try to find one based on the supplied data directory path.
`sample_names` | `str` collection | `None` | List of sample names to process. If provided, only those samples specified will be processed. Otherwise all samples found in the sample sheet will be processed.

### Methpype Command Line Interface (CLI)

Methpype provides a command line interface (CLI) so the package can be used directly in bash/batchfile scripts as part of building your custom processing pipeline.

All invocations of the MethPype CLI will provide contextual help, supplying the possible arguments and/or options available based on the invoked command. If you specify verbose logging the package will emit log output of DEBUG levels and above.

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

The MethPype cli provides two top-level commands:

- `process` to process methylation data
- `sample_sheet` to find/read a sample sheet and output its contents

#### `process`

Process the methylation data for a group of samples listed in a single sample sheet.

If you do not provide the file path for the project's sample_sheet the module will try to find one based on the supplied data directory path.
You must supply either the name of the array being processed or the file path for the array's manifest file. If you only specify the array type, the array's manifest file will be downloaded from a Life Epigenetics archive.

```Shell
>>> python -m methylprep process

usage: methylprep idat [-h] -d DATA_DIR [-a {custom,450k,epic,epic+}]
                     [-m MANIFEST] [-s SAMPLE_SHEET]
                     [--sample_name [SAMPLE_NAME [SAMPLE_NAME ...]]]
                     [--export]

Process Illumina IDAT files

optional arguments:
  -h, --help            show this help message and exit
  -d, --data_dir        Base directory of the sample sheet and associated IDAT
                        files
  -a, --array_type      Type of array being processed
                        Choices: {custom,450k,epic,epic+}
  -m, --manifest        File path of the array manifest file
  -s, --sample_sheet    File path of the sample sheet
  --sample_name         Sample(s) to process
  --export              Export data to csv
```

#### `sample_sheet`

Find and parse the sample sheet in a given directory and emit the details of each sample. This is not required for actually processing data.

```Shell
>>> python -m methylprep sample_sheet

usage: methylprep sample_sheet [-h] -d DATA_DIR

Process Illumina sample sheet file

optional arguments:
  -h, --help            show this help message and exit
  -d, --data_dir        Base directory of the sample sheet and associated IDAT
                        files
```


---

## Low-Level Processing

These are some functions that you can use within methylprep. `run_pipeline` calls them for you as needed.

#### `get_sample_sheet`

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

#### `get_manifest`

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

#### `get_raw_datasets`

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


