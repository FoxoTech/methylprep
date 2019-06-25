# MethPype Package

The MethPype package contains both high-level APIs for processing data from local files and low-level functionality allowing you to customize the flow of data and how it is processed.

---

## High-Level Processing

The primary MethPype API provides methods for the most common data processing and file retrieval functionality.

### `run_pipeline`

Run the complete methylation processing pipeline for the given project directory, optionally exporting the results to file.

Returns: A collection of DataContainer objects for each processed sample

```python
from methpype import run_pipeline

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

### `get_sample_sheet`

Find and parse the sample sheet for the provided project directory path.

Returns: A SampleSheet object containing the parsed sample information from the project's sample sheet file

```python
from methpype import get_sample_sheet

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
from methpype import get_manifest

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
from methpype import get_raw_datasets

raw_datasets = get_raw_datasets(sample_sheet, sample_names=None)
```

Argument | Type | Default | Description
--- | --- | --- | ---
`sample_sheet` | `SampleSheet` | - | A SampleSheet instance from a valid project sample sheet file.
`sample_names` | `str` collection | `None` | List of sample names to process. If provided, only those samples specified will be processed. Otherwise all samples found in the sample sheet will be processed.

---

## Low-Level Processing
