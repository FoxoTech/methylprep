# Command Line Interface (CLI)

MethPype provides a command line interface (CLI) so the package can be used directly without the need of another custom library.

---

## Help and Logging

All invocations of the MethPype CLI will provide contextual help, supplying the possible arguments and/or options available based on the invoked command. If you specify verbose logging the package will emit log output of DEBUG levels and above.

```Shell
>>> python -m methpype

usage: methpype [-h] [-v] {process,sample_sheet} ...

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

## Commands

The MethPype cli provides two top-level commands:

- `process` to process methylation data
- `sample_sheet` to find/read a sample sheet and output its contents

### `process`

Process the methylation data for a group of samples listed in a single sample sheet.

If you do not provide the file path for the project's sample_sheet the module will try to find one based on the supplied data directory path.
You must supply either the name of the array being processed or the file path for the array's manifest file. If you only specify the array type, the array's manifest file will be downloaded from a Life Epigenetics archive.

```Shell
>>> python -m methpype process

usage: methpype idat [-h] -d DATA_DIR [-a {custom,450k,epic,epic+}]
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

### `sample_sheet`

Find and parse the sample sheet in a given directory and emit the details of each sample.

```Shell
>>> python -m methpype sample_sheet

usage: methpype sample_sheet [-h] -d DATA_DIR

Process Illumina sample sheet file

optional arguments:
  -h, --help            show this help message and exit
  -d, --data_dir        Base directory of the sample sheet and associated IDAT
                        files
```
