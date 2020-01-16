# Command Line Interface (CLI)

`methylprep` provides a command line interface (CLI) so the python package can be used directly from windows/cmd or linux/bash or macos/terminal.
---

## Help and Logging

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

## Commands

The methylprep cli provides three top-level commands:

- `process` to process methylation data
- `download` provide a GEO accession number like `-i GSE48472` and it will download and process the raw data for you
- `sample_sheet` to find/read a sample sheet and output its contents

### `process`

Process the methylation data for a group of samples listed in a single sample sheet.

If you do not provide the file path for the project's sample_sheet the module will try to find one based on the supplied data directory path.
You must supply either the name of the array being processed or the file path for the array's manifest file. If you only specify the array type, the array's manifest file will be downloaded from a Life Epigenetics archive.

```Shell
>>> python -m methylprep process

usage: methylprep process [-h] -d DATA_DIR [-a {custom,27k,450k,epic,epic+}]
                          [-m MANIFEST] [-s SAMPLE_SHEET] [--no_sample_sheet]
                          [-n [SAMPLE_NAME [SAMPLE_NAME ...]]] [-b] [-v]
                          [--batch_size BATCH_SIZE] [-u] [-e] [-x]
                          [-i {float64,float32,float16}]

Process Illumina IDAT files, producing NOOB, beta-value, or m_value corrected
scores per probe per sample

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
                        sample names too. If you need to add more meta data
                        into the sample_sheet, look at the create sample_sheet
                        CLI option.
  -n [SAMPLE_NAME [SAMPLE_NAME ...]], --sample_name [SAMPLE_NAME [SAMPLE_NAME ...]]
                        Sample(s) to process. You can pass multiple sample
                        names with multiple -n params.
  -b, --betas           If passed, output returns a dataframe of beta values
                        for samples x probes. Local file beta_values.npy is
                        also created.
  -v, --m_value         If passed, output returns a dataframe of M-values for
                        samples x probes. Local file m_values.npy is also
                        created.
  --batch_size BATCH_SIZE
                        If specified, samples will be processed and saved in
                        batches no greater than the specified batch size
  -u, --uncorrected     If specified, processed csv will contain two
                        additional columns (meth and unmeth) that have not
                        been NOOB corrected.
  -e, --no_export       Default is to export data to csv in same folder where
                        IDAT file resides. Pass in --no_export to suppress
                        this.
  -x, --no_meta_export  Default is to convert the sample sheet into a pickled
                        DataFrame, recognized in methylcheck and methylize.
                        Pass in --no_meta_export to suppress this.
  -i {float64,float32,float16}, --bit {float64,float32,float16}
                        Change the processed beta or m_value data_type output
                        from float64 to float16 or float32, to save disk
                        space.
```

### `download`

Lookup the given GEO dataset (by accession number), download it, along with the meta data stored in MINiML file format, process the idats, and save everything to disk.

```Shell
>>> python -m methylprep download

usage: methylprep download [-h] -d DATA_DIR [-i ID] [-l LIST] [-o]
                           [-b BATCH_SIZE] [-c]

Download and process a public dataset, either from GEO or ArrayExpress

optional arguments:
  -h, --help            show this help message and exit
  -d DATA_DIR, --data_dir DATA_DIR
                        Directory to download series to
  -i ID, --id ID        Unique ID of the series (either GEO or ArrayExpress
                        ID)
  -l LIST, --list LIST  List of series IDs (can be either GEO or ArrayExpress)
  -o, --dict_only       If passed, will only create dictionaries and not
                        process any samples
  -b BATCH_SIZE, --batch_size BATCH_SIZE
                        Number of samples to process at a time, 100 by default
  -c, --no_clean        Leave processing and raw data files in folders. By
                        default, these files are removed during processing.
```


### `sample_sheet`

Find and parse the sample sheet in a given directory and emit the details of each sample.

```Shell
>>> python -m methylprep sample_sheet

usage: methylprep sample_sheet [-h] -d DATA_DIR [-c] [-o OUTPUT_FILE]
                               [-t SAMPLE_TYPE] [-s SAMPLE_SUB_TYPE]

Create an Illumina sample sheet file from idat filenames and user-defined meta
data, or parse an existing sample sheet.

optional arguments:
  -h, --help            show this help message and exit
  -d DATA_DIR, --data_dir DATA_DIR
                        Base directory of the sample sheet and associated IDAT
                        files.
  -c, --create          If specified, this creates a sample sheet from idats
                        instead of parsing an existing sample sheet. The
                        output file will be called "samplesheet.csv".
  -o OUTPUT_FILE, --output_file OUTPUT_FILE
                        If creating a sample sheet, you can provide an
                        optional output filename (CSV).
  -t SAMPLE_TYPE, --sample_type SAMPLE_TYPE
                        Create sample sheet: Adds a "Sample_Type" column and
                        labels all samples in this sheet with this type. If
                        you have a batch of samples that have multiple types,
                        you must create multiple samplesheets and pass in
                        sample names and types to use this, or create your
                        sample sheet manually.
  -s SAMPLE_SUB_TYPE, --sample_sub_type SAMPLE_SUB_TYPE
                        Create sample sheet: Adds a "Sample_Sub_Type" column
                        and labels all samples in this sheet with this type.
                        If you have a batch of samples that have multiple
                        types, you must create multiple samplesheets and pass
                        in sample names and types to use this, or create your
                        sample sheet manually.
```
