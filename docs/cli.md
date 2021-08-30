# using methylprep from the command line

The most common use case is processing `.idat` files on a computer within a command line interface. This can also be done in a Jupyter notebook, but large data sets take hours to run and Jupyter will take longer to run these than command line.

![processing pipeline](https://github.com/FoxoTech/methylprep/blob/master/docs/methylprep-processing-pipeline.png?raw=true)

## getting help from command line

`methylprep` provides a command line interface (CLI) so the package can be used directly in bash/batchfile or windows/cmd scripts as part of building your custom processing pipeline.

All invocations of the `methylprep` CLI will provide contextual help, supplying the possible arguments and/or options available based on the invoked command. If you specify verbose logging the package will emit log output of DEBUG levels and above.

For more information about `methylprep` as a package:

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

For more information about any of `methylprep`'s functions, simply type the name of the function:
```Shell
>>> python -m methylprep process
[Error]:
the following arguments are required: -d/--data_dir

usage: methylprep process [-h] -d DATA_DIR
                          [--array_type {custom,27k,450k,epic,epic+,mouse}]
                          [-m MANIFEST] [-s SAMPLE_SHEET] [--no_sample_sheet]
                          [-n [SAMPLE_NAME [SAMPLE_NAME ...]]] [-b] [-v]
                          [--batch_size BATCH_SIZE] [-u] [-e] [-x]
                          [-i {float64,float32,float16}] [-c] [--poobah]
                          [--export_poobah] [--minfi] [--no_quality_mask] [-a]

Process Illumina IDAT files, producing NOOB, beta-value, or m_value corrected
scores per probe per sample

optional arguments:
  -h, --help            show this help message and exit
  -d DATA_DIR, --data_dir DATA_DIR
                        Base directory of the sample sheet and associated IDAT
                        files. If IDAT files are in nested directories, this
                        will discover them.
  --array_type {custom,27k,450k,epic,epic+,mouse}
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
                        names like this: `python -m methylprep process -d .
                        --all --no_sample_sheet -n Sample_1 Sample_2 Sample_3`
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
  -c, --save_control    If specified, saves an additional "control_probes.pkl"
                        file that contains Control and SNP-I probe data in the
                        data_dir.
  --poobah              By default, any beta-values or m-values output will
                        contain all probes. If True, those probes that fail
                        the p-value signal:noise detection are replaced with
                        NaNs in dataframes in beta_values and m_value output.
  --export_poobah       If specified, exports a pickled dataframe of the
                        poobah p-values per sample.
  --minfi               If specified, processing uses legacy parameters based
                        on minfi. By default, v1.4.0 and higher mimics sesame
                        output.
  --no_quality_mask     If specified, processing to RETAIN all probes that
                        would otherwise be excluded using the quality_mask
                        sketchy-probe list from sesame. --minfi processing
                        does not use a quality_mask.
  -a, --all             If specified, saves everything: (beta_values.pkl,
                        m_value.pkl, control_probes.pkl, CSVs for each sample,
                        uncluding uncorrected raw values, and meta data, and
                        poobah_values.pkl). And removes failed probes using
                        sesame pOOBah method from these files. This overrides
                        individual CLI settings.

```

---

## `process`

```shell
python -m methylprep -v process -d <filepath> --all
```

`-d` (data file path) is the only required option. `-v` (short for `--verbose`) specifies verbose logging. And the `--all` option tells `methylprep process` to save output for ALL of the associated processing steps:<br>

- beta_values.pkl
- poobah_values.pkl
- control_probes.pkl
- m_values.pkl
- noob_meth_values.pkl
- noob_unmeth_values.pkl
- meth_values.pkl
- unmeth_values.pkl
- sample_sheet_meta_data.pkl

By default, the output is usually: 
- beta_values.pkl
- noob_meth_values.pkl
- noob_unmeth_values.pkl
- control_probes.pkl

The default settings are designed to match the output of `R`'s `sesame` processing. Prior to `methylprep v1.4.0`, the defaults matched `minfi`'s output.

Here are some high level options:

Argument | Type | Default | Description
--- | --- | --- | ---
`data_dir` | `str`, `Path` | *Required* | Base directory of the sample sheet and associated IDAT files
`minfi` | `bool` | `False` | Changes many settings to match `minfi` output. Default is `sesame`.

Use these options to specify file locations and array type:

Argument | Type | Default | Description
--- | --- | --- | ---
`array_type` | `str` | `None` | Code of the array type being processed. Possible values are `custom`, `27k`, `450k`, `epic`, and `epic+`. If not provided, the package will attempt to determine the array type based on the number of probes in the raw data. If the batch contains samples from different array types, this may not work. Our data `download` function attempts to split different arrays into separate batches for processing to accommodate this.
`sample_name` | `str` to list | `None` | List of sample names to process, in the CLI format of `-n sample1 sample2 sample3 etc`. If provided, only those samples specified will be processed. Otherwise all samples found in the sample sheet will be processed.
`manifest_filepath` | `str`, `Path` | `None` | File path for the array's manifest file. If not provided, this file will be downloaded from a Life Epigenetics archive.
`no_sample_sheet` | `bool` | `None` | pass in "--no_sample_sheet" from command line to trigger sample sheet auto-generation. Sample names will be based on idat filenames. Useful for public GEO data sets that lack sample sheets.
`sample_sheet_filepath` | `str`, `Path` | `None` | File path of the project's sample sheet. If not provided, the package will try to find one based on the supplied data directory path.

Use these options to specify what gets saved from processing, and how it gets saved:

Argument | Type | Default | Description
--- | --- | --- | ---
`no_export` | `bool` | `False` | Add to prevent saving the processed samples to CSV files.
`no_meta_export` | `bool` | `False` | Add to prevent saving the meta data to pickle files.
`betas` | `bool` | `False` | Add flag to output a pickled dataframe of beta values of sample probe values.
`m_value` | `bool` | `False` | Add flag to output a pickled dataframe of m_values of samples probe values.
`uncorrected` | `bool` | `False` | Saves raw florescence intensities in CSV and pickle output.
`save_control` | `bool` | `False` | Add to save control probe data. Required for some `methylcheck` QC functions.
`export_poobah` | `bool` | `False` | Include probe p-values in output files.
`bit` | `str` | `float32` | Specify data precision, and file size of output files (float16, float32, or float64)
`batch_size` | `int` | `None` | Optional: splits the batch into smaller sized sets for processing. Useful when processing hundreds of samples that can't fit into memory. This approach is also used by the package to process batches that come from different array types.
`poobah` | `bool` | `True` | calculates probe detection p-values and filters failed probes from pickled output files, and includes this data in a column in CSV files.

`data_dir` is the one required parameter. If you do not provide the file path for the project's sample_sheet CSV, it will find one based on the supplied data directory path. It will also auto detect the array type and download the corresponding manifest file for you.
<br>

---

## other commands

The `methylprep` CLI provides these top-level commands, which make it easier to use GEO datasets:

- `sample_sheet` will find/read/validate/create a sample sheet for a data set, or display its contents (given the directory of the data).
    This is also part of `process` and can be applied using the `--no_sample_sheet` flag.
- `download` download and process public data sets in NIH GEO or ArrayExpress collections. Provide the public Accession ID and `methylprep` will handle the rest.
- `beta_bake` combines `download`, `meta_data`, and file format conversion functions to produce a package
    that can be processed (with `process`) or loaded with `methylcheck.load` for analysis.
- `alert` scan GEO database and construct a CSV / dataframe of sample meta data and phenotypes for all studies matching a keyword
- `composite` download a bunch of datasets from a list of GEO ids, process them all, and combine into a large dataset
- `meta_data` will download just the meta data for a GEO dataset (using the MINiML file from the GEO database) and convert it to a samplesheet CSV

### `sample_sheet`

Find and parse the sample sheet in a given directory and emit the details of each sample. This is not required for actually processing data. ```methylprep``` will automatically create a sample sheet as part of the `process` or `download` pipelines.

optional arguments:

  Argument | Type | Description
  --- | --- | ---   
  -h, --help || show this help message and exit
  -d, --data_dir | `str` | Base directory of the sample sheet and associated IDAT files
  -c, --create | `bool` | If specified, this creates a sample sheet from idats instead of parsing an existing sample sheet. The output file will be called "samplesheet.csv".
  -o OUTPUT_FILE, --output_file OUTPUT_FILE | `str` | If creating a sample sheet, you can provide an optional output filename (CSV).                        


### `download`

There are thousands of publically accessible DNA methylation data sets available via the GEO (US NCBI NIH) https://www.ncbi.nlm.nih.gov/geo/ and ArrayExpress (UK) https://www.ebi.ac.uk/arrayexpress/ websites. This function makes it easy to import them and build a reference library of methylation data.

The CLI includes a `download` option. Supply the GEO ID or ArrayExpress ID and it will locate the files, download the idats, process them, and build a dataframe of the associated meta data. This dataframe format is compatible with `methylcheck` and `methylize`.

Argument | Type | Default | Description
--- | --- | --- | ---  
-h, --help |-|-|       show this help message and exit
  -d , --data_dir | `str` | [required path] | path to where the data series will be saved. Folder must exist already.
  -i ID, --id ID | `str` | [required ID] |The dataset's reference ID (Starts with `GSE` for GEO or `E-MTAB-` for ArrayExpress)
  -l LIST, --list LIST | `multiple strings` | optional | List of series IDs (can be either GEO or ArrayExpress), for partial downloading
  -o, --dict_only | `True` | pass flag only | If passed, will only create dictionaries and not process any samples
  -b BATCH_SIZE, --batch_size BATCH_SIZE | `int` | optional | Number of samples to process at a time, 100 by default.

When processing large batches of raw `.idat` files, specify `--batch_size` to break the processing up into smaller batches so the computer's memory won't overload. This is off by default when using `process` but is ON when using `download` and set to batch_size of 100. Set to 0 to force processing everything as one batch. The output files will be split into multiple files afterwards, and you can recomine them using `methylcheck.load`.

### `beta_bake`

`beta_bake` is a function intended for combining data that are not processed in exactly the same way. If IDAT files are present, it will download them for you to run `process` on. If there are no idats, but there is uncorrected methylated/unmethylated data, it will download that instead. If there is no unprocessed data, it will parse processed beta values for you. 

This is intended for creating datasets that sacrifice some data quality in exchange for size. For example, using a machine learning algorithm on 10,000 noisy samples can yield better results than using that algorithm on a more curated set of 1,000 samples. ML algorithms can be trained to read through the noise and benefit from more data to train on. 

Note: less than half of the GEO datasets include raw idat files! Most of the data on GEO has already been processed into beta values. This is a part of why `beta_bake` is so useful. 

Argument | Type | Default | Description
--- | --- | --- | ---  
-h, --help |||        show this help message and exit
-i ID, --id ID | `str` || GEO_ID of the dataset to download
  -d DATA_DIR, --data_dir DATA_DIR | `str` || Folder where series data will appear.
  -v, --verbose |||        if specified, this will turn on more verbose processing messages.
  -s, --save_source |||     if specified, this will retain .idat and/or -tbl-1.txt files used to generate beta_values dataframe pkl files.
  -b BUCKET, --bucket BUCKET | `str` || AWS S3 bucket where downloaded files are stored
  -e EFS, --efs EFS | `str`||    AWS elastic file system name, for lambda or AWS batch processing
  -p PROCESSED_BUCKET, --processed_bucket PROCESSED_BUCKET |`str`|| AWS S3 bucket where final files are saved
  -n, --no_clean|||        If specified, this LEAVES processing and raw data files in temporary folders. By default, these files are removed during processing, and useful files moved to data_dir.


```Shell
>>> python -m methylprep beta_bake -i GSE74013 -d GSE74013
INFO:methylprep.download.miniml:Downloading GSE74013_family.xml.tgz
INFO:methylprep.download.miniml:MINiML file does not provide (Sentrix_id_R00C00) for 24/24 samples.
INFO:methylprep.download.miniml:Final samplesheet contains 24 rows and 9 columns
```
Output file containing a beta_values pickled dataframe: GSE74013_beta_values.pkl.gz<br>

Output file containing meta data: GSE74013_GPL13534_meta_data.pkl.gz

### `composite`

A tool to build a data set from a list of public datasets. This function basically just loops `download` through the provided list of datasets. 

optional arguments:

Argument | Type | Description
--- | --- | ---  
  -h, --help || show this help message and exit
  -l LIST, --list LIST | `str`, filepath | A text file containg several GEO/ArrayExpress series ids. One ID per line in file. Note: The GEO Accession Viewer lets you export search results in this format.
  -d DATA_DIR, --data_dir DATA_DIR | `str`, filepath | Folder where to save data (and read the ID list file).
  -c, --control | `bool` | If flagged, this will only save samples that have the word "control" in their meta data.
  -k KEYWORD --keyword KEYWORD | `str` | Only retain samples that include this keyword (e.g. blood) somewhere in their meta data.
  -e, --export | `bool` | If passed, saves raw processing file data for each sample. (unlike meth-process, this is off by default)
  -b, --betas  | `bool` | If passed, output returns a dataframe of beta values for samples x probes. Local file beta_values.npy is also created.
  -m, --m_value  | `bool` | If passed, output returns a dataframe of M-values for samples x probes. Local file m_values.npy is also created.

### `alert`

Function to check for new datasets on GEO and update a csv each time it is run. Usable as a weekly cron command line function. Saves data to a local csv to compare with old datasets in <pattern>_meta.csv. Saves the dates of each dataset from GEO; calculates any new ones as new rows. updates csv.

optional arguments:

   Argument | Type | Description
   --- | --- | ---  
   keyword | `str` | Specify a word or phrase to narrow the search, such as "spleen blood".

### `meta_data`

Provides a more feature-rich meta data parser for public MINiML (formatted) GEO datasets. You can use `meta_data` to identify 'controls' or samples containing a specific keyword (e.g. blood,
tumor, etc) and remove any samples from sheet that lack these criteria, and delete the associated idats that don't have these keywords. After, run `process` on the rest, saving time. You can effectively ignore the parts of datasets that you don't need based on the associated meta data.

optional arguments:

Argument | Type | Description
--- | --- | ---
  -h, --help  ||   show this help message and exit
  -i ID, --id ID |`str`|   Unique ID of the series (the GEO GSExxxx ID)
  -d DATA_DIR, --data_dir DATA_DIR |`str` or `path`| Directory to search for MINiML file.
  -c, --control  |`str`| [experimental]: If flagged, this will look at the sample sheet and only save samples that appear to be"controls".
  -k KEYWORD, --keyword KEYWORD |`str`| [experimental]: Retain samples that include this keyword (e.g. blood, case insensitive) somewhere in samplesheet values.
  -s, --sync_idats | `bool` |      [experimental]: If flagged, this will scan the `data_dir` and remove all idat files that are not in the filtered samplesheet, so they won't be processed.
  -o, --dont_download  | `bool` | By default, this will first look at the local filepath (--data-dir) for `GSE..._family.xml` files. IF this is specified, it wont later look online to download the file. Sometimes a series has multiple files and it is easier to download, extract, and point this parser to each file instead.
