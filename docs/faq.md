# Tips and FAQs for the methylation suite

## loading processed data

If you have a small data set (under 200 samples), and did not use the `batch_size` option:

```python
import pandas as pd
data = pd.read_pickle('{path_to_file}/beta_values.pkl')
```

Otherwise, loading hundreds of samples can be slow. But there's a helper function that loads a bunch of smaller
batches into one dataframe, if they were processed with methylprep's `batch_size` option:

```python
import methylize
data = methylize.load(path_to_file)
```
That will load 1600 samples in 8 seconds, compared to taking minutes if you `read_pickle` them and `concat` them using pandas.

### Processing Notes

1. DataFrames are the format that methylcheck/methylize functions expect. For saving and loading, the probe names belong in columns and the sample names belong in the index.
But when processing data, some functions auto-transpose this to decrease processing time.
2. By default, `methylprep process` also creates a file called `sample_sheet_meta_data.pkl` from various data sources
   - if there is a GEO series MINiML xml file, it reads this data preferentially
   - if there is a samplesheet, it will convert this to a meta DataFrame
   - if these sources are missing required data, such as Sentrix_Position and Sentrix_ID, it will look for idat files and
   read these from the filenames.

## idat filenames
   - There are two acceptable formats:
     - `<GSM_ID>_<Sentrix_ID>_<Sentrix_Position>_<Red|Grn>.idat<.gz>`
     - `<Sentrix_ID>_<Sentrix_Position>_<Red|Grn>.idat<.gz>`
     - Methylprep will convert `.gz` files to `.idat` uncompressed files when processing.   
   - methylprep does not recognize the older 27k array filename format:
     `<GSM_ID>_<Sentrix_ID>_<Sentrix_Position>_<SOME_LETTER>.idat`

## Why didn't the `methylprep download` function work for GEO dataset GSEnnn?
A significant number of GEO datasets do not store their data in a consistent format. Here are some reasons a GEO dataset fails to download:

1. `idat` filename format is off (missing R00C00 position) The `download` function ONLY supports `.idat` files.
2. no raw idats in zip, only processed data
3. The meta data in MINiML format xml file is incomplete. (There are ways to make it work without meta data)
4. the NIH GEO FTP server is down (yes, we've experienced this whilst testing too)
5. `idat` files in dataset have varying number of probes. (If a dataset combines results from two array types (EPIC and 450k), it can sometimes split the data into two sample sheets and you can process each folder separately.)

In most cases where `download` fails, there will be processed data available in some other format.
Use `methylprep beta_bake` to download and convert this data instead.

In other cases where `beta_bake` also fails, we've made an effort to giving you clear error messages, detailing why. IF you detect a pattern that could be parsed with our code, we're happy to support additional
file formats for methylation data found on GEO.

## What if `methylprep download` and `methylprep beta_bake` fails. How do I use the data anyway?
1. Download the raw idat zipfile manually to a folder.
2. Uncompress it.
3. Confirm that there are `.idat` files present.
4. If the files end in `.idat.gz`, gunzip them first. (In a mac/linux bash window you can navigate to the folder and type `gunzip *` to uncompress all those files. On a PC, use some software like `7zip` or `winzip`.)
5. IF the filenames don't include Sentrix IDs and Sentrix array positions, like `<GSM ID>_<some long number>_<R01C01>_<Red or Grn>.idat`, you'll need to manually edit the samplesheet to match the files.
6. Run `methylprep process` on the folder, possibly with the `--no_sample_sheet` option. It should work, but you won't have any of the sample meta data bundled for you for analysis with datasets in nonstandard formats.

## How to process only part of a GEO dataset

The `methylprep meta_data` command line interface (CLI) option allows you to specify which samples to process using keyword matching against the GEO meta data. You can use this before `methylprep download` and `methylprep process` to create smaller data sets, faster.

You can also specify sample names individually from a large set like this: `methylprep process --sample_name Bob -n Suzy -n Doug`. This will reduce processing time.

### Examples

#### (1) Grab a samplesheet for a GEO data set

```python
python -m methylprep meta_data -i GSE125105
```
Yields three files on local disk:

- GSE125105_family.xml (the original GEO meta data file in MINiML format)
- GSE125105_GPL13534_samplesheet.csv (used for processing)
- GSE125105_GPL13534_meta_data.pkl (used in analysis to describe samples)

| GSM_ID     | Sample_Name                | Sentrix_ID | Sentrix_Position | source              | diagnosis | age | Sex | tissue      | cellcount-cd8t | cellcount-cd4t | cellcount-nk | cellcount-bcell | cellcount-mono | cellcount-gran | description                |
|------------|----------------------------|------------|------------------|---------------------|-----------|-----|-----|-------------|----------------|----------------|--------------|-----------------|----------------|----------------|----------------------------|
| GSM3562834 | genomic DNA from sample291 | 3999840035 | R01C01           | control_whole blood | control   | 73  | F   | whole blood | 0.07679        | 0.09099        | 0.06041      | 0.08542         | 0.09072        | 0.62266        | whole blood control sample |
| GSM3562835 | genomic DNA from sample612 | 3999840035 | R01C02           | case_whole blood    | case      | 32  | M   | whole blood | 0.05544        | 0.07946        | 0.0159       | 0.09557         | 0.05515        | 0.72663        | whole blood case sample    |
| GSM3562836 | genomic DNA from sample611 | 3999840035 | R02C01           | case_whole blood    | case      | 51  | F   | whole blood | 0.08279        | 0.22216        | 0.03107      | 0.0769          | 0.07915        | 0.54165        | whole blood case sample    |
| GSM3562837 | genomic DNA from sample375 | 3999840035 | R02C02           | case_whole blood    | case      | 30  | M   | whole blood | 0.03779        | 0.07368        | 0.00385      | 0.07548         | 0.0891         | 0.74809        | whole blood case sample    |

#### (2) Grab just the samplesheet and create a csv and dataframe pkl of it.

```python
python -m methylprep -v meta_data -i GSE84727 -d GSE84727
```

Verbose CLI output
```python
INFO:methylprep.download.miniml:Downloading GSE84727_family.xml
INFO:methylprep.download.miniml:Downloaded GSE84727_family.xml                                                   
INFO:methylprep.download.miniml:Unpacking GSE84727_family.xml
INFO:methylprep.download.miniml:Downloaded and unpacked GSE84727
INFO:methylprep.download.miniml:MINiML file does not provide `methylprep_name` (sentrix_id_R00C00) for 847/847 samples.
INFO:methylprep.download.miniml:dropped Sentrix_ID (empty)
INFO:methylprep.download.miniml:dropped Sentrix_Position (empty)
INFO:methylprep.download.miniml:source == description; dropping source
INFO:methylprep.download.miniml:dropped `platform` ({'GPL13534'})
INFO:methylprep.download.miniml:title == Sample_Name; dropping title
INFO:methylprep.download.miniml:Final samplesheet contains 847 rows and 7 columns
INFO:methylprep.download.miniml:['GSM_ID', 'Sample_Name', 'sentrixids', 'Sex', 'age', 'disease_status', 'description']
```
Resulting samplesheet

|GSM_ID	|Sample_Name	|sentrixids	|Sex	|age	|disease_status	|description|
|-------|---------------|-------|--------|-------|-----------|--------------|
|GSM2250273	|33262604	|3998567027_R01C01	|M	|47.3	|1	|control blood|
|GSM2250274	|33261623	|3998567027_R02C01	|M	|60.4	|1	|control blood|
|GSM2250275	|33262614	|3998567027_R04C01	|M	|30.1	|1	|control blood|


You'll notice that one column `sentrixids` is misnamed. It should be split into `Sentrix_Position` and `Sentrix_ID` columns for processing to work on this GEO series.
You can edit the csv and fix that prior to running the pipeline with `methylprep process`. If you don't, you'll get an error saying, "methylprep could not find the samplesheet." This error is caused by the researcher putting an arbitrary name into the GSE84727_family.xml `MINiML` meta data file:

```
      <Characteristics tag="sentrixids">
3998567027_R02C02
      </Characteristics>
```
If you use the `methylprep download` option by itself, it can generally avoid this type of XML parsing error, but it will download everything. Doing analysis on just part of a dataset requires some debugging like this.


#### (3) Samplesheet with only "normal" samples
```python
python -m methylprep -v meta_data -i GSE52270 -d GSE52270 -k normal

```

`-k` is shorthand for `--keyword`. The resulting sample sheet only includes samples that include the keyword `normal`

|GSM_ID |	Sample_Name	|source	|disease state	|description|
|-------|---------------|-------|---------------|-----------|
|GSM1278809	|Colon 61	|Large Intestine	|normal	|1011N|
|GSM1278812	|Colon 64	|Large Intestine	|normal	|1082N|
|GSM1278823	|Colon 75	|Large Intestine	|normal	|1184N|
|GSM1278825	|White matter 77	|Central Nervous System	|normal	|12_03|



#### (4) Generate filtered samplesheet with only control samples from blood

```
python -m methylprep  meta_data -i GSE125105 -d GSE125105 --control -k blood
```

| GSM_ID     | Sample_Name                | Sentrix_ID | Sentrix_Position | source              | diagnosis | age | Sex | tissue      | cellcount-cd8t | cellcount-cd4t | cellcount-nk | cellcount-bcell | cellcount-mono | cellcount-gran | description                |
|------------|----------------------------|------------|------------------|---------------------|-----------|-----|-----|-------------|----------------|----------------|--------------|-----------------|----------------|----------------|----------------------------|
| GSM3562834 | genomic DNA from sample291 | 3999840035 | R01C01           | control_whole blood | control   | 73  | F   | whole blood | 0.07679        | 0.09099        | 0.06041      | 0.08542         | 0.09072        | 0.62266        | whole blood control sample |
| GSM3562839 | genomic DNA from sample176 | 3999840035 | R03C02           | control_whole blood | control   | 43  | M   | whole blood | 0.06946        | 0.12989        | 0.04703      | 0.09808         | 0.14105        | 0.54662        | whole blood control sample |
| GSM3562842 | genomic DNA from sample161 | 3999840035 | R05C01           | control_whole blood | control   | 44  | F   | whole blood | 0.10986        | 0.13565        | 0.07657      | 0.09125         | 0.12521        | 0.49317        | whole blood control sample |
| GSM3562846 | genomic DNA from sample270 | 3999840037 | R01C01           | control_whole blood | control   | 64  | M   | whole blood | 0.11508        | 0.14116        | 0.0679       | 0.09415         | 0.15311        | 0.47829        | whole blood control sample |
| GSM3562855 | genomic DNA from sample162 | 3999840037 | R06C01           | control_whole blood | control   | 65  | M   | whole blood | 0.01668        | 0.14318        | 0.1096       | 0.05545         | 0.09695        | 0.60283        | whole blood control sample |

This only retains 211 of the 699 samples in a samplesheet. Next, you download the `.idat` files with with `methylprep download` and then remove the `idat` files you won't need like this:

```python
python -m methylprep -v download -i GSE125105 -d GSE125105
python -m methylprep -v meta_data -i GSE125105 -d GSE125105 --sync_idats --control -k blood
```

And then process the 6.1GB file using this samplesheet, like this:
```python
python -m methylprep -v process -d GSE125105 --betas --m_value -e
```
Those options will create two big files. One is a dataframe of beta_values for each sample. The other, m_values for each sample (kind of the same thing, but sometimes you want m_values). the `-e` or `--no_export` option will suppress the function creating files of probe values for each sample, as these are not needed by most `methylize` and `methylcheck` functions. There is also a `--save_uncorrected` option that prevents any sort of background and NOOB signal enhancement during processing. Uncorrected files are needed for a few analysis functions, namely `p-value probe detection`.

In general, partial-dataset processing fails because the meta data for a GEO dataset is incomplete. Either the array positions are missing, or misnamed. Careful checking can allow one to fix this and build a large data set from multiple GEO datasets.

#### (4b) Another condensed example of downloading GEO data and only processing control samples

```python
python -m methylprep -v download -i GSE130030 -d GSE130030
# next, remove the treatment samples using `-c` and remove extra idats with `-s`
python -m methylprep -v meta_data -i GSE130030 -d GSE130030 --control -s
# finally, process it
python -m methylprep -v process -d GSE130030 --betas --m_value --no_export
```

This creates two files, `beta_values.pkl` and `GSE130030_GPL13534_meta_data.pkl`, that you can work with in `methylize` like this:

Navigate to the `GSE130030` folder created by `methylrep`, and start a python interpreter:
```python
import methylize
data,meta = methylize.load_both()
INFO:methylize.helpers:Found several meta_data files; using: GSE130030_GPL13534_meta_data.pkl)
INFO:methylize.helpers:loaded data (485512, 14) from 1 pickled files (0.159s)
INFO:methylize.helpers:meta.Sample_IDs match data.index (OK)
```

Or if you are running in a notebook, specify the full path:

```python
import methylize
data,meta = methylize.load_both('<path_to...>/GSE105018')
```

## Why won't `methylprep composite` parse some GEO data sets' meta data?

Here are some examples of logical, but unexpected ways data can be stored in the MINiML file format:

The composite expects "control" to be in one of the rows in a spreadsheet. Instead, the authors have recoded "control" as a number, "1" in the column header name. Our parser just isn't smart enough to read that.

```
['GSM_ID', 'Sample_Name', 'diseasestatus (1=control, 2=scz patient)', 'source', 'gender', 'sample type', 'plate', 'sentrix barcode', 'sentrix position', 'well id', 'age',  'used_in_analysis', 'description']
```
Here, instead of having the same names for data for each sample, they have split the smoking status into a bunch of columns, and not provided values for every sample. (smoking_evernever and smoke_free_years don't add up to 95.) Fixing this requires putting in null values for each incomplete column in a sample sheet.

```
ValueError - array lengths vary in sample meta data: [('GSM_ID', 95), ('Sample_Name', 95),
('smoking_evernever', 52), ('smoke_free_years', 30), ('Sentrix_ID', 95), ('Sentrix_Position', 95), ('source', 95), ('gender', 95), ('slide', 95), ('array', 95), ('array_name', 95), ('sample_group', 52),  ('smoking_5years', 52), ('ms_case_control', 52), ('sample_year', 52), ('age_sampling', 52), ('py', 52),  ('description', 95), ('Sample_ID', 95)]
```
