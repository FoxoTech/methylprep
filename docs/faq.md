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
     
