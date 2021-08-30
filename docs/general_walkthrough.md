# General Walkthrough

Here we will provide a few examples of how to use various `methylprep` functions. We'll focus on examples using the CLI, as that is the recommended interface for `methylprep`, but there's also a section near the end of the tutorial that demonstrates how to run `methlyprep` from an IDE.

Note: sample sheets are recommended but not necessary to run `methylprep`. We will cover a few ways to work with a data set without a sample sheet. 

## Set Up

If you haven't done so already, run this command from your terminal to install methylprep:
```shell
>>> pip install methylprep
```
## Downloading from GEO
The first step in this tutorial will be using ```methylprep``` to acquire a dataset from GEO. ```methylcheck``` and `methylize`  work best when users include a sample sheet along with the IDAT files. The meta data contained in sample sheets is useful when running QC or analyses. 

When downloading from a GEO dataset, `download` will attempt to find and download the associated sample sheet. If there is none, ```methylprep``` will automatically generate one. Users may make one of their own as well. The [Illumina sample sheet](https://support.illumina.com/downloads/infinium-methylationepic-sample-sheet.html) is the standard format. 

Note: if you already know you'd like to work with pre-processed GEO data (like beta values), check the [beta_bake](docs/special_cases.md#using-beta_bake-for-preprocessed-data) section of the special cases tutorial. 

### Our GEO dataset
For our tutorial, we will download GEO data from [this experiment](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE147391) where researchers examined differential methylation among high and low grade gliomas (an aggressive form of brain tumors). This is a relatively small dataset (16 samples) of surgical resections processed using the Infinium MethylationEPIC assay.

The authors analyzed methylation data from 6 II-grade, 2 III-grade and 8 IV-grade gliomas from individual surgical resections. Stage I and II gliomas are defined as low grade, while stage III and IV are defined as high grade. 

We will run the `download` command with the `-v` option to get more information about what processes are running. This should be run from the command line. 

```shell
>>> python -m methylprep -v download -i GSE147391 -d <filepath>
```
 Where `<filepath>` is the directory you would like the processed data and output files to be stored.

```methylprep``` will search the GEO database and find the dataset we indicated. It will also find the file with associated metadata and parse through it, indicating any extra fields in the sample sheet. In our case, there are 5 additional fields, which will come in handy in later analysis: 

```shell
INFO:methylprep.processing.pipeline:Found 5 additional fields in sample_sheet: source | histological diagnosis --> histological_diagnosis | gender | description | Sample_ID
```

```methylprep``` will begin downloading, unpacking, and processing the IDAT files after downloading the file manifest and the parsing sample sheet. The automatic processing can be turned off with the `-o` option. 

Processing time depends on the type of array (450k arrays process faster than EPIC arrays) and on the size of the dataset. Datasets with more than 100 samples will, by default, be chunked into batches of 100 for more efficient processing. Batch size is adjustable with the `--batch_size` argument. 

After the files are processed, we're ready to load the files into `methylcheck` for QC. See `methylcheck` documentation for instructions and more details.

## Processing Your Own Data
It is often the case that users have their own idat files that they would like to process.  Instead of using the ```download``` command, we will use the ```process``` command. This is the main workhorse of the ```methylprep``` package, and includes extra options for customization that the ```download``` command lacks.

Users should note that the sample sheet is optional, but the ***manifest file*** is not. Make sure there is a manifest file included with the IDAT data, *especially* if you are analyzing data from a custom array, as ```methlyprep``` will use the manifest file to determine what probes are included in the custom array, which control probes are present, etc. 

Note: If users are interested in processing their own files, the folder containing the IDATs, manifest, and sample sheet needs to be unzipped before being processed. 

Once again, we'll run this from the command line:

```shell
>>> python -m methylprep process -d <filepath> --all
```

`<filepath>` specifies where the manifest and IDAT files (and sample sheet, if any) are stored. 

The `--all` option at the end tells ```methylprep``` to save output for ALL of the associated processing steps.<br>

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

These files will be saved in the same directory as the data. ```methylprep``` will also create two folders to store processed .csv files. If users are interested in the probes that failed poobah, they are included in the .csv file outputs.

The default output is all that is needed to run qc reports from ```methylcheck```. However, other useful outputs like the poobal_values.pkl can optionally be included in the qc functions for greater details on sample quality.  

## Using methylprep from a Jupyter Notebook
`methylprep` also offers a scikit-learn style interface for users to run within jupyter notebooks or a similar IDE. We recommend using the CLI for `methylprep`--IDEs tend to process more slowly, especially for large batches--but users are able to get most of the package's functionality from within an IDE.

## `methylprep.run_pipeline`
This command is the IDE equivalent of `methylprep process`. The output, if left to default, is a set of data containers. Alternatively, users may specify `betas=True` or `m_values=True` to get a dataframe of their chosen values. We also recommend setting `export=True` to get output files saved as well. That way you can easily load the data in the future instead of running this command every time you open the notebook. The only required argument for this function is the directory where the raw data is stored. 

`run_pipeline` will take almost all of the same arguments that `process` will take from the CLI.

```python
from methylprep import run_pipeline
from pathlib import Path
filepath = Path('~/tutorial/GPL21145')

data_containers = run_pipeline(filepath, export=True)
```

For more information on optional arguments, exports, etc, run `help(methylprep.run_pipeline)`


## `methylprep.make_pipeline`
Users may specify ['all'] for both the 'steps' and 'exports' arguments to get the same functionality as the CLI code `methlyprep process -d <filepath> --all`. Alternatively, if users are interested in customizing their pipeline, they may list the steps they'd like to include and leave out any they would like to skip. 

```python
import methylprep
from pathlib import Path
filepath = Path('~/tutorial')

# choose any combination of:
# ['infer_channel_switch', 'poobah', 'quality_mask', 'noob', 'dye_bias']
mysteps = ['all']

# choose any combination of:
# ['csv', 'poobah', 'meth', 'unmeth', 'noob_meth', 'noob_unmeth', 'sample_sheet_meta_data', 'mouse', 'control']
myexports = ['all'] 

# make_pipeline will auto-run the pipeline 
# make_pipeline will also take the same kwargs as run_pipeline
methylprep.make_pipeline(data_dir=filepath, steps=mysteps, exports=myexports, estimator=None)
```

For more information on this function, run `help(methylprep.make_pipeline)`