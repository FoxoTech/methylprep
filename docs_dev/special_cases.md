# welcome to the specialized functions of `methylprep` tutorial!

We cover the most high level use cases in our general walkthrough. However, there are some functions available in `methylprep` for less common (more specialized) use cases that we'll cover here. 

## Building a composite dataset with `alert` and `composite`

`methylprep` includes a few functions that can be used in a pipeline to merge datasets from GEO into one larger dataset. Say, for example, that you are particularly interested in investigating data from glioblastoma (an aggressive form of cancer that can occur in the brain or spine). 

You may use the `alert` command to compile a csv of GEO datasets that mention the keyword glioblastoma. The csv will contain, at most, 500 rows (a query limit set by GEO). 

```shell
>>> python -m methylprep alert -k "glioblastoma"
```

After looking through the .csv file and eliminating any datasets that aren't useful, you can pull the column of GSE id's and run `methylprep composite` from the CLI.  `composite` will run the `download` command in a loop given a text file with a list of GSE id's. 

From our geo alert glioblastoma .csv file, here are 4 datasets that we might pull into one larger composite dataset. One thing to be aware of is to make sure that the platform types are all the same. Otherwise you may be combining 450k array data with EPIC array data. 

Here is a peek at the txt file with the list of GSEs. Make sure that the format is such that there is just one GSE id per line. 

```shell
>>> cat glioblastoma.txt 
GSE161175
GSE122808
GSE143842
GSE143843
```

Now we can run the `composite` command and it'll download/process data from each GSE in the list. You need to use `-l` to specify the .txt file that contains the GSEs and `-d` to specify the directory where to save the data. 

```shell
 >>> python -m methylprep composite -l ~/tutorial/glioblastoma.txt -d ~/tutorial
```

The `meta_data` function has a similar application for building a composite dataset when using the `--control` option.  See below for the walkthrough of how to use that command to build a dataset of healthy control samples instead of experimental/disease case samples. 

## using the `meta_data` function

### Creating a meta_data.pkl file
GEO data will always be accompanied by a MINiML file (an .xml file that contains information about the platform, series, and samples). `methylprep` will handle downloading and parsing this as part of the `download` or `process` functions, with the output stored in the meta_data.pkl file. So if you're going to run either of those, you don't need to worry about running `meta_data` separately. 

However, if you are interested in creating a meta_data.pkl file without downloading/processing data, the `meta_data` function will come in handy. 

```python
>>> python -m methylprep meta_data -i GSE147391 -d ~/tutorial
INFO:methylprep.download.miniml:found 16 idat files; updated meta_data.
INFO:methylprep.download.miniml:Final samplesheet contains 16 rows and 9 columns
```

### Filtering samples with `meta_data`
`meta_data` also has a few experimental commands that may be useful in filtering unwanted samples and reducing the time needed to process datasets. For example, consider the glioma dataset from [GSE147391](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE147391) (the dataset discussed in the general walkthrough). This dataset contains methylation profiles of high and low grade gliomas from 16 patients (8 high grade gliomas and 8 low grade gliomas). If we wanted to ignore the high grade gliomas and focus on the low grade gliomas, we could use `meta_data`'s `--keyword` (`-k` for short) command to search for samples that are indicated as "grade II" in the sample sheet and only retain those values. 

```shell
>>> python -m methylprep meta_data -i GSE147391 -d ~/tutorial -k "grade II"
INFO:methylprep.download.miniml:found 16 idat files; updated meta_data.
INFO:methylprep.download.miniml:Filtering keyword: Retained 8/16 samples matching `grade II``.
INFO:methylprep.download.miniml:Final samplesheet contains 8 rows and 9 columns
```

Now our meta_data.pkl file only contains information on the low grade gliomas, and we can exclude the other samples from our future analyses. If we haven't already used `methylprep process` on our data, we can also include the `--sample_name` or `-n` flag with the list of the samples we want to run (the low grade gliomas) to save time processing the data. 


### Building a composite dataset using `meta_data`
To better demonstrate the use case for `meta_data's --control` command, we will work with a new dataset.  [GSE163970](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE163970) is a dataset that examines differential methylation in 253 patients with Paget's disease of bone (PDB) as compared to 280 controls. 

We might be able to cut down processing time for this large dataset if we are just interested in the control data and want to exclude the PDB samples. One example of a scenario in which this would be helpful is building a large composite dataset that includes only healthy control samples from multiple GEO datasets. 

Before downloading any data from GEO, we'll run this:

```shell
>>> python -m methylprep meta_data --control -i GSE163970          
```

`methylprep` will look through our computer and realize that we don't have any data or metadata files downloaded. It will then pull the MINiML file from the GEO database and search it for samples with the word "control" or "ctrl" in the metadata.

```shell
INFO:methylprep.download.miniml:Downloading GSE163970_family.xml.tgz
INFO:methylprep.download.miniml:Filtering controls: Retained 260/492 samples matching ['control', 'ctrl'].
INFO:methylprep.download.miniml:Final samplesheet contains 260 rows and 13 columns
```
As we noted above, these functions may not work on every dataset. `methylprep` does its best to identify controls based on keywords like 'ctrl' or 'control', but controls may be labelled differently in various datasets. 


## `beta_bake`
`beta_bake` walkthrough coming soon!


<br><br>
---

