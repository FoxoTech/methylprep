# welcome to the specialized functions of `methylprep` tutorial!

We cover the most high level use cases in our general walkthrough. However, there are functions available in `methylprep` for less common use cases that we'll cover here. 

## building a composite dataset
alert -> meta_data --controls/--keyword -> composite

walk through pipeline to build a custom composite dataset using these functions

## beta_bake
walk through beta_bake use case


<br><br>
---

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
`meta_data` also has a few experimental commands that may be useful in filtering unwanted samples and reducing the time needed to process datasets. For example, consider our glioma dataset. If we wanted to ignore the high grade gliomas and focus on the low grade gliomas, we could use `meta_data`'s `--keyword` (`-k` for short) command to search for samples that are indicated as "grade II" in the sample sheet and only retain those values. 

```shell
>>> python -m methylprep meta_data -i GSE147391 -d ~/tutorial -k "grade II"
INFO:methylprep.download.miniml:found 16 idat files; updated meta_data.
INFO:methylprep.download.miniml:Filtering keyword: Retained 8/16 samples matching `grade II``.
INFO:methylprep.download.miniml:Final samplesheet contains 8 rows and 9 columns
```

Now our meta_data.pkl file only contains information on the low grade gliomas, and we can exclude the other samples from our future analyses. If we haven't already used `methylprep process` on our data, we can also include the `--sample_name` or `-n` flag with the list of the samples we want to run (the low grade gliomas) to save time processing the data. 


### Building a composite dataset using `meta_data --control`
To demonstrate the use case for `meta_data's --control` command, we need to work with a new dataset.  [GSE163970](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE163970) is a dataset that examines differential methylation in 253 patients with Paget's disease of bone (PDB) as compared to 280 controls. These samples were run on 450k arrays, but since the dataset is so large, it will still take a long time to process.

We might be able to cut down that processing time if we are interested in the control data and want to exclude the PDB samples. One example of a scenario in which this would be helpful is building a large composite dataset that includes only healthy control samples from multiple GEO datasets. 

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
As we noted above, these functions may not work on every dataset. Our code does its best to identify controls based on keywords, but controls may be labelled differently in various datasets. 