# Public GEO datasets

`methylprep` provides methods to use public data in a variety of formats.
- idat
- processed tab delimited (`txt`)
- processed `csv`
- processed `xlsx`
- pickled dataframes (`pkl`) created using `methylprep process` or `methylprep.run_pipeline`
  - format: dataframe should have probe names as columns or rows, and sample probe values in the other dimension.
  - meta data: as a dataframe with values for samples, so long as one of those characteristics, the sample name, matches the Sentrix_Position sample name that is the default output of Illumina arrays.
  - Note: if you use `methylprep.load_both` and pass in the folder location of your methylprep processed data, it will construct the meta data frame for you. Ultimately, it is reading the GEO `miniml` format (XML) file for the public data set, or your samplesheet if you provided one.


### download from GEO

```Shell
(base) $ python -m methylprep download -i GSE122126 -d GEO/GSE122126
INFO:methylprep.download.geo:Downloading GSE122126_family.xml
GSE122126:   3%|█▉                                                            | 12.3M/407M [00:07<05:57, 1.10Mb/s]

INFO:methylprep.download.geo:Downloaded GSE122126_family.xml
INFO:methylprep.download.geo:Unpacking GSE122126_family.xml
GSE122126:   7%|████▎                                                          | 121M/1.77G [01:24<42:48, 644kb/s]

```

If you choose a dataset that lacks raw idat files, it will warn you.

```Shell
(base) $ python -m methylprep download -i GSE123211 -d GEO/GSE123211
ERROR:methylprep.download.process_data:[!] Geo data set GSE123211 probably does NOT contain usable raw data (in .idat format). Not downloading.
ERROR:methylprep.download.process_data:Series failed to download successfully.
```

If you want to use the author's processed data instead of reprocessing it yourself,
download the `.gz` file using a web browser, then `gunzip` it to create a `txt | pkl | xlsx | csv` file,
and then load that using `methylprep.read_geo`.

### loading processed GEO data

```python
import methylprep
import methylcheck
from pathlib import Path

df = methylprep.read_geo(Path('~/Downloads', 'GSE115278_Matrix_processed.txt'))
# or
df = methylprep.read_geo(Path('~/Downloads', 'GSE111165_data_processed_detection_p_val_EPIC.csv'))
methylcheck.beta_density_plot(df)
```

![Fig.19](tutorial_figs/fig19.png)
