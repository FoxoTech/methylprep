# Adding New Array Type
### (based on mouse array notes from March/April 2020)
### this is a developer's note on how to add a new Illumina array type, when it becomes available.

## MANIFEST STRUCTURE
(for each, normals include first control probe at end)

- EPIC
  - controls at end, with one header row
  - 635
- EPIC+
  - controls at end, with NO header row
  - 635
- mouse
  - controls after normals, with NO header row
  - 635
- 450k
  - controls at end, with one header row
  - normals: all but controls
  - 850 controls
- mouse:
  - 577 probes are missing (expected 0)
  - meth has 0, unmeth has 577 missing
  - testing: these are all IG/IR type. update processing code to treat thse are type I when calculating unmeth values.
  - (UNMETHYLATED_PROBE_SUBSETS looks ok. treated as 'I'.)
  - NOTE: I realized there are no type 'I' probes in this manifest. Converted/renamed IR and IG to 'I'.

## REQUIRED EDITS TO MOUSE MANIFEST

- update methylprep with mouse manifest (csv)
- replace columns using textedit first
- replace NA with '' everywhere in csv
- renaming
  - M == AddressA_ID
  - U == AddressB_ID
- use python on full dataframe of manifest csv  

```
man = man.rename(columns={'U':'AddressA_ID', 'M':'AddressB_ID', 'Probe_ID':'IlmnID', 'DESIGN':'Infinium_Design_Type', 'COLOR_CHANNEL':'Color_Channel', 'AlleleA_Probe_Sequence': 'AlleleA_ProbeSeq', 'AlleleB_Probe_Sequence': 'AlleleB_ProbeSeq'})
```

- methylprep expects the cgxxxxxx part of IlmnID to be in a 'Name' column in manifest.
  - add Name col (pull out part of IlmnID)

B1 V2 solution
```
pattern = re.compile('([0-9a-zA-Z]+(_\d+)?)')
man['Name'] = man.apply(lambda x: re.match(pattern, x['Probe_ID']).groups()[0], axis='columns')
```

B3 V2 solution (names were more complex, with extra underscores)
I noticed that all non-control probe names either had a `cg` or a '_F' or a '_R' in them, before the illumina-junk I wanted to cut.
```
names = []
for val in df['IlmnID']:
    if re.match('(cg\d+)', val):
        name = re.match('cg\d+', val).group()
        names.append(name)
    elif '_F' in val:
        names.append(val.split('_F')[0])
    elif '_R' in val:
        names.append(val.split('_R')[0])
    else:
        names.append(val)
print(len(names), df.shape)
print(names[1020:1060])
df['Name'] = pd.Series(names)
# add and reorder 'Name' as 2nd column, to match order of old manifest
df.insert(1, 'Name', pd.Series(names))
df
```

- Mouse Array had to new types of probes.  Normal probes are [cg, ch] and mouse-specific were [rs mu] types.
- manifest missing
  - missing genome build
  - missing CHR
  - missing MAPINFO
  - missing Strand
- When saving a pandas datafram as a CSV, remember to remove the index and order columns: IlmnID must be FIRST.
  - `man.to_csv('LEGX_B1_manifest_mouse_v1_min.csv', index=False, index_label=False)`
  - for reordering columns to put IlmnID first: https://stackoverflow.com/questions/37070759/preserving-column-order-in-the-pandas-to-csv-method
- gzip and copy to hidden manifest folder, '~/.methylprep_manifest_files'
  - also copy to the S3 `array_manifests` bucket: https://s3.console.aws.amazon.com/s3/buckets/array-manifest-files/?region=us-east-1
  - and to the pipeline's private copy of same file (dunno why it exists, but it does)

- MISSING control probes.
  - found them. they start with ctl_.... but methylprep expected other names.
  - `ctrl = man.loc[man['IlmnID'].str.contains('ctl')]`

```
names = [name[:2] for name in list(man['IlmnID'])]
types = Counter(names)
types.most_common(50)
[('cg', 266086), ('rp', 4344), ('ch', 2746), ('ct', 635), ('rs', 536), ('mu', 43)]
```

- breakdown of probe categories
  - cg + ch are "normal"
  - 'mu' = multi-index
  - 'rp' = reapeated probes
  - 'ct = ctl (controls) 635, like EPIC
  - 'rs' = SNPs
- manifest requirements
  - normal (cg,ch) must be at start of manifest file.
  - controls must be AFTER last normal probe in file (update: they can be anywhere after normals)
  - checked: rs SNPs will load ANYWHERE in file. Loaded using regex.
  - (based on these rules, I concluded it must be cg,ch,ctl, ... rs order)
- check total probe count (for code edits later)
  - num_probes = 268832 (cp+ch in manifest)
- to test under the hood:
  - `m = Manifest(ArrayType.ILLUMINA_MOUSE, 'LEGX_B1_manifest_mouse_v1_min.csv.gz')`

- BECAUSE control probes only use the first four columns in all other manifests (and the column names don't matter), this manifest needs to use a specific order for first four:
  - IlmnID
  - AddressA_ID --> Address_ID
  - Color_Channel --> Color
  - Extended_Type --> (2nd half of Probe_ID, after removing the 'ctl_' prefix)
  - [just copied the epic manifest control section into this manifest, because it should be the same, and already formatted to work]
- `man.to_csv('LEGX_B1_manifest_mouse_v1_min.csv', index=False, index_label=False)`

- manifest AddressB_ID does not need to be int. floats work too.
- renamed 'IG' and 'IR' to type 'I', as methylprep doesn't recognize 'IG' and 'IR'.

```
man['Infinium_Design_Type'] = man['Infinium_Design_Type'].replace({'IG':'I', 'IR':'I'})
man.to_csv('LEGX_B1_manifest_mouse_v1_min.csv', index=False, index_label=False)
```

# Manifest: ensure same order for all columns, because dtypes() of first four columns matter.
# this also ensures prev col names are all there and spelled the same way
# to reorder DF columns...
```
df = df[['IlmnID','Name','AddressA_ID','AlleleA_ProbeSeq','AddressB_ID','Infinium_Design_Type',
'Color_Channel','Genome_Build','CHR','MAPINFO','Strand','AlleleB_ProbeSeq', 'Probe_Type']]
```
BUG whilst Reading manifest file: LEGX_B3_manifest_mouse_v2_min.csv:
TypeError: Cannot cast array from dtype('O') to dtype('float64') according to the rule 'safe'
ValueError: could not convert string to float: 'DNP(20K)'

After reordered the new column order to match the old column order precisely, it worked.


## list of code edits by file in methylprep for mouse support
- manifest.py
  - download_default bypassed mouse b/c this manifest was stored locally already. but would require S3 upload to test.
  - read_probes -- nrows load was off by one for ALL manifest types. Not sure how this affected everything past, but tried to resolve.
- arrays.py
  - ArrayType (ILLUMINA_MOUSE)
  - detect array - from_probe_count
  - probe_count?? run process on one mouse idat file to get it from the error code
    - `Unknown array type: (315639 probes detected)`
	- num_probes (274390 rows in manifest)
	- num_controls (from wanding's powerpoint)
	- num_snps (from Bret's screen share 'sn'? 251 'rs')
- probes.py (from_manifest_values)
  - IR, IG types (I added support for these, even though IG/IR was changed to I in manifest. but future manifests may need this.)
- raw_dataset.py
  - filter_oob_probes (temp added log messages to debug) but nothing changed codewise.
  - oob_green is empty DF (filter_oob_probes was not working, because some of the failed probes, I think).
- SampleDataContainer now has a mouse_probes variable. this gets exported.
  - calls postprocess.consolidate_mouse_probes
  - using container.mouse_meth / .mouse_unmeth values
- added to postprocess.py
  - this will separate mouse probes and mouse snps after processing in a postprocessing step.
  - splits mouse from normal probes
    - confirms normal probe values are only cg/ch probes
  - added to run_pipeline a step export these mouse probes, if they exist in SampleDataContainer
    - consolidate_mouse_probes()
- additional edits I did, then removed later

```
DID, then UNDID:
x added MethDataSet.mouse_meth/unmeth methods
x added ProbeType.MOUSE_ONE / MOUSE_TWO
	(like SNP_ONE and SNP_TWO)
x updated ProbeType.from_manifest_values()
	[note: pain later to split mu from rp]
	will return MOUSE_ONE or MOUSE_TWO by
	reading probe name patterns.
x added to probes METHYLATED_MOUSE_PROBES
	UNMETHYLATED_MOUSE_PROBES
x add these to models/__init__.py
x add this to meth_dataset.py imports
x added SampleDataContainer to use this methdataset
```

## testing

```
python -m methylprep -v process -d idats_BETA_LEGX_mm10_RND1 --betas --no_sample_sheet
```
GOT batch to run, except for 577 IG/IR probes (listed as "I" now in manifest version)
(later: these are probably failed probes at the synthesis level)

2nd result: works, but now we have 900-ish missing probes. All the mu and rp probes do appear to work in the processed CSV now because they're treated as "normal" probes (fix: keep separate next) easiest way to keep separate is to undo my num_probes to stop at last cg/ch probe. trying: separate them out in postprocess step. and confirm normals are only cg/ch

#### other bugs discovered:
  - 27k controls: non-existent but arrays.py says 140 exist. not in manifest.
  - BUG?? (check to see if SNPs are being loaded twice (along with cg+ch)

#### data flow through methylprep

```
from_manifest_values() <<- read_probes() <<- manifest.data_frame['probe_type'] <<- get_probe_details() <<-- _get_subset_means() <<- probes.py SUBSETS defined (here snps are separated from non-snp probes.)

so SampleDataContainer
	RawDataSet
	  MethDataset
        METHYLATED_PROBE_SUBSETS
```
