# Release History

## v1.5.0
- MAJOR refactor/overhaul of all the internal classes. This was necessary to fully support the mouse array.
  - new SigSet class object that mirror's sesame's SigSet and SigDF object.
    - Combines idats, manifest, and sample sheet into one object that is inherited by SampleDataContainer
  - RawDataset, MethylationDataset, ProbeSubtype all deprecated and replaced by SigSet
  - SamleDataContainer class is now basically the SigSet plus all pipeline processing settings
  - new mouse manifest
  


## v1.4.7
- mouse manifest updated to conform with illumina Genome Studio / sesame probe naming convention.
  - mouse_probes.pkl now includes different probe types. Previously, if a probe type was 'mu' (multi)
  or 'rp' (repeat) or IlmnID started with 'uk' (unknown?), it was moved to experimental mouse_probes.pkl.
  This was about 6300 probes.
  Now, all 'Multi' and 'Random' probes are moved and stored in mouse_probes.pkl, about 30,000.
  - mouse manifest has a 'design' column with tons of human-readable notes on different probe origins,
  including analogous EPIC human-mapped probes.

## v1.4.6
- pipeline CSV output will now include meth, unmeth, beta, and m-values for all probes, including failed probes.
    version 1.4.0 to 1.4.5 was replacing these values with NaN if a probe was filtered by the quality_mask.
    Pickled beta, M-value, noob_meth, noob_unmeth output files will continue to exclude (e.g. show NaN) probes that failed poobah_pval or quality_mask.

## v1.4.5
- fixed qualityMask for epic+

## v1.4.4
- faster circleci testing
- mouse probes have duplicate names, breaking dye-bias step, so it will fallback to linear-dye when duplicates are present
- added more mouse array test coverage

## v1.4.0
- now uses sesame's infer_type_I_channel function to detect and correct probe switching, if sesame=True
- uses sesame's nonlinear dye bias correction function, if sesame=True
    instead of the previous linear-dye-correction in the NOOB function.
- as part of the run_pipeline(sesame=True) default ON settings, it will apply sesame's "quality_mask"
    that automatically removes probes that are unreliable from all data.
- reads more IDAT raw data (run_info, probe nbeads, probe standard deviation)
    - idat.py IdatDataset has new kwargs, including bit='float16' option to cut file/memory usage in half
    by clipping max intensity at 32127 (which cuts off ~0.01% of probes)
- processing will mirror sesame more closely now, instead of minfi (to revert, use sesame=False in run_pipeline)
- adds sesame quality_mask, which auto-hides known set of sketchy probes.
- internal objects updated so that values align in every stage of processing
    (i.e. if you apply the sesame quality mask, the output files and the SampleDataContainer will exclude those probes)
- make_pipeline provides a scikit-learn style interface, as alternative to run_pipeline

## v1.3.3
- ensures methylprep output matches sesame output
- order of probes in CSVs, pickles, and SampleDataContainer doesn't match
- fixes bug where are few probes had negative meth/unmeth values because of int16 limits.
    Now it uses unsigned int16 data type and unit tests confirm no negative values appear.

## v1.3.2
- updated support for Illumina mouse array
- summarized processing warnings at end, to make tqdm progress bar cleaner

## v1.3.1
- run_pipeline() has 50% shorter processing time due to user-submitted changes
- idats can be processed while gzipped (.idat.gz) and saved this way using --no_uncompress flag
- 'download' function manages FTP connection better
- improved unit testing: download and process_series
- run_pipeline() function has two new optional parameters
    - poobah_decimals: if you want more than the default 3 decimals saved in poobah_values.pkl and _processed.csv files-, then specify a higher limit.
    - poobah_sig: default significance level for excluding probes is 0.05. You can set it to something
    else in the 0.1 to 0.000001 range, if you want.

## v1.3.0
- improved methylprep's run_pipeline() process to use less memory and avoid memory usage spikes
- files created are also smaller too (because they use float32 or int16 instead of 64 bit data)
- output files like beta_values.pkl are automatically consolidated at end of pipeline.
    batch_size will split large batches into multiple files during processing to save memory,
    but the output will be merged at the end.
- improved support for Illumina Mouse Arrays.
- methylprep.download has a new all-encompassing pipeline that will read GEO data sets and convert
    any data file type into a pickle of beta_values, whether from idats or processed matrix files.
