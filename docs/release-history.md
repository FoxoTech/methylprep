# Release History

## v1.7.1
- fixed edge case bug where infer_type_I_probes fails because there are no type-I-red probes
- fixed 16-bit overflow bug reading IDATs with fluorescence > 32766 units created negative values; Reason: pandas version >
    1.2.5 can't coerce the uint16 numpy array into float32, but casting it as int64 first resolves the negative fluorescence errors. 

## v1.7.0
- added support for `parquet` file formats (as an alternative to pickle that is readable
    by other languages)
  - run_pipeline (beta, M, noob, raw_meth/unmeth, samplesheet_data_frame) as parquet
- minor fixes to GEO download processing

## v1.6.2
- Minor bug fixes

## v1.6.1
- samplesheet: ensures csv and meta_data pickle match
- better error message when multiple samplesheet csvs are detected, and more stringent detection parameters
- updated CI/CD with github actions, faster testing, dropped flaky unit tests
- updated documentation

## v1.6.0
- qualityMask: All versions of 1.5.x used a different convention for probes to drop or keep. In version 1.6.0 and above,
methylprep reverts to the previous convention used in v1.4.6 and below: For the `quality_mask` column in the processed output
CSV and the in-memory SampleDataContainer data_frame, 1.0 means keep the probe (good) and 0.0 means drop. (This reverses
a previous convention and fixes bugs where the 0.0 were sometimes NaN.)

## v1.5.9
- Fixed bug: the SampleDataContainer returned by run_pipeline did not have 1.0 and 0 for quality_mask. Matches the CSV export now.

## v1.5.8
- Fixed bug in CSV export; quality_mask (probes to be excluded from pickles) were NaN instead of 1.0

## v1.5.7
- Merely a maintenance release to deal with some dependency conflicts in the python science stack.
    Pandas version 1.3+ now supported.

## v1.5.6
- completely rewritten and updated documentation, including a more extensive tutorial
- updated all manifests to include probe-to-locus mapping for two genome builds
  - Uses OLD_ in front of 4 genomic columns to denote the older genome build in each respective array/organism
  - mouse manifest has mm39 (newer) and mm10 (older) genome assemblies.
  - human hg38 is a corrected and improved version of (OLD) hg19.
- 27k array is no longer supported. (It turns out it was never supported, but we just never unit-tested it. And
    because there are no control probes and no type-II probes on the 27k first generation array, it would be a lot
    of work to support it, and nobody has ever asked about it.)
- removed read_geo and detect_header_pattern from methylprep, because it is maintained in methylcheck and imported
    into methylprep now.
- new beta_bake CLI option `-z / --compress` will put all downloaded files into a zipfile. This used to be the default
    behavior, and now it is optional.
- fixed minor bug where malformed custom user-supplied manifest files will return a useful error message.
- better processing error detection with a `check_for_probe_loss` functions that warns if probes are dropped

## v1.5.5
- Fixed Reading IDATs progress bar in 'process'
- `download` now uses HTTPS requests to fetch GSExxxxxx_RAW.TAR data files instead of FTP, because this is way more reliable and avoids errors.
- Better error messages when trying to download and process GEO datasets containing multiple array types. The process splits these samples into separate folders based on the meta data, and tries to run each array type separately.
- Improved documentation everywhere
- Improved support for GEO series_matrix.txt.gz files and _family.xml -tbl-1.txt files
- Fixed bug where quality_mask was removing SNP (rs) probes from CSV or SampleDataContainer output.
- Removed detect_header_pattern and read_geo from methylprep. These are maintained in methylcheck and
  not called anywhere. Instead, methylprep tries to import methylcheck in the Download functions that need
  to read GEO data files. But if a user has not installed methylcheck, they will only see an error if they
  run specific commands in Download that require it.

## v1.5.2
- Bug fix: added 'Extended_Type' into control_probes.pkl output. Required by methylcheck.plot_controls().
- Minor bug fixes and improved unit test coverage.
- Fixed bug where `process --minfi` was not working with `--all`. Added more test coverage for CLI.
- updated read_geo to handle more edge cases
- deprecated some never-used functions.
  - instead of methylprep.files.idat.RunInfo use IdatDataset( verbose=True )


## v1.5.0, v1.5.1
- MAJOR refactor/overhaul of all the internal classes. This was necessary to fully support the mouse array.
- new SigSet class object that mirror's sesame's SigSet and SigDF object.
- Combines idats, manifest, and sample sheet into one object that is inherited by SampleDataContainer
- RawDataset, MethylationDataset, ProbeSubtype all deprecated and replaced by SigSet
- SampleDataContainer class is now basically the SigSet plus all pipeline processing settings
- new mouse manifest covers all probes and matches sesame's output
- Processing will work even if a batch of IDATs have differing probe counts for same array_type, though those
differing probes in question may not be saved.
- unit tests confirm that methylprep, sesame, and minfi beta values output match to within 1% of each other now. Note that the intermediate stages of processing (after NOOB and after DYE) do not match
with sesame in this version. Can be +/- 100 intensity units, likely due to differences in order of
steps and/or oob/mask probes used.

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

## Older versions exist on pypi, but no changelog
