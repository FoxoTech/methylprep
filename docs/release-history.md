# Release History


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
