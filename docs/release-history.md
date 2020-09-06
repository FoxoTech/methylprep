# Release History

## v1.3.0
- improved methylprep's run_pipeline() process to use less memory and avoid memory usage spikes
- files created are also smaller too (because they use float32 or int16 instead of 64 bit data)
- output files like beta_values.pkl are automatically consolidated at end of pipeline.
    batch_size will split large batches into multiple files during processing to save memory,
    but the output will be merged at the end.
- improved support for Illumina Mouse Arrays.
- methylprep.download has a new all-encompassing pipeline that will read GEO data sets and convert
    any data file type into a pickle of beta_values, whether from idats or processed matrix files.
    
