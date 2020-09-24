from pathlib import Path
import pandas as pd
import numpy as np
import re
from collections import Counter
import logging

LOGGER = logging.getLogger(__name__)
LOGGER.setLevel( logging.INFO )

__all__ = ['read_geo', 'detect_header_pattern']

''' circular imports problems --- https://stackabuse.com/python-circular-imports/
try:
    # first: try to map the canonical version here
    import methylprep
    read_geo = methylprep.read_geo
except ImportError as error:
    # if user doesn't have methylprep installed for the canonical version of this function, import this copy below
'''

def read_geo_v1(filepath, verbose=False, debug=False):
    """Use to load preprocessed GEO data into methylcheck. Attempts to find the sample beta/M_values
    in the CSV/TXT/XLSX file and turn it into a clean dataframe, with probe ids in the index/rows.
    VERSION 1.0 (deprecated June 2020 for v3, called "read_geo")

    - reads a downloaded file, either in csv, xlsx, pickle, txt
    - looks for /d_RxxCxx patterned headings and an probe index
    - sets index in df to probes
    - sets columns to sample names
    - forces probe values to be floats, if strings/mixed
    - if filename has 'intensit' or 'signal' in it, this converts to betas and saves
      even if filename doesn't match, if columns have Methylated in them, it will convert and save
    - detect multi-line headers and adjusts dataframe columns accordingly
    - returns the usable dataframe

TODO:
    - handle files with .Signal_A and .Signal_B instead of Meth/Unmeth
    - handle processed files with sample_XX

notes:
    this makes inferences based on strings in the filename, and based on the column names.
    """
    this = Path(filepath)

    if '.csv' in this.suffixes:
        raw = pd.read_csv(this)
    elif '.xlsx' in this.suffixes:
        raw = pd.read_excel(this)
    elif '.pkl' in this.suffixes:
        raw = pd.read_pickle(this)
        return raw
    elif '.txt' in this.suffixes:
        raw = pd.read_csv(this, sep='\t')
        if raw.shape[1] == 1: # pandas doesn't handle \r\n two char line terminators, but seems to handle windows default if unspecified.
            raw = pd.read_csv(this, sep='\t', lineterminator='\r') # leaves \n in values of first column, but loads
            # lineterminator='\r\n')
            # or use codecs first to load and parse text file before dataframing...
    else:
        LOGGER.error(f'ERROR: this file type (){this.suffix}) is not supported')
        return

    # next, see if betas are present of do we need to calculate them?
    test = raw.iloc[0:100]
    unmeth = False
    if 'intensit' in str(this.name).lower() or 'signal' in str(this.name).lower(): # signal intensities
        unmeth = True # need to calculate beta from unmeth/meth columns
        LOGGER.info('Expecting raw meth/unmeth probe data')
    else:
        #meth_pattern_v1 = re.compile(r'.*[_ \.]Methylated[_ \.]', re.I)
        meth_pattern = re.compile(r'.*[_ \.]?(Un)?methylated[_ \.]?', re.I)
        meth_cols = len([col for col in test.columns if re.match(meth_pattern, col)])
        if meth_cols > 0:
            unmeth = True
            # this should work below, so that even if betas are present, it will use betas first, then fall back to meth/unmeth

    def calculate_beta_value(methylated_series, unmethylated_series, offset=100):
        """ borrowed from methylprep.processing.postprocess.py """
        methylated = max(methylated_series, 0)
        unmethylated = max(unmethylated_series, 0)
        total_intensity = methylated + unmethylated + offset
        intensity_ratio = methylated / total_intensity
        return intensity_ratio

    # look for probe names in values (of first 100 responses)
    index_name = None
    multiline_header = False
    sample_pattern = re.compile(r'\w?\d+_R\d{2}C\d{2}$') # $ ensures column ends with the regex part
    sample_pattern_loose = re.compile(r'\w?\d+_R\d{2}C\d{2}.*beta', re.I)
    probe_pattern = re.compile(r'(cg|rs|ch\.\d+\.|ch\.X\.|ch\.Y\.)\d+')
    samples = []
    for col in test.columns:
        probes = [i for i in test[col] if type(i) == str and re.match(probe_pattern,i)] #re.match('cg\d+',i)]
        if len(probes) == len(test):
            index_name = col
            if verbose:
                LOGGER.info(f"Found probe names in `{col}` column and setting as index.")
        elif len(probes)/len(test) > 0.8:
            index_name = col
            multiline_header = True
            break

        if re.match(sample_pattern, col):
            samples.append(col)

    if multiline_header: # start over with new column names
        try:
            start_index = len(test) - len(probes) - 1
            # recast without header, starting at row before first probe
            new_column_names = pd.Series(list(raw.iloc[start_index])).replace(np.nan, 'No label')
            probe_list = raw[index_name].iloc[start_index + 1:]
            probe_list = probe_list.rename(raw[index_name].iloc[start_index + 1])
            bad_probe_list = [probe for probe in probe_list if not re.match(probe_pattern, probe)] # not probe.startswith('cg')]
            if bad_probe_list != []:
                LOGGER.error(f'ERROR reading probes with multiline header: {bad_probe_list[:200]}')
                return
            raw = raw.iloc[start_index + 1:]
            raw.columns = new_column_names
            test = raw.iloc[0:100]
            samples = []
            for col in test.columns:
                if re.match(sample_pattern, col):
                    samples.append(col)
            # raw has changed.
            out_df = pd.DataFrame(index=probe_list)
        except Exception as e:
            LOGGER.error("ERROR: Unable to parse the multi-line header in this file. If you manually edit the file headers to ensure the sample intensities unclude 'Methylated' and 'Unmethylated' in column names, it might work on retry.")
            return
    else:
        out_df = pd.DataFrame(index=raw[index_name]) # only used with unmethylated data sets

    if samples == []:
        # in some cases there are multiple columns matching sample_ids, and we want the 'beta' one
        for col in test.columns:
            if re.match(sample_pattern_loose, col):
                samples.append(col)

        # or we need TWO columns per sample and we calculate 'beta'.
        if samples == [] and unmeth:
            unmeth_samples = []
            meth_samples = []
            #unmeth_pattern_v1 = re.compile(r'.*[_ \.]Unmethylated[_ \.].*', re.I)
            #meth_pattern_v1 = re.compile(r'.*[_ \.]Methylated[_ \.].*', re.I)
            unmeth_pattern = re.compile(r'.*[_ \.]?Unmethylated[_ \.]?', re.I)
            meth_pattern = re.compile(r'.*[_ \.]?(?<!Un)Methylated[_ \.]?', re.I)
            for col in test.columns:
                if re.match(unmeth_pattern, col):
                    unmeth_samples.append(col)
                    if debug:
                        LOGGER.info(col)
                if re.match(meth_pattern, col):
                    meth_samples.append(col)
                    if debug:
                        LOGGER.info(col)
            if unmeth_samples != [] and meth_samples != [] and len(unmeth_samples) == len(meth_samples):
                # next: just need to match these up. they should be same if we drop the meth/unmeth part
                if verbose:
                    LOGGER.info(f"{len(unmeth_samples)} Samples with Methylated/Unmethylated probes intensities found. Calculating Beta Values.")
                linked = []
                for col in unmeth_samples:
                    test_name = col.replace('Unmethylated','Methylated')
                    if test_name in meth_samples:
                        linked.append([col, test_name])
                # Here, we calculate betas for full raw data frame
                for col_u, col_m in linked:
                    col_name = col_u.replace('Unmethylated','').replace('Signal','').strip()
                    unmeth_series = raw[col_u]
                    meth_series = raw[col_m]
                    betas = calculate_beta_value(meth_series, unmeth_series)
                    try:
                        out_df[col_name] = betas
                        samples.append(col_name)
                    except Exception as e:
                        LOGGER.error(f"ERROR {col_name} {len(betas)} {out_df.shape} {e}")
            elif unmeth:
                LOGGER.info(f"File appears to contain probe intensities, but the column names don't match up for samples, so can't calculate beta values.")

        if samples != [] and verbose and not unmeth:
            LOGGER.info(f"Found {len(samples)} samples on second pass, apparently beta values with a non-standard sample_naming convention.")
        elif samples != [] and verbose and unmeth:
            pass
        elif samples == []:
            # no samples matched, so show the columns instead
            LOGGER.info(f"No samples found. Here are some column names:")
            LOGGER.info(list(test.columns)[:20])
            return

    if index_name == None:
        LOGGER.error("Error: probe names not found in any columns")
        return
    if unmeth and samples != [] and out_df.shape[1] > 1:
        # column names are being merged and remapped here as betas
        df = out_df # index is already set
    elif multiline_header:
        df = raw.loc[:, samples]
        df.index = probe_list
    else:
        df = raw[[index_name] + samples]
        df = df.set_index(index_name)

    # finally, force probe values to be floats
    num_converted = 0
    for col in df.columns:
        if df[col].dtype.kind != 'f' and df[col].dtype.kind == 'O':
            # convert string to float
            try:
                #df[col] = df[col].astype('float16')
                # use THIS when mix of numbers and strings
                df[col] = pd.to_numeric(df[col], errors='coerce')
                num_converted += 1
            except:
                LOGGER.error('error')
                df = df.drop(columns=[col])
    if verbose:
        if num_converted > 0:
            LOGGER.info(f"Converted {num_converted} samples from string to float16.")
        LOGGER.info(f"Found {len(samples)} samples and dropped {len(raw.columns) - len(samples)} meta data columns.")
    return df


def read_geo(filepath, verbose=False, debug=False, as_beta=True, column_pattern=None, test_only=False, rename_probe_column=True, decimals=3):
    """Use to load preprocessed GEO data into methylcheck. Attempts to find the sample beta/M_values
    in the CSV/TXT/XLSX file and turn it into a clean dataframe, with probe ids in the index/rows.
    Version 3 (introduced June 2020)

    - reads a downloaded file, either in csv, xlsx, pickle, txt
    - looks for /d_RxxCxx patterned headings and an probe index
    - sets index in df to probes
    - sets columns to sample names
    - forces probe values to be floats, if strings/mixed
    - if filename has 'intensit' or 'signal' in it, this converts to betas and saves
      even if filename doesn't match, if columns have Methylated in them, it will convert and save
    - detect multi-line headers and adjusts dataframe columns accordingly
    - returns the usable dataframe

    as_beta == True -- converts meth/unmeth into a df of sample betas.
    column_pattern=None (Sample21 | Sample_21 | Sample 21) -- some string of characters that precedes the number part of each sample in the columns of the file to be ingested.

FIXED:
    [x] handle files with .Signal_A and .Signal_B instead of Meth/Unmeth
    [x] BUG: can't parse matrix_... files if uses underscores instead of spaces around sample numbers, or where sampleXXX has no separator.
    [x] handle processed files with sample_XX
    [x] returns IlmnID as index/probe column, unless 'rename_probe_column' == False
    [x] pass in sample_column names from header parser so that logic is in one place
        (makes the output much larger, so add kwarg to exclude this)
    [x] demicals (default 3) -- round all probe beta/intensity/p values returned to this number of decimal places.
    [x] bug: can only recognize beta samples if 'sample' in column name, or sentrix_id pattern matches columns.
        need to expand this to handle arbitrary sample naming styles (limited to one column per sample patterns)
TODO:
    [-] BUG: meth_unmeth_pval works `as_beta` but not returning full data yet
    [-] multiline header not working with all files yet.


notes:
    this makes inferences based on strings in the filename, and based on the column names.
    """

    #1. load and read the file, regardless of type and CSV delimiter choice.
    this = Path(filepath)
    kwargs = {'nrows':200} if test_only else {}

    raw = pd_load(filepath, **kwargs)
    if not isinstance(raw, pd.DataFrame):
        LOGGER.error(f"Did not detect a file: {type(raw)} aborting")
        return
    if debug: LOGGER.info(f"{filepath} loaded.")

    def calculate_beta_value(methylated_series, unmethylated_series, offset=100):
        """ borrowed from methylprep.processing.postprocess.py """
        methylated = np.clip(methylated_series, 0, None)
        unmethylated = np.clip(unmethylated_series, 0, None)

        total_intensity = methylated + unmethylated + offset
        intensity_ratio = methylated / total_intensity
        return intensity_ratio

    #if as_beta:
    def convert_meth_to_beta(raw, meta, rename_probe_column=True, decimals=3, verbose=False):
        #if (meta['column_pattern'] == 'meth_unmeth_pval' and
        #    meta['sample_columns_meth_unmeth_count'] > 0 and
        #    'sample_numbers_range' != []):
        # beta = methylated / total intensity (meth + unmeth + 100)
        # use sample_column_names and sample_numbers_list to find each set of columns
        out_df = pd.DataFrame(index=raw[meta['probe_column_name']])
        probe_name_msg = meta['probe_column_name']
        if rename_probe_column:
            out_df.index.name = 'IlmnID' # standard for methyl-suite, though ID_REF is common in GEO.
            probe_name_msg = f"{meta['probe_column_name']} --> IlmnID"

        if not meta.get('sample_names') and not meta['sample_names'].get('meth'):
            # when "sample" not in columns.
            LOGGER.error("ERROR: missing sample_names")
            pass
        if meta['sample_numbers_range'] and meta['sample_numbers_range'][1] <= meta['total_samples']:
            # if non-sequential sample numbers, cannot rely on this range.
            sample_min, sample_max = meta['sample_numbers_range']
        else:
            sample_min = 1
            sample_max = meta['total_samples']

        # JSON returns 1 to N (as human names) but index is 0 to N.
        for sample_number in range(sample_min -1, sample_max):
            # here need to get the corresponding parts of each sample triplet of columns.
            try:
                col_m = meta['sample_names']['meth'][sample_number]
                col_u = meta['sample_names']['unmeth'][sample_number]
                #col_pval = meta['sample_columns']['pval']
            except Exception as e:
                LOGGER.error(f"ERROR {e} {meta['sample_names']['meth'][sample_number]} {list(meta['sample_names']['unmeth'].columns)[sample_number]}")
                continue

            unmeth_series = raw[col_u]
            meth_series = raw[col_m]
            betas = calculate_beta_value(meth_series, unmeth_series)

            # try to lookup and retain any unique part of column names here
            col_name = f"Sample_{sample_number + 1}"
            if isinstance(meta['sample_names_stems'], list):
                # assume order of sample_names_stems matches order of sample names
                try:
                    this_stem = meta['sample_names_stems'][sample_number]
                    if this_stem not in out_df.columns:
                        col_name = this_stem
                except IndexError:
                    if debug: LOGGER.error(f"ERROR: unable to assign Sample {sample_number} using original column stem.")

                """ This code was trying to match samples up using regex/difflib but it was unreliable.
                this_stem = [stem for stem in meta['sample_names_stems'] if stem in col_m]
                if len(this_stem) > 0:
                    if len(this_stem) > 1:
                        # difflib to separatet Sample 1 vs Sample 11
                        best_match = difflib.get_close_matches(col_m, this_stem, 1)
                        if best_match != [] and best_match[0] not in out_df.columns:
                            col_name = best_match[0]
                        elif best_match != []: # ensure unique in out_df, even if not a perfect transfer of labels.
                            col_name = f"{best_match[0]}_{sample_number}"
                        else:
                            if debug: LOGGER.info(f"WARNING: multiple similar sample names detected but none were a close match: {col_m} : {this_stem}")
                            col_name = f"Sample_{sample_number}"
                    elif len(this_stem) == 1 and this_stem[0] not in out_df.columns:
                        # only one match, and is unique.
                        col_name = this_stem[0]
                    else: # only one match, but already in out_df, so can't reuse.
                        col_name = f"{this_stem[0]}_Sample_{sample_number}"
                else:
                    col_name = f"Sample_{sample_number}"

            else:
                col_name = f"Sample_{sample_number}"
                """

            try:
                out_df[col_name] = betas
                sample_number += 1
            except Exception as e:
                LOGGER.error(f"ERROR {col_name} {len(betas)} {out_df.shape} {e}")
        out_df = out_df.round(decimals)

        beta_value_range = True if all([all(out_df[col_name].between(0,1)) == True for col_name in out_df.columns]) else False
        intensity_value_range = True if all([all(out_df[col_name].between(0,1000000)) == True for col_name in out_df.columns]) else False
        try:
            value_mean = round(sum(out_df.mean(axis=0))/len(out_df.columns),2)
        except:
            value_mean = 0
        if verbose:
            LOGGER.info(f"Returning {len(out_df.columns)} samples (mean: {value_mean}). Structure appears to be {meta['columns_per_sample']} columns per sample; raw data column pattern was '{meta['column_pattern']}'; probes in rows; and {probe_name_msg} as the probe names.")
        return out_df

    #2. test file structure
    # next, see if betas are present of do we need to calculate them?
    test = raw.iloc[0:100]
    meta = detect_header_pattern(test, filepath, return_sample_column_names=True)
    if meta['multiline_header'] and meta['multiline_header_rows'] > 0:
        if verbose: print("Reloading raw data, excluding header.")
        old_non_blank_col_count = len(raw.loc[:, ~raw.columns.str.contains('^Unnamed:')].columns) #~ before test counts non-blank columns
        kwargs['skiprows'] = meta['multiline_header_rows']
        raw = pd_load(filepath, **kwargs)
        new_non_blank_col_count = len(raw.loc[:, ~raw.columns.str.contains('^Unnamed:')].columns)
        if verbose: print(f"After ignoring multiline header: {old_non_blank_col_count} old columns are now {new_non_blank_col_count} new named columns.")
    if debug:
        LOGGER.info(f"file shape: {raw.shape}")
        concise = meta.copy()
        concise.pop('sample_names')
        from pprint import pprint
        pprint(concise)
        del concise

    #3. use header_meta to parse and return data as dataframe
    """
{'all_sample_columns': True,
 'column_pattern': 'meth_unmeth_pval',
 'columns_per_sample': 3,
 'fraction_sample pval meth unmeth signal intensity': (1.0, 0.33, 0.33, 0.33, 0.66, 0.0),
 'has_beta_values': False,
 'has_meth_unmeth_values': True,
 'has_p_values': True,
 'multiline_header': False,
 'multiline_rows': 0,
 'one_sample_beta': False,
 'probe_column_name': 'ID_REF',
 'sample_columns_meth_unmeth_count': 222,
 'sample_numbers_range': [1, 111], # or None
 'sequential_numbers': True,
 'total_samples': 111}"""

    if (meta['column_pattern'] == 'sample_beta_sentrix_id'):
        sample_columns = meta['sample_names']
        out_df = pd.DataFrame(data=raw[sample_columns])
        out_df.index=raw[meta['probe_column_name']]
        out_df = out_df.round(decimals) #{sample:demicals for sample in sample_columns})

        probe_name_msg = meta['probe_column_name']
        if rename_probe_column:
            out_df.index.name = 'IlmnID' # standard for methyl-suite, though ID_REF is common in GEO.
            probe_name_msg = f"{meta['probe_column_name']} --> IlmnID"

        beta_value_range = True if all([all(out_df[col_name].between(0,1)) == True for col_name in out_df.columns]) else False
        intensity_value_range = True if all([all(out_df[col_name].between(0,1000000)) == True for col_name in out_df.columns]) else False
        try:
            value_mean = round(sum(out_df.mean(axis=0))/len(out_df.columns),2)
        except:
            value_mean = 0
        if verbose and meta['columns_per_sample'] != 1 and beta_value_range:
            LOGGER.info(f"Returning raw data. Structure appears to be {meta['columns_per_sample']} columns per sample; numbered samples: {numbered_samples}; column_pattern: {column_pattern}; probes in rows; and {probe_name_msg} as the probe names.")
        elif verbose and beta_value_range:
            LOGGER.info(f"Returning raw data. Appears to be sample beta values in columns (mean: {value_mean}), probes in rows, and {probe_name_msg} as the probe names.")
        elif verbose and not beta_value_range and intensity_value_range and value_mean > 10:
            LOGGER.info(f"Returning raw data. Appears to be sample fluorescence intensity values in columns (mean: {int(value_mean)}), probes in rows, and {probe_name_msg} as the probe names.")
        elif verbose and not beta_value_range and intensity_value_range and value_mean > 10:
            LOGGER.info(f"Returning raw data of UNKNOWN type. Not beta values or fluorescence intensity values in columns. Mean probe value was {value_mean}. Probes are in rows, and {probe_name_msg} as the probe names.")
        return out_df

    elif ((meta['column_pattern'] == 'sample_beta_numbered') or

        (meta['all_sample_columns'] and meta['has_beta_values']) or
        (meta['sequential_numbers'] and meta['sample_columns_meth_unmeth_count'] == 0) or
        (meta['sample_numbers_range'] and meta['has_p_values'] is False) or
        (meta['all_sample_columns'] and meta['sample_numbers_range'] and meta['columns_per_sample'] == 1)
       ):
        sample_columns = [column for column in test.columns if 'sample' in column.lower()]
        out_df = pd.DataFrame(data=raw[sample_columns])
        out_df.index=raw[meta['probe_column_name']]
        out_df = out_df.round(decimals) #{sample:demicals for sample in sample_columns})
        probe_name_msg = meta['probe_column_name']
        if rename_probe_column:
            out_df.index.name = 'IlmnID' # standard for methyl-suite, though ID_REF is common in GEO.
            probe_name_msg = f"{meta['probe_column_name']} --> IlmnID"
        #if debug: LOGGER.info(f"DEBUG: out_df.columns: {out_df.columns} --- out_df.index {out_df.index.name} out_df.shape {out_df.shape}")

        beta_value_range = True if all([all(out_df[col_name].between(0,1)) == True for col_name in out_df.columns]) else False
        intensity_value_range = True if all([all(out_df[col_name].between(0,1000000)) == True for col_name in out_df.columns]) else False
        try:
            value_mean = round(sum(out_df.mean(axis=0))/len(out_df.columns),2)
        except:
            value_mean = 0
        if verbose and meta['columns_per_sample'] != 1 and beta_value_range:
            LOGGER.info(f"Returning raw data. Structure appears to be {meta['columns_per_sample']} columns per sample; {len(out_df.columns)} numbered samples; column_pattern: {column_pattern}; probes in rows; and {probe_name_msg} as the probe names.")
        elif verbose and beta_value_range:
            LOGGER.info(f"Returning raw data. Appears to contain {len(out_df.columns)} numbered samples; beta values in columns (mean: {value_mean}), probes in rows, and {probe_name_msg} as the probe names.")
        elif verbose and not beta_value_range and intensity_value_range and value_mean > 10:
            LOGGER.info(f"Returning raw data. Appears to be sample fluorescence intensity values in columns ({len(out_df.columns)} samples with mean: {int(value_mean)}), probes in rows, and {probe_name_msg} as the probe names.")
        elif verbose and not beta_value_range and intensity_value_range and value_mean > 10:
            LOGGER.info(f"Returning raw data of UNKNOWN type. Not beta values or fluorescence intensity values in columns. {len(out_df.columns)} samples with a mean probe value of {value_mean}. Probes are in rows, and {probe_name_msg} as the probe names.")

        return out_df

    elif meta['column_pattern'] == 'beta_intensity_pval':
        if verbose: LOGGER.info("returning raw data without processing. Column pattern was 'beta_intensity_pval'.")
        return raw
    elif meta['column_pattern'] == 'beta_pval_a_b':
        if as_beta and meta['sample_names'] != None:
            # only returning the beta part. Signal_A and Signal_B are the meth/unmeth parts.
            out_df = pd.DataFrame(data=raw[meta['sample_names']])
            out_df.index=raw[meta['probe_column_name']]
            out_df = out_df.round(decimals) #{sample:demicals for sample in sample_columns})
            probe_name_msg = meta['probe_column_name']
            if rename_probe_column:
                out_df.index.name = 'IlmnID' # standard for methyl-suite, though ID_REF is common in GEO.
                probe_name_msg = f"{meta['probe_column_name']} --> IlmnID"
            try:
                value_mean = round(sum(out_df.mean(axis=0))/len(out_df.columns),2)
            except:
                value_mean = 0
            if verbose:
                LOGGER.info(f"Returning beta values for {len(out_df.columns)} samples in columns (mean: {value_mean}), probes in rows, and {probe_name_msg} as the probe names.")
            return out_df
        else:
            if verbose: LOGGER.info("returning raw data without processing. Column pattern was 'beta_pval_a_b'.")
            return raw

    if meta['column_pattern'] == 'meth_unmeth':
        LOGGER.info("*** column_pattern was 'meth_unmeth' this has litterally never happened before. ***")
    if meta['column_pattern'] in ('meth_unmeth','meth_unmeth_pval') and as_beta:
        if debug: LOGGER.info("Converting meth and unmeth intensities to beta values.")
        return convert_meth_to_beta(raw, meta, rename_probe_column=rename_probe_column, verbose=verbose)
    if meta['column_pattern'] in ('meth_unmeth','meth_unmeth_pval') and not as_beta:
        if verbose: LOGGER.info(f"Returning raw data without processing. Column pattern was {meta['column_pattern']}.")
        return raw

    if meta['column_pattern'] is None:
        if debug: LOGGER.info("Returning raw data without processing. No file header column pattern was detected.")
        return raw
    raise Exception("Unable to identify file structure")


def detect_header_pattern(test, filename, return_sample_column_names=False):
    """test is a dataframe with first 100 rows of the data set, and all columns.
    makes all the assumptions easier to read in one place.

    betas
    non-normalized
    matrix_processed
    matrix_signal
    series_matrix
    methylated_signal_intensities and unmethylated_signal_intensities
    _family

    TODO: GSM12345-tbl-1.txt type files (in _family.tar.gz packages) are possible, but needs more work.
    TODO: combining two files with meth/unmeth values

    - numbered samples handled differently from sample_ids in columns
    - won't detect columns with no separators in strings
    """
    if test.shape[0] != 100:
        raise ValueError("test dataset must be exactly 100 rows")
    if test.shape[1] == 1:
        raise ValueError("this dataset has only one sample. it is likely that the columns were not parsed correctly.")

    seps = [' ', '_', '.', '-'] # for parsing columns. also try without any separators
    index_names = ['IlmnID', 'ID_REF', 'illumina_id']

    # sample patterns
    sample_pattern = re.compile(r'\w?\d+_R\d{2}C\d{2}$') # $ ensures column ends with the regex part
    sample_pattern_loose = re.compile(r'\w?\d+_R\d{2}C\d{2}.*beta', re.I)
    samplelike_pattern = re.compile(r'.*(?:\w?\d+_R\d{2}C\d{2}|sample).*', re.I)
    probe_pattern = re.compile(r'(cg|rs|ch\.\d+\.|ch\.X\.|ch\.Y\.)\d+')
    #pval_pattern = re.compile(r'(.*)(?:\bPval\b|\.Pval|_Pval_|-Pval|Pval|\bDetection\bPval|_Detection_Pval|\._Detection\.Pval||\._Detection\bPval\b).*', re.I)
    pval_pattern = re.compile(r'(.*)(?:\bPval\b|\.Pval|_Pval_|-Pval|Pval).*', re.I)
    meth_pattern = re.compile(r'(.*)(?:\bmeth|\.meth|_meth|-meth|(?<!un)meth\w+).*', re.I)
    unmeth_pattern = re.compile(r'(.*)(?:\bunmeth|\.unmeth|_unmeth|-unmeth|unmeth\w+).*', re.I)
    intensity_pattern = re.compile(r'(.*)(?:\bintensity\b|\.intensity|_intensity|-intensity|intensity).*', re.I)
    signal_pattern = re.compile(r'(.*)(?:\bsignal\b|\.signal|_signal|-signal|signal).*', re.I)
    betalike_pattern = re.compile(r'(.*)(?:\Wbeta\W|\Wavg\Wbeta\W|\Wavg_beta|_avg_beta|\Wavgbeta).*', re.I)
    # use the last part to find any sample identifiers; then re.sub() replace them. then pass rest into the beta_converter function.
    residual_pattern = re.compile(r'(.*)(?:[\b\._-]Pval[\b\._-]|[\b\._-]Detection[\b\._-]Pval|[\b\._-]meth|(?<!un)meth\w+|\bunmeth|\.unmeth|_unmeth|-unmeth|unmeth\w+|\bintensity\b|\.intensity|_intensity|-intensity|intensity).*', re.I)
    file_gsm_txt_in_family_tar = re.compile(r'gsm\d+-tbl-1.txt') # all matches are lower() -- not used yet
    #meth_unmeth_pattern = re.compile(r'.*[_ \.]?(Un)?methylated[_ \.]?', re.I)

    # filename patterns
    this = str(Path(filename).name).lower()
    has_beta_values = True if ('matrix' in this and 'signal' not in this) else False #exceptions to this rule found.
    # for meth_unmeth, need to calculate beta from columns
    has_meth_unmeth_values = True if 'intensit' in this or 'signal' in this else False
    # no rule for this yet
    has_p_values = None
    # GSM1280914-tbl-1.txt contains probe names, one sample beta val column, and an optional p-val column with no header.
    one_sample_beta = True if this.startswith('gsm') and re.search(this, file_gsm_txt_in_family_tar) else False

    # column parts
    sample_columns = [column for column in test.columns if 'sample' in column.lower()]
    # -- next: test first 100 rows and calculate the faction of each column that matches the probe_pattern. list of floats.
    fraction_probelike = [len([row for row in test[column] if isinstance(row, str) and re.match(probe_pattern,row)])/len(test) for column in test.columns]
    fraction_pvallike = round(sum([(True if re.search(pval_pattern,column) else False) for column in test.columns])/len(list(test.columns)),2)
    fraction_meth = round(sum([(True if re.search(meth_pattern,column) else False) for column in test.columns])/len(list(test.columns)),2)
    fraction_unmeth = round(sum([(True if re.search(unmeth_pattern,column) else False) for column in test.columns])/len(list(test.columns)),2)
    fraction_intensity = round(sum([(True if re.search(intensity_pattern,column) else False) for column in test.columns])/len(list(test.columns)),2)
    fraction_signal = round(sum([(True if re.search(signal_pattern,column) else False) for column in test.columns])/len(list(test.columns)),2)
    fraction_beta = round(sum([(True if re.search(betalike_pattern,column) else False) for column in test.columns])/len(list(test.columns)),2)
    fraction_samplelike = round(sum([(True if re.search(samplelike_pattern,column) else False) for column in test.columns])/len(list(test.columns)),2)
    probe_column_position = fraction_probelike.index( max(fraction_probelike) )
    probe_column_name = list(test.columns)[probe_column_position] # used for index_name
    # grab residual column name parts to add to output column names (tricky)
    sample_column_residuals = [re.search(residual_pattern, column).groups()[0] for column in test.columns if re.search(residual_pattern, column)]
    # a trick to deal with 3-column samples that have similar names; retain any repeated names, if repeats exist.
    if sample_column_residuals != [] and Counter(sample_column_residuals).most_common(1)[0][1] > 1:
        sample_column_residuals = [column.strip() for column,freq in Counter(sample_column_residuals).most_common() if freq > 1]

    ## DETECT column repeating structure.
    # how many columns PER SAMPLE in structure?
    # meth or (un)meth appear anywhere in column name along with 'sample' anywhere
    sample_columns_meth_unmeth_count = sum([(True if (re.search(meth_pattern,column) or re.search(unmeth_pattern,column)) else False) for column in test.columns])
    #sum([('meth' in column.lower().replace('sample', '')) for column in sample_columns])

    # extract any "first" numbers found in column parts, where parts are strings separated any any logical separators.
    sample_id_columns = []
    sample_numbers_list = []
    sample_numbers_range = None
    sample_count = 0
    sequential_numbers = None # starting from 1 to xx
    avg_sample_repeats = 1
    for sep in seps:
        try:
            sample_numbers_list = [[part for part in column.split(sep) if re.match(r'\d+', part)][0] for column in sample_columns]
            if len(sample_numbers_list) > 0:
                sorted_sample_numbers_list = sorted([int(j) for j in list(set(sample_numbers_list))])
                sample_repeats = list(Counter([int(j) for j in list(sample_numbers_list)]).values())
                avg_sample_repeats = int(round(sum(sample_repeats)/len(sample_repeats)))
                #sequential_numbers = True if all([(i+1 == j) for i,j in list(enumerate(sorted_sample_numbers_list))]) else False
                sequential_numbers = True if all([(i+1 == j) for i,j in list(enumerate(sorted_sample_numbers_list))]) else False
                sample_count = len(sample_numbers_list)
                sample_numbers_range = [min(sorted_sample_numbers_list), max(sorted_sample_numbers_list)]
                break
        except Exception as e:
            continue

    # detect if some part of columns are named like sample_ids
    if sample_numbers_list == []:
        for sep in seps:
            sample_id_columns = [column for column in test.columns if any([re.search(sample_pattern, part) for part in column.split(sep)])]
        sample_count = len(sample_id_columns)

    # tests / data attributes: pval, multiline, probes, sample/nonsample columns
    if max(fraction_probelike) < 0.75:
        LOGGER.warning(f"WARNING: Unable to identify the column with probe names ({max(fraction_probelike)})")
    multiline_rows = 0
    if 0.75 <= max(fraction_probelike) < 1.00:
        multiline_rows = int(round(100*max(fraction_probelike)))
        header_rows = 100 - multiline_rows
        LOGGER.info(f"Multiline header detected with {header_rows} rows.")
    multiline_header = True if multiline_rows > 0 else False
    all_sample_columns = all([('sample' in column.lower() or 'ID_REF' in column) for column in list(set(test.columns))])
    has_p_values = True if fraction_pvallike > 0 else False

    # overall pattern logic
    column_pattern = None
    if (sample_numbers_list != [] or sample_id_columns != []):
        columns_per_sample = int(round(sample_count/len(set(sample_numbers_list or sample_id_columns))))
    elif (0.31 <= fraction_meth <= 0.34 and 0.31 <= fraction_unmeth <= 0.34 and 0.31 <= fraction_pvallike <= 0.34):
        columns_per_sample = 3
    else:
        columns_per_sample = 1

    if sample_numbers_list != [] and columns_per_sample != avg_sample_repeats:
        LOGGER.warning('WARNING: inconsistent columns per sample')
    if (sample_count > 0 and
        sample_numbers_list != [] and
        sample_columns_meth_unmeth_count == 0 and
        0.9 < columns_per_sample < 1.1):
        column_pattern = 'sample_beta_numbered'
        total_samples = len(sample_columns)

    elif (sample_count > 0 and
          sample_id_columns != [] and
          sample_columns_meth_unmeth_count == 0 and
          0.9 < columns_per_sample < 1.1):
        column_pattern = 'sample_beta_sentrix_id'
        total_samples = len(sample_columns)

    elif (sample_columns_meth_unmeth_count > 0 and
          sample_count > 0 and
          sample_columns_meth_unmeth_count > 0 and
          1.8 < len(sample_numbers_list) / sample_columns_meth_unmeth_count < 2.2):
        column_pattern = 'meth_unmeth'
        total_samples = int(round(sample_columns_meth_unmeth_count/2))

    elif (sample_columns_meth_unmeth_count > 0 and
          sample_count > 0 and
          sample_columns_meth_unmeth_count > 0 and
          1.4 < len(sample_numbers_list)/sample_columns_meth_unmeth_count < 1.6):
        column_pattern = 'meth_unmeth_pval'
        total_samples = int(round(sample_columns_meth_unmeth_count/2))

    elif (0.03 <= fraction_pvallike <= 0.36 and
          0.03 <= fraction_meth <= 0.36 and
          0.03 <= fraction_unmeth <= 0.36):
        column_pattern = 'meth_unmeth_pval'
        total_samples = int(round(len(test.columns)*fraction_pvallike))

    elif (0.30 <= fraction_pvallike <= 0.36 and sample_columns_meth_unmeth_count == 0):
        column_pattern = 'beta_intensity_pval'
        total_samples = int(round(len(test.columns)*fraction_pvallike))

    elif (0.45 <= fraction_signal <= 0.55 and
          0.22 <= fraction_pvallike <= 0.26):
        column_pattern = 'beta_pval_a_b' #SignalA ... SignalB
        total_samples = int(round(len(test.columns)*fraction_pvallike))

    else:
        column_pattern = None
        total_samples = 0

    # return_sample_column_names
    # -- depends on structure.
    # -- sample_columns will be a flat list for "sample"-like ones.
    # -- Or a dict of lists if meth_unmeth_pval.
    if return_sample_column_names:
        if fraction_meth > 0 and fraction_unmeth > 0 and column_pattern in ('meth_unmeth_pval', 'meth_unmeth'):
                sample_columns = {
                                  'meth': [column for column in test.columns if re.search(meth_pattern,column)],
                                  'unmeth': [column for column in test.columns if re.search(unmeth_pattern,column)],
                                  'pval': [column for column in test.columns if re.search(pval_pattern,column)],
                                 }
        elif sample_columns == [] and sample_id_columns != []:
            sample_columns = sample_id_columns
        elif sample_columns == [] and column_pattern == 'beta_intensity_pval':
            sample_columns = test.columns
        elif column_pattern == 'beta_pval_a_b': # and sample_columns == []
            if 0.2 <= fraction_beta <= 0.28:
                sample_columns = [column for column in test.columns if re.search(betalike_pattern,column)]
            else:
                LOGGER.warning("WARNING: test columns for beta_pval_a_b don't look right.")
                sample_columns = test.columns
    else:
        sample_columns = None
        sample_columns_residuals = None

    return {'has_beta_values': has_beta_values, #  or (beta_assumed and columns_per_sample == 1),
            'has_meth_unmeth_values': has_meth_unmeth_values,
            'has_p_values': has_p_values,
            'one_sample_beta': one_sample_beta,
            'sample_columns_meth_unmeth_count': sample_columns_meth_unmeth_count,
            'all_sample_columns': all_sample_columns,
            'multiline_header': multiline_header,
            'multiline_header_rows': header_rows,
            'column_pattern': column_pattern,
            'columns_per_sample': columns_per_sample,
            'sequential_numbers': sequential_numbers,
            'sample_numbers_range': sample_numbers_range,
            'total_samples': total_samples,
            'fraction_sample pval meth unmeth signal intensity beta': (fraction_samplelike, fraction_pvallike, fraction_meth, fraction_unmeth, fraction_signal, fraction_intensity, fraction_beta),
            'probe_column_name': probe_column_name,
            'sample_names': sample_columns,
            'sample_names_stems': sample_column_residuals, # unique parts of sample names
           }

def pd_load(filepath, **kwargs):
    """ helper function that reliably loads any GEO file, by testing for the separator that gives the most columns """
    this = Path(filepath)

    if this.suffix not in ('.xlsx', '.pkl'):
        # first, check that we're getting the max cols
        test_csv = pd.read_csv(this, nrows=100, skiprows=kwargs.get('skiprows',0))
        test_t = pd.read_csv(this, sep='\t', nrows=100, skiprows=kwargs.get('skiprows',0))
        test_space = pd.read_csv(this, sep=r',\s+', quoting=csv.QUOTE_ALL, engine='python', skiprows=kwargs.get('skiprows',0))

        params = [
            {'method':'auto', 'cols': test_csv.shape[1], 'kwargs': {}},
            {'method':'tab', 'cols': test_t.shape[1], 'kwargs': {'sep':'\t'}},
            {'method':'quoted', 'cols': test_space.shape[1], 'kwargs': {'sep':r',\s+', 'quoting':csv.QUOTE_ALL, 'engine': 'python'}},
        ]
        best_params = sorted([(parts['method'], parts['cols'], parts['kwargs']) for parts in params], key= lambda param_tuple: param_tuple[1], reverse=True)
        kwargs.update(best_params[0][2])

    if '.csv' in this.suffixes:
        raw = pd.read_csv(this, **kwargs)
    elif '.xlsx' in this.suffixes:
        raw = pd.read_excel(this, **kwargs)
    elif '.pkl' in this.suffixes:
        raw = pd.read_pickle(this, **kwargs)
        #return raw
    elif '.txt' in this.suffixes:
        try:
            raw = pd.read_csv(this, **kwargs) # includes '\t' if testing worked best with tabs
        except ParserError as e:
            if debug: print(f"{e}: look at file and deleted few header rows first.")
            return
        if raw.shape[1] == 1: # pandas doesn't handle \r\n two char line terminators, but seems to handle windows default if unspecified.
            raw = pd.read_csv(this, sep='\t', lineterminator='\r', **kwargs) # leaves \n in values of first column, but loads
            # lineterminator='\r\n')
            # or use codecs first to load and parse text file before dataframing...
    else:
        print(f'ERROR: this file type ({this.suffix}) is not supported')
        return
    return raw
