from pathlib import Path
import pandas as pd
import numpy as np
import re

def read_geo(filepath, verbose=False, debug=False):
    """Use to load preprocessed GEO data into methylcheck. Attempts to find the sample beta/M_values
    in the CSV/TXT/XLSX file and turn it into a clean dataframe, with probe ids in the index/rows.

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

    if debug=True: does nothing.
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
        print(f'ERROR: this file type (){this.suffix}) is not supported')
        return

    # next, see if betas are present of do we need to calculate them?
    test = raw.iloc[0:100]
    unmeth = False
    if 'intensit' in str(this.name).lower() or 'signal' in str(this.name).lower(): # signal intensities
        unmeth = True # need to calculate beta from unmeth/meth columns
        print('Expecting raw meth/unmeth probe data')
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
                print(f"Found probe names in `{col}` column and setting as index.")
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
                print(f'ERROR reading probes with multiline header: {bad_probe_list[:200]}')
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
            print("ERROR: Unable to parse the multi-line header in this file. If you manually edit the file headers to ensure the sample intensities unclude 'Methylated' and 'Unmethylated' in column names, it might work on retry.")
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
            vectorized_beta_func = np.vectorize(calculate_beta_value)
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
                        print(col)
                if re.match(meth_pattern, col):
                    meth_samples.append(col)
                    if debug:
                        print(col)
            if debug and len(unmeth_samples) != len(meth_samples):
                print(f"DEBUG: meth/unmeth don't match!")
            if unmeth_samples != [] and meth_samples != [] and len(unmeth_samples) == len(meth_samples):
                # next: just need to match these up. they should be same if we drop the meth/unmeth part
                if verbose:
                    print(f"{len(unmeth_samples)} Samples with Methylated/Unmethylated probes intensities found. Calculating Beta Values.")
                linked = []
                for col in unmeth_samples:
                    # test_name looks for samples that are the same, except for the 'Un' part
                    test_name = col.replace('Unmethylated','Methylated')
                    print(col, test_name, test_name in meth_samples)
                    if test_name in meth_samples:
                        linked.append([col, test_name])
                        if debug:
                            print(linked[-1])
                # Here, we calculate betas for full raw data frame
                if verbose and linked != []:
                    print(f"Found {len(linked)} paired columns that appear to be linked meth/unmeth data for samples.")
                for col_u, col_m in linked:
                    col_name = col_u.replace('Unmethylated','').replace('Signal','').strip()
                    unmeth_series = raw[col_u]
                    meth_series = raw[col_m]
                    betas = vectorized_beta_func(meth_series, unmeth_series)
                    try:
                        out_df[col_name] = betas
                        samples.append(col_name)
                    except Exception as e:
                        print('ERROR', col_name, len(betas), out_df.shape, e)
            elif unmeth:
                print(f"File appears to contain probe intensities, but the column names don't match up for samples, so can't calculate beta values.")

        if samples != [] and verbose and not unmeth:
            print(f"Found {len(samples)} samples on second pass, apparently beta values with a non-standard sample_naming convention.")
        elif samples != [] and verbose and unmeth:
            pass
        elif samples == []:
            # no samples matched, so show the columns instead
            print(f"No samples found. Here are some column names:")
            print(list(test.columns)[:20])
            return

    if index_name == None:
        print("Error: probe names not found in any columns")
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
                print('error')
                df = df.drop(columns=[col])
    if verbose:
        if num_converted > 0:
            print(f"Converted {num_converted} samples from string to float16.")
        print(f"Found {len(samples)} samples and dropped {len(raw.columns) - len(samples)} meta data columns.")
    return df
