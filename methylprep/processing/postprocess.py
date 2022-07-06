# Lib
import numpy as np
import pandas as pd
import os
import pickle
from pathlib import Path
import logging
# app
from ..utils import is_file_like
#from ..utils.progress_bar import * # context tqdm

os.environ['NUMEXPR_MAX_THREADS'] = "8" # suppresses warning


__all__ = ['calculate_beta_value', 'calculate_m_value', 'consolidate_values_for_sheet', 'consolidate_mouse_probes']

LOGGER = logging.getLogger(__name__)

def calculate_beta_value(methylated_noob, unmethylated_noob, offset=100):
    """ the ratio of (methylated_intensity / total_intensity)
    where total_intensity is (meth + unmeth + 100) -- to give a score in range of 0 to 1.0.
    minfi offset is 100 and sesame (default) offset is zero."""
    methylated = np.clip(methylated_noob, 1, None)
    unmethylated = np.clip(unmethylated_noob, 1, None)

    total_intensity = methylated + unmethylated + offset
    with np.errstate(all='raise'):
        intensity_ratio = np.true_divide(methylated, total_intensity)
    return intensity_ratio


def calculate_m_value(methylated_noob, unmethylated_noob, offset=0):
    """ the log(base 2) (1+meth / 1+unmeth) intensities (with a min clip intensity of 1 to avoid divide-by-zero-errors, like sesame)"""
    methylated = np.clip(methylated_noob, 1, None) + offset
    unmethylated = np.clip(unmethylated_noob, 1, None) + offset

    with np.errstate(all='raise'):
        intensity_ratio = np.true_divide(methylated, unmethylated)
    return np.log2(intensity_ratio)


def calculate_copy_number(methylated_noob, unmethylated_noob, offset=None):
    """ the log(base 2) of the combined (meth + unmeth AKA green and red) intensities """
    # Note: offset is a kwarg to match other calculate functions above, but there is no offset used in this function.
    total_intensity = methylated_noob + unmethylated_noob
    copy_number = np.log2(total_intensity)
    return copy_number


def consolidate_values_for_sheet(data_containers, postprocess_func_colname='beta_value', bit='float32', poobah=True, poobah_sig=0.05, exclude_rs=True):
    """ Transforms results into a single dataframe with all of the function values,
    with probe names in rows, and sample beta values for probes in columns.

    Input:
        data_containers -- the output of run_pipeline() is this, a list of data_containers.
        (a list of processed SampleDataContainer objects)

    Arguments for postprocess_func_colname:
        calculate_beta_value --> 'beta_value'
        calculate_m_value --> 'm_value'
        calculate_copy_number --> 'cm_value'

    note: these functions are hard-coded in pipeline.py as part of process_all() step.
    note: if run_pipeline included 'sesame' option, then quality mask is automatically applied to all pickle outputs, and saved as column in processed CSV.

    Options:
        bit (float16, float32, float64) -- change the default data type from float32
            to another type to save disk space. float16 works fine, but might not be compatible
            with all numnpy/pandas functions, or with outside packages, so float32 is default.
            This is specified from methylprep process command line.
        poobah
            If true, filters by the poobah_pval column. (beta m_val pass True in for this.)
        data_container.quality_mask (True/False)
            If 'quality_mask' is present in df, True filters these probes from pickle output.
        exclude_rs
            as of v1.5.0 SigSet keeps snp ('rs') probes with other probe types (if qualityMask is false); need to separate them here
            before exporting to file."""
    poobah_column = 'poobah_pval'
    quality_mask = 'quality_mask'
    for idx,sample in enumerate(data_containers):
        sample_id = f"{sample.sample.sentrix_id}_{sample.sample.sentrix_position}"

        if poobah == True and poobah_column in sample._SampleDataContainer__data_frame.columns:
            # remove all failed probes by replacing with NaN before building DF.
            sample._SampleDataContainer__data_frame.loc[sample._SampleDataContainer__data_frame[poobah_column] >= poobah_sig, postprocess_func_colname] = np.nan
        elif poobah == True and poobah_column not in sample._SampleDataContainer__data_frame.columns:
            LOGGER.warning('DEBUG: missing poobah')

        if sample.quality_mask == True and quality_mask in sample._SampleDataContainer__data_frame.columns:
            # blank there probes where quality_mask == 0
            sample._SampleDataContainer__data_frame.loc[sample._SampleDataContainer__data_frame[quality_mask] == 0, postprocess_func_colname] = np.nan

        this_sample_values = sample._SampleDataContainer__data_frame[postprocess_func_colname]

        if exclude_rs: # dropping rows before exporting
            mask_snps = (sample._SampleDataContainer__data_frame.index.str.startswith('rs'))
            this_sample_values = this_sample_values.loc[ ~mask_snps ]

        if idx == 0:
            merged = pd.DataFrame(this_sample_values, columns=[postprocess_func_colname])
            merged.rename(columns={postprocess_func_colname: sample_id}, inplace=True)
            continue
        merged = pd.concat([merged, this_sample_values], axis=1)
        merged.rename(columns={postprocess_func_colname: sample_id}, inplace=True)
    if bit != 'float32' and bit in ('float64','float16'):
        merged = merged.astype(bit)
    return merged


def one_sample_control_snp(container):
    """Creates the control_probes.pkl dataframe for export
Unlike all the other postprocessing functions, this uses a lot of SampleDataContainer objects that get removed to save memory,
so calling this immediately after preprocessing a sample so whole batches of data are not kept in memory.

Input: one SampleDataContainer object.

Returns:
    a pickled dataframe with non-CpG probe values.
    Control: red and green intensity
    SNP: beta values, based on uncorrected meth and unmeth intensity.
    Where
        0-0.25 ~~ homogyzous-recessive
        0.25--0.75 ~~ heterozygous
        0.75--1.00 ~~ homozygous-dominant
    for each of 50-60 SNP locii on the genome.
    methylcheck can plot these to genotype samples.

Notes:
    Each data container has a ctrl_red and ctrl_green dataframe.
    This shuffles them into a single dictionary of dataframes with keys as Sample_ID matching the meta_data.
    Where Sample_ID is:
        {container.sentrix_id_sentrix_position}
    and each dataframe contains both the ctrl_red and ctrl_green data. Column names will specify which channel
    (red or green) the data applies to.

    snp values are stored in container.snp_methylated.data_frame
    and container.snp_unmethylated.data_frame

    methylcheck.plot_controls requires these columns in output: ['Control_Type','Color','Extended_Type']
        copy from container.ctl_man

    v1.5.0+ saves either RAW or NOOB meth/unmeth SNP values, depending on whether NOOB/dye was processed.
    """
    RED = container.ctrl_red.rename(columns={'mean_value': 'Mean_Value_Red'})[['Mean_Value_Red']]
    GREEN = container.ctrl_green.rename(columns={'mean_value': 'Mean_Value_Green'})[['Mean_Value_Green']]
    CONTROL = RED.join(GREEN)
    CONTROL = pd.concat([CONTROL, container.ctl_man], axis=1)

    meth_col = 'Meth' if container.do_noob == False else 'noob_Meth'
    unmeth_col = 'Unmeth' if container.do_noob == False else 'noob_Unmeth'
    SNP = container.snp_methylated.rename(columns={meth_col: 'snp_meth'})
    SNP_UNMETH = container.snp_unmethylated.rename(columns={unmeth_col: 'snp_unmeth'})[['snp_unmeth']]
    SNP = SNP.join(SNP_UNMETH)
    # below (snp-->beta) is analogous to:
    # SampleDataContainer._postprocess(input_dataframe, calculate_beta_value, 'beta_value')
    # except that it doesn't use the predefined noob columns.
    vectorized_func = np.vectorize(calculate_beta_value)
    SNP['snp_beta'] = vectorized_func(
        SNP['snp_meth'].values,
        SNP['snp_unmeth'].values,
    )
    SNP = SNP[['snp_beta','snp_meth','snp_unmeth']]

    # space saving, but catches NaN bugs without breaking.
    if SNP[['snp_meth','snp_unmeth']].isna().sum().sum() == 0:
        SNP['snp_meth'] = SNP['snp_meth'].apply(pd.to_numeric, downcast='unsigned')
        SNP['snp_unmeth'] = SNP['snp_unmeth'].apply(pd.to_numeric, downcast='unsigned')
    if CONTROL[['Mean_Value_Green','Mean_Value_Red']].isna().sum().sum() == 0:
        CONTROL['Mean_Value_Green'] = CONTROL['Mean_Value_Green'].apply(pd.to_numeric, downcast='unsigned')
        CONTROL['Mean_Value_Red'] = CONTROL['Mean_Value_Red'].apply(pd.to_numeric, downcast='unsigned')

    CONTROL = CONTROL.join(SNP, how='outer').round({'snp_beta':3})
    # finally, copy 'design' col from manifest, if exists
    if 'design' in container.man.columns:
        probe_designs = container.man[['design']]
        CONTROL = CONTROL.join(probe_designs, how='left')
    return CONTROL


def consolidate_mouse_probes(data_containers, filename_or_fileobj, file_format, object_name='mouse_data_frame', poobah_column='poobah_pval', poobah_sig=0.05):
    """ these probes have 'Multi'|'Random' in `design` col of mouse manifest. used to populate 'mouse_probes.pkl'.
    pre v1.4.6: ILLUMINA_MOUSE specific probes (starting with 'rp' for repeat sequence or 'mu' for murine, 'uk' for unknown-experimental)
    stored as data_container.mouse_data_frame.

    saves as a dataframe just like controls:
        a dict of dataframes like processed.csv format, but only mouse probes.
        keys are sample_ids -- values are dataframes"""
    out = dict()
    for idx,sample in enumerate(data_containers):
        sample_id = f"{sample.sample.sentrix_id}_{sample.sample.sentrix_position}"
        data_frame = getattr(sample, object_name)
        data_frame = data_frame.round({'noob_meth':0, 'noob_unmeth':0, 'm_value':3, 'beta_value':3, 'cm_value':3,
            'meth':0, 'unmeth':0, 'poobah_pval':3})
        out[sample_id] = data_frame

    if is_file_like(filename_or_fileobj): # and file_format == 'pickle':
        pickle.dump(out, filename_or_fileobj)
    #elif is_file_like(filename_or_fileobj) and file_format == 'parquet':

    else: #except TypeError: # File must have a write attribute
        with open(filename_or_fileobj, 'wb') as f:
            pickle.dump(out, f)
    return

def merge_batches(num_batches, data_dir, filepattern, file_format):
    """for each of the output pickle file types,
    this will merge the _1, _2, ..._X batches into a single file in data_dir.
    """
    dfs = []
    suffix = 'parquet' if file_format == 'parquet' else 'pkl'
    for num in range(num_batches):
        try:
            filename = f"{filepattern}_{num+1}.{suffix}"
            part = Path(data_dir, filename)
            if part.exists() and file_format == 'parquet':
                dfs.append( pd.read_parquet(part) )
            elif part.exists():
                dfs.append( pd.read_pickle(part) )
        except Exception as e:
            LOGGER.error(f'error merging batch {num} of {filepattern}')
    if dfs: # pipeline passes in all filenames, but not all exist
        dfs = pd.concat(dfs, axis='columns', join='inner')
        outfile_name = Path(data_dir, f"{filepattern}.{suffix}")
        LOGGER.info(f"{filepattern}: {dfs.shape}")
        func = dfs.to_parquet if file_format == 'parquet' else dfs.to_pickle
        func(str(outfile_name))
        del dfs # save memory.

        # confirm file saved ok.
        if not Path(outfile_name).exists():
            print("error saving consolidated file: {outfile_name}; use methylcheck.load() to merge the parts")
            return

    # now delete the parts
    for num in range(num_batches):
        filename = f"{filepattern}_{num+1}.{suffix}"
        part = Path(data_dir, filename)
        if part.exists():
            part.unlink() # delete it
