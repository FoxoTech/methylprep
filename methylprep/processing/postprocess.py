# Lib
import numpy as np
import pandas as pd
import os
import pickle
# app
from ..utils import is_file_like

os.environ['NUMEXPR_MAX_THREADS'] = "8" # suppresses warning


__all__ = ['calculate_beta_value', 'calculate_m_value', 'consolidate_values_for_sheet',
    'consolidate_control_snp', 'consolidate_mouse_probes']


def calculate_beta_value(methylated_noob, unmethylated_noob, offset=100):
    """ the ratio of (methylated_intensity / total_intensity)
    where total_intensity is (meth + unmeth + 100) -- to give a score in range of 0 to 1.0"""
    methylated = max(methylated_noob, 0)
    unmethylated = max(unmethylated_noob, 0)

    total_intensity = methylated + unmethylated + offset
    with np.errstate(all='raise'):
        intensity_ratio = np.true_divide(methylated, total_intensity)
    return intensity_ratio


def calculate_m_value(methylated_noob, unmethylated_noob, offset=1):
    """ the log(base 2) (1+meth / (1+unmeth_ intensities (with an offset to avoid divide-by-zero-errors)"""
    methylated = methylated_noob + offset
    unmethylated = unmethylated_noob + offset

    with np.errstate(all='raise'):
        intensity_ratio = np.true_divide(methylated, unmethylated)
    return np.log2(intensity_ratio)


def calculate_copy_number(methylated_noob, unmethylated_noob):
    """ the log(base 2) of the combined (meth + unmeth AKA green and red) intensities """
    total_intensity = methylated_noob + unmethylated_noob
    copy_number = np.log2(total_intensity)
    return copy_number


def consolidate_values_for_sheet(data_containers, postprocess_func_colname='beta_value', bit='float32', poobah=False):
    """ with a data_containers (list of processed SampleDataContainer objects),
    this will transform results into a single dataframe with all of the function values,
    with probe names in rows, and sample beta values for probes in columns.

    Input:
        data_containers -- the output of run_pipeline() is this, a list of data_containers.

    Arguments for postprocess_func_colname:
        calculate_beta_value --> 'beta_value'
        calculate_m_value --> 'm_value'
        calculate_copy_number --> 'cm_value'

    note: these functions are hard-coded in pipeline.py as part of process_all() step.

    Options:
        bit (float16, float32, float64) -- change the default data type from float32
            to another type to save disk space. float16 works fine, but might not be compatible
            with all numnpy/pandas functions, or with outside packages, so float32 is default.
            This is specified from methylprep process command line.

        poobah
            If true, filters by the poobah_pval column. (beta m_val pass True in for this.)
        """
    poobah_column = 'poobah_pval'
    pval_cutoff = 0.05
    for idx,sample in enumerate(data_containers):
        sample_id = f"{sample.sample.sentrix_id}_{sample.sample.sentrix_position}"

        if poobah == True and poobah_column in sample._SampleDataContainer__data_frame.columns:
            # remove all failed probes by replacing with NaN before building DF.
            sample._SampleDataContainer__data_frame.loc[sample._SampleDataContainer__data_frame[poobah_column] >= pval_cutoff, postprocess_func_colname] = np.nan
        elif poobah == True and poobah_column not in sample._SampleDataContainer__data_frame.columns:
            print('DEBUG: missing poobah')

        this_sample_values = sample._SampleDataContainer__data_frame[postprocess_func_colname]
        if idx == 0:
            merged = pd.DataFrame(this_sample_values, columns=[postprocess_func_colname])
            merged.rename(columns={postprocess_func_colname: sample_id}, inplace=True)
            continue
        merged = pd.concat([merged, this_sample_values], axis=1)
        merged.rename(columns={postprocess_func_colname: sample_id}, inplace=True)
    if bit != 'float32' and bit in ('float64','float16'):
        merged = merged.astype(bit)
    return merged


def consolidate_control_snp(data_containers, filename_or_fileobj):
    """saves a pickled dataframe with non-CpG probe values.
Returns:
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
        {sample.sentrix_id_sample.sentrix_position}
    and each dataframe contains both the ctrl_red and ctrl_green data. Column names will specify which channel
    (red or green) the data applies to.

    snp values are stored in container.snp_methylated.data_frame
    and container.snp_unmethylated.data_frame
    """
    out = {}
    for idx,sample in enumerate(data_containers):
        sample_id = f"{sample.sample.sentrix_id}_{sample.sample.sentrix_position}"
        RED = sample.ctrl_red.rename(columns={'mean_value': 'Mean_Value_Red'})
        GREEN = sample.ctrl_green.rename(columns={'mean_value': 'Mean_Value_Green'})
        GREEN = GREEN.drop(['Control_Type', 'Color', 'Extended_Type'], axis='columns')

        SNP = sample.snp_methylated.data_frame.rename(columns={'mean_value': 'snp_meth'})
        SNP_UNMETH = sample.snp_unmethylated.data_frame.rename(columns={'mean_value': 'snp_unmeth'})
        SNP_UNMETH = SNP_UNMETH.loc[:, ['snp_unmeth']]
        SNP = pd.merge(SNP, SNP_UNMETH, left_index=True, right_index=True, how='outer')
        # below (snp-->beta) is analogous to:
        # SampleDataContainer._postprocess(input_dataframe, calculate_beta_value, 'beta_value')
        # except that it doesn't use the predefined noob columns.
        vectorized_func = np.vectorize(calculate_beta_value)
        SNP['snp_beta'] = vectorized_func(
            SNP['snp_meth'].values,
            SNP['snp_unmeth'].values,
        )
        SNP = SNP[['snp_beta','snp_meth','snp_unmeth']]
        SNP = SNP.astype({
            'snp_meth': 'int32',
            'snp_unmeth': 'int32'})

        merged = pd.merge(RED, GREEN, left_index=True, right_index=True, how='outer')
        merged = merged.astype({
            'Mean_Value_Green': 'int32',
            'Mean_Value_Red': 'int32'})
        merged = pd.merge(merged, SNP, left_index=True, right_index=True, how='outer')
        merged = merged.round({'snp_beta':3})
        out[sample_id] = merged

    if is_file_like(filename_or_fileobj):
        pickle.dump(out, filename_or_fileobj)
    else:
        #except TypeError: # File must have a write attribute
        with open(filename_or_fileobj, 'wb') as f:
            pickle.dump(out, f)
    return


def consolidate_mouse_probes(data_containers, filename_or_fileobj, object_name='mouse_data_frame'):
    """ ILLUMINA_MOUSE specific probes (starting with 'rp' for repeat sequence or 'mu' for murine)
    stored as data_container.mouse_data_frame.

    saves as a dataframe just like controls:
        a dict of dataframes like processed.csv format, but only mouse probes.
        keys are sample_ids -- values are dataframes"""
    poobah_column = 'poobah_pval'
    pval_cutoff = 0.05
    out = dict()
    for idx,sample in enumerate(data_containers):
        sample_id = f"{sample.sample.sentrix_id}_{sample.sample.sentrix_position}"
        out[sample_id] = getattr(sample, object_name)

    if is_file_like(filename_or_fileobj):
        pickle.dump(out, filename_or_fileobj)
    else:
        #except TypeError: # File must have a write attribute
        with open(filename_or_fileobj, 'wb') as f:
            pickle.dump(out, f)
    return
