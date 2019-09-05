# Lib
import numpy as np
import pandas as pd


__all__ = ['calculate_beta_value', 'calculate_m_value', 'consolidate_values_for_sheet']


def calculate_beta_value(methylated_noob, unmethylated_noob, offset=100):
    methylated = max(methylated_noob, 0)
    unmethylated = max(unmethylated_noob, 0)

    total_intensity = methylated + unmethylated + offset
    intensity_ratio = methylated / total_intensity
    return intensity_ratio


def calculate_m_value(methylated_noob, unmethylated_noob, offset=1):
    methylated = methylated_noob + offset
    unmethylated = unmethylated_noob + offset

    intensity_ratio = np.true_divide(methylated, unmethylated)
    return np.log2(intensity_ratio)


def calculate_copy_number(methylated_noob, unmethylated_noob):
    total_intensity = methylated_noob + unmethylated_noob
    copy_number = np.log2(total_intensity)
    return copy_number


def consolidate_values_for_sheet(data_containers, postprocess_func_colname='beta_value'):
    """ with a data_containers (list of processed SampleDataContainer objects),
    this will transform results into a single dataframe with all of the function values,
    with probe names in rows, and sample beta values for probes in columns.

    Input:
        data_containers -- the output of run_pipeline() is this, a list of data_containers.

    Arguments for postprocess_func_colname:
        calculate_beta_value --> 'beta_value'
        calculate_m_value --> 'm_value'
        calculate_copy_number --> 'cm_value'

    note: these are hard-coded in pipeline.py as part of process_all() step.
    """
    for idx,sample in enumerate(data_containers):
        sample_id = "{0}_{1}".format(sample.sample.sentrix_id,sample.sample.sentrix_position)
        this_sample_values = sample._SampleDataContainer__data_frame[postprocess_func_colname]
        if idx == 0:
            merged = pd.DataFrame(this_sample_values, columns=[postprocess_func_colname])
            merged.rename(columns={postprocess_func_colname: sample_id}, inplace=True)
            continue
        merged = pd.concat([merged, this_sample_values], axis=1)
        merged.rename(columns={postprocess_func_colname: sample_id}, inplace=True)
    return merged
