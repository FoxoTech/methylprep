# Lib
import numpy as np
import pandas as pd


__all__ = ['calculate_beta_value', 'calculate_m_value', 'consolidate_values_for_sheet']


def calculate_beta_value(methylated_noob, unmethylated_noob, offset=100):
    methylated = max(methylated_noob, 0)
    unmethylated = max(unmethylated_noob, 0)

    total_intensity = methylated + unmethylated + offset
    with np.errstate(all='raise'):
        intensity_ratio = methylated / total_intensity
    return intensity_ratio


def calculate_m_value(methylated_noob, unmethylated_noob, offset=1):
    methylated = methylated_noob + offset
    unmethylated = unmethylated_noob + offset

    with np.errstate(all='raise'):
        intensity_ratio = np.true_divide(methylated, unmethylated)
    return np.log2(intensity_ratio)


def calculate_copy_number(methylated_noob, unmethylated_noob):
    total_intensity = methylated_noob + unmethylated_noob
    copy_number = np.log2(total_intensity)
    return copy_number


def consolidate_values_for_sheet(data_containers, postprocess_func_colname='beta_value', bit='float64'):
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
        bit (float16, float32, float64) -- change the default data type from float64
            to another type to save disk space. float16 works fine, but might not be compatible
            with all numnpy/pandas functions, or with outside packages, so float64 is default.
            This is specified from methylprep process command line."""
    for idx,sample in enumerate(data_containers):
        sample_id = f"{sample.sample.sentrix_id}_{sample.sample.sentrix_position}"
        this_sample_values = sample._SampleDataContainer__data_frame[postprocess_func_colname]
        if idx == 0:
            merged = pd.DataFrame(this_sample_values, columns=[postprocess_func_colname])
            merged.rename(columns={postprocess_func_colname: sample_id}, inplace=True)
            continue
        merged = pd.concat([merged, this_sample_values], axis=1)
        merged.rename(columns={postprocess_func_colname: sample_id}, inplace=True)
    if bit != 'float64' and bit in ('float32','float16'):
        merged = merged.astype(bit)
    return merged


def consolidate_control_snp(data_containers, control_filename):
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
    import pickle
    with open(control_filename, 'wb') as f:
        pickle.dump(out, f)
    return
