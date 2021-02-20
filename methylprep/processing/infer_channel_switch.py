# Lib
import logging
import numpy as np
import pandas as pd
# App

__all__ = ['infer_type_I_probes']

LOGGER = logging.getLogger(__name__)

def infer_type_I_probes(container, debug=False):
    """ Adapted from sesaame from https://github.com/zwdzwd/sesame/blob/RELEASE_3_12/R/channel_inference.R.
    -- pass in a SampleDataContainer
    -- runs in SampleDataContainer.__init__ this BEFORE qualityMask step, so NaNs are not present
    -- changes raw_data idat probe_means
    -- runs on raw_dataset, before meth-dataset is created, so @IR property doesn't exist yet; but get_infer has this"""
    # this first step combines all I-red and I-green channel intensities, so IG+oobG and IR+oobR.
    channels = container.raw_dataset.get_infer_channel_probes(container.manifest)
    green_I_channel = channels['green']
    red_I_channel = channels['red']
    ## NAN probes occurs when manifest is not complete
    ## If there are NA in the probe intensity, exclude these probes.
    n_red = channels['IR'].count() # == container.IR.count() ## nrow(IR(sset)) | here, I am using count() to ignore missing probes
    green_I_channel = green_I_channel.dropna()
    red_I_channel = red_I_channel.dropna()

    # get the higher of each channel per probe (thus there are 4 values per probe compared here; red meth, red unmeth, green meth, green unmeth)
    red_max = red_I_channel.max(axis=1)
    green_max = green_I_channel.max(axis=1)
    red_idx = (red_max > green_max) # TRUE mask; FALSE means the channel will be swapped

    # min_ib: take the lower of the channels and calculate quantile score,
    # then exclude if lower than the value where 95% of values would be above this range
    min_ib = pd.DataFrame(np.minimum(red_I_channel.min(axis=1), green_I_channel.min(axis=1))) .quantile(0.95, axis=0)
    # min_ib is ONE number, the low-cutoff intensity. == 644 in sesame testing
    min_ib = min_ib.values[0]
    # now compare the higher of each channel and confirm it is always greater than the min_ib
    big_idx = (np.maximum(red_max, green_max) > min_ib) # a TRUE mask, probes that are OK
    percent_probes_ok = 100* sum( np.maximum(red_max, green_max) > min_ib ) / len(red_max)
    count_probes_to_swap = sum(np.maximum(red_max, green_max) <= min_ib)
    # goal here: create a mask with TRUE/FALSE for every probe that is swapped or not. Then update data.
    # print(f"big_idx: {big_idx}")
    # if test fails, swap the value that corresponds to RED at the original index, using red_max0
    # red_idx = np.where(big_idx, red_idx, False) # another TRUE mask
    R2R = sum(np.where(red_I_channel.index.isin(channels['IR'].index) & red_idx & big_idx, True, False))
    G2G = sum(np.where(green_I_channel.index.isin(channels['IG'].index) & ~red_idx & big_idx, True, False))
    R2G = sum(np.where(red_I_channel.index.isin(channels['IR'].index) & ~red_idx & big_idx, True, False))
    G2R = sum(np.where(green_I_channel.index.isin(channels['IG'].index) & red_idx & big_idx, True, False))
    FailedR = sum(np.where(red_I_channel.index.isin(channels['IR'].index) & ~big_idx, True, False))
    FailedG = sum(np.where(green_I_channel.index.isin(channels['IG'].index) & ~big_idx, True, False))

    if debug:
        print(f"min_ib: {min_ib}, %swapped: {round(100-percent_probes_ok,3)} ({count_probes_to_swap})")
        print('R2R', R2R, 'G2G', G2G)
        print('R2G', R2G, 'G2R', G2R)
        print('FailedR', FailedR, 'FailedG', FailedG)

    # finally, actually swap these probe values in the container and return nothing.
    R2G_mask = np.where(red_I_channel.index.isin(channels['IR'].index) & ~red_idx & big_idx, True, False)
    G2R_mask = np.where(green_I_channel.index.isin(channels['IG'].index) & red_idx & big_idx, True, False)

    # this runs EARLY in processing, so modifying red_idat and green_idat directly.
    lookup = channels['lookup'] # index are cpg probe names;
    # -- illumina_ids are in 'AddressA_ID'(red) and 'AddressB_ID'(green) for in-band matching
    lookupIG = dict(zip(lookup.index,lookup['AddressB_ID']))
    lookupIR = dict(zip(lookup.index,lookup['AddressA_ID']))

    # swap probe values and save a probe list for R2G and G2R idat probes
    R2G_illumina_ids = [lookupIR[i] for i in red_I_channel.index[R2G_mask]]
    G2R_illumina_ids = [lookupIG[i] for i in green_I_channel.index[G2R_mask]]
    mask = container.raw_dataset.red_idat.probe_means.index.isin(R2G_illumina_ids + G2R_illumina_ids)
    pre_red = container.raw_dataset.red_idat.probe_means.copy()
    pre_green = container.raw_dataset.green_idat.probe_means.copy()
    post_red = container.raw_dataset.red_idat.probe_means.copy()
    post_green = container.raw_dataset.green_idat.probe_means.copy()

    # green --> red
    post_red.loc[mask, 'mean_value'] = pre_green.loc[mask, 'mean_value']
    # original red --> green
    post_green.loc[mask, 'mean_value'] = pre_red.loc[mask, 'mean_value']

    container.red_switched = list(red_I_channel.index[R2G_mask])
    container.green_switched = list(green_I_channel.index[G2R_mask])
    #print(f"switched {len(container.red_switched)} red and {len(container.green_switched)} green probes")
    container.raw_dataset.red_idat.probe_means = post_red
    container.raw_dataset.green_idat.probe_means = post_green
    return
