# Lib
import logging
import numpy as np
import pandas as pd
# App
from ..models import ProbeType, Channel

__all__ = ['infer_type_I_probes']

LOGGER = logging.getLogger(__name__)

def infer_type_I_probes(container, debug=False):
    """ Adapted from sesaame from https://github.com/zwdzwd/sesame/blob/RELEASE_3_12/R/channel_inference.R.
    -- pass in a SampleDataContainer
    -- runs in SampleDataContainer.__init__ this BEFORE qualityMask step, so NaNs are not present
    -- changes raw_data idat probe_means
    -- runs on raw_dataset, before meth-dataset is created, so @IR property doesn't exist yet; but get_infer has this"""
    # this first step combines all I-red and I-green channel intensities, so IG+oobG and IR+oobR.
    channels = get_infer_channel_probes(container.manifest, container.green_idat, container.red_idat, debug=debug)
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
        if len(red_max) == 0:
            print('No probes were swapped because there are no type-I-ref probes detected!')
        else:
            percent_probes_ok = 100* sum( np.maximum(red_max, green_max) > min_ib ) / len(red_max)
            print(f"min_ib: {min_ib}, %swapped: {round(100-percent_probes_ok,3)} ({count_probes_to_swap})")
        print('R2R', R2R, 'G2G', G2G)
        print('R2G', R2G, 'G2R', G2R)
        print('FailedR', FailedR, 'FailedG', FailedG)

    # finally, actually swap these probe values in the container and return nothing.
    R2G_mask = np.where(red_I_channel.index.isin(channels['IR'].index) & ~red_idx & big_idx, True, False)
    G2R_mask = np.where(green_I_channel.index.isin(channels['IG'].index) & red_idx & big_idx, True, False)

    # this runs EARLY in processing, so modifying red_idat and green_idat directly.
    lookup = channels['lookup'] # index are cpg probe names;
    # Filter lookup to probes that need to be swapped and get a list of all the addresses that belong to them
    R2G_lookup = lookup[lookup.index.isin(red_I_channel.index[R2G_mask])]
    R2G_illumina_ids = R2G_lookup["AddressA_ID"].to_list() + R2G_lookup["AddressB_ID"].to_list()

    G2R_lookup = lookup[lookup.index.isin(green_I_channel.index[G2R_mask])]
    G2R_illumina_ids = G2R_lookup["AddressA_ID"].to_list() + G2R_lookup["AddressB_ID"].to_list()
    # swap probe values
    mask = container.red_idat.probe_means.index.isin(R2G_illumina_ids + G2R_illumina_ids)
    pre_red = container.red_idat.probe_means.copy()
    pre_green = container.green_idat.probe_means.copy()
    post_red = container.red_idat.probe_means.copy()
    post_green = container.green_idat.probe_means.copy()

    # green --> red
    post_red.loc[mask, 'mean_value'] = pre_green.loc[mask, 'mean_value']
    # original red --> green
    post_green.loc[mask, 'mean_value'] = pre_red.loc[mask, 'mean_value']

    container.red_switched = list(red_I_channel.index[R2G_mask])
    container.green_switched = list(green_I_channel.index[G2R_mask])
    #print(f"switched {len(container.red_switched)} red and {len(container.green_switched)} green probes")
    container.red_idat.probe_means = post_red
    container.green_idat.probe_means = post_green
    return



def get_infer_channel_probes(manifest, green_idat, red_idat, debug=False):
    """ like filter_oob_probes, but returns two dataframes for green and red channels with meth and unmeth columns
    effectively criss-crosses the red-oob channels and appends to green, and appends green-oob to red
    returns a dict with 'green' and 'red' channel probes

    THIS runs before processing in SampleDataContainer, so that infer_type_I_probes() can modify the IDAT probe_means directly.    """
    probe_details_IR = manifest.get_probe_details(
        probe_type=ProbeType.ONE,
        channel=Channel.RED,
    )
    probe_details_IG = manifest.get_probe_details(
        probe_type=ProbeType.ONE,
        channel=Channel.GREEN,
    )
    # need: illumina_id in index, (green)'meth', (red)'unmeth'
    idat_meth_unmeth = (
        green_idat.probe_means
        .rename(columns={'mean_value':'meth'})
        .sort_index()
        .merge(
        red_idat.probe_means.rename(columns={'mean_value':'unmeth'}).sort_index(),
        sort=False,
        left_index=True,
        right_index=True)
        )

    # OOB PROBE values are IR(unmeth) and IG(meth); I'll replace IR(meth) and IG(unmeth) below
    # RED channel; uses AddressA_ID for oob IR(unmeth)
    oobR = probe_details_IG.merge(
        idat_meth_unmeth,
        how='inner',
        left_on='AddressA_ID',
        right_index=True,
        suffixes=(False, False),
    )
    oobR = oobR[['meth','unmeth']].sort_index()
    # GREEN channel; AddressB_ID for oob IG(meth)
    oobG = probe_details_IR.merge(
        idat_meth_unmeth,
        how='inner',
        left_on='AddressB_ID',
        right_index=True,
        suffixes=(False, False),
    )
    oobG = oobG[['meth','unmeth']].sort_index()

    # IN BAND probes should be IR(meth) and IG(unmeth)
    # NOTE: below uses same idat DF as before, but the probe_details from manifest are swapped.
    red_in_band = probe_details_IR.merge(
        idat_meth_unmeth,
        how='inner',
        left_on='AddressA_ID',
        right_index=True,
        suffixes=(False, False),
    )
    red_in_band = red_in_band[['meth','unmeth']].sort_index()
    green_in_band = probe_details_IG.merge(
        idat_meth_unmeth,
        how='inner',
        left_on='AddressB_ID',
        right_index=True,
        suffixes=(False, False),
    )
    green_in_band = green_in_band[['meth','unmeth']].sort_index()

    ## HACK: I can't read/get idats to match sesame exactly, so moving columns around to match
    # - swap oob-green[unmeth] with red[meth]
    # - swap oob-red[meth] with green[unmeth]
    oobG_unmeth = oobG['unmeth'].copy()
    oobR_meth = oobR['meth'].copy()
    oobG['unmeth'] = red_in_band['meth'].copy()
    oobR['meth'] = green_in_band['unmeth'].copy()
    red_in_band['meth'] = oobG_unmeth
    green_in_band['unmeth'] = oobR_meth

    # next, add the green-in-band to oobG and red-in-band to oobR
    oobG_IG = oobG.append(green_in_band).sort_index()
    oobR_IR = oobR.append(red_in_band).sort_index()

    # channel swap requires a way to update idats with illumina_ids
    lookupIR = probe_details_IR.merge(
        idat_meth_unmeth,
        how='inner',
        left_on='AddressA_ID',
        right_index=True,
        suffixes=(False, False),
    )[['AddressA_ID','AddressB_ID']]
    lookupIG = probe_details_IG.merge(
        idat_meth_unmeth,
        how='inner',
        left_on='AddressB_ID',
        right_index=True,
        suffixes=(False, False),
    )[['AddressA_ID','AddressB_ID']]
    lookup = lookupIG.append(lookupIR).sort_index()

    if debug:
        return {'green': oobG_IG, 'red': oobR_IR, 'oobG': oobG, 'oobR':oobR, 'IG': green_in_band, 'IR': red_in_band, 'lookup': lookup}
    return {'green': oobG_IG, 'red': oobR_IR, 'IR': red_in_band, 'IG': green_in_band, 'lookup': lookup}
