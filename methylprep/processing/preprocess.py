# Normal-exponential using out-of-band probes
# normex: negative control probes
# noob: ‘out-of-band’ Infinium I probes

# Lib
import logging
import numpy as np
import pandas as pd
from statsmodels import robust
from scipy.stats import norm, lognorm
# App
from ..models import ControlType, ArrayType
from ..models.sketchy_probes import qualityMask450, qualityMaskEPIC, qualityMaskEPICPLUS, qualityMaskmouse


__all__ = ['preprocess_noob']


LOGGER = logging.getLogger(__name__)


def preprocess_noob(container, offset=15, pval_probes_df=None, quality_mask_df=None, nonlinear_dye_correction=True, debug=False, unit_test_oob=False): # v1.4.5+
    """ NOOB pythonized copy of https://github.com/zwdzwd/sesame/blob/master/R/background_correction.R
    - The function takes a SigSet and returns a modified SigSet with the background subtracted.
    - Background is modelled in a normal distribution and true signal in an exponential distribution.
    - The Norm-Exp deconvolution is parameterized using Out-Of-Band (oob) probes.
    - includes snps, but not control probes yet
    - output should replace the container instead of returning debug dataframes
    - II RED and II GREEN both have data, but manifest doesn't have a way to track this, so function tracks it.
    - keep IlmnID as index for meth/unmeth snps, and convert fg_green

    if nonlinear_dye_correction=True, this uses a sesame method in place of minfi method, in a later step.
    if unit_test_oob==True, returns the intermediate data instead of updating the SigSet/SampleDataContainer.
    """
    if debug:
        print(f"DEBUG NOOB {debug} nonlinear_dye_correction={nonlinear_dye_correction}, pval_probes_df={pval_probes_df.shape if isinstance(pval_probes_df,pd.DataFrame) else 'None'}, quality_mask_df={quality_mask_df.shape if isinstance(quality_mask_df,pd.DataFrame) else 'None'}")
    # stack- need one long list of values, regardless of Meth/Uneth
    ibG = pd.concat([
        container.ibG.reset_index().rename(columns={'Meth': 'mean_value'}).assign(used='M'),
        container.ibG.reset_index().rename(columns={'Unmeth': 'mean_value'}).assign(used='U')
    ])
    ibG = ibG[ ~ibG['mean_value'].isna() ].drop(columns=['Meth','Unmeth'])

    ibR = pd.concat([
        container.ibR.reset_index().rename(columns={'Meth': 'mean_value'}).assign(used='M'), #.drop(columns=['Meth','Unmeth']),
        container.ibR.reset_index().rename(columns={'Unmeth': 'mean_value'}).assign(used='U') #.drop(columns=['Meth','Unmeth'])
    ])
    ibR = ibR[ ~ibR['mean_value'].isna() ].drop(columns=['Meth','Unmeth'])

    # out-of-band is Green-Unmeth and Red-Meth
    # exclude failing probes
    pval = pval_probes_df.loc[ pval_probes_df['poobah_pval'] > container.poobah_sig ].index if isinstance(pval_probes_df, pd.DataFrame) else []
    qmask = quality_mask_df.loc[ quality_mask_df['quality_mask'] == 0 ].index if isinstance(quality_mask_df, pd.DataFrame) else []
    # the ignored errors here should only be from probes that are both pval failures and qmask failures.
    Rmeth = list(container.oobR['Meth'].drop(index=pval, errors='ignore').drop(index=qmask, errors='ignore'))
    Runmeth = list(container.oobR['Unmeth'].drop(index=pval, errors='ignore').drop(index=qmask, errors='ignore'))
    oobR = pd.DataFrame( Rmeth + Runmeth, columns=['mean_value'])
    Gmeth = list(container.oobG['Meth'].drop(index=pval, errors='ignore').drop(index=qmask, errors='ignore'))
    Gunmeth = list(container.oobG['Unmeth'].drop(index=pval, errors='ignore').drop(index=qmask, errors='ignore'))
    oobG = pd.DataFrame( Gmeth + Gunmeth, columns=['mean_value'])
    # minfi test
    # ref fg_green = 442614 | vs ibG 442672 = 396374 + 46240
    # ref fg_red = 528410 | vs ibR 528482 = 439279 + 89131
    # ref oob_green = 178374
    # ref oob_red = 92578
    #oobR = pd.DataFrame( data={'mean_value': container.oobR['Meth']})
    #oobG = pd.DataFrame( data={'mean_value': container.oobG['Unmeth']})
    #print(f" oobR {oobR.shape} oobG {oobG.shape}")
    #import pdb;pdb.set_trace()

    debug_warnings = ""
    if oobR['mean_value'].isna().sum() > 0:
        debug_warnings += f" NOOB: oobG had {oobG['mean_value'].isna().sum()} NaNs"
        oobR = oobR.dropna()
    if oobG['mean_value'].isna().sum() > 0:
        debug_warnings += f" NOOB: oobG had {oobG['mean_value'].isna().sum()} NaNs"
        oobG = oobG.dropna()
    if ibG['mean_value'].isna().sum() > 0 or ibR['mean_value'].isna().sum() > 0:
        raise ValueError("ibG or ibR is missing probe intensities. need to filter them out.")

    if debug:
        print(f"ibG {len(ibG)} ibR {len(ibR)} oobG {len(oobG)} oobR {len(oobR)} | {debug_warnings}")

    # set minimum intensity to 1
    ibG_affected = len(ibG.loc[ ibG['mean_value'] < 1 ].index)
    ibR_affected = len(ibR.loc[ ibR['mean_value'] < 1 ].index)
    ibG.loc[ ibG['mean_value'] < 1, 'mean_value'] = 1
    ibR.loc[ ibR['mean_value'] < 1, 'mean_value'] = 1
    oobG_affected = len(oobG[ oobG['mean_value'] < 1])
    oobR_affected = len(oobR[ oobR['mean_value'] < 1])
    oobG.loc[ oobG.mean_value < 1, 'mean_value'] = 1
    oobR.loc[ oobR.mean_value < 1, 'mean_value'] = 1
    if debug:
        if ibR_affected > 0 or ibR_affected > 0:
            print(f"ib: Set {ibR_affected} red and {ibG_affected} green to 1.0 ({len(ibR[ ibR['mean_value'] == 1 ].index)}, {len(ibG[ ibG['mean_value'] == 1 ].index)})")
        if oobG_affected > 0 or oobR_affected > 0:
            print(f"oob: Set {oobR_affected} red and {oobG_affected} green to 1.0 ({len(oobR[ oobR['mean_value'] == 1 ].index)}, {len(oobG[ oobG['mean_value'] == 1 ].index)})")

    # do background correction in each channel; returns "normalized in-band signal"
    ibG_nl, params_green = normexp_bg_corrected(ibG, oobG, offset, sample_name=container.sample.name)
    ibR_nl, params_red   = normexp_bg_corrected(ibR, oobR, offset, sample_name=container.sample.name)
    noob_green = ibG_nl.round({'bg_corrected':0})
    noob_red = ibR_nl.round({'bg_corrected':0})

    if unit_test_oob:
        return {
            'oobR': oobR,
            'oobG': oobG,
            'noob_green': noob_green,
            'noob_red': noob_red,
        }

    # by default, this last step is omitted for sesame
    if nonlinear_dye_correction == True:
        # update() expects noob_red/green to have IlmnIDs in index, and contain bg_corrected for ALL probes.
        container.update_probe_means(noob_green, noob_red)
    elif nonlinear_dye_correction == False:
        # this "linear" method may be anologous to the ratio quantile normalization described in Nature: https://www.nature.com/articles/s41598-020-72664-6
        normexp_bg_correct_control(container.ctrl_green, params_green)
        normexp_bg_correct_control(container.ctrl_red, params_red)
        mask_green = container.ctrl_green['Control_Type'].isin(ControlType.normalization_green())
        mask_red = container.ctrl_red['Control_Type'].isin(ControlType.normalization_red())
        avg_green = container.ctrl_green[mask_green]['bg_corrected'].mean()
        avg_red = container.ctrl_red[mask_red]['bg_corrected'].mean()
        rg_ratios = avg_red / avg_green
        red_factor = 1 / rg_ratios
        container.update_probe_means(noob_green, noob_red, red_factor)
        container._SigSet__minfi_noob = True
    elif nonlinear_dye_correction is None:
        if debug:
            LOGGER.info("skipping linear/nonlinear dye-bias correction step")
        # skips the minfi-linear step and won't trigger the sesame nonlinear dye bias step downstream, if you REALLY want it uncorrected. Mostly for debugging / benchmarking.
        container.update_probe_means(noob_green, noob_red)


class BackgroundCorrectionParams():
    """ used in apply_bg_correction """
    __slots__ = (
        'bg_mean',
        'bg_mad',
        'mean_signal',
        'offset',
    )

    def __init__(self, bg_mean, bg_mad, mean_signal, offset):
        # note: default offset was 15. In v1.3.3 (Jan 2020) I kept 15, after finding this made results match sesame's NOOB output exactly, if dye step ommitted.
        # offset is specified in the preprocess_noob() function.
        self.bg_mean = bg_mean
        self.bg_mad = bg_mad
        self.mean_signal = mean_signal
        self.offset = offset


def normexp_bg_corrected(fg_probes, ctrl_probes, offset, sample_name=None):
    """ analogous to sesame's backgroundCorrectionNoobCh1 """
    fg_means = fg_probes['mean_value']
    if fg_means.min() == fg_means.max():
        LOGGER.error(f"{sample_name}: min and max intensity are same. Sample probably bad.")
        params = BackgroundCorrectionParams(bg_mean=1.0, bg_mad=1.0, mean_signal=1.0, offset=15)
        fg_probes['bg_corrected'] = 1.0
        return fg_probes, params
    fg_mean, _fg_mad = huber(fg_means)
    bg_mean, bg_mad = huber(ctrl_probes['mean_value'])
    mean_signal = np.maximum(fg_mean - bg_mean, 10) # "alpha" in sesame function

    params = BackgroundCorrectionParams(bg_mean, bg_mad, mean_signal, offset)

    corrected_signals = apply_bg_correction(fg_means, params)
    fg_probes['bg_corrected'] = corrected_signals
    fg_probes['bg_corrected'] = fg_probes['bg_corrected'].round(1)
    return fg_probes, params


def normexp_bg_correct_control(control_probes, params):
    """Function for getting xcs controls for preprocessNoob"""
    control_means = control_probes['mean_value']
    corrected_signals = apply_bg_correction(control_means, params)
    control_probes['bg_corrected'] = corrected_signals
    return control_probes


def apply_bg_correction(mean_values, params):
    """ this function won't work with float16 in practice (underflow). limits use to float32 """
    if not isinstance(params, BackgroundCorrectionParams):
        raise ValueError('params is not a BackgroundCorrectionParams instance')

    np.seterr(under='ignore') # 'raise to explore fixing underflow warning here'

    bg_mean = params.bg_mean #mu
    bg_mad = params.bg_mad #sigma
    mean_signal = params.mean_signal #alpha
    offset = params.offset

    mu_sf = mean_values - bg_mean - (bg_mad ** 2) / mean_signal

    #try:
    #    signal_part_one = mu_sf + (bg_mad ** 2)
    #    signal_part_two = np.exp(norm(mu_sf, bg_mad).logpdf(0) - norm(mu_sf, bg_mad).logsf(0))
    #    signal = signal_part_one * signal_part_two
    #except:
    #    print(signal_part_one, norm(mu_sf, bg_mad).logpdf(0),  norm(mu_sf, bg_mad).logsf(0))
    # norm is from scipy.stats
    signal = mu_sf + (bg_mad ** 2) * np.exp(norm(mu_sf, bg_mad).logpdf(0) - norm(mu_sf, bg_mad).logsf(0))

    """ COMPARE with sesame:
    signal <- mu.sf + sigma2 * exp(
        dnorm(0, mean = mu.sf, sd = sigma, log = TRUE) -
            pnorm(
                0, mean = mu.sf, sd = sigma,
                lower.tail = FALSE, log.p = TRUE))
    """

    # sesame: "Limit of numerical accuracy reached with very low intensity or very high background:
    # setting adjusted intensities to small value"
    signal = np.maximum(signal, 1e-6)
    true_signal = signal + offset
    return true_signal


def huber(vector):
    """Huber function. Designed to mirror MASS huber function in R

    Parameters
    ----------
    vector: list
        list of float values

    Returns
    -------
    local_median: float
        calculated mu value
    mad_scale: float
        calculated s value
    """
    num_values = len(vector)
    positive_factor = 1.5
    convergence_tol = 1.0e-6
    mad_scale = robust.mad(vector)
    local_median = np.median(vector)
    init_local_median = local_median

    if not (local_median or mad_scale):
        return local_median, mad_scale

    while True:
        yy = np.minimum(
            np.maximum(
                local_median - positive_factor * mad_scale,
                vector,
            ),
            local_median + positive_factor * mad_scale,
        )

        init_local_median = sum(yy) / num_values

        if abs(local_median - init_local_median) < convergence_tol * mad_scale:
            return local_median, mad_scale

        local_median = init_local_median


def _apply_sesame_quality_mask(data_container):
    """ adapted from sesame's qualityMask function, which is applied just after poobah
    to remove probes Wanding thinks are sketchy.
    OUTPUT: this pandas DataFrame will have NaNs for probes to be excluded and 0.0 for probes to be retained. NaNs converted to 1.0 in final processing output.

    SESAME:
        masked <- sesameDataGet(paste0(sset@platform, '.probeInfo'))$mask
        to use TCGA masking, only applies to HM450

    """
    if data_container.array_type not in (
        # ArrayType.ILLUMINA_27K,
        ArrayType.ILLUMINA_450K,
        ArrayType.ILLUMINA_EPIC,
        ArrayType.ILLUMINA_EPIC_PLUS,
        ArrayType.ILLUMINA_MOUSE):
        LOGGER.info(f"Quality masking is not supported for {data_container.array_type}.")
        return
    # load set of probes to remove from local file
    if data_container.array_type == ArrayType.ILLUMINA_450K:
        probes = qualityMask450
    elif data_container.array_type == ArrayType.ILLUMINA_EPIC:
        probes = qualityMaskEPIC
    elif data_container.array_type == ArrayType.ILLUMINA_EPIC_PLUS:
        # this is a bit of a hack; probe names don't match epic, so I'm temporarily renaming, then filtering, then reverting.
        probes = qualityMaskEPICPLUS
    elif data_container.array_type == ArrayType.ILLUMINA_MOUSE:
        probes = qualityMaskmouse

    # v1.6+: the 1.0s are good probes and the 0.0 are probes to be excluded.
    cgs = pd.DataFrame( np.zeros((len(data_container.man.index), 1)), index=data_container.man.index, columns=['quality_mask'])
    cgs['quality_mask'] = 1.0
    snps = pd.DataFrame( np.zeros((len(data_container.snp_man.index), 1)), index=data_container.snp_man.index, columns=['quality_mask'])
    snps['quality_mask'] = 1.0
    df = pd.concat([cgs, snps])
    df.loc[df.index.isin(probes), 'quality_mask'] = 0
    #LOGGER.info(f"DEBUG quality_mask: {df.shape}, {df['quality_mask'].value_counts()} from {probes.shape} probes")
    return df


""" ##### DEPRECATED (<v1.5.0) #####
def _old_reprocess_noob_sesame_v144(container, offset=15, debug=False):
    ''' NOOB pythonized copy of https://github.com/zwdzwd/sesame/blob/master/R/background_correction.R

    - The function takes a SigSet and returns a modified SigSet with that background subtracted.
    - Background is modelled in a normal distribution and true signal in an exponential distribution.
    - The Norm-Exp deconvolution is parameterized using Out-Of-Band (oob) probes.
    - includes snps, but not control probes yet
    - output should replace the container instead of returning debug dataframes
    - II RED and II GREEN both have data, but manifest doesn't have a way to track this, so function tracks it.
    '''
    # get in-band red and green channel probe means
    #ibR <- c(IR(sset), II(sset)[,'U'])    # in-band red signal = IR_meth + IR_unmeth + II[unmeth]
    #ibG <- c(IG(sset), II(sset)[,'M'])    # in-band green signal = IG_meth + IG_unmeth + II[meth]
    # cols: mean_value, IlmnID, probe_type (I,II); index: illumina_id
    #CHECKED: AddressA or AddressB for each probe subtype matches probes.py
    raw = container.snp_methylated.data_frame
    snp_IR_meth = (raw[(raw['Infinium_Design_Type'] == 'I') & (raw['Color_Channel'] == 'Red')][['mean_value','AddressB_ID']]
                   .reset_index().rename(columns={'AddressB_ID':'illumina_id'}).set_index('illumina_id'))
    snp_IR_meth['Channel'] = 'Red'
    snp_IG_meth = (raw[(raw['Infinium_Design_Type'] == 'I') & (raw['Color_Channel'] == 'Grn')][['mean_value','AddressB_ID']]
                   .reset_index().rename(columns={'AddressB_ID':'illumina_id'}).set_index('illumina_id'))
    snp_IG_meth['Channel'] = 'Grn'
    snp_II_meth = (raw[(raw['Infinium_Design_Type'] == 'II')][['mean_value','AddressA_ID']]
                   .reset_index().rename(columns={'AddressA_ID':'illumina_id'}).set_index('illumina_id'))
    snp_II_meth['Channel'] = 'Grn'
    raw = container.snp_unmethylated.data_frame
    snp_IR_unmeth = (raw[(raw['Infinium_Design_Type'] == 'I') & (raw['Color_Channel'] == 'Red')][['mean_value','AddressA_ID']]
                   .reset_index().rename(columns={'AddressA_ID':'illumina_id'}).set_index('illumina_id'))
    snp_IR_unmeth['Channel'] = 'Red'
    snp_IG_unmeth = (raw[(raw['Infinium_Design_Type'] == 'I') & (raw['Color_Channel'] == 'Grn')][['mean_value','AddressA_ID']]
                   .reset_index().rename(columns={'AddressA_ID':'illumina_id'}).set_index('illumina_id'))
    snp_IG_unmeth['Channel'] = 'Grn'
    snp_II_unmeth = (raw[(raw['Infinium_Design_Type'] == 'II')][['mean_value','AddressA_ID']]
                   .reset_index().rename(columns={'AddressA_ID':'illumina_id'}).set_index('illumina_id'))
    snp_II_unmeth['Channel'] = 'Red'
    if debug:
        print('snp probes:', snp_IR_meth.shape, snp_IG_unmeth.shape, snp_II_meth.shape, snp_II_unmeth.shape)
    #--> copy over snps, but first get snps with illumina_id in index
    # swap index on all snps from IlmnID to illumina_id

    ## note: 350076 II + 89203 IR + 46298 IG = 485577 (including rs probes, but excl controls)
    ibG = container.fg_green # --> self.raw_dataset.get_fg_values(self.manifest, Channel.GREEN)
    ibG['Channel'] = 'Grn'
    ibG.index.name = 'illumina_id'
    ibR = container.fg_red # --> self.raw_dataset.get_fg_values(self.manifest, Channel.RED)
    ibR['Channel'] = 'Red'
    ibR.index.name = 'illumina_id'
    # to match sesame, extra probes are IR_unmeth and IG_unmeth in ibR red and ibG green, respectively.
    ibG = pd.concat([ibG,
                      snp_IG_meth,
                      snp_IG_unmeth,
                      snp_II_meth
                     ], sort=True).drop('probe_type', axis=1)
    # sort=True, because column order varies
    ibR = pd.concat([ibR,
                      snp_IR_meth,
                      snp_IR_unmeth,
                      snp_II_unmeth
                     ], sort=True).drop('probe_type', axis=1)
    if debug:
        print('in-bound Green:', ibG.shape) # green IG is AddressB, (meth) according to PROBE_SUBSETS
        print('in-bound Red:', ibR.shape) # red IR is AddressA (unmeth) according to PROBE_SUBSETS
        ### at this point, ibG ibR probe counts match sesame EXACTLY

    # set minimum intensity to 1
    ibR_affected = len(ibR.loc[ ibR['mean_value'] < 1 ].index)
    ibG_affected = len(ibG.loc[ ibG['mean_value'] < 1 ].index)
    ibR.loc[ ibR['mean_value'] < 1, 'mean_value'] = 1
    ibG.loc[ ibG['mean_value'] < 1, 'mean_value'] = 1
    if debug:
        print(f"IB: Set {ibR_affected} red and {ibG_affected} green to 1.0 ({len(ibR[ ibR['mean_value'] == 1 ].index)}, {len(ibG[ ibG['mean_value'] == 1 ].index)})")

    red_dupes = len(ibR.index)-len(ibR.drop_duplicates().index)
    grn_dupes = len(ibG.index)-len(ibG.drop_duplicates().index)
    if debug and (red_dupes or grn_dupes):
        print(f"duplicate probes: {red_dupes} red and {grn_dupes} green")

    ref = container.manifest.data_frame # [['Infinium_Design_Type','Color_Channel']]
    # using a copy .oobG and .oobR here; does not update the idat or other source data probe_means
    # adopted from raw_dataset.filter_oob_probes here
    oobR = (container.oobR.merge(container.manifest.data_frame[['AddressB_ID']],
                how='left',
                left_index=True,
                right_index=True)
            .reset_index()
            .rename(columns={'AddressB_ID':'illumina_id', 'Unnamed: 0': 'IlmnID'})
            .set_index('illumina_id')
           )
    oobR = pd.DataFrame(list(oobR['meth']) + list(oobR['unmeth']), columns=['mean_value'])
    oobG = (container.oobG.merge(container.manifest.data_frame[['AddressA_ID']],
                how='left',
                left_index=True,
                right_index=True)
            .reset_index()
            .rename(columns={'AddressA_ID':'illumina_id', 'Unnamed: 0': 'IlmnID'})
            .set_index('illumina_id')
           )
    oobG = pd.DataFrame(list(oobG['meth']) + list(oobG['unmeth']), columns=['mean_value'])

    oobG_affected = len(oobG[ oobG['mean_value'] < 1])
    oobG.loc[ oobG.mean_value < 1, 'mean_value'] = 1
    oobR_affected = len(oobR[ oobR['mean_value'] < 1])
    oobR.loc[ oobR.mean_value < 1, 'mean_value'] = 1

    # here: do bg_subtract AND normalization step here ...
    ## do background correction in each channel; returns "normalized in-band signal"
    ibR_nl, params_red   = normexp_bg_corrected(ibR, oobR, offset, sample_name=container.sample.name)
    #<- .backgroundCorrectionNoobCh1(ibR, oobR(sset), ctl(sset)$R, getBackgroundR(sset, bgR), offset=offset)
    ibG_nl, params_green = normexp_bg_corrected(ibG, oobG, offset, sample_name=container.sample.name)
    # <- .backgroundCorrectionNoobCh1(ibG, oobG(sset), ctl(sset)$G, getBackgroundG(sset, bgG), offset=offset)
    ibG_nl = ibG_nl.round({'bg_corrected':0})
    ibR_nl = ibR_nl.round({'bg_corrected':0})
    #print('ibG_nl', ibG_nl.shape)
    #print('ibR_nl', ibR_nl.shape)
    noob_green = ibG_nl
    noob_red = ibR_nl
    if debug:
        print(f"OOB: Set {oobR_affected} red and {oobG_affected} green to 1.0; shapes: {oobG.shape}, {oobR.shape}")
        print(f"noob_red with Grn: {len(noob_red[noob_red['Channel'] == 'Grn'])} noob_green with Red: {len(noob_green[noob_green['Channel'] == 'Red'])}")
        ref_IG = ref[(ref['Color_Channel']=='Grn') & (ref['Infinium_Design_Type']=='I')]
        ref_IR = ref[(ref['Color_Channel']=='Red') & (ref['Infinium_Design_Type']=='I')]
        ref_II = ref[ref['Infinium_Design_Type']=='II'] # II channel is NaN, but BOTH channels have data
        print(f"from manifest: ref_IG {ref_IG.shape} ref_IR {ref_IR.shape} ref_II {ref_II.shape}")

    # Combine and return red (IG + IR + II_unmeth) and green (IG + IR + II_meth)
    # ibR_nl has IlmnID and illumina_id (index); ref has IlmnID as index
    # ref_meth/ref_unmeth from probes.py
    ref_meth = pd.concat([
            ref[(ref['Color_Channel'].isna()) & (ref['Infinium_Design_Type']=='II')]['AddressA_ID'].reset_index().rename(columns={'AddressA_ID':'illumina_id'}),
            ref[(ref['Color_Channel']=='Grn') & (ref['Infinium_Design_Type']== 'I')]['AddressB_ID'].reset_index().rename(columns={'AddressB_ID':'illumina_id'}),
            ref[(ref['Color_Channel']=='Red') & (ref['Infinium_Design_Type']== 'I')]['AddressB_ID'].reset_index().rename(columns={'AddressB_ID':'illumina_id'}),
                             ]) #.set_index('illumina_id') # .drop('illumina_id', axis=1)
    ref_unmeth = pd.concat([
            ref[(ref['Color_Channel'].isna()) & (ref['Infinium_Design_Type']=='II')]['AddressA_ID'].reset_index().rename(columns={'AddressA_ID':'illumina_id'}),
            ref[(ref['Color_Channel']=='Grn') & (ref['Infinium_Design_Type']== 'I')]['AddressA_ID'].reset_index().rename(columns={'AddressA_ID':'illumina_id'}),
            ref[(ref['Color_Channel']=='Red') & (ref['Infinium_Design_Type']== 'I')]['AddressA_ID'].reset_index().rename(columns={'AddressA_ID':'illumina_id'}),
                             ]) #.set_index('illumina_id') # .drop('illumina_id', axis=1)
    noob_meth_G = noob_green[noob_green.index.isin(ref_meth['illumina_id'])]
    noob_unmeth_G = noob_green[noob_green.index.isin(ref_unmeth['illumina_id'])]
    noob_meth_R = noob_red[noob_red.index.isin(ref_meth['illumina_id'])]
    noob_unmeth_R = noob_red[noob_red.index.isin(ref_unmeth['illumina_id'])]
    noob_meth_dupes = pd.concat([noob_meth_G, noob_meth_R])
    noob_unmeth_dupes = pd.concat([noob_unmeth_G, noob_unmeth_R])
    # CONFIRMED: this dedupe method below matches sesame's output exactly for noob_meth
    noob_meth = (noob_meth_dupes[~noob_meth_dupes.index.duplicated(keep='first')]
                 .set_index('IlmnID')
                 .sort_index()
                 .rename(columns={'bg_corrected':'meth'})
                )
    # conveniently, the FIRST value of each duplicate probe appears to be the one we want for both meth/unmeth R/G channels
    noob_unmeth = (noob_unmeth_dupes[~noob_unmeth_dupes.index.duplicated(keep='first')]
                   .set_index('IlmnID')
                   .sort_index()
                   .rename(columns={'bg_corrected':'unmeth'})
                  )

    # update II, IG, IR, oobR, oobG, ctrl_red, ctrl_green
    # --> --> probes.py subsets concatenate these:
    # fg_green
    #   GREEN + AddressA + II
    #   GREEN + AddressA + IG
    #   GREEN + AddressB + IG
    # oob_green
    #   RED   + AddressA + IR
    # fg_red
    #   RED   + AddressA + II
    #   RED   + AddressA + IR
    #   RED   + AddressB + IR
    # oob_red
    #   GREEN + AddressB + IG
    #
    # methylated
    #   GREEN + AddressA + II
    #   GREEN + AddressB + I
    #   RED   + AddressB + I
    # unmethylated
    #   RED   + AddressA + II
    #   GREEN + AddressA + I
    #   RED   + AddressA + I
    # RETROFITTING BELOW -- may not work, as sesame works with noob_meth / noob_unmeth instead

    try:
        container.methylated.set_bg_corrected(noob_green, noob_red)
        container.unmethylated.set_bg_corrected(noob_green, noob_red)
        container.methylated.set_noob(1.0)
        container.unmethylated.set_noob(1.0)
    except ValueError as e:
        print(e)
        if debug:
            LOGGER.warning("could not update container methylated / unmethylated noob values, because preprocess_sesame_noob has already run once.")

    # output df should have sample meth or unmeth in a column with sample name and IlmnID as index. 485512 rows
    if debug:
        return {
            'noob_meth': noob_meth,
            'noob_unmeth': noob_unmeth,
            'oobR': oobR,
            'oobG': oobG,
            'noob_green': noob_green,
            'noob_red': noob_red,
            'dupe_meth': noob_meth_dupes,
            'dupe_unmeth': noob_unmeth_dupes,
        }
    return # noob_meth, noob_unmeth
"""
