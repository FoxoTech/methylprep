import sys
import numpy as np
import pandas as pd
from pathlib import Path
#import matplotlib.pyplot as plt

# App
import methylprep

#patching
try:
    # python 3.4+ should use builtin unittest.mock not mock package
    from unittest.mock import patch
except ImportError:
    from mock import patch

#LOCAL = Path('/Volumes/LEGX/GSE69852/idats/') # run_sesame.R drops files in here.
LOCAL = Path('docs/example_data/GSE69852/')


def sesame_convert(filename, drop_rs=True):
    """ generic converter, sesame to mprep dataframe, assuming II,IG,IR were rbind() into one df then CSV'ed."""
    df = pd.read_csv(filename).rename(columns={'Unnamed: 0':'IlmnID'}).set_index('IlmnID')
    # separate into meth and unmeth dataframes
    meth_columns = [c for c in df.columns if '.M' in c]
    unmeth_columns = [c for c in df.columns if '.U' in c]
    meth_df = df[meth_columns]
    unmeth_df = df[unmeth_columns]
    # drop the leading X and trailing .M or .U from sample names here
    meth_sample_dict = dict(zip(meth_df.columns, [c[1:-2] for c in meth_df.columns]))
    unmeth_sample_dict = dict(zip(unmeth_df.columns, [c[1:-2] for c in unmeth_df.columns]))
    sesame_meth = meth_df.rename(columns=meth_sample_dict)
    sesame_unmeth = unmeth_df.rename(columns=unmeth_sample_dict)
    # sort probes and columns and force floats
    sesame_meth = sesame_meth.sort_index().reindex(sorted(sesame_meth.columns), axis=1).astype(float)
    sesame_unmeth = sesame_unmeth.sort_index().reindex(sorted(sesame_unmeth.columns), axis=1).astype(float)
    # drop SNPs, so same as methylprep
    if drop_rs:
        sesame_meth = sesame_meth[~sesame_meth.index.str.startswith('rs')] # dropping SNP probes, which are included in sesame output
        sesame_unmeth = sesame_unmeth[~sesame_unmeth.index.str.startswith('rs')]
    # dimensions
    print(sesame_meth.shape,sesame_unmeth.shape)
    return sesame_meth, sesame_unmeth


def test_compare_methylprep_sesame__raw_oob_noob_dye():
    """ this test confirms that these values match closely with sesame:
    - raw meth, raw unmeth
    - oobG, oobR and preprocess_noob_sesame oobG/oobR
    - noob values
    - noob + dye corected values

    Notes:
    - it is clear that 0.001% of raw meth/unmeth values did not match for THIS test sample (GSE69852),
        corresponding to 77 meth and 308 unmeth cg probes. They are WAY off.
        These are 308 R2G and 77 G2R infer channel switched probes. They must not be switched in the sesame ref data.
    - the final dye meth values vary in the range of -3 to +4,
        and dye unmeth values varied in range of -15 to +12, but MEAN was <1."""
    # does not include rs probes, because .methylated DF does not include them for comparison.

    #containers = methylprep.run_pipeline(LOCAL, sesame=True, save_uncorrected=True, low_memory=False, debug=False)
    samp = 1 # '9247377085_R04C02'
    steps = ['infer_channel_switch', 'quality_mask', 'noob', 'dye_bias'] # no poobah 'poobah'
    containers = methylprep.processing.pipeline.make_pipeline(LOCAL, steps=steps, exports=None, estimator=None, low_memory=False, debug=False)

    # will run noob and dye and sesame==True
    s_meth, s_unmeth           = sesame_convert(Path(LOCAL,'sesame_raw.csv'), drop_rs=True)
    #s_noob_meth, s_noob_unmeth = sesame_convert(Path(LOCAL,'sesame_noob.csv'), drop_rs=True) # bg_sub + noob applied to raw intensities

    # after adding detectmask before noob, the structure changed! -- below is for this
    s_noob_df = pd.read_csv(Path(LOCAL,'sesame_noob.csv')).set_index('Unnamed: 0').sort_index()
    s_noob_df = s_noob_df[~s_noob_df.index.str.startswith('rs')]
    s_noob_df.index.name = 'IlmnID'
    s_noob_meth = s_noob_df[['M']].rename(columns={'M':'9247377085_R04C02'})
    s_noob_unmeth = s_noob_df[['U']].rename(columns={'U':'9247377085_R04C02'})


    m_meth = containers[samp].methylated.data_frame.sort_index()[['mean_value']]
    m_unmeth = containers[samp].unmethylated.data_frame.sort_index()[['mean_value']]
    #(s_meth['9247377085_R04C02'] - m_meth['mean_value']).plot.hist(bins=100)
    raw_meth_mean_diff = (s_meth['9247377085_R04C02'] - m_meth['mean_value']).mean()
    raw_unmeth_mean_diff = (s_unmeth['9247377085_R04C02'] - m_unmeth['mean_value']).mean()
    raw_meth_isclose = sum(np.isclose( s_meth['9247377085_R04C02'], m_meth['mean_value'], atol=1.0))/len(s_meth)
    raw_unmeth_isclose = sum(np.isclose(s_unmeth['9247377085_R04C02'], m_unmeth['mean_value'], atol=1.0))/len(s_unmeth)
    if not (-1 < raw_meth_mean_diff < 1) or not (-1 < raw_unmeth_mean_diff < 1) or raw_meth_isclose < 0.99 or raw_unmeth_isclose < 0.99:
        raise AssertionError(f"raw meth/unmeth values differ between sesame and methylprep (expect exact match on 99% or more): meth: ({raw_meth_mean_diff} {raw_meth_isclose}) unmeth: ({raw_unmeth_mean_diff} {raw_unmeth_isclose})")
    print("preprocess_sesame() raw meth/unmeth match; mean_diff +/- 1.0 unit {raw_meth_isclose} {raw_unmeth_isclose}")

    # note that methylated data_frame noob will match sesame, but SDC will be dye-corrected noob.
    noob_meth = containers[samp].methylated.data_frame[['bg_corrected']].sort_index()
    noob_unmeth = containers[samp].unmethylated.data_frame[['bg_corrected']].sort_index()
    #meth_matches = all(data['noob_meth'][['meth']].rename(columns={'meth':'9247377085_R04C02'}) == s_noob_meth[['9247377085_R04C02']])
    meth_matches = all(noob_meth[['bg_corrected']].rename(columns={'bg_corrected':'9247377085_R04C02'}) == s_noob_meth[['9247377085_R04C02']])
    print( 'noob meth matches:', meth_matches )
    #unmeth_matches = all(data['noob_unmeth'][['unmeth']].rename(columns={'unmeth':'9247377085_R04C02'}) == s_noob_unmeth[['9247377085_R04C02']])
    unmeth_matches = all(noob_unmeth[['bg_corrected']].rename(columns={'bg_corrected':'9247377085_R04C02'}) == s_noob_unmeth[['9247377085_R04C02']])
    print( 'noob unmeth matches:', unmeth_matches )
    print("preprocess_sesame() noob output matches sesame exactly.")

    data = methylprep.processing.preprocess.preprocess_noob_sesame(containers[samp], debug=True)
    ref = containers[samp].manifest.data_frame

    s_oobG = pd.read_csv(Path(LOCAL,'sesame_oobG.csv')).set_index('Unnamed: 0')
    s_oobG = (s_oobG.merge(ref[['AddressA_ID']], how='left', left_index=True, right_index=True)
            .reset_index().rename(columns={'AddressA_ID':'illumina_id', 'Unnamed: 0': 'IlmnID'})
            .set_index('illumina_id').sort_values('IlmnID')
           )
    s_oobG = pd.DataFrame(list(s_oobG['M']) + list(s_oobG['U']), columns=['mean_value'])
    oobG_matches = sorted(s_oobG) == sorted(data['oobG'])

    s_oobR = pd.read_csv(Path(LOCAL,'sesame_oobR.csv')).set_index('Unnamed: 0')
    s_oobR = (s_oobR.merge(ref[['AddressB_ID']], how='left', left_index=True, right_index=True)
            .reset_index().rename(columns={'AddressB_ID':'illumina_id', 'Unnamed: 0': 'IlmnID'})
            .set_index('illumina_id')
           )
    s_oobR = pd.DataFrame(list(s_oobR['U']) + list(s_oobR['M']), columns=['mean_value'])
    oobR_matches = sorted(s_oobR) == sorted(data['oobR'])

    print('oobG_matches', oobG_matches)
    if oobG_matches is False or oobR_matches is False or meth_matches is False or unmeth_matches is False:
        raise AssertionError("oob or noob values don't match")

    # LAST, compare dye results (from SDC.noob_meth and .noob_unmeth) vs sesame loaded dye+noob values
    sesame_dye = pd.read_csv(Path(LOCAL, 'sesame_noob_dye.csv')).set_index('Unnamed: 0').sort_index()
    # 9247377085_R04C02
    s_dye_meth = sesame_dye[['M']][~sesame_dye.index.str.startswith('rs')].rename(columns={'M':'meth'})
    m_dye_meth = containers[samp].methylated.data_frame.sort_index()[['noob']].rename(columns={'noob':'meth'})
    #(s_dye_meth - m_dye_meth).plot.hist(bins=500,xlim=(-50,50))
    #plt.show()
    dye_meth_diff = float((s_dye_meth - m_dye_meth).mean())

    s_dye_unmeth = sesame_dye[['U']][~sesame_dye.index.str.startswith('rs')].rename(columns={'U':'unmeth'})
    m_dye_unmeth = containers[samp].unmethylated.data_frame.sort_index()[['noob']].rename(columns={'noob':'unmeth'})
    #(a - b).plot.hist(bins=500,xlim=(-50,50), ylim=(0,500000))
    #plt.show()
    dye_unmeth_diff = float((s_dye_unmeth - m_dye_unmeth).mean())
    if dye_meth_diff > 1.0 or dye_unmeth_diff > 1.1: # manual testing found a mean diff of ~0.5 and ~1.03 for each.
        #print(f"WARNING noob + dye corrected values don't match in sesame vs methylprep: meth {dye_meth_diff} unmeth {dye_unmeth_diff} (MEAN)")
        raise AssertionError(f"noob + dye corrected values don't match in sesame vs methylprep: meth {dye_meth_diff} unmeth {dye_unmeth_diff} (MEAN)")

    m_betas = (pd.DataFrame(containers[samp]._SampleDataContainer__data_frame['beta_value'])
               .sort_index()
               .rename(columns={'beta_value':'9247377085_R04C02'})
              )
    s_betas = (pd.read_csv(Path(LOCAL, 'sesame_noob_dye_betas.csv'))
               .set_index('Unnamed: 0')
               .rename(columns=(lambda x: x[1:]))
               .sort_index()
              )
    s_betas = s_betas[~s_betas.index.str.startswith('rs')]
    #(s_betas[['9247377085_R04C02']] - m_betas).plot.hist(bins=300, xlim=(-0.01,0.01))
    #plt.show()
    beta_diff = float((s_betas[['9247377085_R04C02']] - m_betas).mean())
    if not (-0.001 < beta_diff < 0.001):
        # got -0.001037 beta_diff on 3/30/2021, so barely failed, with a noob+dye diff of -4 or +8.
        raise AssertionError(f"betas values differ (sesame/methylprep) NOOB + DYE + BETA mean diff: {beta_diff}")


def test_open_sesame_betas_vs_methylprep():
    """ simplest test: does the openSesame beta output match the run_pipeline beta output? """
    LOCAL = Path('docs/example_data/GSE69852/')
    m_betas= methylprep.run_pipeline(LOCAL, betas=True, sample_name=['FetalLiver1'])
    m_betas = m_betas.sort_index()
    # columns={'beta_value':'9247377085_R04C02'}
    s_betas = (pd.read_csv(Path(LOCAL, 'sesame_open_betas.csv'))
               .set_index('Unnamed: 0')
               .rename(columns=(lambda x: x[1:]))
               .sort_index()
              )
    s_betas.index.name = 'IlmnID'
    s_betas = s_betas[~s_betas.index.str.startswith('rs')]
    s_betas = s_betas[['247377085_R04C02']]
    s_betas = s_betas.rename(columns={'247377085_R04C02':'9247377085_R04C02'}) # opensesame had wrong sentrix ID, missing 9, for some reason
    print((m_betas - s_betas).mean())
    #(m_betas - s_betas).plot.hist(bins=200)
    assert abs( float((m_betas - s_betas).mean()) ) < 0.002 # actual is 0.001023


""" how I created the sesame reference datasets

library(BiocManager)
library(sesame)
library(preprocessCore)
in_dir = paste0('/Users/mmaxmeister/methylprep/docs/example_data/GSE69852/')

ssets <- lapply(searchIDATprefixes(in_dir), readIDATpair)
ssets.noob <- lapply(ssets, noob)
ssets <- lapply(ssets, detectionMask) # poobah
ssets.dye <- lapply(ssets.noob, dyeBiasCorrTypeINorm)

ssets.noob.df = data.frame(rbind( IG(ssets.noob[[1]]), IR(ssets.noob[[1]]), II(ssets.noob[[1]]) ))
write.csv(file= file.path(in_dir, 'sesame_noob.csv'), x=ssets.noob.df, row.names=TRUE)

ssets.dye.df = data.frame(rbind( IG(ssets.dye[[1]]), IR(ssets.dye[[1]]), II(ssets.dye[[1]]) ))
write.csv(file= file.path(in_dir, 'sesame_noob_dye.csv'), x=ssets.dye.df, row.names=TRUE)

ssets.betas = lapply(ssets.dye, getBetas)
write.csv(file= file.path(in_dir, 'sesame_noob_dye_betas.csv'), x=ssets.betas, row.names=TRUE)

"""
