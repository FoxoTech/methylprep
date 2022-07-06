try:
    from methylcheck import mean_beta_compare
except ImportError:
    mean_beta_compare = None
import types
from scipy import stats
from statsmodels.distributions.empirical_distribution import ECDF
import pandas as pd
import numpy as np


def _pval_sesame_preprocess(data_container, combine_neg=True):
    """Performs p-value detection of low signal/noise probes. This ONE SAMPLE version uses meth/unmeth before it is contructed into a _SampleDataContainer__data_frame.
    - returns a dataframe of probes and their detected p-value levels.
    - this will be saved to the csv output, so it can be used to drop probes at later step.
    - output: index are probes (IlmnID or illumina_id); one column [poobah_pval] contains the sample p-values.
    - called by pipeline CLI --poobah option.
    - confirmed that this version produces identical results to the pre-v1.5.0 version on 2021-06-16
    """
    if data_container.debug == True: 
        print("DEBUG: running 1.6.2 poobah method instead")
        # 2021-03-22 assumed 'mean_value' for red and green MEANT meth and unmeth (OOBS), respectively.
        # oob[G/R]['Unmeth'] is the out of band signal from the probe capturing the unmethylated state
        # oob[G/R]['Meth'] is the out of band signal from the probe capturing the methylated state
        funcG = ECDF(data_container.oobG['Unmeth'].values)
        funcR = ECDF(data_container.oobR['Meth'].values)
    else:
        bgG = (
            list(data_container.oobG['Unmeth'].values)
            + list(data_container.oobG['Meth'].values)
        )
        bgR = (
            list(data_container.oobR['Unmeth'].values)
            + list(data_container.oobR['Meth'].values)
        )
        # SeSAMe by default includes negative controls as a part of the background green and red intensities
        if combine_neg:
            # Add green signal from negative controls to background green intensities
            dfG = data_container.ctrl_green
            dfG = dfG[dfG['Control_Type']=='NEGATIVE']
            dfG = dfG[['Extended_Type','mean_value']]
            bgG += list(dfG['mean_value'].values)
            # Add reg signal from negative controls to background red intensitites
            dfR = data_container.ctrl_red
            dfR = dfR[dfR['Control_Type']=='NEGATIVE']
            dfR = dfR[['Extended_Type','mean_value']]
            bgR += list(dfR['mean_value'].values)
        # Build empirical cumulative distribution functions
        funcG = ECDF(bgG)
        funcR = ECDF(bgR)

    # Apply function of background red intensity to red probes and background green intensity to green probes
    pIR = pd.DataFrame(
        index=data_container.IR.index,
        data=1-np.maximum(funcR(data_container.IR['Meth']), funcR(data_container.IR['Unmeth'])),
        columns=['poobah_pval'])
    pIG = pd.DataFrame(
        index=data_container.IG.index,
        data=1-np.maximum(funcG(data_container.IG['Meth']), funcG(data_container.IG['Unmeth'])),
        columns=['poobah_pval'])
    pII = pd.DataFrame(
        index=data_container.II.index,
        data=1-np.maximum(funcG(data_container.II['Meth']), funcR(data_container.II['Unmeth'])),
        columns=['poobah_pval'])
    # pval output: index is IlmnID; and threre's one column, 'poobah_pval' with p-values
    pval = pd.concat([pIR,pIG,pII])
    return pval

def _pval_neg_ecdf(data_container):
    dfR = data_container.ctrl_red
    dfR = dfR[dfR['Control_Type']=='NEGATIVE']
    dfR = dfR[['Extended_Type','mean_value']]
    funcR = ECDF(dfR['mean_value'].values)
    dfG = data_container.ctrl_green
    dfG = dfG[dfG['Control_Type']=='NEGATIVE']
    dfG = dfG[['Extended_Type','mean_value']]
    funcG = ECDF(dfG['mean_value'].values)
    pIR = pd.DataFrame(
        index=data_container.IR.index,
        data=1-np.maximum(funcR(data_container.IR['Meth']), funcR(data_container.IR['Unmeth'])),
        columns=['pNegECDF_pval'])
    pIG = pd.DataFrame(
        index=data_container.IG.index,
        data=1-np.maximum(funcG(data_container.IG['Meth']), funcG(data_container.IG['Unmeth'])),
        columns=['pNegECDF_pval'])
    pII = pd.DataFrame(
        index=data_container.II.index,
        data=1-np.maximum(funcG(data_container.II['Meth']), funcR(data_container.II['Unmeth'])),
        columns=['pNegECDF_pval'])
    # concat and sort
    pval = pd.concat([pIR, pIG, pII])
    return pval



""" DEPRECATED FUNCTIONS (<v1.5.0)
def detect_probes(data_containers, method='sesame', save=False, silent=True):
    '''
About:
    a wrapper for the p-value probe detection methods. Tries to check inputs and rationalize them
    with methyl-suite's standard data objects.

Inputs:
    a list of sample data_containers. The dataframe must have a 'IlMnID' for the index of probe names.
    And probe names should be `cgxxxxxx` format to work with other methylcheck functions

    To create this, use:

    data_containers = methylprep.run_pipeline(data_dir,
            save_uncorrected=True,
            betas=True)
    (if there are more than 200 samples, you'll need to load the data from disk instead, as nothing will be returned)

    method:
        sesame -- use the sesame ported algorithm
        minfi -- use the minfi ported algorithm

Checks:
    data_containers must have 'meth' and 'unmeth' columns (uncorrected data)
    the values for these columns should be between 0 and 10000s
    (beta: 0 to 1; m_value: -neg to pos+ range)

Returns:
    dataframe of p-value filtered probes
    '''
    if not ('unmeth' in data_containers[0]._SampleDataContainer__data_frame.columns and
            'meth' in data_containers[0]._SampleDataContainer__data_frame.columns):
            raise ValueError("Provide a list of data_containers that includes uncorrected data (with 'meth' and 'unmeth' columns, using the 'save_uncorrected' option in run_pipeline)")
    if method == 'minfi':
        pval = pval_minfi(data_containers)
    else:
        pval = pval_sesame(data_containers)

    if silent == False and type(mean_beta_compare) is types.FunctionType:
        # plot it
        # df1 and df2 are probe X sample_id matrices
        mean_beta_compare(df1, df2, save=save, verbose=False, silent=silent)
    return pval


def _pval_minfi(data_containers):
    # negative control p-value
    # Pull M and U values
    meth = pd.DataFrame(data_containers[0]._SampleDataContainer__data_frame.index)
    unmeth = pd.DataFrame(data_containers[0]._SampleDataContainer__data_frame.index)

    for i,c in enumerate(data_containers):
        sample = data_containers[i].sample
        m = c._SampleDataContainer__data_frame.rename(columns={'meth':sample})
        u = c._SampleDataContainer__data_frame.rename(columns={'unmeth':sample})
        meth = pd.merge(left=meth,right=m[sample],left_on='IlmnID',right_on='IlmnID',)
        unmeth = pd.merge(left=unmeth,right=u[sample],left_on='IlmnID',right_on='IlmnID')

    # Create empty dataframes for red and green negative controls
    negctlsR = pd.DataFrame(data_containers[0].ctrl_red['Extended_Type'])
    negctlsG = pd.DataFrame(data_containers[0].ctrl_green['Extended_Type'])

    # Fill red and green dataframes
    for i,c in enumerate(data_containers):
        sample = str(data_containers[i].sample)
        dfR = c.ctrl_red
        dfR = dfR[dfR['Control_Type']=='NEGATIVE']
        dfR = dfR[['Extended_Type','mean_value']].rename(columns={'mean_value':sample})
        dfG = c.ctrl_green
        dfG = dfG[dfG['Control_Type']=='NEGATIVE']
        dfG = dfG[['Extended_Type','mean_value']].rename(columns={'mean_value':sample})
        negctlsR = pd.merge(left=negctlsR,right=dfR,on='Extended_Type')
        negctlsG = pd.merge(left=negctlsG,right=dfG,on='Extended_Type')

    # Reset index on dataframes
    negctlsG = negctlsG.set_index('Extended_Type')
    negctlsR = negctlsR.set_index('Extended_Type')

    # Get M and U values for IG, IR and II

    # first pull out sections of manifest (will be used to identify which probes belong to each IG, IR, II)
    manifest = data_containers[0].manifest.data_frame[['Infinium_Design_Type','Color_Channel']]
    IG = manifest[(manifest['Color_Channel']=='Grn') & (manifest['Infinium_Design_Type']=='I')]
    IR = manifest[(manifest['Color_Channel']=='Red') & (manifest['Infinium_Design_Type']=='I')]
    II = manifest[manifest['Infinium_Design_Type']=='II']

    # second merge with meth and unmeth dataframes
    IG_meth = pd.merge(left=IG,right=meth,on='IlmnID').drop(columns=['Infinium_Design_Type','Color_Channel']).set_index('IlmnID')
    IG_unmeth = pd.merge(left=IG,right=unmeth,on='IlmnID').drop(columns=['Infinium_Design_Type','Color_Channel']).set_index('IlmnID')
    IR_meth = pd.merge(left=IR,right=meth,on='IlmnID').drop(columns=['Infinium_Design_Type','Color_Channel']).set_index('IlmnID')
    IR_unmeth = pd.merge(left=IR,right=unmeth,on='IlmnID').drop(columns=['Infinium_Design_Type','Color_Channel']).set_index('IlmnID')
    II_meth = pd.merge(left=II,right=meth,on='IlmnID').drop(columns=['Infinium_Design_Type','Color_Channel']).set_index('IlmnID')
    II_unmeth = pd.merge(left=II,right=unmeth,on='IlmnID').drop(columns=['Infinium_Design_Type','Color_Channel']).set_index('IlmnID')

    # Calcuate parameters
    sdG = stats.median_absolute_deviation(negctlsG)
    muG = np.median(negctlsG,axis=0)
    sdR = stats.median_absolute_deviation(negctlsR)
    muR = np.median(negctlsR,axis=0)

    # calculate p values for type 1 Red
    pIR = pd.DataFrame(index=IR_meth.index,
                   data=1 - stats.norm.cdf(IR_meth+IR_unmeth,2*muR,2*sdR),
                   columns=IR_meth.columns)

    # calculate p values for type 1 Green
    pIG = pd.DataFrame(index=IG_meth.index,
                   data=1 - stats.norm.cdf(IG_meth+IG_unmeth,2*muG,2*sdG),
                   columns=IG_meth.columns)

    # calculat4e p values for type II
    pII = pd.DataFrame(index=II_meth.index,
                  data=1-stats.norm.cdf(II_meth+II_unmeth,muR+muG,sdR+sdG),
                  columns=II_meth.columns)
    # concat and sort
    pval = pd.concat([pIR, pIG, pII])
    pval = pval.sort_values(by='IlmnID')
    return pval

def _pval_sesame(data_containers):
    # pOOHBah, called using ___ in sesame
    # Pull M and U values
    meth = pd.DataFrame(data_containers[0]._SampleDataContainer__data_frame.index)
    unmeth = pd.DataFrame(data_containers[0]._SampleDataContainer__data_frame.index)

    for i,c in enumerate(data_containers):
        sample = data_containers[i].sample
        m = c._SampleDataContainer__data_frame.rename(columns={'meth':sample})
        u = c._SampleDataContainer__data_frame.rename(columns={'unmeth':sample})
        meth = pd.merge(left=meth, right=m[sample], left_on='IlmnID', right_on='IlmnID',)
        unmeth = pd.merge(left=unmeth, right=u[sample], left_on='IlmnID', right_on='IlmnID')

    # Separate M and U values for IG, IR and II
    # first pull out sections of manifest (will be used to identify which probes belong to each IG, IR, II)
    manifest = data_containers[0].manifest.data_frame[['Infinium_Design_Type','Color_Channel']]
    IG = manifest[(manifest['Color_Channel']=='Grn') & (manifest['Infinium_Design_Type']=='I')]
    IR = manifest[(manifest['Color_Channel']=='Red') & (manifest['Infinium_Design_Type']=='I')]
    II = manifest[manifest['Infinium_Design_Type']=='II']

    # second merge with meth and unmeth dataframes
    IG_meth = pd.merge(left=IG,right=meth,on='IlmnID').drop(columns=['Infinium_Design_Type','Color_Channel']).set_index('IlmnID')
    IG_unmeth = pd.merge(left=IG,right=unmeth,on='IlmnID').drop(columns=['Infinium_Design_Type','Color_Channel']).set_index('IlmnID')
    IR_meth = pd.merge(left=IR,right=meth,on='IlmnID').drop(columns=['Infinium_Design_Type','Color_Channel']).set_index('IlmnID')
    IR_unmeth = pd.merge(left=IR,right=unmeth,on='IlmnID').drop(columns=['Infinium_Design_Type','Color_Channel']).set_index('IlmnID')
    II_meth = pd.merge(left=II,right=meth,on='IlmnID').drop(columns=['Infinium_Design_Type','Color_Channel']).set_index('IlmnID')
    II_unmeth = pd.merge(left=II,right=unmeth,on='IlmnID').drop(columns=['Infinium_Design_Type','Color_Channel']).set_index('IlmnID')

    pval = pd.DataFrame(data=manifest.index, columns=['IlmnID'])
    for i,c in enumerate(data_containers):
        # 2021-03-22 assumed 'mean_value' for red and green MEANT meth and unmeth (OOBS), respectively.
        funcG = ECDF(data_containers[i].oobG['unmeth'].values)
        funcR = ECDF(data_containers[i].oobR['meth'].values)
        sample = data_containers[i].sample
        pIR = pd.DataFrame(index=IR_meth.index,data=1-np.maximum(funcR(IR_meth[sample]), funcR(IR_unmeth[sample])),columns=[sample])
        pIG = pd.DataFrame(index=IG_meth.index,data=1-np.maximum(funcG(IG_meth[sample]), funcG(IG_unmeth[sample])),columns=[sample])
        pII = pd.DataFrame(index=II_meth.index,data=1-np.maximum(funcG(II_meth[sample]), funcR(II_unmeth[sample])),columns=[sample])
        p = pd.concat([pIR,pIG,pII]).reset_index()
        pval = pd.merge(pval,p)
    return pval
"""
