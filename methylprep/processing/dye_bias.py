# Lib
import logging
import math
import numpy as np
import pandas as pd
import scipy
# App
import methylprep

__all__ = ['nonlinear_dye_bias_correction']

LOGGER = logging.getLogger(__name__)


def get_ranks(x):
    """ get_ranks - Python version of the C function get_ranks() --- part of qnorm_using_target """
    i = 0
    n = len(x)
    rank = np.zeros(n)
    while (i < n):
        j = i
        while ((j < n - 1) and (x[j]['data'] == x[j + 1]['data'])):
            j += 1

        if (i != j):
            for k in range(i, j+1):
                rank[k] = (i + j + 2) / 2.0
        else:
            rank[i] = i + 1
        i = j + 1

    return rank

def qnorm_using_target(data, target):
    """ using_target - Python version of the C function using_target() ;
    data and target in must be ndarray like np.transpose(np.array([IR1]))"""
    nrows = data.shape[0]
    ncols = data.shape[1]
    targetrows = target.shape[0]
    float_eps = np.finfo(np.float32).eps

    if nrows != targetrows:
        raise NotImplementedError('Data and target are different lengths')
    else:
        for j in range(ncols):
            non_na = 0
            dimat = []
            for i in range(nrows):
                if ~np.isnan(data[i, j]):
                    dataitem = {'data': data[i, j], 'rank': i}
                    dimat.append(dataitem)
                    non_na += 1

            dimat = sorted(dimat, key=lambda k: k['data'])
            ranks = get_ranks(dimat)
            if non_na == nrows:
                for i in range(nrows):
                    ind = dimat[i]['rank']
                    if (ranks[i] - math.floor(ranks[i])) > 0.4:
                        data[ind, j] = 0.5*(target[int(math.floor(ranks[i])-1)] + target[int(math.floor(ranks[i]))])
                    else:
                        data[ind, j] = target[int(math.floor(ranks[i])-1)]
            else:
                for i in range(non_na):
                    samplepercentile = float(ranks[i] - 1)/float(non_na - 1)
                    target_ind_double = 1.0 + (float(targetrows) - 1.0) * samplepercentile
                    target_ind_double_floor = math.floor(target_ind_double + 4*float_eps)
                    target_ind_double = target_ind_double - target_ind_double_floor
                    if (math.fabs(target_ind_double) <= 4*float_eps):
                        target_ind_double = 0.0

                    if (target_ind_double  == 0.0):
                        target_ind = int(math.floor(target_ind_double_floor + 0.5))
                        ind = dimat[i]['rank']
                        data[ind, j] = target[target_ind-1]
                    elif (target_ind_double == 1.0):
                        target_ind = int(math.floor(target_ind_double_floor + 1.5))
                        ind = dimat[i]['rank']
                        data[ind, j] = target[target_ind-1]
                    else:
                        target_ind = int(math.floor(target_ind_double_floor + 0.5))
                        ind = dimat[i]['rank']
                        if ((target_ind < targetrows) and (target_ind > 0)):
                             data[ind, j] = (1.0- target_ind_double)*target[target_ind-1] + target_ind_double*target[target_ind]
                        elif (target_ind >= targetrows):
                             data[ind, j] = target[targetrows-1]
                        else:
                            data[ind, j] = target[0]
    # assuming I only need to return a single column here
    return np.transpose(data)[0]


def nonlinear_dye_bias_correction(container, debug=False):
    """ transforms Red and Green probe intensities to better align with each other.
    - equivalent to sesame's dyeBiasCorrTypeINorm function
    - function order: read-idats --> poobah --> noob --> dye-bias (before beta or m-values calculated)
    - overwrites container probe values with new ones
        - container.methylated.data_frame | container.unmethylated.data_frame
            - update [mean_value | bg_corrected | noob] in each
            - update green_idat and red_idat probe_values?
        - container._SampleDataContainer__data_frame
            - container.II
            - container.IG
            - container.IR

    - does not change SNPs or control probes
    - SampleDataContainer will pass in noob or raw values, depending on `do_noob`
        ... but columns will always be named noob_...

    """
    container._SigSet__dye_bias_corrected = False # sets to true when successful; otherwise, will run linear correction if sample fails.

    if not isinstance(container, methylprep.processing.SampleDataContainer):
        raise TypeError("You must provide a sample data container object.")
    if debug:
        import matplotlib.pyplot as plt # not required by package for normal users

    # get the IG & IR probes that pass the pvalue qualityMask; drops failed probes
    if 'poobah_pval' in container._SampleDataContainer__data_frame.columns:
        mask = (container._SampleDataContainer__data_frame['poobah_pval'] < container.poobah_sig)
        if mask.index.duplicated().sum() > 0:
            # equivalent to len(mask.index) > len(set(mask.index))
            LOGGER.info("Duplicate probe names found; switching to linear-dye correction.")
            mask = None
            print(f'DEBUG dupes IR: {container.IR.index.duplicated().sum()} IG: {container.IG.index.duplicated().sum()}')
            return container
    else:
        mask = None # fetches everything

    # dye-correct NOOB or RAW intensities, depending on preprocessing flags here.
    columns = {'noob_Meth':'Meth','noob_Unmeth':'Unmeth'} if container.do_noob == True else {'Meth':'Meth','Unmeth':'Unmeth'}
    drop_columns = ['Meth', 'Unmeth', 'poobah_pval', 'used', 'AddressA_ID', 'AddressB_ID'] if container.do_noob == True else ['noob_Meth', 'noob_Unmeth', 'poobah_pval', 'used', 'AddressA_ID', 'AddressB_ID']
    if container.pval is False:
        drop_columns.remove('poobah_pval')

    if isinstance(mask,pd.Series):
        sub_mask = mask[mask.index.isin(container.IG.index)]
        IG0 = container.IG.join(sub_mask, how='inner')
        IG0 = IG0[IG0['poobah_pval'] == True].drop(columns=drop_columns).rename(columns=columns).sort_index()

        sub_mask = mask[mask.index.isin(container.IR.index)]
        IR0 = container.IR.join(sub_mask, how='inner')
        IR0 = IR0[IR0['poobah_pval'] == True].drop(columns=drop_columns).rename(columns=columns).sort_index()
    else:
        IG0 = container.IG.copy().drop(columns=drop_columns).rename(columns=columns).sort_index()
        IR0 = container.IR.copy().drop(columns=drop_columns).rename(columns=columns).sort_index()

    # IG/IR includes snps now
    #IR0 = pd.concat( [IR0, container.snp_IR.rename(columns={'meth':'noob_meth', 'unmeth':'noob_unmeth'})] ).sort_index()
    #IG0 = pd.concat( [IG0, container.snp_IG.rename(columns={'meth':'noob_meth', 'unmeth':'noob_unmeth'})] ).sort_index()

    if debug:
        pd.options.mode.chained_assignment = 'raise' # only needed during debug
        if isinstance(mask, pd.Series):
            print(f"pval mask probes passing overall: {round(100*mask.sum()/len(mask),2)}%")
            print(f"Usable probes; IG0: {len(IG0)} of {len(container.IG)} ({round(100*len(IG0)/len(container.IG),2)}%) IR0: {len(IR0)} of {len(container.IR)} ({round(100*len(IR0)/len(container.IR),2)}%)")

    maxIG = np.nanmax(IG0); minIG = np.nanmin(IG0)
    maxIR = np.nanmax(IR0); minIR = np.nanmin(IR0)
    if maxIG <= 0 or maxIR <= 0 or minIG <= 0 or minIR <= 0:
        LOGGER.error(f"{container.sample.name} one of (maxIG,maxIR,minIG,minIR) was zero; cannot run dye-bias correction")
        return container

    # make Meth + Unmeth a long sorted list of probe values, drop index
    IR1 = sorted(IR0['Meth'].tolist() + IR0['Unmeth'].tolist())
    IG1 = sorted(IG0['Unmeth'].tolist() + IG0['Meth'].tolist())

    def same_N_interpol(ig,target):
        """ stretches  IG to IR's number of points, and visa versa, using linear interpolation., before feeding into qnorm
        - in: two lists (IG1, IR1)
        - sesame's inputs were IR1, target=IG0."""
        ig_length = len(ig); target_length = len(target)
        # make some linear functions for ir and ig
        #ig_interp = scipy.interpolate.interp1d(np.arange(np.array(ig).size), np.array(ig), fill_value="extrapolate")
        target_interp = scipy.interpolate.interp1d(np.arange(np.array(target).size), np.array(target), fill_value="extrapolate")
        # re-cast as new length of other data set (num= <new size>)
        #ig_stretch = ig_interp(np.linspace(0, ig_length, num=target_length))
        ig_stretch = target_interp(np.linspace(0, target_length, num=ig_length))
        return np.array(ig_stretch) #, np.array(ir_stretch)

    IG_stretch = np.sort(same_N_interpol(IR1, IG1))
    IR_stretch = np.sort(same_N_interpol(IG1, IR1))

    if len(IG1) != len(IR_stretch):
        raise ValueError("wrong length")
    #if debug:
    #    print(f"{len(IG1)} {len(IR1)} --> {IG_stretch.shape[0]} {IR_stretch.shape[0]}")
    #    print(f"IG1 {len(IG1)} | IR1 {len(IR1)} --stretched--> IG {len(IG_stretch)} | IR {len(IR_stretch)}")
    IR2 = qnorm_using_target(np.transpose(np.array([IR1])), np.transpose(np.array([IG_stretch])))
    IG2 = qnorm_using_target(np.transpose(np.array([IG1])), np.transpose(np.array([IR_stretch])))

    #IRmid = (IR1 + IR2) / 2.0 # avg of IR_meth and qnorm -- interpolated IG_unmeth values
    IRmid = (np.array(IR1) + IR2) / 2.0
    maxIRmid = max(IRmid)
    minIRmid = min(IRmid)

    IGmid = (np.array(IG1) + IG2) / 2.0
    maxIGmid = max(IGmid)
    minIGmid = min(IGmid)

    def fit_func_red(data):
        """ handles II probes that are out of IG, IR range by transforming all proportionally. """
        #data = data.copy()
        _IR1 = np.array(IR1)
        _IRmid = np.array(IRmid)
        #insupp = (minIR <= data <= maxIR) and not data.isna() # filter of data within allowed range
        #insupp = data >= minId & data <= maxIR & ~data.isna()
        insupp = (( data >= np.nanmin( minIR ) ) & ( data <= np.nanmax( maxIR ) ) & ~data.isna())
        #oversupp = data > maxIR and not data.isna() # a filter of data over max
        oversupp = (data > maxIR) & ~data.isna()
        #undersupp = data < minIR and not data.isna() # a filter of data under min
        undersupp = (data < minIR) & ~data.isna()
        # approx() Returns a list of points which linearly interpolate given data points, or a function performing the linear (or constant) interpolation.
        # transform the in-support probes proportionally (probes where function is defined)
        # data[insupp] <- approx(x=IG1, y=IGmid, xout=data[insupp], ties=mean)$y
        mask = ~np.isnan(_IR1) & ~np.isnan(_IRmid)
        yinterp = np.interp(x=data.loc[insupp], xp=_IR1[mask], fp=_IRmid[mask], period=None, left=None, right=None)
        data.loc[insupp] = yinterp
        data.loc[oversupp] = data.loc[oversupp] - maxIR + maxIRmid
        data.loc[undersupp] = data.loc[undersupp] * (minIRmid / minIR)
        """
        slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(x=_IR1[mask], y=_IRmid[mask])
        #print(slope, intercept, r_value, p_value, std_err)
        # transform insupp using the function, where (x) is insupp prev values
        data.loc[insupp] = slope*data.loc[insupp] + intercept
        #data[oversupp] = data[oversupp] - maxIR + maxIRmid
        data.loc[oversupp] = data.loc[oversupp] - maxIR + maxIRmid
        # transform low-intensity values, avoiding negative or zero values
        data.loc[undersupp] = minIRmid / minIR * data.loc[undersupp]
        #data = data.round(1) # or Int64 with NA
        """
        return data

    def fit_func_green(data):
        #data = data.copy()
        _IG1 = np.array(IG1)
        _IGmid = np.array(IGmid)
        insupp = ( data >= np.nanmin( minIG ) ) & ( data <= np.nanmax( maxIG ) ) & ~data.isna()
        #oversupp = data > maxIR and not data.isna() # a filter of data over max
        oversupp = (data > maxIG) & ~data.isna()
        #undersupp = data < minIR and not data.isna() # a filter of data under min
        undersupp = (data < minIG) & ~data.isna()
        # approx() Returns a list of points which linearly interpolate given data points, or a function performing the linear (or constant) interpolation.
        # transform the in-support probes proportionally (probes where function is defined)
        # data[insupp] <- approx(x=IG1, y=IGmid, xout=data[insupp], ties=mean)$y
        mask = ~np.isnan(_IG1) & ~np.isnan(_IGmid)
        yinterp = np.interp(x=data.loc[insupp], xp=_IG1[mask], fp=_IGmid[mask], period=None, left=None, right=None)
        data.loc[insupp] = yinterp
        data.loc[oversupp] = data.loc[oversupp] - maxIG + maxIGmid
        data.loc[undersupp] = data.loc[undersupp] * (minIGmid / minIG)

        """ previous working method
        slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(x=_IG1[mask], y=_IGmid[mask])
        #print(slope, intercept, r_value, p_value, std_err)
        # transform insupp using the function, where (x) is insupp prev values
        data.loc[insupp] = slope*data.loc[insupp] + intercept
        data.loc[oversupp] = data.loc[oversupp] - maxIG + maxIGmid
        # transform low-intensity values, avoiding negative or zero values
        data.loc[undersupp] = minIGmid / minIG * data.loc[undersupp]
        #data = data.round(1)
        """
        return data

    meth = 'noob_Meth' if container.do_noob else 'Meth'
    unmeth = 'noob_Unmeth' if container.do_noob else 'Unmeth'

    transformed_II_meth = fit_func_green(container.II[meth].astype('float32').copy()).round()
    transformed_II_unmeth = fit_func_red(container.II[unmeth].astype('float32').copy()).round()
    transformed_IR_meth = fit_func_red(container.IR[meth].astype('float32').copy()).round()
    transformed_IR_unmeth = fit_func_red(container.IR[unmeth].astype('float32').copy()).round()
    transformed_IG_meth = fit_func_green(container.IG[meth].astype('float32').copy()).round()
    transformed_IG_unmeth = fit_func_green(container.IG[unmeth].astype('float32').copy()).round()
    #oobR = fit_func_red(container.oobR[meth].astype('float32').copy()) # 2021-03-22 assumed 'mean_value' for red and green MEANT meth and unmeth (OOBS), respectively.
    #oobG = fit_func_green(container.oobG[unmeth].astype('float32').copy()) # v1.5.0+ uses noob version now, if available.

    if len(container.ctrl_red) == 0 or len(container.ctrl_green) == 0:
        pass # not correcting these if missing; sesame had this caveat too
    else:
        # THIS IS NOT SAVED BELOW... yet.
        ctrl_red = fit_func_red(container.ctrl_red['mean_value'].astype('float32').copy()).round()
        ctrl_green = fit_func_red(container.ctrl_green['mean_value'].astype('float32').copy()).round()

    if debug:
        pass
        """
        fig,ax = plt.subplots(3, 2, figsize=(12,8))
        scipy.stats.probplot(container.II[meth], plot=ax[0,1])
        scipy.stats.probplot(transformed_II_meth, plot=ax[0,1])
        ax[0,1].set_title("transformed_II_meth")
        ax[0,1].get_lines()[0].set_markerfacecolor('r')
        ax[0,1].get_lines()[0].set_color('r')
        scipy.stats.probplot(container.II[unmeth], plot=ax[0,0])
        scipy.stats.probplot(transformed_II_unmeth, plot=ax[0,0])
        ax[0,0].get_lines()[0].set_markerfacecolor('g')
        ax[0,0].get_lines()[0].set_color('g')
        ax[0,0].set_title("transformed_II_unmeth")
        scipy.stats.probplot(container.IR[meth], plot=ax[1,0])
        scipy.stats.probplot(transformed_IR_meth, plot=ax[1,0])
        ax[1,0].set_title("transformed_IR_meth")
        ax[1,0].get_lines()[0].set_markerfacecolor('r')
        ax[1,0].get_lines()[0].set_color('r')
        scipy.stats.probplot(container.IG[unmeth], plot=ax[1,1])
        scipy.stats.probplot(transformed_IG_unmeth, plot=ax[1,1])
        ax[1,1].set_title("transformed_IG_unmeth")
        ax[1,1].get_lines()[0].set_markerfacecolor('g')
        ax[1,1].get_lines()[0].set_color('g')
        scipy.stats.probplot(container.IR[unmeth], plot=ax[2,0])
        scipy.stats.probplot(transformed_IR_unmeth, plot=ax[2,0])
        ax[2,0].set_title("transformed_IR_unmeth")
        ax[2,0].get_lines()[0].set_markerfacecolor('r')
        ax[2,0].get_lines()[0].set_color('r')
        scipy.stats.probplot(container.IG[meth], plot=ax[2,1])
        scipy.stats.probplot(transformed_IG_meth, plot=ax[2,1])
        ax[2,1].set_title("transformed_IG_meth")
        ax[2,1].get_lines()[0].set_markerfacecolor('g')
        ax[2,1].get_lines()[0].set_color('g')

        plt.show()
        transformed_II_meth.plot.hist(bins=100, alpha=0.5, legend=True)
        transformed_II_unmeth.plot.hist(bins=100, alpha=0.5, legend=True)
        transformed_IR_meth.plot.hist(bins=100, alpha=0.5, legend=True)
        transformed_IG_unmeth.plot.hist(bins=100, alpha=0.5, legend=True)
        transformed_IR_unmeth.plot.hist(bins=100, alpha=0.5, legend=True)
        transformed_IG_meth.plot.hist(bins=100, alpha=0.5, legend=True)
        plt.show()
        (container.II[meth] - transformed_II_meth).plot.hist(bins=100, alpha=0.5, legend=True)
        (container.II[unmeth] - transformed_II_unmeth).plot.hist(bins=100, alpha=0.5, legend=True)
        (container.IR[meth] - transformed_IR_meth).plot.hist(bins=100, alpha=0.5, legend=True)
        (container.IG[unmeth] - transformed_IG_unmeth).plot.hist(bins=100, alpha=0.5, legend=True)
        (container.IG[meth] - transformed_IG_meth).plot.hist(bins=100, alpha=0.5, legend=True)
        (container.IR[unmeth] - transformed_IR_unmeth).plot.hist(bins=100, alpha=0.5, legend=True)
        (container.IR[meth] - transformed_IR_meth).plot.hist(bins=100, alpha=0.5, legend=True)
        (container.IG[unmeth] - transformed_IG_unmeth).plot.hist(bins=100, alpha=0.5, legend=True)
        plt.show()
        """

    # IG, IR, II, oobG, oobR must be updated. -- if input was raw_IG/IR/II the output column is still called noob_meth in DF internally.
    # [mean_value | bg_corrected | noob] -- only noob updated
    # updates work if indexes match and column names match
    noob = 'noob' if container.do_noob else 'Meth'
    transformed_II_meth.name = noob
    transformed_IG_meth.name = noob
    transformed_IR_meth.name = noob
    #oobR.name = noob
    container.methylated.update(transformed_II_meth)
    container.methylated.update(transformed_IG_meth)
    container.methylated.update(transformed_IR_meth)
    container.II.update(transformed_II_meth)
    container.IG.update(transformed_IG_meth)
    container.IR.update(transformed_IR_meth)
    #container.oobR.update(oobR)
    container._SampleDataContainer__data_frame['noob_meth'] = container.methylated[noob].round()

    noob = 'noob' if container.do_noob else 'Unmeth'
    transformed_II_unmeth.name = noob
    transformed_IG_unmeth.name = noob
    transformed_IR_unmeth.name = noob
    #oobG.name = noob
    container.unmethylated.update(transformed_II_unmeth)
    container.unmethylated.update(transformed_IG_unmeth)
    container.unmethylated.update(transformed_IR_unmeth)
    container.II.update(transformed_II_unmeth)
    container.IG.update(transformed_IR_unmeth)
    container.IR.update(transformed_IR_unmeth)
    #container.oobG.update(oobG)
    container._SampleDataContainer__data_frame['noob_unmeth'] = container.unmethylated[noob].round()

    container.check_for_probe_loss(f"dye_bias - {noob}") # looking for probes that got dropped by accident.

    # CONTROLS are pulled directly from manifest; not updated
    container.ctrl_green = container.ctrl_green.assign(noob=ctrl_green)
    container.ctrl_red = container.ctrl_red.assign(noob=ctrl_red)

    container._SigSet__dye_bias_corrected = True

    if debug:
        return {'IIu':transformed_II_unmeth,
               'IIm': transformed_II_meth,
               'IR': transformed_IR_meth,
               'IRu': transformed_IR_unmeth,
               'IG': transformed_IG_unmeth,
               'IGm': transformed_IG_meth,
               'ctrlR': ctrl_red,
               'ctrlG': ctrl_green,
               #'oobG': oobG,
               #'oobR': oobR,
               'IG1': IG1,
               'IR1': IR1,
               'IGmid': IGmid,
               'IRmid': IRmid,
               'IG2': IG2,
               'IR2': IR2,
               'IG_stretch': IG_stretch,
               'IR_stretch': IR_stretch}
    return
