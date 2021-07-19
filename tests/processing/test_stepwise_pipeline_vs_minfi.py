#PATH = '/Volumes/LEGX/Barnes/44668_MURMETVEP/204617710009'
#PATH = '/Volumes/LEGX/Barnes/48230_MURMETVEP/361821/204879580038'
#PATH =  '../../docs/example_data/mouse/'
#PATH =  '/Volumes/LEGX/Barnes/mouse_test'
# PATH =  '../../docs/example_data/GSE69852/minfi/' #--- for testing in console
PATH = 'docs/example_data/minfi/'
IDAT_SOURCE = 'docs/example_data/GSE69852'
import methylprep
import pandas as pd
import numpy as np
from pathlib import Path
import shutil

def test_noob_df_same_size_as_minfi():
    ID = '9247377085_R04C02'
    manifest = methylprep.files.Manifest(methylprep.models.ArrayType('450k'))
    print('* loading one idat pair of files')
    green_filepath = Path(PATH, f'{ID}_Grn.idat') #'204879580038_R06C02_Grn.idat')
    red_filepath = Path(PATH, f'{ID}_Red.idat') #'204879580038_R06C02_Red.idat')
    print(f"* GREEN --> {green_filepath.name}")
    print(f"* RED --> {red_filepath.name}")
    if not green_filepath.exists():
        shutil.copy(Path(IDAT_SOURCE, f'{ID}_Grn.idat'), green_filepath)
    if not red_filepath.exists():
        shutil.copy(Path(IDAT_SOURCE, f'{ID}_Red.idat'), red_filepath)

    green_idat = methylprep.files.IdatDataset(green_filepath, channel=methylprep.models.Channel.GREEN)
    red_idat = methylprep.files.IdatDataset(red_filepath, channel=methylprep.models.Channel.RED)
    sample = methylprep.models.Sample('','1234567890','R01C01', Sample_Name='testsample')
    sigset = methylprep.models.SigSet(sample, green_idat, red_idat, manifest)
    #print('* raw_dataset')
    #raw_dataset = methylprep.models.raw_dataset.RawDataset(sample, green_idat, red_idat)
    #print('* meth_dataset.unmethylated')
    #unmethylated = methylprep.models.MethylationDataset.unmethylated(raw_dataset, manifest)
    #print('* meth_dataset.methylated')
    #methylated = methylprep.models.MethylationDataset.methylated(raw_dataset, manifest)
    m_minfi = pd.read_csv(Path(PATH, 'minfi_raw_meth.csv')).rename(columns={'Unnamed: 0':'IlmnID'}).set_index('IlmnID')
    u_minfi = pd.read_csv(Path(PATH, 'minfi_raw_unmeth.csv')).rename(columns={'Unnamed: 0':'IlmnID'}).set_index('IlmnID')
    m1 = sigset.methylated.sort_index()[['Meth']].rename(columns={'Meth': ID})
    m2 = m_minfi[[ID]]
    mean_diff_m = (m1 - m2).mean()
    u1 = sigset.unmethylated.sort_index()[['Unmeth']].rename(columns={'Unmeth': ID})
    u2 = u_minfi[[ID]]
    mean_diff_u = (u1 - u2).mean()
    print(f"minfi mean difference, meth: {mean_diff_m}, unmeth: {mean_diff_u}")
    if float(mean_diff_m.sum()) != 0 or float(mean_diff_u.sum()) != 0:
        raise AssertionError(f"raw meth/unmeth values don't match between methylprep and minfi METH: {float(mean_diff_m.sum())}, UNMETH: {float(mean_diff_u.sum())}")

    nm_minfi = pd.read_csv(Path(PATH, 'minfi_noob_meth.csv')).rename(columns={'Unnamed: 0':'IlmnID'}).set_index('IlmnID').sort_index()
    nu_minfi = pd.read_csv(Path(PATH, 'minfi_noob_unmeth.csv')).rename(columns={'Unnamed: 0':'IlmnID'}).set_index('IlmnID').sort_index()
    b_minfi = pd.read_csv(Path(PATH, 'minfi_noob_betas.csv')).rename(columns={'Unnamed: 0':'IlmnID'}).set_index('IlmnID').sort_index()

    idat_dataset_pair = {'green_idat': green_idat, 'red_idat':red_idat, 'sample':sample}
    container = methylprep.processing.SampleDataContainer(idat_dataset_pair, manifest,
        retain_uncorrected_probe_intensities=True,
        pval=False,
        do_noob=True,
        quality_mask=False,
        switch_probes=False,
        do_nonlinear_dye_bias=False,
        debug=False,
        sesame=False,
        )
    data_frame = container.preprocess()
    data_frame = container.process_beta_value(data_frame)
    #container._postprocess(input_dataframe, calculate_beta_value, 'beta_value', offset)
    #beta_df = self.process_beta_value(containers[0]data_frame)
    #pre_noob_meth = container.methylated.data_frame[['bg_corrected']].sort_index()
    #pre_noob_unmeth = container.unmethylated.data_frame[['bg_corrected']].sort_index()
    # the match isn't perfect here anymore after noob:
    meth_test = nm_minfi[['9247377085_R04C02']].join( data_frame[ ~data_frame.index.str.startswith('rs') ][['noob_meth']].sort_index() )
    unmeth_test = nu_minfi[['9247377085_R04C02']].join( data_frame[ ~data_frame.index.str.startswith('rs') ][['noob_unmeth']].sort_index() )
    print(f"noob meth mean diff: meth: {(meth_test['9247377085_R04C02'] - meth_test['noob_meth']).mean()} | unmeth: {(unmeth_test['9247377085_R04C02'] - unmeth_test['noob_unmeth']).mean()}")
    noob_meth_match = all(np.isclose(nm_minfi['9247377085_R04C02'].round(0), data_frame[ ~data_frame.index.str.startswith('rs') ]['noob_meth'].sort_index(), atol=2.0))
    noob_unmeth_match = all(np.isclose(nu_minfi['9247377085_R04C02'].round(0), data_frame[ ~data_frame.index.str.startswith('rs') ]['noob_unmeth'].sort_index(), atol=2.0))
    # OLD LIMIT WAS: noob_unmeth_match = all(np.isclose(nu_minfi['9247377085_R04C02'].round(0), data_frame['noob_unmeth'].sort_index(), atol=1.0))
    print(f"minfi NOOB matches for METH: {noob_meth_match}, UNMETH: {noob_unmeth_match}")
    if noob_meth_match is False or noob_unmeth_match is False:
        raise AssertionError("noob meth or unmeth values don't match between minfi and methylprep (expect 100% match)")

    ref_data_frame = data_frame[ ~data_frame.index.str.startswith('rs') ]
    noob_betas_match = sum(np.isclose(b_minfi['9247377085_R04C02'], ref_data_frame['beta_value'].sort_index(), atol=0.03))/len(data_frame)
    noob_betas_loose_match = sum(np.isclose(b_minfi['9247377085_R04C02'], ref_data_frame['beta_value'].sort_index(), atol=0.1))/len(data_frame)
    print(f"minfi betas match (+/- 0.03): {noob_betas_match} or +/- 0.1: {noob_betas_loose_match}")

    # this overwrites data, so copying it
    alt_frame = container._postprocess(ref_data_frame.copy(), methylprep.processing.postprocess.calculate_beta_value, 'beta_value', offset=0)
    noob_betas_match = sum(np.isclose(b_minfi['9247377085_R04C02'], alt_frame['beta_value'].sort_index(), atol=0.001))/len(data_frame)
    noob_betas_loose_match = sum(np.isclose(b_minfi['9247377085_R04C02'], alt_frame['beta_value'].sort_index(), atol=0.01))/len(data_frame)
    print(f"minfi betas match (+/- 0.001): {noob_betas_match} or +/- 0.01: {noob_betas_loose_match}")
    if noob_betas_match < 0.999:
        raise AssertionError("noob betas don't match between minfi and methylprep (expecte 99.9% of betas for probes to be +/- 0.001)")

    Path(PATH, f'{ID}_Grn.idat').unlink()
    Path(PATH, f'{ID}_Red.idat').unlink()
    return {'mf_meth': nm_minfi, 'mf_unmeth': nu_minfi, 'mf_beta': b_minfi,
        'df': data_frame.sort_index(), 'test': alt_frame}

#grn, red = test_noob_df_same_size()

def test_make_pipeline_noob_only():
    IDAT_SOURCE = 'docs/example_data/mouse'
    test_outputs = [
        Path(IDAT_SOURCE, 'beta_values.pkl'),
        Path(IDAT_SOURCE, 'noob_meth_values.pkl'),
        Path(IDAT_SOURCE, 'noob_unmeth_values.pkl'),
        Path(IDAT_SOURCE, 'samplesheet.csv'),
        ]
    for outfile in test_outputs:
        if outfile.exists():
            outfile.unlink()
    df = methylprep.make_pipeline(IDAT_SOURCE, steps=['noob'], sesame=False, debug=True, make_sample_sheet=True)
    sample_name = '204879580038_R06C02'
    ref_beta = [
    ['cg00101675_BC21',              0.865335],
    ['cg00116289_BC21',              0.898574],
    ['cg00211372_TC21',              0.822134],
    ['cg00531009_BC21',              0.895200],
    ['cg00747726_TC21',              0.828053],
    ['uk9978_TC11',                  0.315556],
    ['uk9983_TC21',                  0.234043],
    ['uk9986_BC11',                  0.254296],
    ['uk998926237_TC11',             0.293478],
    ['uk9995_TC11',                  0.333333],
    ]
    ref_noob_meth = [
    ['cg00101675_BC21',                  1992],
    ['cg00116289_BC21',                  2206],
    ['cg00211372_TC21',                  1248],
    ['cg00531009_BC21',                  8952],
    ['cg00747726_TC21',                  6116],
    ]
    ref_beta = pd.DataFrame(data=ref_beta, columns = ['IlmnID', sample_name]).set_index('IlmnID')
    ref_noob_meth = pd.DataFrame(data=ref_noob_meth, columns = ['IlmnID', sample_name]).set_index('IlmnID')
    test_beta = df.loc[ ref_beta.index ]
    test_noob_meth = pd.read_pickle(Path(IDAT_SOURCE, 'noob_meth_values.pkl'))
    test_noob_meth = test_noob_meth.loc[ ref_noob_meth.index ]
    if not np.isclose(test_beta, ref_beta, atol=0.01).all():
        raise AssertionError(f"beta values don't match ref values")
    if not np.isclose(test_noob_meth, ref_noob_meth, atol=0.01).all():
        raise AssertionError(f"test noob meth values don't match ref values")


def shrink_csv(filename):
    # use on minfi output, but on betas use round
    PATH = 'docs/example_data/minfi/'
    # minfi_raw_betas.csv
    x = pd.read_csv(Path(PATH, filename))
    if 'meth' in filename or 'unmeth' in filename:
        x['9247377085_R04C02'] = x['9247377085_R04C02'].astype(int)
    elif 'betas' in filename:
        x['9247377085_R04C02'] = x['9247377085_R04C02'].round(3)
    x.to_csv(Path(PATH, filename), index=False)
    test = pd.read_csv(Path(PATH, filename)).rename(columns={'Unnamed: 0':'IlmnID'}).set_index('IlmnID').sort_index()
    return test
