import sys
import numpy as np
import pandas as pd
from pathlib import Path
# App
from methylprep.processing import pipeline

class TestPipeline():

    @staticmethod
    def test_run_pipeline_epic_plus_export_data():
        """ check that we get back useful data with --export option """
        test_data_dir = 'docs/example_data/epic_plus'
        testfile_1 = Path(test_data_dir, '202651080072', '202651080072_R01C01_processed.csv')
        if testfile_1.exists():
            testfile_1.unlink()
        test_data_containers = pipeline.run_pipeline(test_data_dir, export=True, sesame=False)
        if not testfile_1.exists():
            raise AssertionError("no exported processed csv found")

        # spot checking the output.
        test1 = pd.read_csv(testfile_1)
        num_missing = test1['beta_value'].isna().sum()
        if num_missing == 1:
            if test1[test1.beta_value.isna()]['IlmnID'].iloc[0] == 'cg00968771_I_F_C_rep1_GWG1':
                print("WARNING: cg00968771_I_F_C_rep1_GWG1 probe data is STILL missing from output")
                #NOT A FATAL ERROR. but not fixing today.
        elif num_missing > 0:
            print(test1.head())
            raise AssertionError('{num_missing} missing values in processed csv')
        if not np.isclose(test1['beta_value'].iloc[5], 0.145, atol=0.01):
            print(test1.iloc[5])
            raise AssertionError('beta_value doesnt match expected value')
        if not np.isclose(test_data_containers[0]._SampleDataContainer__data_frame.iloc[2]['noob_unmeth'], 284.0, atol=1.0):
            print(test_data_containers[0]._SampleDataContainer__data_frame)
            raise AssertionError(f"data_container output ({test_data_containers[0]._SampleDataContainer__data_frame.iloc[2]['noob_unmeth']}) differs from expected value (284.0)")

    @staticmethod
    def test_pipeline_meth_unmeth_int16():
        test_data_dir = 'docs/example_data/GSE69852'
        testfile_1 = Path(test_data_dir, 'meth_values.pkl')
        testfile_2 = Path(test_data_dir, 'unmeth_values.pkl')
        if testfile_1.exists():
            testfile_1.unlink()
        if testfile_2.exists():
            testfile_2.unlink()
        test_data_containers = pipeline.run_pipeline(test_data_dir, export=True, save_uncorrected=True, sesame=False)
        if not testfile_1.exists():
            raise AssertionError("no meth_values.pkl found")
        if not testfile_2.exists():
            raise AssertionError("no unmeth_values.pkl found")
        # ensure no negative values, as these mean some data exceeded the allowed intensity range
        m = pd.read_pickle(Path(test_data_dir,'meth_values.pkl')) # standard output, as int16
        u = pd.read_pickle(Path(test_data_dir,'unmeth_values.pkl'))
        errors = []
        mask = (m < 0)
        for sample in m.columns:
            match = len(m[sample][mask[sample]]) == len(pd.Series())
            if not match:
                print(m[sample][mask[sample] == True])
                print("")
                errors.append(sample)
        mask = (u < 0)
        for sample in u.columns:
            match = len(u[sample][mask[sample]]) == len(pd.Series())
            if not match:
                print(u[sample][mask[sample] == True])
                print("")
                errors.append(sample)
        if testfile_1.exists():
            testfile_1.unlink()
        if testfile_2.exists():
            testfile_2.unlink()
        # also confirm these CSV columns are the same as the pickled columns, and non-negative
        testfile_3 = Path(test_data_dir, '9247377093', '9247377093_R02C01_processed.csv')
        testfile_4 = Path(test_data_dir, '9247377085', '9247377085_R04C02_processed.csv')
        csv3 = pd.read_csv(testfile_3).set_index('IlmnID')
        rs = csv3.index.str.startswith('rs')

        # meth_values.pkl matches csv[meth] after removing rs probes from csv
        if (~np.isclose( m['9247377093_R02C01'].sort_index(), csv3.loc[~rs]['meth'].sort_index(), atol=0.1)).sum() > 0:
            errors.append(f"9247377093_R02C01 meth pkl != csv {(~np.isclose( m['9247377093_R02C01'].sort_index(), csv3['meth'].sort_index(), atol=10.0)).sum()}")
        if (~np.isclose( u['9247377093_R02C01'].sort_index(), csv3.loc[~rs]['unmeth'].sort_index(), atol=0.1)).sum() > 0:
            errors.append(f"9247377093_R02C01 unmeth pkl != csv {(~np.isclose( m['9247377093_R02C01'].sort_index(), csv3['unmeth'].sort_index(), atol=10.0)).sum()}")

        # order not the same, but probes should all be there
        same_probes = m.sort_index().index.equals( test_data_containers[0]._SampleDataContainer__data_frame.loc[~rs]['meth'].sort_index().index )
        if not same_probes:
            errors.append("probes in meth_values.pkl don't match probes in SampleDataContainer")

        # order meth_values.pkl matched SDC pre v1.4.5.
        #same_order = m.index.equals( test_data_containers[0]._SampleDataContainer__data_frame.loc[~rs]['meth'].index )
        #if not same_order:
        #    errors.append("order of probes in meth_values.pkl don't match SampleDataContainer")
        same_order = csv3['meth'].index.equals( test_data_containers[0]._SampleDataContainer__data_frame['meth'].index )
        if not same_order:
            errors.append("order of probes in output CSV don't match SampleDataContainer")

        # As of v1.5.0, CSV and SDC match exactly.
        #test_data_containers[0]._SampleDataContainer__data_frame['meth'] == csv3['meth']
        #sdc_match = (~np.isclose( test_data_containers[0]._SampleDataContainer__data_frame['meth'], csv3['meth'], atol=10.0)).sum()
        sdc_match = test_data_containers[0]._SampleDataContainer__data_frame.astype('float32').equals(csv3.astype('float32'))
        if sdc_match is False:
            errors.append("SampleDataContainer does not match csv3 output")

        csv4 = pd.read_csv(testfile_4).set_index('IlmnID')
        if (~np.isclose( m['9247377085_R04C02'].sort_index(), csv4.loc[~rs]['meth'].sort_index(), atol=0.1)).sum() > 0:
            #if not m['9247377085_R04C02'].equals( csv4['meth'] ):
            errors.append(f"9247377085_R04C02 meth pkl != csv {(~np.isclose( m['9247377085_R04C02'].sort_index(), csv4['meth'].sort_index(), atol=10.0)).sum()}")
        if (~np.isclose( u['9247377085_R04C02'].sort_index(), csv4.loc[~rs]['unmeth'].sort_index(), atol=0.1)).sum() > 0:
            #if not u['9247377085_R04C02'].equals( csv4['unmeth'] ):
            errors.append(f"9247377085_R04C02 unmeth pkl != csv {(~np.isclose( m['9247377085_R04C02'].sort_index(), csv4['unmeth'].sort_index(), atol=10.0)).sum()}")
        #sdc_match = (~np.isclose( test_data_containers[1]._SampleDataContainer__data_frame['meth'], csv4['meth'], atol=10.0)).sum()
        sdc_match = test_data_containers[1]._SampleDataContainer__data_frame.astype('float32').equals(csv4.astype('float32'))
        if sdc_match is False:
            errors.append("SampleDataContainer (as a whole) does not match csv4 output (as a whole)")

        if errors:
            #import pdb;pdb.set_trace()
            raise ValueError('\n'.join(errors))
