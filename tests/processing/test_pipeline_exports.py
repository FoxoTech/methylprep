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
        probes = ['cg02713524_II_R_C_rep1_EPIC', 'cg18835320_II_R_C_rep1_EPIC', 'cg05205094_II_F_C_rep1_EPIC', 'cg23230554_I_R_C_rep1_EPIC', 'cg05369274_I_F_C_rep1_EPIC', 'cg06844455_II_F_C_rep1_EPIC', 'cg07506599_I_R_C_rep1_EPIC', 'cg01482676_II_R_C_rep1_EPIC', 'cg15418368_II_R_C_rep1_EPIC', 'cg25991066_II_R_C_rep1_EPIC']
        test1 = test1.set_index('IlmnID').loc[probes]

        if num_missing > 0:
            raise ValueError(f"{num_missing} EPIC+ beta values are missing from beta_value of processed CSV file.")
            # in older versions, the EPIC+ manifest had one probe mistake.
            #if test1[test1.beta_value.isna()]['IlmnID'].iloc[0] == 'cg00968771_I_F_C_rep1_GWG1':
            #    print("WARNING: cg00968771_I_F_C_rep1_GWG1 probe data is STILL missing from output")
            #    #NOT A FATAL ERROR. but not fixing today.
        test1_ref = [
            ['cg00000029_II_F_C_rep1_EPIC',     1180.0,        274.0,       0.759,    2.107],
            ['cg00000103_II_F_C_rep1_EPIC',     1363.0,        681.0,       0.636,    1.001],
            ['cg00000109_II_F_C_rep1_EPIC',     2427.0,        284.0,       0.863,    3.095],
            ['cg00000155_II_F_C_rep1_EPIC',     4392.0,        212.0,       0.934,    4.373],
            ['cg00000158_II_F_C_rep1_EPIC',     5742.0,        211.0,       0.949,    4.766],
            ['rs6626309_I_N_C_rep1_EPIC',        189.0,       2376.0,       0.071,   -3.652],
            ['rs6991394_I_F_C_rep1_GWG1',        261.0,       3956.0,       0.060,   -3.922],
            ['rs6991394_I_N_C_rep1_EPIC',        262.0,       5099.0,       0.048,   -4.283],
            ['rs9292570_I_F_C_rep1_GWG1',       7714.0,        208.0,       0.962,    5.213],
            ['rs9292570_I_N_C_rep1_EPIC',       9152.0,        254.0,       0.963,    5.171]
        ]
        test1_ref = pd.DataFrame(test1_ref, columns=['IlmnID', 'noob_meth','noob_unmeth','beta_value','m_value']).set_index('IlmnID')
        test1 = pd.read_csv(testfile_1).set_index('IlmnID').loc[test1_ref.index]
        if not test1_ref.equals(test1):
            raise AssertionError('CSV output doesnt match expected values for a selection of probes: {list(test1_ref.index)}')
        probes = ['cg00000158_II_F_C_rep1_EPIC','cg00000029_II_F_C_rep1_EPIC']
        if not np.isclose(test_data_containers[0]._SampleDataContainer__data_frame.loc[probes]['noob_unmeth'], [211, 274], atol=1.0).all():
            #np.isclose(test_data_containers[0]._SampleDataContainer__data_frame.iloc[2]['noob_unmeth'], 284.0, atol=1.0):
            print(test_data_containers[0]._SampleDataContainer__data_frame)
            raise AssertionError(f"data_container output ({test_data_containers[0]._SampleDataContainer__data_frame.iloc[2]['noob_unmeth']}) differs from expected value (284.0)")
        if not np.isclose(test_data_containers[0]._SampleDataContainer__data_frame.loc[probes]['noob_meth'], [5742, 1180], atol=1.0).all():
            #np.isclose(test_data_containers[0]._SampleDataContainer__data_frame.iloc[2]['noob_unmeth'], 284.0, atol=1.0):
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
