import os
import sys
import numpy as np
import pandas as pd
from pathlib import Path
# App
from methylprep.processing import pipeline
#patching
import unittest
try:
    # python 3.4+ should use builtin unittest.mock not mock package
    from unittest.mock import patch
except ImportError:
    from mock import patch


class TestPipeline():
    columns = { # the atol tolerance for each csv column returned in SampleDataContainer DF
        'meth':1,
        'noob_meth':1,
        'unmeth':5,
        'noob_unmeth':5,
        'poobah_pval':0.01,
        'quality_mask':0.01,
        'beta_value':0.01,
        'm_value':0.01}

    def test_pipeline_cli_minfi(self):
        test_data_dir = 'docs/example_data/GSE69852'
        test_outputs = [
            Path(test_data_dir, 'control_probes.pkl'),
            Path(test_data_dir, 'beta_values.pkl'),
            Path(test_data_dir, 'm_values.pkl'),
            Path(test_data_dir, 'meth_values.pkl'),
            Path(test_data_dir, 'unmeth_values.pkl'),
            Path(test_data_dir, 'noob_meth_values.pkl'),
            Path(test_data_dir, 'noob_unmeth_values.pkl'),
            Path(test_data_dir, 'sample_sheet_meta_data.pkl'),
            Path(test_data_dir, '9247377085', '9247377085_R04C02_processed.csv'),
            Path(test_data_dir, '9247377093', '9247377093_R02C01_processed.csv'),
            ]
        for outfile in test_outputs:
            if outfile.exists():
                outfile.unlink()

        exit_status = os.system(f'python -m methylprep process -d {test_data_dir} --all --minfi')
        if exit_status != 0:
            for outfile in test_outputs:
                if outfile.exists():
                    outfile.unlink()
            raise AssertionError("methylprep process CLI failed with error(s)")

        testfile_1 = Path(test_data_dir, '9247377085', '9247377085_R04C02_processed.csv')
        test1 = pd.read_csv(testfile_1).set_index('IlmnID')
        if test1['beta_value'].isna().sum() > 0:
            print(test1.head())
            raise AssertionError('missing beta values in processed csv')

        testfile_2 = Path(test_data_dir, '9247377093', '9247377093_R02C01_processed.csv')
        test2 = pd.read_csv(testfile_2)
        if test2['beta_value'].isna().sum() > 0:
            print(test2.head())
            raise AssertionError('missing beta values in processed csv')


        test1_ref_v151 = [ # bug with --all was applying nonlinear dye correction with --minfi here. fixed in v1.5.2.
            ['cg00000029',   3325.0,     2242.0,  1062.0,      823.0, 0.001, 0.0, 0.731,    1.446],
            ['cg00000108',  10266.0,     7892.0,  1035.0,      788.0,     0, 0.0, 0.909,    3.324],
            ['cg00000109',   5055.0,     3527.0,   839.0,      534.0,     0, 0.0, 0.869,    2.724],
            ['cg00000165',    879.0,      350.0,  2556.0,     3021.0, 0.006, 0.0, 0.104,   -3.110],
            ['cg00000236',   5161.0,     3612.0,   604.0,      279.0,     0, 0.0, 0.928,    3.694],
            ['rs9363764',    4267.0,     2942.0,  1978.0,     2138.0, 0.001, 0.0, 0.579,    0.461],
            ['rs939290',     5724.0,     4040.0,  2908.0,     3555.0,     0, 0.0, 0.532,    0.185],
        ]
        # note that poobah values get a little worse when you apply minfi, linear dye bias correction instead of sesame.
        test1_ref_v162 = [
            ['cg00000029',   3325.0,     2798.0,  1062.0,     1340.0,  0.002, 1.0, 0.660,    1.062],
            ['cg00000108',  10266.0,     9739.0,  1035.0,     1288.0,      0, 1.0, 0.875,    2.919],
            ['cg00000109',   5055.0,     4528.0,   839.0,      914.0,  0.001, 1.0, 0.817,    2.309],
            ['cg00000165',    879.0,      373.0,  2556.0,     4196.0,  0.006, 1.0, 0.080,   -3.492],
            ['cg00000236',   5161.0,     4634.0,   604.0,      505.0,  0.001, 1.0, 0.885,    3.198],
            ['rs9363764',    4267.0,     3740.0,  1978.0,     3091.0,  0.001, 1.0, 0.540,    0.275],
            ['rs939290',     5724.0,     5197.0,  2908.0,     4869.0,  0.001, 1.0, 0.511,    0.094],
        ]
        test1_ref = [ # v1.6.3 -- Mauro @ illumina made a few changes to improve NOOB/poobah
            ['cg00000029',   3325.0,     2798.0,  1062.0,   1338.0,  0.003,  1.0,   0.661,   1.064],
            ['cg00000108',  10266.0,     9739.0,  1035.0,   1288.0,  0.000,  1.0,   0.875,   2.919],
            ['cg00000109',   5055.0,     4528.0,   839.0,    914.0,  0.001,  1.0,   0.817,   2.309],
            ['cg00000165',    879.0,      374.0,  2556.0,   4194.0,  0.005,  1.0,   0.080,  -3.487],
            ['cg00000236',   5161.0,     4634.0,   604.0,    505.0,  0.001,  1.0,   0.885,   3.198],
            ['rs9363764',    4267.0,     3740.0,  1978.0,   3089.0,  0.002,  1.0,   0.540,   0.276],
            ['rs939290',     5724.0,     5197.0,  2908.0,   4867.0,  0.001,  1.0,   0.511,   0.095],
        ]
        test1_ref = pd.DataFrame(data=test1_ref, columns=['IlmnID', 'meth', 'noob_meth',  'unmeth', 'noob_unmeth',  'poobah_pval', 'quality_mask', 'beta_value', 'm_value']).set_index('IlmnID')
        test1_sub = test1.loc[ test1.index.isin(test1_ref.index) ]

        for column,tol in self.columns.items():
            sub = test1_sub[column]
            ref = test1_ref[column]
            print(f"testing {column} at {tol}")
            if not np.isclose(sub, ref, atol=tol).all():
                print("calculted results:")
                print(test1_sub[['meth', 'noob_meth',  'unmeth', 'noob_unmeth']])
                print(test1_sub[['poobah_pval', 'quality_mask', 'beta_value', 'm_value']])
                print("ref results:")
                print(test1_ref[['meth', 'noob_meth',  'unmeth', 'noob_unmeth']])
                print(test1_ref[['poobah_pval', 'quality_mask', 'beta_value', 'm_value']])
                raise AssertionError(f"test_data_containers[1]._SampleDataContainer__data_frame doesn't match ref values")
        total_nas = test1['beta_value'].isna().sum()
        if total_nas > 0:
            print(f'found {total_nas} missing beta_values (N/A or inf) in sample')
            raise AssertionError()

        for outfile in test_outputs:
            if not outfile.exists():
                raise FileNotFoundError(f"Expected {outfile.name} to be generated by run_pipeline() but it was missing.")
            else:
                print('+', outfile)
                outfile.unlink()


    @staticmethod
    def test_run_pipeline_all():
        """ check that we get back useful data.
        check that output files exist, then remove them.

        - combined with test_run_pipeline_export_data_450k()
        - checks that we get back useful data with --export option
        """
        test_data_dir = 'docs/example_data/GSE69852'
        test_outputs = [
            Path(test_data_dir, 'control_probes.pkl'),
            # Path(test_data_dir, 'beta_values.pkl'), -- cannot create this file AND return SDCs
            # Path(test_data_dir, 'm_values.pkl'),
            Path(test_data_dir, 'meth_values.pkl'),
            Path(test_data_dir, 'unmeth_values.pkl'),
            Path(test_data_dir, 'noob_meth_values.pkl'),
            Path(test_data_dir, 'noob_unmeth_values.pkl'),
            Path(test_data_dir, 'sample_sheet_meta_data.pkl'),
            Path(test_data_dir, '9247377085', '9247377085_R04C02_processed.csv'),
            Path(test_data_dir, '9247377093', '9247377093_R02C01_processed.csv'),
            ]
        for outfile in test_outputs:
            if outfile.exists():
                outfile.unlink()
        # betas=True, m_value=True returns betas instead
        test_data_containers = pipeline.run_pipeline(test_data_dir, export=True, save_uncorrected=True, save_control=True,  batch_size=None, sesame=False)

        testfile_1 = Path(test_data_dir, '9247377085', '9247377085_R04C02_processed.csv')
        test1 = pd.read_csv(testfile_1).set_index('IlmnID')
        if test1['beta_value'].isna().sum() > 0:
            print(test1.head())
            raise AssertionError('missing beta values in processed csv')
        testfile_2 = Path(test_data_dir, '9247377093', '9247377093_R02C01_processed.csv')
        test2 = pd.read_csv(testfile_2)
        if test2['beta_value'].isna().sum() > 0:
            print(test2.head())
            raise AssertionError('missing beta values in processed csv')
        test1_ref = [
            ['cg00000029',   3325.0,     2785.0,  1062.0,       1323.0,       0.662,    1.074],
            ['cg00000108',  10266.0,     9726.0,  1035.0,       1271.0,       0.876,    2.936],
            ['cg00000109',   5055.0,     4515.0,   839.0,        898.0,       0.819,    2.330],
            ['cg00000165',    879.0,      366.0,  2556.0,       4182.0,       0.079,   -3.510],
            ['cg00000236',   5161.0,     4621.0,   604.0,        498.0,       0.885,    3.214],
            ['rs9363764',    4267.0,     3727.0,  1978.0,       3076.0,       0.540,    0.277],
            ['rs939290',     5724.0,     5184.0,  2908.0,       4856.0,       0.511,    0.094],
        ]
        test1_ref = pd.DataFrame(data=test1_ref, columns=['IlmnID', 'meth', 'noob_meth',  'unmeth', 'noob_unmeth', 'beta_value', 'm_value']).set_index('IlmnID')
        test1_sub = test1.loc[ test1.index.isin(test1_ref.index) ]

        # spot checking the output.
        if not np.isclose(test1_sub, test1_ref, atol=0.01).all():
            raise AssertionError(f"test_data_containers[1]._SampleDataContainer__data_frame doesn't match ref values")
        # spot checking the output.
        total_nas = test_data_containers[0]._SampleDataContainer__data_frame['beta_value'].isna().sum()
        if total_nas > 0:
            print(f'found {total_nas} missing beta_values (N/A or inf) in sample')
            raise AssertionError()

        for outfile in test_outputs:
            if not outfile.exists():
                raise FileNotFoundError(f"Expected {outfile.name} to be generated by run_pipeline() but it was missing.")
            else:
                print('+', outfile)
                outfile.unlink()

    def test_minfi(self):
        test_data_dir = 'docs/example_data/GSE69852'
        test_outputs = [
            Path(test_data_dir, 'control_probes.pkl'),
            #Path(test_data_dir, 'beta_values.pkl'),
            #Path(test_data_dir, 'm_values.pkl'),
            #Path(test_data_dir, 'meth_values.pkl'),
            #Path(test_data_dir, 'unmeth_values.pkl'),
            Path(test_data_dir, 'noob_meth_values.pkl'),
            Path(test_data_dir, 'noob_unmeth_values.pkl'),
            Path(test_data_dir, 'sample_sheet_meta_data.pkl'),
            Path(test_data_dir, '9247377085', '9247377085_R04C02_processed.csv'),
            Path(test_data_dir, '9247377093', '9247377093_R02C01_processed.csv'),
            ]
        for outfile in test_outputs:
            if outfile.exists():
                outfile.unlink()
        #testargs = ["__program__", '-v', 'process', '-d', test_data_dir, '--no_export', '--sample_name', 'AdultLiver1', 'FetalLiver1', '--minfi', '--debug']
        #with patch.object(sys, 'argv', testargs): # these inputs feed the function, not CLI app
        #testargs = ["__program__", '-d', test_data_dir, '--no_export', '--sample_name', 'AdultLiver1', 'FetalLiver1', '--betas', '--all']
        #with patch.object(sys, 'argv', testargs):
        test_data_containers = pipeline.run_pipeline(test_data_dir, sesame=False, export=True, save_control=True)
        test2_ref = [ #minfi, without poobah and quality_mask
        # DEBUG: sesame False switch False noob True poobah False mask False, dye False == minfi mode
            ['cg00000029',     1288.0,       1969.0,    0.383676, -0.612330],
            ['cg00000108',     4631.0,       5237.0,    0.464587, -0.177417],
            ['cg00000109',     4615.0,        309.0,    0.918591,  3.900652],
            ['cg00000165',      597.0,       3518.0,    0.141637, -2.558953],
            ['cg00000236',     4267.0,        651.0,    0.850339,  2.712493],
            ['rs9363764',      3642.0,       2797.0,    0.556966,  0.380851],
            ['rs939290',        380.0,       7564.0,    0.047240, -4.315078],
            ['rs951295',       6391.0,       5675.0,    0.525316,  0.171421],
            ['rs966367',        126.0,       5467.0,    0.022132, -5.439254],
            ['rs9839873',      5485.0,        189.0,    0.949948,  4.859034],
        ]
        test2_ref = pd.DataFrame(data=test2_ref, columns=['IlmnID', 'noob_meth',  'noob_unmeth',  'beta_value',   'm_value']).set_index('IlmnID').astype('float32')
        test2 = test_data_containers[0]._SampleDataContainer__data_frame
        test2_sub = test2.loc[ test2.index.isin(test2_ref.index) ].astype('float32')
        if not np.isclose(test2_sub, test2_ref, atol=0.01).all():
            raise AssertionError()

        for outfile in test_outputs:
            if not outfile.exists():
                raise FileNotFoundError(f"Expected {outfile.name} to be generated by run_pipeline() but it was missing.")
            else:
                print('+', outfile)
                outfile.unlink()

    def test_run_pipeline_with_create_sample_sheet_plus(self):
        test_data_dir = 'docs/example_data/epic_plus'
        test_data_containers = pipeline.run_pipeline(test_data_dir, export=False, sample_name=['Sample_1'],
            meta_data_frame=False, make_sample_sheet=True, sesame=False)
        test1 = test_data_containers[0]._SampleDataContainer__data_frame
        #test2 = test_data_containers[1]._SampleDataContainer__data_frame
        test1_ref = [
            ['cg00000029_II_F_C_rep1_EPIC',    1180.0,        274.0,    0.759331,  2.106539],
            ['cg00000103_II_F_C_rep1_EPIC',    1363.0,        681.0,    0.635728,  1.001059],
            ['cg00000109_II_F_C_rep1_EPIC',    2427.0,        284.0,    0.863394,  3.095211],
            ['cg00000155_II_F_C_rep1_EPIC',    4392.0,        212.0,    0.933673,  4.372742],
            ['cg00000158_II_F_C_rep1_EPIC',    5742.0,        211.0,    0.948620,  4.766238],
            ['rs951295_I_N_C_rep1_EPIC',        210.0,       7981.0,    0.025329, -5.248108],
            ['rs966367_II_F_C_rep1_GWG1',      1587.0,       1569.0,    0.487408,  0.016457],
            ['rs966367_II_N_C_rep1_EPIC',      2419.0,       2774.0,    0.457019, -0.197557],
            ['rs9839873_II_F_C_rep1_GWG1',     3430.0,        173.0,    0.926276,  4.309365],
            ['rs9839873_II_N_C_rep1_EPIC',     3558.0,        170.0,    0.929467,  4.387460],
        ]
        test1_ref = pd.DataFrame(data=test1_ref, columns=['IlmnID', 'noob_meth',  'noob_unmeth',  'beta_value',   'm_value']).set_index('IlmnID').astype('float32')
        # spot checking the output.
        test1_sub = test1.loc[ test1.index.isin(test1_ref.index) ].astype('float32')
        if not np.isclose(test1_sub, test1_ref, atol=0.01).all():
            # test1_ref.equals( test1_sub ) will fail, and pd.testing.assert_frame_equal(test1_sub, test1_ref) fails because of rounding at 7th decimal place.
            raise AssertionError("data container values don't match")

    def test_run_pipeline_sesame_defaults(self):
        """ check that we get back useful data.
        checks SDC, CSV outputs, and pickles after sesame=True processing
        check that output files exist, then remove them.
        """
        test_data_dir = 'docs/example_data/GSE69852'
        test_outputs = [
            Path(test_data_dir, 'control_probes.pkl'),
            Path(test_data_dir, 'beta_values.pkl'),
            Path(test_data_dir, 'm_values.pkl'),
            Path(test_data_dir, 'meth_values.pkl'),
            Path(test_data_dir, 'unmeth_values.pkl'),
            Path(test_data_dir, 'noob_meth_values.pkl'),
            Path(test_data_dir, 'noob_unmeth_values.pkl'),
            Path(test_data_dir, 'sample_sheet_meta_data.pkl'),
            Path(test_data_dir, '9247377085', '9247377085_R04C02_processed.csv'),
            Path(test_data_dir, '9247377093', '9247377093_R02C01_processed.csv'),
            ]
        for outfile in test_outputs:
            if outfile.exists():
                outfile.unlink()

        test_data_containers = pipeline.run_pipeline(test_data_dir, sesame=True, export=True)
        test_probes = ['cg00063477', 'cg00121626', 'cg00223952', 'cg27614706', 'cg27619353', 'cg27620176', 'cg27647370', 'cg27652464']
        # for version 1.4.0
        minfi_reference_data = [
            ['cg00035864',     2040.0,       4480.0,    0.308157, -1.134930],
            ['cg00061679',     5946.0,       5276.0,    0.525172,  0.172475],
            ['cg00063477',     5759.0,        315.0,    0.932783,  4.192395],
            ['cg00121626',     3811.0,       7636.0,    0.330042, -1.002648],
            ['cg00223952',      277.0,      12107.0,    0.022188, -5.449811],
            ['cg27614706',     5831.0,        265.0,    0.941091,  4.459679],
            ['cg27619353',     7466.0,      14894.0,    0.332413, -0.996324],
            ['cg27620176',    11753.0,        222.0,    0.973333,  5.726326],
            ['cg27647370',    15752.0,       2212.0,    0.872011,  2.832112],
            ['cg27652464',      656.0,      15224.0,    0.041051, -4.536508],
        ]
        minfi_ref = pd.DataFrame(minfi_reference_data, columns=['IlmnID','noob_meth','noob_unmeth','beta_value','m_value']).set_index('IlmnID')
        NaN = np.nan # this matches '9247377093_R02C01'
        reference_data_old_noob = [
            ['cg00063477',     4107.0,        172.0,           1.0,       0.960,    4.578],
            ['cg00121626',     3542.0,       3397.0,           1.0,       0.510,    0.060],
            ['cg00223952',        NaN,          NaN,           NaN,         NaN,      NaN],
            ['cg27614706',        NaN,          NaN,           NaN,         NaN,      NaN],
            ['cg27619353',     2226.0,       9714.0,           1.0,       0.186,   -2.126],
            ['cg27620176',     6057.0,         94.0,           1.0,       0.985,    6.010],
            ['cg27647370',     8897.0,        167.0,           1.0,       0.982,    5.735],
            ['cg27652464',      398.0,       8832.0,           1.0,       0.043,   -4.472],
        ]
        reference_data = [ #CSV file
            ['cg00063477',     4115.0,        172.0,           1.0,       0.960,    4.580],
            ['cg00121626',     3552.0,       3381.0,           1.0,       0.512,    0.071],
            ['cg00223952',      420.0,       7058.0,           0.0,       0.056,   -4.071],
            ['cg27614706',     3612.0,         90.0,           0.0,       0.976,    5.327],
            ['cg27619353',     2204.0,       9713.0,           1.0,       0.185,   -2.140],
            ['cg27620176',     6052.0,         94.0,           1.0,       0.985,    6.010],
            ['cg27647370',     8895.0,        167.0,           1.0,       0.982,    5.735],
            ['cg27652464',      396.0,       8829.0,           1.0,       0.043,   -4.479],
        ]
        reference_container_data = [
            ['cg00063477',     4115.0,        172.0,           1.0,       0.960,    4.580],
            ['cg00121626',     3552.0,       3381.0,           1.0,       0.512,    0.071],
            ['cg00223952',        NaN,          NaN,           NaN,       0.056,   -4.071],
            ['cg27614706',        NaN,          NaN,           NaN,       0.976,    5.327],
            ['cg27619353',     2204.0,       9713.0,           1.0,       0.185,   -2.140],
            ['cg27620176',     6052.0,         94.0,           1.0,       0.985,    6.010],
            ['cg27647370',     8895.0,        167.0,           1.0,       0.982,    5.735],
            ['cg27652464',      396.0,       8829.0,           1.0,       0.043,   -4.479],
        ]
        ref_noobfix_data = [ #CSV file; pre v1.5.0 NOOB wasn't excluding poobah/qualityMask failed probes from oobG/oobR
            ['cg00063477',     4125.0,        171.0,           1.0,       0.960,    4.592],
            ['cg00121626',     3562.0,       3391.0,           1.0,       0.512,    0.071],
            ['cg00223952',      428.0,       7068.0,           0.0,       0.057,   -4.046],
            ['cg27614706',     3622.0,         87.0,           0.0,       0.977,    5.363],
            ['cg27619353',     2214.0,       9722.0,           1.0,       0.185,   -2.135],
            ['cg27620176',     6062.0,         90.0,           1.0,       0.985,    6.074],
            ['cg27647370',     8905.0,        166.0,           1.0,       0.982,    5.745],
            ['cg27652464',      404.0,       8840.0,           1.0,       0.043,   -4.452],
        ]
        ref_v163_noobfix_data = [ #CSV file; pre v1.5.0 NOOB wasn't excluding poobah/qualityMask failed probes from oobG/oobR
            ['cg00063477',     4118.0,        172.0,           1.0,       0.960,    4.581],
            ['cg00121626',     3554.0,       3402.0,           1.0,       0.512,    0.063],
            ['cg00223952',      426.0,       7070.0,           0.0,       0.057,   -4.046],
            ['cg27614706',     3633.0,         88.0,           0.0,       0.977,    5.368],
            ['cg27619353',     2228.0,       9723.0,           1.0,       0.185,   -2.126],
            ['cg27620176',     6066.0,         90.0,           1.0,       0.985,    6.075],
            ['cg27647370',     8906.0,        167.0,           1.0,       0.982,    5.737],
            ['cg27652464',      406.0,       8841.0,           1.0,       0.043,   -4.445],
        ]
        ref_v165_csv_data = [ #CSV file; pre v1.5.0 NOOB wasn't excluding poobah/qualityMask failed probes from oobG/oobR
            ['cg00063477',     4120.0,        172.0,           1.0,       0.960,    4.582],
            ['cg00121626',     3554.0,       3402.0,           1.0,       0.511,    0.063],
            ['cg00223952',      426.0,       7068.0,           0.0,       0.057,   -4.052],
            ['cg27614706',     3633.0,         88.0,           0.0,       0.976,    5.368],
            ['cg27619353',     2230.0,       9722.0,           1.0,       0.187,   -2.124],
            ['cg27620176',     6064.0,         90.0,           1.0,       0.985,    6.074],
            ['cg27647370',     8905.0,        167.0,           1.0,       0.982,    5.737],
            ['cg27652464',      406.0,       8840.0,           1.0,       0.044,   -4.444],
        ]

        ref_noobfix_container_data = [ #CSV file; pre v1.5.0 NOOB wasn't excluding poobah/qualityMask failed probes from oobG/oobR
            ['cg00063477',     4125.0,        171.0,           1.0,       0.960,    4.592],
            ['cg00121626',     3562.0,       3391.0,           1.0,       0.512,    0.071],
            ['cg00223952',        NaN,          NaN,           NaN,       0.057,   -4.046],
            ['cg27614706',        NaN,          NaN,           NaN,       0.977,    5.363],
            ['cg27619353',     2214.0,       9722.0,           1.0,       0.185,   -2.135],
            ['cg27620176',     6062.0,         90.0,           1.0,       0.985,    6.074],
            ['cg27647370',     8905.0,        166.0,           1.0,       0.982,    5.745],
            ['cg27652464',      404.0,       8840.0,           1.0,       0.043,   -4.452],
        ]
        ref_v163_container_data = [
            ['cg00063477',     4118.0,        172.0,           1.0,       0.960,    4.581],
            ['cg00121626',     3554.0,       3402.0,           1.0,       0.512,    0.063],
            ['cg00223952',        NaN,          NaN,           NaN,       0.057,   -4.053],
            ['cg27614706',        NaN,          NaN,           NaN,       0.977,    5.368],
            ['cg27619353',     2228.0,       9723.0,           1.0,       0.185,   -2.126],
            ['cg27620176',     6066.0,         90.0,           1.0,       0.985,    6.075],
            ['cg27647370',     8906.0,        167.0,           1.0,       0.982,    5.737],
            ['cg27652464',      406.0,       8841.0,           1.0,       0.043,   -4.445],
        ]
        ref_v165_container_data = [
            ['cg00063477',     4120.0,        172.0,           1.0,       0.960,    4.582],
            ['cg00121626',     3554.0,       3402.0,           1.0,       0.511,    0.063],
            ['cg00223952',        NaN,          NaN,           NaN,       0.057,   -4.053],
            ['cg27614706',        NaN,          NaN,           NaN,       0.977,    5.368],
            ['cg27619353',     2230.0,       9722.0,           1.0,       0.187,   -2.124],
            ['cg27620176',     6064.0,         90.0,           1.0,       0.985,    6.074],
            ['cg27647370',     8905.0,        167.0,           1.0,       0.982,    5.737],
            ['cg27652464',      406.0,       8840.0,           1.0,       0.044,   -4.444],
        ]

        container_ref = pd.DataFrame(ref_v165_container_data, columns=['IlmnID','noob_meth','noob_unmeth','quality_mask','beta_value','m_value']).set_index('IlmnID')
        # checking outputs.
        idata = test_data_containers[0]._SampleDataContainer__data_frame.index
        iref = container_ref.index
        subdata = test_data_containers[0]._SampleDataContainer__data_frame[idata.isin(iref)]
        # print('subdata', subdata) -- each time things change, need to update slight shifts in values
        meth = all(np.isclose(subdata[['noob_meth']], container_ref[['noob_meth']], atol=self.columns['noob_meth'], equal_nan=True))
        unmeth = all(np.isclose(subdata[['noob_unmeth']], container_ref[['noob_unmeth']], atol=self.columns['noob_unmeth'], equal_nan=True))
        beta = all(np.isclose(subdata[['beta_value']], container_ref[['beta_value']], atol=self.columns['beta_value'], equal_nan=True))
        m = all(np.isclose(subdata[['m_value']], container_ref[['m_value']], atol=self.columns['m_value'], equal_nan=True))

        if meth is False:
            raise AssertionError(f"container meth values don't match in data container:\n{subdata[['noob_meth']]}\n{container_ref[['noob_meth']]}")
        if unmeth is False:
            raise AssertionError(f"container unmeth values don't match in data container:\n{subdata[['noob_unmeth']]}\n{container_ref[['noob_unmeth']]}")
        if beta is False:
            raise AssertionError(f"container beta values don't match in data container:\n{subdata[['beta_value']]}\n{container_ref[['beta_value']]}")
        if m is False:
            raise AssertionError(f"container m values don't match in data container:\n{subdata[['m_value']]}\n{container_ref[['m_value']]}")

        csv_ref = pd.DataFrame(ref_v165_csv_data, columns=['IlmnID','noob_meth','noob_unmeth','quality_mask','beta_value','m_value']).set_index('IlmnID')
        csv_ref = csv_ref[ csv_ref.index.isin(test_probes)]
        csv_data = pd.read_csv(Path(test_data_dir, '9247377093', '9247377093_R02C01_processed.csv')).set_index('IlmnID')
        csv_data = csv_data[ csv_data.index.isin(test_probes)]
        csv_meth = all(np.isclose(csv_data[['noob_meth']], csv_ref[['noob_meth']], atol=self.columns['noob_meth'], equal_nan=True))
        csv_unmeth = all(np.isclose(csv_data[['noob_unmeth']], csv_ref[['noob_unmeth']], atol=self.columns['noob_unmeth'], equal_nan=True))
        csv_beta = all(np.isclose(csv_data[['beta_value']], csv_ref[['beta_value']], atol=self.columns['beta_value'], equal_nan=True))
        csv_m = all(np.isclose(csv_data[['m_value']], csv_ref[['m_value']], atol=self.columns['m_value'], equal_nan=True))

        if csv_meth is False:
            raise AssertionError(f"csv meth values don't match in data container:\n{csv_data[['noob_meth']]}\n{csv_ref[['noob_meth']]}")
        if csv_unmeth is False:
            raise AssertionError(f"csv unmeth values don't match in data container:\n{csv_data[['noob_unmeth']]}\n{csv_ref[['noob_unmeth']]}")
        if csv_beta is False:
            raise AssertionError(f"csv beta values don't match in data container:\n{csv_data[['beta_value']]}\n{csv_ref[['beta_value']]}")
        if csv_m is False:
            raise AssertionError(f"csv m values don't match in data container:\n{csv_data[['m_value']]}\n{csv_ref[['m_value']]}")

        #beta = pd.read_pickle(Path(test_data_dir, 'beta_values.pkl'))
        noob_meth = pd.read_pickle(Path(test_data_dir, 'noob_meth_values.pkl'))
        noob_unmeth = pd.read_pickle(Path(test_data_dir, 'noob_unmeth_values.pkl'))
        ref_meth0 = [ # pre v1.5.0
            ['cg00000029',                   2231],
            ['cg00000108',                   7880],
            ['cg00000109',                   3516],
            ['cg00000165',                    344],
            ['cg00000236',                   3601],
        ]
        ref_meth1 = [ # v1.5.0+ with noob-poobah-quality_mask fix
            ['cg00000029',                   2242],
            ['cg00000108',                   7892],
            ['cg00000109',                   3527],
            ['cg00000165',                    350],
            ['cg00000236',                   3612],
        ]
        ref_meth2 = [ # v1.6.3
            ['cg00000029',                   2235],
            ['cg00000108',                   7890],
            ['cg00000109',                   3520],
            ['cg00000165',                    350],
            ['cg00000236',                   3604],
        ]
        ref_meth = [ # v1.6.5+
            ['cg00000029',                   2235],
            ['cg00000108',                   7892],
            ['cg00000109',                   3522],
            ['cg00000165',                    350],
            ['cg00000236',                   3605],
        ]

        ref_meth = pd.DataFrame(ref_meth, columns = ['IlmnID', '9247377085_R04C02']).set_index('IlmnID')
        test_noob_meth = noob_meth['9247377085_R04C02'][noob_meth.index.isin(ref_meth.index)]
        meth = all(np.isclose(test_noob_meth, ref_meth['9247377085_R04C02'], atol=self.columns["noob_meth"], equal_nan=True))
        if meth is False:
            raise AssertionError(f"meth values don't match in pickle (ref: {ref_meth} vs {test_noob_meth})")

        test_data_dir = 'docs/example_data/GSE69852'
        test_outputs = [
            Path(test_data_dir, 'control_probes.pkl'),
            Path(test_data_dir, 'beta_values.pkl'),
            Path(test_data_dir, 'm_values.pkl'),
            Path(test_data_dir, 'meth_values.pkl'),
            Path(test_data_dir, 'unmeth_values.pkl'),
            Path(test_data_dir, 'noob_meth_values.pkl'),
            Path(test_data_dir, 'noob_unmeth_values.pkl'),
            Path(test_data_dir, 'sample_sheet_meta_data.pkl'),
            Path(test_data_dir, '9247377085', '9247377085_R04C02_processed.csv'),
            Path(test_data_dir, '9247377093', '9247377093_R02C01_processed.csv'),
            ]
        for outfile in test_outputs:
            if outfile.exists():
                outfile.unlink()

    @staticmethod
    def test_process_mouse():
        """ catches anything serious broken about mouse array processing
        in v1.4.4 / v0.7.4 I expect this to use the linear dye fallback within the sesame method, because of dupe probe names."""
        PATH = 'docs/example_data/mouse'
        ID = '204879580038_R06C02'
        print('* loading mouse manifest')
        import methylprep
        manifest = methylprep.files.Manifest(methylprep.models.ArrayType('mouse'))
        print('* loading one idat pair of files')
        green_filepath = Path(PATH, f'{ID}_Grn.idat') #'204879580038_R06C02_Grn.idat')
        red_filepath = Path(PATH, f'{ID}_Red.idat') #'204879580038_R06C02_Red.idat')
        print(f"* GREEN --> {green_filepath.name}")
        print(f"* RED --> {red_filepath.name}")
        if not (green_filepath.exists() and green_filepath.is_file()):
            raise FileNotFoundError("mouse test data missing")
        if not (red_filepath.exists() and red_filepath.is_file()):
            raise FileNotFoundError("mouse test data missing")
        files_to_remove = ['samplesheet.csv', 'control_probes.pkl', 'mouse_probes.pkl',
            'sample_sheet_meta_data.pkl', 'noob_meth_values.pkl', 'noob_unmeth_values.pkl']
        for _file in files_to_remove:
            if Path(PATH, _file).is_file():
                Path(PATH, _file).unlink()
        data = methylprep.run_pipeline(PATH, make_sample_sheet=True)
        df = data[0]._SampleDataContainer__data_frame
        #print( np.isclose(list(df['beta_value'][:3]), [0.905712,  0.841185,  0.129731]) )
        #assert np.isclose(list(df['beta_value'][:3]), [0.905712,  0.841185,  0.129731]).all() == True
        for _file in files_to_remove:
            if Path(PATH, _file).is_file():
                Path(PATH, _file).unlink()

        # compare with 10000 randomly selected probes processed in sesame
        # SEE docs/debug_notebooks/sesame_mouse.R
        # sesame data came from openSesame() without any other kwargs
        # test = ses.loc[ sesame_mask].loc[ ~ses['beta_value'].isna() ].sample(10000)
        # 10000 is a sample of all probes that are found in methylprep betas output | excludes mouse probes and rs probes.
        sesame = pd.read_csv(Path(PATH,'open_sesame_mouse_betas_subdata.csv')).set_index('IlmnID')
        # pd.read_csv(Path(PATH,f'sesame_{ID}_beta_subset.csv')).set_index('IlmnID')
        # because sesame's output is ALL probes, I need to filter to just the overlapping ones
        sesame_mask = set(sesame.index) & set(df.index)
        df_test = df[ df.index.isin(sesame_mask) ]
        diff = (df_test[['beta_value']] - sesame)
        # diff.hist(bins=300, range=[-0.1, 0.1])
        # import matplotlib.pyplot as plt;plt.show()
        if diff.mean()[0] > 0.01: # actual is 0.0350574
            raise AssertionError(f"sesame betas are too different from methylprep's for mouse data. Mean beta diff: {df_test.mean()}")


class UnitTestCase(unittest.TestCase):
    def test_pipeline_wrong_sample_name_fails(self):
        LOCAL = Path('docs/example_data/GSE69852/')
        with self.assertRaises(SystemExit):
            pipeline.run_pipeline(LOCAL, betas=True, sample_name=['blahblah_wrong_sample_name'])
