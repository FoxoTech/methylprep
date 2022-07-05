import os
import sys
import numpy as np
import pandas as pd
from pathlib import Path
# App
import methylprep
#patching
import unittest
try:
    # python 3.4+ should use builtin unittest.mock not mock package
    from unittest.mock import patch
except ImportError:
    from mock import patch


class TestPipeline():

    @staticmethod
    def test_pipeline_two_samples():
        """ pass in --sample_name with 2 samples -- from fake command line args """
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
            #Path(test_data_dir, '9247377085', '9247377085_R04C02_processed.csv'),
            #Path(test_data_dir, '9247377093', '9247377093_R02C01_processed.csv'),
            ]
        for outfile in test_outputs:
            if outfile.exists():
                outfile.unlink()
        # CLI not working in unit tests
        #testargs = ["__program__", '-d', test_data_dir, '--no_export', '--sample_name', 'AdultLiver1', 'FetalLiver1', '--minfi', '--betas']
        #with patch.object(sys, 'argv', testargs):
        test_data_containers = methylprep.run_pipeline(test_data_dir, export=False, sample_name=['AdultLiver1', 'FetalLiver1'])
        # spot checking the output.
        test1 = test_data_containers[0]._SampleDataContainer__data_frame
        test1_v162_ref = [
            ['cg00000029',     1144.0,       1369.0,        0.003,           1.0,    0.455233, -0.259035],
            ['cg00000108',     3841.0,       3926.0,        0.000,           1.0,    0.494528, -0.031578],
            ['cg00000109',     3828.0,        190.0,        0.000,           1.0,    0.952713,  4.332520],
            ['cg00000165',      555.0,       2580.0,        0.006,           1.0,    0.177033, -2.216811],
            ['cg00000236',     3545.0,        428.0,        0.000,           1.0,    0.892273,  3.050103],
            ['rs9363764',      3038.0,       2009.0,        0.001,           1.0,    0.601942,  0.596644],
            ['rs939290',        360.0,       5615.0,        0.001,           1.0,    0.060251, -3.963217],
            ['rs951295',       5317.0,       4706.0,        0.000,           1.0,    0.530480,  0.176111],
            ['rs966367',        120.0,       4104.0,        0.002,           1.0,    0.028409, -5.095924],
            ['rs9839873',      4545.0,        112.0,        0.000,           1.0,    0.975950,  5.342710]
        ]
        test1_v163_ref = [
            ['cg00000029',     1137.0,       1380.0,        0.003,           1.0,    0.451, -0.279],
            ['cg00000108',     3834.0,       3935.0,        0.000,           1.0,    0.493, -0.037],
            ['cg00000109',     3821.0,        190.0,        0.000,           1.0,    0.953,  4.330],
            ['cg00000165',      551.0,       2593.0,        0.006,           1.0,    0.175, -2.234],
            ['cg00000236',     3538.0,        431.0,        0.000,           1.0,    0.891,  3.037],
            ['rs9363764',      3030.0,       2021.0,        0.001,           1.0,    0.599,  0.584],
            ['rs939290',        359.0,       5619.0,        0.001,           1.0,    0.060, -3.968],
            ['rs951295',       5312.0,       4700.0,        0.000,           1.0,    0.530,  0.177],
            ['rs966367',        120.0,       4112.0,        0.002,           1.0,    0.028, -5.099],
            ['rs9839873',      4540.0,        112.0,        0.000,           1.0,    0.976,  5.341]
        ]
        """ FULL 6-digit values for 1.6.3
        0.451728 -0.279436
        0.493500 -0.037513
        0.952630  4.329879
        0.175254 -2.234498
        0.891408  3.037174
        0.599881  0.584248
        0.060054 -3.968258
        0.530563  0.176594
        0.028355 -5.098734
        0.975924  5.341122
        """
        test1_v165_ref = [
            ['cg00000029',     1136.0,       1382.0,        0.003,           1.0,    0.451, -0.283],
            ['cg00000108',     3834.0,       3934.0,        0.000,           1.0,    0.493, -0.037],
            ['cg00000109',     3822.0,        190.0,        0.000,           1.0,    0.953,  4.330],
            ['cg00000165',      551.0,       2594.0,        0.006,           1.0,    0.175, -2.235],
            ['cg00000236',     3538.0,        432.0,        0.000,           1.0,    0.891,  3.034],
            ['rs9363764',      3030.0,       2024.0,        0.001,           1.0,    0.600,  0.582],
            ['rs939290',        359.0,       5618.0,        0.001,           1.0,    0.060, -3.968],
            ['rs951295',       5314.0,       4702.0,        0.000,           1.0,    0.530,  0.177],
            ['rs966367',        120.0,       4112.0,        0.002,           1.0,    0.028, -5.099],
            ['rs9839873',      4542.0,        112.0,        0.000,           1.0,    0.976,  5.341]
        ]
        test1_ref = pd.DataFrame(data=test1_v165_ref, columns=['IlmnID', 'noob_meth',  'noob_unmeth', 'poobah_pval', 'quality_mask', 'beta_value',   'm_value']).set_index('IlmnID').astype('float32')
        test1_sub = test1.loc[ test1.index.isin(test1_ref.index) ].astype('float32')

        if not np.isclose(test1_sub, test1_ref, atol=0.01).all():
            print('REF',test1_ref)
            print('SUB',test1_sub)
            raise AssertionError()

        for outfile in test_outputs:
            if not outfile.exists():
                raise FileNotFoundError(f"Expected {outfile.name} to be generated by run_pipeline() but it was missing.")
            else:
                print('+', outfile)
                outfile.unlink()
