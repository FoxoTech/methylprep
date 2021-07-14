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
        test_data_containers = pipeline.run_pipeline(test_data_dir, export=False, sample_name=['AdultLiver1', 'FetalLiver1'])
        # spot checking the output.
        test1 = test_data_containers[0]._SampleDataContainer__data_frame
        test1_ref = [
            ['cg00000029',     1144.0,       1369.0,        0.003,           0.0,    0.455233, -0.259035],
            ['cg00000108',     3841.0,       3926.0,        0.000,           0.0,    0.494528, -0.031578],
            ['cg00000109',     3828.0,        190.0,        0.000,           0.0,    0.952713,  4.332520],
            ['cg00000165',      555.0,       2580.0,        0.006,           0.0,    0.177033, -2.216811],
            ['cg00000236',     3545.0,        428.0,        0.000,           0.0,    0.892273,  3.050103],
            ['rs9363764',      3038.0,       2009.0,        0.001,           0.0,    0.601942,  0.596644],
            ['rs939290',        360.0,       5615.0,        0.001,           0.0,    0.060251, -3.963217],
            ['rs951295',       5317.0,       4706.0,        0.000,           0.0,    0.530480,  0.176111],
            ['rs966367',        120.0,       4104.0,        0.002,           0.0,    0.028409, -5.095924],
            ['rs9839873',      4545.0,        112.0,        0.000,           0.0,    0.975950,  5.342710]
        ]
        test1_ref = pd.DataFrame(data=test1_ref, columns=['IlmnID', 'noob_meth',  'noob_unmeth', 'poobah_pval', 'quality_mask', 'beta_value',   'm_value']).set_index('IlmnID').astype('float32')
        test1_sub = test1.loc[ test1.index.isin(test1_ref.index) ].astype('float32')
        if not np.isclose(test1_sub, test1_ref, atol=0.01).all():
            raise AssertionError()

        for outfile in test_outputs:
            if not outfile.exists():
                raise FileNotFoundError(f"Expected {outfile.name} to be generated by run_pipeline() but it was missing.")
            else:
                print('+', outfile)
                outfile.unlink()
