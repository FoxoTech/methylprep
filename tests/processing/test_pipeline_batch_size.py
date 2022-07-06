from methylprep.processing import pipeline
import pandas as pd
from pathlib import Path
import numpy as np
import unittest

class TestBatchSize(unittest.TestCase):

    """ TOO SLOW for CI - and redundant with test_pipeline
    @staticmethod
    def test_no_batch_size():
        test_data_dir = 'docs/example_data/GSE69852'
        df = pipeline.run_pipeline(test_data_dir, export=True)
        if not np.isclose(df[0].unmethylated.data_frame.iloc[0]['noob'], 4119.633578946326):
            raise AssertionError()
        if not Path('docs/example_data/GSE69852/9247377093/9247377093_R02C01_processed.csv').is_file():
            raise AssertionError()

    @staticmethod
    def test_no_batch_size_betas():
        test_data_dir = 'docs/example_data/GSE69852'
        betas = pipeline.run_pipeline(test_data_dir, betas=True)
        if not np.isclose(betas.iloc[0]['9247377093_R02C01'], 0.23623395577166542):
            raise AssertionError()
        if not Path('beta_values.pkl').is_file():
            raise AssertionError()

    @staticmethod
    def test_no_batch_size_m_values():
        test_data_dir = 'docs/example_data/GSE69852'
        m_values = pipeline.run_pipeline(test_data_dir, m_value=True)
        if not np.isclose(m_values.iloc[0]['9247377093_R02C01'], -1.6575579144514734):
            raise AssertionError()
        if not Path('m_values.pkl').is_file():
            raise AssertionError()
    """

    def test_pipeline_sample_name_must_be_list(self):
        test_data_dir = 'docs/example_data/GSE69852'
        with self.assertRaises(SystemExit):
            df = pipeline.run_pipeline(test_data_dir, sample_name='AdultLiver1')

    def test_with_batch_size(self):
        test_data_dir = 'docs/example_data/GSE69852'
        df = pipeline.run_pipeline(test_data_dir, export=True, batch_size=1, sample_name=['AdultLiver1'])

        ref_v145 = [
            ['cg00063477',     4115.0,        172.0,        0.000,           0.0,       0.960,    4.580],
            ['cg00121626',     3552.0,       3381.0,        0.000,           0.0,       0.512,    0.071],
            ['cg27619353',     2204.0,       9713.0,        0.000,           0.0,       0.185,   -2.140],
            ['cg27620176',     6052.0,         94.0,        0.001,           0.0,       0.985,    6.009],
        ]
        ref_v162 = [
            ['cg00063477',     4125.0,        171.0,        0.000,           0.0,       0.960,    4.592],
            ['cg00121626',     3562.0,       3391.0,        0.000,           0.0,       0.512,    0.071],
            ['cg27619353',     2214.0,       9722.0,        0.000,           0.0,       0.185,   -2.135],
            ['cg27620176',     6062.0,         90.0,        0.001,           0.0,       0.985,    6.074],
        ]
        ref_v163 = [
            ['cg00063477',     4118.0,        172.0,        0.001,           1.0,       0.960,    4.581],
            ['cg00121626',     3554.0,       3402.0,        0.001,           1.0,       0.511,    0.063],
            ['cg27619353',     2228.0,       9723.0,        0.000,           1.0,       0.186,   -2.126],
            ['cg27620176',     6066.0,         90.0,        0.001,           1.0,       0.985,    6.075],
        ]
        ref_v165 = [
            ['cg00063477',     4120.0,        172.0,        0.001,           1.0,       0.960,    4.582],
            ['cg00121626',     3554.0,       3402.0,        0.001,           1.0,       0.511,    0.063],
            ['cg27619353',     2230.0,       9722.0,        0.000,           1.0,       0.186,   -2.124],
            ['cg27620176',     6064.0,         90.0,        0.001,           1.0,       0.985,    6.074],
        ]
        ref_data = pd.DataFrame(ref_v165, columns=['IlmnID','noob_meth','noob_unmeth','poobah_pval','quality_mask','beta_value','m_value']).set_index('IlmnID')
        data = df[0]._SampleDataContainer__data_frame.loc[ref_data.index]
        #if not np.isclose(df[0].unmethylated.data_frame.iloc[0]['noob'], 4119.633578946326):
        #if not np.isclose(df[0]._SampleDataContainer__data_frame.iloc[0]['beta_value'], 0.236):
        if not np.isclose(data, ref_data, atol=1.0).all():
            print("SDC mismatch")
            print('REF',ref_data)
            print('data',data)
            raise AssertionError()
        if not Path('docs/example_data/GSE69852/9247377093/9247377093_R02C01_processed.csv').is_file():
            raise AssertionError()
        csv_data = pd.read_csv(Path(test_data_dir, '9247377093', '9247377093_R02C01_processed.csv')).set_index('IlmnID')
        test_csv = csv_data.loc[ref_data.index]
        if not np.isclose(test_csv, ref_data, atol=1.0).all():
            print("csv mismatch")
            print(test_csv)
            print(ref_data)
            raise AssertionError()
        print(f"TEST OUTPUT FILES: {list(Path(test_data_dir).rglob('*'))}")
        for file in Path(test_data_dir).rglob('*.pkl'):
            file.unlink()


    def test_batch_size_betas(self):
        test_data_dir = 'docs/example_data/GSE69852'
        betas = pipeline.run_pipeline(test_data_dir, betas=True, batch_size=1)
        ref_v145 = [
            ['cg00063477',    0.959879,           0.961307],
            ['cg00121626',    0.512332,           0.351993],
            ['cg27619353',    0.184946,           0.358009],
            ['cg27620176',    0.984706,           0.983877],
        ]
        ref_162 = [ # v1.5.0 to 1.6.2
            ['cg00063477',    0.960196,           0.961819],
            ['cg00121626',    0.512297,           0.352334],
            ['cg27619353',    0.185489,           0.358258],
            ['cg27620176',    0.985371,           0.984264],
        ]
        ref_v163 = [
            ['cg00063477',    0.959907,           0.961307],
            ['cg00121626',    0.510926,           0.351993],
            ['cg27619353',    0.186428,           0.358009],
            ['cg27620176',    0.985380,           0.983877],
        ]

        ref_data = pd.DataFrame(ref_v163, columns=['IlmnID', '9247377093_R02C01','9247377085_R04C02']).set_index('IlmnID')
        print('ref',ref_data)
        test_betas = betas.loc[ref_data.index]
        if not np.isclose(ref_data, test_betas, atol=0.01).all():
            print(f"Betas returned don't match:\n{ref_data}\ntest:\ntest_betas")
            # raise AssertionError(f"Betas returned don't match:\n{ref_data}\ntest:\ntest_betas")
        #if not np.isclose(betas.iloc[0]['9247377093_R02C01'], 0.23624517): #0.23623395577166542):
        #    raise AssertionError(f"{betas.iloc[0]['9247377093_R02C01']} != 0.23623395577166542")
        if not Path(test_data_dir, 'beta_values.pkl').is_file():
            raise AssertionError(f"beta_values.pkl not found")
        betas_pkl = pd.read_pickle(Path(test_data_dir, 'beta_values.pkl'))['9247377093_R02C01']
        if not np.isclose(betas_pkl.loc[['cg00063477', 'cg00121626','cg27619353','cg27620176']], ref_data['9247377093_R02C01'], atol=0.01).all():
            print(betas_pkl.loc[['cg00063477', 'cg00121626','cg27619353','cg27620176']])
            print(ref_data['9247377093_R02C01'])
            raise AssertionError("beta_values.pkl don't match expected test values")
        betas_pkl = pd.read_pickle(Path(test_data_dir, 'beta_values.pkl'))['9247377085_R04C02']
        if not np.isclose(betas_pkl.loc[['cg00063477', 'cg00121626','cg27619353','cg27620176']], ref_data['9247377085_R04C02'], atol=0.01).all():
            print(betas_pkl.loc[['cg00063477', 'cg00121626','cg27619353','cg27620176']])
            print(ref_data['9247377085_R04C02'])
            raise AssertionError("beta_values.pkl don't match expected test values")
        print(f"TEST OUTPUT FILES: {list(Path(test_data_dir).rglob('*'))}")
        for file in Path(test_data_dir).rglob('*.pkl'):
            file.unlink()
