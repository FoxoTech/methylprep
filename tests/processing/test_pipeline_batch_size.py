from methylprep.processing import pipeline
import pandas as pd
from pathlib import Path
import numpy as np

class TestBatchSize():

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

    @staticmethod
    def test_with_batch_size():
        test_data_dir = 'docs/example_data/GSE69852'
        df = pipeline.run_pipeline(test_data_dir, export=True, batch_size=1, sample_name='AdultLiver1')
        #if not np.isclose(df[0].unmethylated.data_frame.iloc[0]['noob'], 4119.633578946326):
        if not np.isclose(df[0]._SampleDataContainer__data_frame.iloc[0]['beta_value'], 0.236):
            raise AssertionError()
        if not np.isclose(df[0]._SampleDataContainer__data_frame.iloc[2]['m_value'], 4.146):
            raise AssertionError()
        if not Path('docs/example_data/GSE69852/9247377093/9247377093_R02C01_processed.csv').is_file():
            raise AssertionError()
        print(f"TEST OUTPUT FILES: {list(Path(test_data_dir).rglob('*'))}")
        for file in Path(test_data_dir).rglob('*.pkl'):
            file.unlink()


    @staticmethod
    def test_batch_size_betas():
        test_data_dir = 'docs/example_data/GSE69852'
        betas = pipeline.run_pipeline(test_data_dir, betas=True, batch_size=1)
        if not np.isclose(betas.iloc[0]['9247377093_R02C01'], 0.23624517): #0.23623395577166542):
            raise AssertionError(f"{betas.iloc[0]['9247377093_R02C01']} != 0.23623395577166542")
        if not Path(test_data_dir, 'beta_values.pkl').is_file():
            raise AssertionError()
        print(f"TEST OUTPUT FILES: {list(Path(test_data_dir).rglob('*'))}")
        for file in Path(test_data_dir).rglob('*.pkl'):
            file.unlink()

    """
    @staticmethod
    def test_batch_size_m_values():
        test_data_dir = 'docs/example_data/GSE69852'
        m_values = pipeline.run_pipeline(test_data_dir, m_value=True, batch_size=1)
        if not np.isclose(m_values.iloc[0]['9247377093_R02C01'], -1.6575579144514734):
            raise AssertionError()
        if not (Path('m_values_1.pkl').is_file() and Path('m_values_2.pkl').is_file()):
            raise AssertionError()
        if not Path('m_values.pkl').is_file():
            raise AssertionError()
    """
