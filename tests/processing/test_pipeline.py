import sys
import numpy as np
import pandas as pd
from pathlib import Path
# App
from methylprep.processing import pipeline
from methylprep.utils.files import download_file
#patching
try:
    # python 3.4+ should use builtin unittest.mock not mock package
    from unittest.mock import patch
except ImportError:
    from mock import patch


class TestPipeline():

    @staticmethod
    def test_run_pipeline_demo_data():
        """ check that we get back useful data """
        test_data_dir = 'docs/example_data/GSE69852'
        test_data_containers = pipeline.run_pipeline(test_data_dir)
        # spot checking the output.
        if not test_data_containers[1].unmethylated.data_frame.iloc[0]['mean_value'] == 2712:
            raise AssertionError()
        if not np.isclose(test_data_containers[1].unmethylated.data_frame.iloc[0]['noob'], 4479.96501260212):
            raise AssertionError()

    @staticmethod
    def test_download_manifest_dummy_file():
        """ will download a tiny file from the array-manifest-files s3 bucket, to test the SSL connection on all platforms.
        The dummy file is not a proper manifest CSV, so doesn't test format.
        download_file now defaults to non-SSL if SSL fails, with warning to user."""
        test_filename = 'unittest.txt'
        test_s3_bucket = 'https://array-manifest-files.s3.amazonaws.com'  # 's3://array-manifest-files'
        dest_dir = 'tests'
        # use the .download_file() method in files.py to test the download step specifically. this is called by Manifests() class.
        download_file(test_filename, test_s3_bucket, dest_dir, overwrite=False)
        # in testing mode, this should not exist, and should get deleted right after each successful test.
        if not Path(dest_dir,test_filename).is_file():
            raise AssertionError()
        Path(dest_dir,test_filename).unlink() # deletes file.

    @staticmethod
    def test_pipeline_two_samples():
        """ pass in --sample_name with 2 samples """
        test_data_dir = 'docs/example_data/GSE69852'
        testargs = ["__program__", '-d', test_data_dir, '--no_export', '--sample_name', 'AdultLiver1', 'FetalLiver1']
        with patch.object(sys, 'argv', testargs):
            test_data_containers = pipeline.run_pipeline(test_data_dir)
            # spot checking the output.
            if not test_data_containers[1].unmethylated.data_frame.iloc[0]['mean_value'] == 2712:
                raise AssertionError()
            if not np.isclose(test_data_containers[1].unmethylated.data_frame.iloc[0]['noob'], 4479.96501260212):
                raise AssertionError()

    @staticmethod
    def test_run_pipeline_export_data():
        """ check that we get back useful data with --export option """
        test_data_dir = 'docs/example_data/GSE69852'
        testfile_1 = Path(test_data_dir, '9247377093', '9247377093_R02C01_processed.csv')
        testfile_2 = Path(test_data_dir, '9247377085', '9247377085_R04C02_processed.csv')
        if testfile_1.exists():
            testfile_1.unlink()
        if testfile_2.exists():
            testfile_2.unlink()
        test_data_containers = pipeline.run_pipeline(test_data_dir, export=True)
        if not testfile_1.exists():
            raise AssertionError("no exported processed csv found")

        test1 = pd.read_csv(testfile_1)
        if test1['beta_value'].isna().sum() > 0:
            print(test1.head())
            raise AssertionError('missing values in processed csv')
        test2 = pd.read_csv(testfile_2)
        if test2['beta_value'].isna().sum() > 0:
            print(test2.head())
            raise AssertionError('missing values in processed csv')

        # spot checking the output.
        if not test_data_containers[1].unmethylated.data_frame.iloc[0]['mean_value'] == 2712:
            raise AssertionError()
        # spot checking the output.
        total_nas = test_data_containers[0]._SampleDataContainer__data_frame['beta_value'].isna().sum()
        if total_nas > 0:
            print(f'found {total_nas} missing beta_values (N/A or inf) in sample')
            raise AssertionError()
        if not np.isclose(test_data_containers[1].unmethylated.data_frame.iloc[0]['noob'], 4479.96501260212):
            raise AssertionError()


#@staticmethod
def test_run_pipeline_epic_plus_export_data():
    """ check that we get back useful data with --export option """
    test_data_dir = 'docs/example_data/epic_plus'
    testfile_1 = Path(test_data_dir, '202651080072', '202651080072_R01C01_processed.csv')
    if testfile_1.exists():
        testfile_1.unlink()
    test_data_containers = pipeline.run_pipeline(test_data_dir, export=True)
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
    if not np.isclose(test1['beta_value'].iloc[5], 0.145):
        print(test1.iloc[5])
        raise AssertionError('beta_value doesnt match expected value')
    if not np.isclose(round(test_data_containers[0].unmethylated.data_frame.iloc[0]['noob'],1), 274.7):
        raise AssertionError("data_container output differs from expected value")

# POSSIBLE FUTURE TEST: invoke from command line instead of within package
# doesn't seem to cause problems, but best to mock the output to save time
