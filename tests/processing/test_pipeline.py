import sys
import numpy as np
# App
from methylprep.processing import pipeline
from methylprep.utils.files import download_file
from pathlib import Path
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
