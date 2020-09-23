import pytest
import sys
import pandas as pd
from pathlib import Path
import os
import shutil
# App
from methylprep.download import process_data
from methylprep.download import convert_miniml
from methylprep.download import run_series
#patching
try:
    # python 3.4+ should use builtin unittest.mock not mock package
    from unittest.mock import patch
except ImportError:
    from mock import patch



class TestProcessData():

    @staticmethod
    def test_download():
        """this is a full integration test of `download` function on smallest data set on GEO
        CLI: python -m methylprep download -i GSE123754 -d GSE123754 --dict_only --no_decompress
        pass in fake command line args.
        verify files exist, then remove them. confirm time takes ~2 (<5) mins"""
        test_geo_id = 'GSE123754'
        test_data_dir = 'docs/example_data/GSE123754'

        exit_status = os.system(f'python -m methylprep download -i {test_geo_id} -d {test_data_dir} --dict_only --no_decompress')
        if exit_status != 0:
            shutil.rmtree(test_data_dir)
            raise AssertionError("methylprep download exited with error(s)")
        # verify files exist
        expected_dirs = [
            Path(test_data_dir),
            Path(test_data_dir, 'GPL13534')
        ]
        expected_files = [
            # Path(test_data_dir, 'GSE123754_RAW.tar'), -- if clean==False, this remains
            Path(test_data_dir, 'GSE123754_family.xml'),
            #Path(test_data_dir, 'GSE123754_family.xml.tgz'), -- if clean==False, this remains
            Path(test_data_dir, 'GPL13534', 'GSE123754_GPL13534_meta_data.pkl'),
            Path(test_data_dir, 'GPL13534', 'GSE123754_GPL13534_samplesheet.csv'),
            Path(test_data_dir, 'GPL13534', 'GSM3510888_9406921128_R06C02_Grn.idat.gz'),
            Path(test_data_dir, 'GPL13534', 'GSM3510888_9406921128_R06C02_Red.idat.gz'),
            Path(test_data_dir, 'GPL13534', 'GSM3510887_9406921128_R04C01_Grn.idat.gz'),
            Path(test_data_dir, 'GPL13534', 'GSM3510887_9406921128_R04C01_Grn.idat.gz'),
        ]
        for dir in expected_dirs:
            if not (dir.exists() and dir.is_dir()):
                #shutil.rmtree(test_data_dir)
                raise AssertionError(f"folder {dir} does not exist after running `download`")
        for file in expected_files:
            if not file.exists():
                #shutil.rmtree(test_data_dir)
                raise AssertionError(f"file {file} does not exist after running `download`")
        print("all tests passed. deleting files now")
        shutil.rmtree(test_data_dir)


    @staticmethod
    def test_run_series():
        """ a simple small GEO dataset to practice downloading. It is 27k so it won't fully parse.
        But obvious errors will be caught with this.
        python -m methylprep -v download -i GSE120341 -d GSE120341
        has 44 idats. GSE17769 has 20 idats."""
        geo_id = 'GSE17769'
        test_data_dir = f'docs/example_data/{geo_id}' # created by this in test environment
        testargs = ["__program__", '-i', geo_id, '-d', test_data_dir]
        with patch.object(sys, 'argv', testargs):
            process_data.run_series(geo_id, test_data_dir)
            # count idats returned
            files_found = list(Path(test_data_dir).rglob('*.idat'))
            if len(files_found) != 20:
                shutil.rmtree(test_data_dir)
                raise AssertionError(f"Not enough idat files appeared after downloading")
            # cleanup
            print(f"{len(files_found)} files found")
            for file in files_found:
                file.unlink()
            for file in Path(test_data_dir).rglob('*.xml'):
                file.unlink()
            #process_data.cleanup(test_data_dir) # removes empty folders
            for folder in Path(test_data_dir).rglob('*'):
                non_empty_dirs = {str(p.parent) for p in Path(folder).rglob('*') if p.is_file()}
                if non_empty_dirs == set():
                    Path(folder).rmdir()
            Path(test_data_dir).rmdir()


    @staticmethod
    def test_convert_miniml():
        geo_id = 'GSE17769'
        test_data_dir = f'docs/example_data/{geo_id}' # created by this in test environment
        testargs = ["__program__", '-i', geo_id, '-d', test_data_dir]
        with patch.object(sys, 'argv', testargs):
            convert_miniml(
                geo_id,
                data_dir=test_data_dir,
                merge=True,
                download_it=True,
                extract_controls=False,
                require_keyword=None,
                sync_idats=True,
                remove_tgz=True,
                verbose=True)
            for file in Path(test_data_dir).rglob('*.xml'):
                file.unlink()
            Path(test_data_dir).rmdir()

    @staticmethod
    def test_convert_miniml_keyword():
        geo_id = 'GSE123754'
        gsm_id = 'GSM3510888'
        test_data_dir = f'docs/example_data/{geo_id}' # created by this in test environment
        testargs = ["__program__", '-i', geo_id, '-d', test_data_dir, '-k', 'SKTR']
        with patch.object(sys, 'argv', testargs):
            convert_miniml(
                geo_id,
                data_dir=test_data_dir,
                merge=True,
                download_it=True,
                extract_controls=False,
                require_keyword='SKTR',
                sync_idats=True,
                remove_tgz=True,
                verbose=True)
            files_found = list(Path(test_data_dir).rglob('*'))
            if len(files_found) != 3:
                raise AssertionError("Did not download all the files.")
            samplesheet = pd.read_csv(Path(test_data_dir,f'{geo_id}_GPL13534_samplesheet.csv'))
            if samplesheet['GSM_ID'][0] != gsm_id:
                raise AssertionError("Samplesheet did not contain the right data.")
            if len(samplesheet['GSM_ID']) != 1:
                raise AssertionError("Keyword filtering failed")
            for file in Path(test_data_dir).rglob('*.xml'):
                file.unlink()
            for file in Path(test_data_dir).rglob('*.pkl'):
                file.unlink()
            for file in Path(test_data_dir).rglob('*.csv'):
                file.unlink()
            Path(test_data_dir).rmdir()

    @staticmethod
    def test_convert_miniml_control():
        geo_id = 'GSE99650'
        gsm_id = 'GSM2649283'
        test_data_dir = f'docs/example_data/{geo_id}' # created by this in test environment
        testargs = ["__program__", '-i', geo_id, '-d', test_data_dir, '-k', 'SKTR']
        with patch.object(sys, 'argv', testargs):
            convert_miniml(
                geo_id,
                data_dir=test_data_dir,
                merge=True,
                download_it=True,
                extract_controls=True,
                require_keyword=None,
                sync_idats=True,
                remove_tgz=True,
                verbose=True)
            files_found = list(Path(test_data_dir).rglob('*'))
            if len(files_found) != 3:
                raise AssertionError("Did not download all the files.")
            samplesheet = pd.read_csv(Path(test_data_dir,f'{geo_id}_GPL13534_samplesheet.csv'))
            if samplesheet['GSM_ID'][0] != gsm_id:
                raise AssertionError("Samplesheet did not contain the right data.")
            if len(samplesheet['GSM_ID']) != 1:
                raise AssertionError("Control filtering failed")
            for file in Path(test_data_dir).rglob('*.xml'):
                file.unlink()
            for file in Path(test_data_dir).rglob('*.pkl'):
                file.unlink()
            for file in Path(test_data_dir).rglob('*.csv'):
                file.unlink()
            for file in Path(test_data_dir).rglob('*.tgz'):
                file.unlink()
            Path(test_data_dir).rmdir()

    @staticmethod
    def _disabled_for_now_test_ae_download():
        """ small dataset. no processing done. """
        ae_id = 'E-MTAB-6331'
        test_data_dir = f'docs/example_data/{ae_id}' # created by this in test environment
        testargs = ["__program__", '-i', ae_id, '-d', test_data_dir]
        with patch.object(sys, 'argv', testargs):
            process_data.run_series(ae_id, test_data_dir)
            # count idats returned
            files_found = list(Path(test_data_dir).rglob('*.idat'))
            if len(files_found) != 12:
                raise AssertionError()
            # cleanup
            for idat in files_found:
                idat.unlink()
            for file in Path(test_data_dir).rglob('*'):
                if file.is_file():
                    file.unlink()
            for folder in Path(test_data_dir).rglob('*'):
                non_empty_dirs = {str(p.parent) for p in Path(folder).rglob('*') if p.is_file()}
                if non_empty_dirs == set():
                    Path(folder).rmdir()
            Path(test_data_dir).rmdir()
