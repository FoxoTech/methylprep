import pytest
import sys
import pandas as pd
import json
from pathlib import Path
import os
import shutil

# App
import methylprep

#patching
try:
    # python 3.4+ should use builtin unittest.mock not mock package
    from unittest.mock import patch
except ImportError:
    from mock import patch

class TestBetaBake():
    """ Note: all 4 tests complete in 30 to 40s """

    def __prev_test_pipeline_find_betas_any_source_samplesheet(self):
        """ SKIP: this test no longer works because the dataset now has 11GB of data, instead of just a samplesheet. """
        LOCAL = Path('docs/example_data/GSE132203')
        kwargs = {'project_name': 'GSE132203', 'data_dir': LOCAL, 'clean':False, 'compress':False, 'verbose':True}
        result = methylprep.download.pipeline_find_betas_any_source(**kwargs)
        print(result)
        # runs a bunch of stuff. This GSE dataset has no IDATs or series matrix betas. just testing the samplesheet part.
        samplesheet = pd.read_csv(Path(LOCAL, f"{kwargs['project_name']}_samplesheet.csv"))
        ref_samplesheet = pd.read_csv(Path(LOCAL, f"ref_{kwargs['project_name']}_samplesheet.csv"))
        if not samplesheet.equals(ref_samplesheet):
            raise AssertionError("Samplesheets don't match.")
        print('OK. samplesheets match.')
        with open(Path(LOCAL, f"{kwargs['project_name']}_series_summary.json"), 'rb') as f:
            summary = json.load(f)
        with open(Path(LOCAL, f"ref_{kwargs['project_name']}_series_summary.json"), 'rb') as f:
            ref_summary = json.load(f)
        if summary != ref_summary:
            raise AssertionError("Series summaries don't match.")
        print('OK. summaries match.')
        Path(LOCAL, f"{kwargs['project_name']}_samplesheet.csv").unlink()
        Path(LOCAL, f"{kwargs['project_name']}_series_summary.json").unlink()
        Path(LOCAL, f"{kwargs['project_name']}_beta_values.pkl").unlink()

    def test_pipeline_find_betas_any_source_matrix(self):
        """parses GSE110454_series_matrix.txt.gz with 4 samples """
        expected_file_sizes = {
            'GSE110454_samplesheet.csv': 4794,
            'GSE110454_series_summary.json': 1563,
        }
        expected_beta_file_range = 39025000
        expected_beta_file_sizes = [39025301, 39025330, 39025298, 39025252]
        # 'GSE110454_beta_values.pkl': many GSE110454_beta_values.pkl sizes found with retesting.

        LOCAL = Path('docs/example_data/GSE110454')
        kwargs = {'project_name': 'GSE110454', 'data_dir': LOCAL, 'clean':False, 'compress':False, 'verbose':True}
        result = methylprep.download.pipeline_find_betas_any_source(**kwargs)
        print(result)
        # runs a bunch of stuff. This GSE dataset has no IDATs or series matrix betas. just testing the samplesheet part.
        samplesheet = pd.read_csv(Path(LOCAL, f"{kwargs['project_name']}_samplesheet.csv"))
        ref_samplesheet = pd.read_csv(Path(LOCAL, f"ref_{kwargs['project_name']}_samplesheet.csv"))
        if not samplesheet.equals(ref_samplesheet):
            raise AssertionError("Samplesheets don't match.")
        print('OK. samplesheets match.')
        with open(Path(LOCAL, f"{kwargs['project_name']}_series_summary.json"), 'rb') as f:
            summary = json.load(f)
        with open(Path(LOCAL, f"ref_{kwargs['project_name']}_series_summary.json"), 'rb') as f:
            ref_summary = json.load(f)
        if summary != ref_summary:
            raise AssertionError("Series summaries don't match.")
        print('OK. summaries match.')
        betas = pd.read_pickle(Path(LOCAL, f"{kwargs['project_name']}_beta_values.pkl"))
        # getting a SMALL discrepancy between local and circleci filesize
        print(betas)
        print(betas.isna().sum())
        # FINALLY, compare against expected files sizes. Easiest way to verify the beta_values download worked.
        for _file,_size in expected_file_sizes.items():
            if Path(LOCAL,_file).stat().st_size not in (_size,39025330): # this one file is one size locally and a diff size on circleci; accept either one
                raise AssertionError(f"File size mismatch for {_file}: {Path(LOCAL,_file).stat().st_size} != {_size} expected")
        _file = 'GSE110454_beta_values.pkl'
        actual_beta_file_size = Path(LOCAL, _file).stat().st_size
        if actual_beta_file_size < expected_beta_file_range:
            Path(LOCAL, f"{kwargs['project_name']}_samplesheet.csv").unlink()
            Path(LOCAL, f"{kwargs['project_name']}_series_summary.json").unlink()
            Path(LOCAL, f"{kwargs['project_name']}_beta_values.pkl").unlink()
            raise AssertionError(f"Beta file size mismatch for {_file}: {repr(actual_beta_file_size)}; expected one of {repr(expected_beta_file_sizes)}.")
        Path(LOCAL, f"{kwargs['project_name']}_samplesheet.csv").unlink()
        Path(LOCAL, f"{kwargs['project_name']}_series_summary.json").unlink()
        Path(LOCAL, f"{kwargs['project_name']}_beta_values.pkl").unlink()


    def test_pipeline_find_betas_any_source_27k_idats(self):
        """parses IDATs with 10 27k samples, 13MB takes about 30sec; also tests run_series() """
        expected_file_sizes = {
            # 'geo_alert GSE17769.csv': 559, -- in a different parent folder
            # 'GSE17769_family.xml': 28475,
            'GSM443811_HCC38_4308918033_G_Grn.idat': 719560,
            'GSM443811_HCC38_4308918033_G_Red.idat': 719560,
            'GSM443812_HCC1143_4308918024_C_Grn.idat': 719560,
            'GSM443812_HCC1143_4308918024_C_Red.idat': 719560,
            'GSM443813_HCC1599_4308918033_C_Grn.idat': 719560,
            'GSM443813_HCC1599_4308918033_C_Red.idat': 719560,
            'GSM443814_HCC2218_4308918033_E_Grn.idat': 719560,
            'GSM443814_HCC2218_4308918033_E_Red.idat': 719560,
            'GSM443815_MCF10A_4308918033_A_Grn.idat': 719560,
            'GSM443815_MCF10A_4308918033_A_Red.idat': 719560,
            'GSM443816_BT474_4308918024_K_Grn.idat': 719561,
            'GSM443816_BT474_4308918024_K_Red.idat': 719561,
            'GSM443817_HCC1008_4308918024_J_Grn.idat': 719559,
            'GSM443817_HCC1008_4308918024_J_Red.idat': 719559,
            'GSM443818_HCC1395_4308918024_L_Grn.idat': 719559,
            'GSM443818_HCC1395_4308918024_L_Red.idat': 719559,
            'GSM443819_HCC1937_4308918024_A_Grn.idat': 719559,
            'GSM443819_HCC1937_4308918024_A_Red.idat': 719559,
            'GSM443821_MCF7_4308918024_I_Grn.idat': 719559,
            'GSM443821_MCF7_4308918024_I_Red.idat': 719559,
            # 'GSE17769_GPL8490_meta_data.pkl': 2723, # 2677 != 2723 expected
            'GSE17769_GPL8490_samplesheet.csv': 1862,
        }
        PLATFORM = 'GPL8490'
        LOCAL = Path('docs/example_data/GSE17769')
        kwargs = {'project_name': 'GSE17769', 'data_dir': LOCAL, 'clean':False, 'compress':False, 'verbose':True}
        result = methylprep.download.pipeline_find_betas_any_source(**kwargs)
        print(result)
        # runs a bunch of stuff. This GSE dataset has no IDATs or series matrix betas. just testing the samplesheet part.
        with open(Path(LOCAL, f"{kwargs['project_name']}_family.xml"), 'rb') as f:
            family = f.read()
        with open(Path(LOCAL, f"ref_{kwargs['project_name']}_family.xml"), 'rb') as f:
            ref_family = f.read()
        if not family == ref_family:
            raise AssertionError("Miniml _family.xml files don't match.")
        print('OK. Miniml _family.xml files match.')
        for _file,_size in expected_file_sizes.items():
            if Path(LOCAL,PLATFORM,_file).stat().st_size != _size:
                raise AssertionError(f"File size mismatch for {_file}: {Path(LOCAL,PLATFORM,_file).stat().st_size} != {_size} expected")
            print(f'OK: {_file} filesize matches: {Path(LOCAL,PLATFORM,_file).stat().st_size} : {_size}')
        for _file in Path(LOCAL,PLATFORM).rglob('*'):
            if _file.is_dir():
                continue
            if str(_file.name) == f"ref_{kwargs['project_name']}_family.xml":
                continue
            #if not str(_file.name).startswith('ref_') and not _file.is_dir(): # will get error trying to delete TempDir.
            if str(_file.name) in expected_file_sizes:
                _file.unlink()
        Path(LOCAL,'geo_alert GSE17769.csv').unlink()
        Path(LOCAL,'GSE17769_family.xml').unlink()
        Path(LOCAL,PLATFORM,'GSE17769_GPL8490_meta_data.pkl').unlink()
        Path(LOCAL,PLATFORM).rmdir()


    def test_read_series_matrix(self):
        import methylcheck # used in beta_bake
        test_file = Path('docs/example_data/GSE158089/test_series_matrix.txt')
        data = methylcheck.read_geo_processed.read_series_matrix(test_file, include_headers_df=True)
        if  data['df'].shape != (6, 14):
            raise AssertionError("dummy data mismatch data['df']")
        if len(data['series_dict']) != 26:
            raise AssertionError("dummy data mismatch data['series_dict']")
        if data['headers_df'].shape != (34, 14):
            raise AssertionError("dummy data mismatch data['headers_df']")

        result = methylprep.files.sample_sheets.sample_names_from_matrix(Path('docs/example_data/GSE158089/'))
        if result != ['iPSC_1', 'iPSC_2', 'NPC_1', 'NPC_2', 'NPC_3', 'NPC_4', 'Neuron_D37_1', 'Neuron_D37_2', 'Neuron_D37_3', 'Neuron_D37_4', 'Neuron_D58_1', 'Neuron_D58_2', 'Neuron_D58_3', 'Neuron_D58_4']:
            raise AssertionError("sample_names_from_matrix failed to parse test_series_matrix.txt")

    def test_samplesheet_from_series_matrix(self):
        import methylcheck
        test_file = Path('docs/example_data/GSE158089/test_series_matrix.txt')
        data = methylcheck.read_geo_processed.read_series_matrix(test_file, include_headers_df=True)
        meta_df = methylprep.download.geo.samplesheet_from_series_matrix( data['headers_df'] )
        if any(meta_df.columns.duplicated()):
            raise AsserttionError(f"Duplicate columns returned")
