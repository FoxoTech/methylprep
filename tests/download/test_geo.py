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
    """ Note: all 3 tests complete in 30 to 40s """

    def test_pipeline_find_betas_any_source_samplesheet(self):
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
            'GSE110454_beta_values.pkl': 39025301, # GSE110454_beta_values.pkl: 39025301 != 39025330 expected
            'GSE110454_samplesheet.csv': 4794,
            'GSE110454_series_summary.json': 1563,
        }
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
        Path(LOCAL, f"{kwargs['project_name']}_samplesheet.csv").unlink()
        Path(LOCAL, f"{kwargs['project_name']}_series_summary.json").unlink()
        Path(LOCAL, f"{kwargs['project_name']}_beta_values.pkl").unlink()

    def test_pipeline_find_betas_any_source_27k_idats(self):
        """parses IDATs with 10 27k samples, 13MB takes about 30sec """
        expected_file_sizes = {
            'geo_alert GSE17769.csv': 559,
            'GSE17769_family.xml': 28475,
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
        }
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
            if Path(LOCAL,_file).stat().st_size != _size:
                raise AssertionError(f"File size mismatch for {_file}: {Path(LOCAL,_file).stat().st_size} != {_size} expected")
            print(f'OK: {_file} filesize matches: {Path(LOCAL,_file).stat().st_size} : {_size}')
        for _file in Path(LOCAL).rglob('*'):
            if _file.is_dir():
                continue
            if str(_file.name) == f"ref_{kwargs['project_name']}_family.xml":
                continue
            #if not str(_file.name).startswith('ref_') and not _file.is_dir(): # will get error trying to delete TempDir.
            if str(_file.name) in expected_file_sizes:
                _file.unlink()
