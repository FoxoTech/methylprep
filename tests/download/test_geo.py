import pytest
import sys
import pandas as pd
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

LOCAL = Path('docs/example_data/GSE132203')

class TestGeo():

    def test_pipeline_find_betas_any_source(self):
        kwargs = {'project_name': 'GSE132203', 'data_dir': 'docs/example_data/GSE132203', 'clean':False, 'compress':False, 'verbose':True}
        result = methylprep.download.pipeline_find_betas_any_source(**kwargs)
        print(result)
        # runs a bunch of stuff. This GSE dataset has no IDATs or series matrix betas. just testing the samplesheet part.
        samplesheet = pd.read_csv(Path(LOCAL, f"{kwargs['project_name']}_samplesheet.csv"))
        ref_samplesheet = pd.read_csv(Path(LOCAL, f"ref_{kwargs['project_name']}_samplesheet.csv"))
        if not samplesheet.equals(ref_samplesheet):
            raise AssertionError("Samplesheets don't match.")
        print('OK. samplesheets match.')
        summary = pd.read_csv(Path(LOCAL, f"{kwargs['project_name']}_series_summary.json"))
        ref_summary = pd.read_csv(Path(LOCAL, f"ref_{kwargs['project_name']}_series_summary.json"))
        if not summary.equals(ref_summary):
            raise AssertionError("Series summaries don't match.")
        print('OK. summaries match.')
        Path(LOCAL, f"{kwargs['project_name']}_samplesheet.csv").unlink()
        Path(LOCAL, f"{kwargs['project_name']}_series_summary.json").unlink()
        Path(LOCAL, f"{kwargs['project_name']}_beta_values.pkl").unlink()
