# Lib
from logging import NullHandler, getLogger
# App
from .files import get_sample_sheet, get_sample_sheet_s3
from .processing import (
    get_manifest,
    get_raw_datasets,
    run_pipeline,
    consolidate_values_for_sheet,
    read_geo,
    )
from .download import run_series, run_series_list, convert_miniml, build_composite_dataset


getLogger(__name__).addHandler(NullHandler())


__all__ = [
    'get_manifest',
    'get_raw_datasets',
    'run_pipeline',
    'get_sample_sheet',
    'consolidate_values_for_sheet',
    'run_series',
    'run_series_list',
    'convert_miniml',
    'read_geo',
    'build_composite_dataset',
]

from pathlib import Path
with open(Path(Path(__file__).parents[1],'VERSION')) as version_file:
    __version__ = version_file.read().strip()
