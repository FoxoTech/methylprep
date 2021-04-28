# Lib
from logging import NullHandler, getLogger
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
# App
from .files import get_sample_sheet, get_sample_sheet_s3
from .processing import (
    get_manifest,
    run_pipeline,
    consolidate_values_for_sheet,
    read_geo,
    detect_header_pattern,
    )
from .download import run_series, run_series_list, convert_miniml, build_composite_dataset
from .models import ArrayType, get_raw_datasets
from .files import Manifest
from .version import __version__

getLogger(__name__).addHandler(NullHandler())

#import numpy as np
#np.seterr(all='raise') -- for debugging overflow / underflow somewhere

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
    'detect_header_pattern',
    'build_composite_dataset',
    'Manifest',
    'ArrayType',
    'make_pipeline',
]
