# Lib
from logging import NullHandler, getLogger
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
# App
from .files import get_sample_sheet, get_sample_sheet_s3
from .processing import (
    run_pipeline,
    make_pipeline,
    consolidate_values_for_sheet
    )
from .download import run_series, run_series_list, convert_miniml, build_composite_dataset
from .models import ArrayType, parse_sample_sheet_into_idat_datasets
from .files import Manifest
from .version import __version__

getLogger(__name__).addHandler(NullHandler())

#import numpy as np
#np.seterr(all='raise') -- for debugging overflow / underflow somewhere

__all__ = [
    'ArrayType',
    'Manifest',
    #'get_manifest',
    #'get_raw_datasets',
    'get_sample_sheet',
    'parse_sample_sheet_into_idat_datasets',
    'consolidate_values_for_sheet',
    'run_series',
    'run_series_list',
    'convert_miniml',
    'build_composite_dataset',
    'run_pipeline',
    'make_pipeline',
]
