# Lib
import logging
import math
import numpy as np
import pandas as pd
import scipy
from pathlib import Path
import shutil
from collections import Counter
# App
from ..models.sigset import parse_sample_sheet_into_idat_datasets
from ..files import find_sample_sheet, create_sample_sheet, SampleSheet

LOGGER = logging.getLogger(__name__)

def check_array_folders(data_dir, verbose=True):
    """Confirms a folder's idats are separated into sub-folders by array_type, but does not move them for you.
    Should deal with GEO multi-array packages."""
    instructions = []
    try:
        sample_sheet_file_path = find_sample_sheet(data_dir, return_all=True)
    except FileNotFoundError as e:
        import pdb;pdb.set_trace()
        create_sample_sheet(data_dir, matrix_file=False, output_file='samplesheet.csv',
            sample_type='', sample_sub_type='')
        sample_sheet_file_path = f"{data_dir}/samplesheet.csv"
        if verbose: LOGGER.info(f"Created {data_dir}/samplesheet.csv")
    # where multiple samplesheets found
    if isinstance(sample_sheet_file_path, list):
        sorted_files = {}
        # check if they are in GEO GPLxxx folders
        folders = ['GPL8490', 'GPL13534', 'GPL21145', 'epic_plus', 'mouse']
        for sample_sheet_file in sample_sheet_file_path:
            if any([array_folder in sample_sheet_file.parts for array_folder in folders]):
                # now verify idats are present. later: do the idats match the samplesheet?
                array_type = [array_folder for array_folder in folders if array_folder in sample_sheet_file.parts][0]
                idats_found = list(Path(sample_sheet_file.parent).rglob('*.idat')) + list(Path(sample_sheet_file.parent).rglob('*.idat.gz'))
                instructions.append(f"For {int(len(idats_found)/2)} {array_type} samples run: `methylprep process -d {sample_sheet_file.parent} --all`")
    return instructions

'''
def __draft_separate_array_files(data_dir, verbose=True):
    """Provide a folder location of idats comprised of more than one array type and this will
    separate the files into separate sub-folders and split the samplesheet (if found) accordingly.
    Then, you can run `methylprep process` on each of these subfolders in tandem it will work."""
    # 0 - find sample_sheet meta data
    # FOR GEO, look for the richer pickle file to read instead!!!
    try:
        sample_sheet_file_path = find_sample_sheet(data_dir)
    except FileNotFoundError as e:
        LOGGER.info(e)
        create_sample_sheet(dir_path, matrix_file=False, output_file='samplesheet.csv',
            sample_type='', sample_sub_type='')
        sample_sheet_file_path = f"{data_dir}/samplesheet.csv"
        if verbose: LOGGER.info(f"Created {data_dir}/samplesheet.csv")
    except Exception as e: # >1 samplesheets found
        # In this messy case, track how many subfolders are present and abort if more than 2,
        # because this is likely being run in a top level folder that SHOULD NOT BE sorted.
        if verbose: LOGGER.info(e)
        idat_files = sample_dir.rglob('*Grn.idat')
        idat_paths_found = Counter()
        for idat_file in idat_files:
            if idat_file.is_file():
                idat_paths_found[idat_file.parent] += 1
        if len(idat_paths_found) > 1:
            raise ValueError(f"Found idat files in {len(idat_paths_found)} different paths, so you must move files manually: {idat_paths_found}")
        create_sample_sheet(dir_path, matrix_file=False, output_file='multi_samplesheet.csv',
            sample_type='', sample_sub_type='')
        sample_sheet_file_path = f"{data_dir}/multi_samplesheet.csv"
        if verbose: LOGGER.info(f"Created {data_dir}/multi_samplesheet.csv")
    sample_sheet = SampleSheet(sample_sheet_file_path, data_dir)
    if verbose: LOGGER.info("Found {sample_sheet_file_path}")
    # 1 - find out which idat files go with which array types
    idat_datasets = parse_sample_sheet_into_idat_datasets(sample_sheet, sample_name=None, from_s3=None, meta_only=False)
    # 2 - split these | green_idat, red_idat, sample, array_type
    ['GPL21145', 'GPL13534', 'GPL23976', 'GPL8490', 'GPL16304', 'GPL18809'] # GPL16304, GPL28271 (Horvath)
    array_folder_names = {
    '27k': 'GPL8490',
    '450k': 'GPL13534',
    'epic': 'GPL21145',
    'epic+': 'epic_plus',
    'mouse': 'mouse',
    }
    instructions = []
    array_types = {idat_dataset['array_type'] for idat_dataset in idat_datasets}
    sorted_files = {}
    for array_type in array_types:
        files_moved = []
        new = Path(data_dir, array_folder_names[array_type])
        if not new.exists():
            new.mkdir()
        else:
            raise OSError("folder {new} already exists, so we can't copy idats into it.")
        # .gz also okay
        for idat in list(Path(data_dir).glob('*.idat')) + list(Path(data_dir).glob('*.idat.gz')):
            filename = shutil.move(idat, Path(new))
            files_moved.append(filename)
        sorted_files[array_type] = files_moved
        if verbose: LOGGER.info(f"Moved {len(files_moved)} idats to {new}")
        instructions.append("For {array_type} ({len(files_moved)/2} samples), run `methylprep process -d {new} --all`")
    for line in instructions:
        print(line)
'''
