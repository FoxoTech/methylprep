# Lib
import logging
from pathlib import Path, PurePath
import pandas as pd
# App
from ..models import Sample
from ..utils import get_file_object, reset_file


__all__ = ['SampleSheet', 'get_sample_sheet']


LOGGER = logging.getLogger(__name__)

REQUIRED_HEADERS = {'Sample_Name', 'Sentrix_ID', 'Sentrix_Position'}


def get_sample_sheet(dir_path, filepath=None):
    """Generates a SampleSheet instance for a given directory of processed data.

    Arguments:
        dir_path {string or path-like} -- Base directory of the sample sheet and associated IDAT files.

    Keyword Arguments:
        filepath {string or path-like} -- path of the sample sheet file if provided, otherwise
            one will try to be found. (default: {None})

    Returns:
        [SampleSheet] -- A SampleSheet instance.
    """
    LOGGER.info('Generating sample sheet')

    if not filepath:
        filepath = find_sample_sheet(dir_path)

    data_dir = PurePath(filepath).parent
    return SampleSheet(filepath, data_dir)


def find_sample_sheet(dir_path):
    """Find sample sheet file for Illumina methylation array

    Arguments:
        dir_path {string or path-like} -- Base directory of the sample sheet and associated IDAT files.

    Raises:
        FileNotFoundError: [description]
        FileNotFoundError: [description]
        Exception: [description]

    Returns:
        [string] -- Path to sample sheet in base directory
    """
    LOGGER.info('Searching for sample_sheet in %s', dir_path)

    sample_dir = Path(dir_path)

    if not sample_dir.is_dir():
        raise FileNotFoundError(f'{dir_path} is not a valid directory path')

    csv_files = sample_dir.glob('*.csv')
    candidates = [
        csv_file
        for csv_file in csv_files
        if SampleSheet.is_sample_sheet(csv_file)
    ]

    num_candidates = len(candidates)

    if num_candidates == 0:
        raise FileNotFoundError()

    if num_candidates > 1:
        raise Exception('Too many sample sheets were found in this directory')

    sample_sheet_file = candidates[0]
    LOGGER.info('Found sample sheet file: %s', sample_sheet_file)
    return sample_sheet_file


class SampleSheet():
    """Validates and parses an Illumina sample sheet file.

    Arguments:
        filepath_or_buffer {file-like} -- the sample sheet file to parse.
        dir_path {string or path-like} -- Base directory of the sample sheet and associated IDAT files.

    Raises:
        ValueError: The sample sheet is not formatted properly or a sample cannot be found.
    """

    __data_frame = None

    def __init__(self, filepath_or_buffer, data_dir):
        self.__samples = []

        self.data_dir = data_dir
        self.headers = []

        with get_file_object(filepath_or_buffer) as sample_sheet_file:
            self.read(sample_sheet_file)

    @staticmethod
    def is_sample_sheet(filepath_or_buffer):
        """Checks if the provided file-like object is a valid sample sheet.

        Arguments:
            filepath_or_buffer {file-like} -- the sample sheet file to parse.

        Returns:
            [boolean] -- Whether the file is a valid sample sheet.
        """
        data_frame = pd.read_csv(filepath_or_buffer, header=None, nrows=10)

        reset_file(filepath_or_buffer)

        for _, row in data_frame.iterrows():
            if REQUIRED_HEADERS.issubset(row.values):
                return True

        return False

    def get_samples(self):
        """Retrieves Sample objects from the processed sample sheet rows,
        building them if necessary."""
        if not self.__samples:
            self.build_samples()
        return self.__samples

    def get_sample(self, sample_name):
        candidates = [
            sample
            for sample in self.__samples
            if sample.name == sample_name
        ]

        num_candidates = len(candidates)
        if num_candidates != 1:
            raise ValueError(f'Expected sample with name {sample_name}. Found {num_candidates}')

        return candidates[0]

    def read(self, sample_sheet_file):
        """Validates and reads a sample sheet file, building a DataFrame from the parsed rows."""
        LOGGER.info('Parsing sample_sheet')

        if not self.is_sample_sheet(sample_sheet_file):
            columns = ', '.join(REQUIRED_HEADERS)
            raise ValueError(f'Cannot find header with values: {columns}')

        available_retries = 25
        cur_line = sample_sheet_file.readline()
        while not cur_line.startswith(b'[Data]'):
            if not available_retries:
                raise ValueError('Sample sheet is invalid. Could not find start of data row.')

            self.headers.append(cur_line.decode())
            cur_line = sample_sheet_file.readline()
            available_retries -= 1

        self.__data_frame = pd.read_csv(
            sample_sheet_file,
            keep_default_na=False,
            skip_blank_lines=True,
            dtype=str,
        )

        reset_file(sample_sheet_file)

    def build_samples(self):
        """Builds Sample objects from the processed sample sheet rows"""

        self.__samples = []

        for _index, row in self.__data_frame.iterrows():
            sentrix_id = row['Sentrix_ID'].strip()
            sentrix_position = row['Sentrix_Position'].strip()

            if not (sentrix_id and sentrix_position):
                continue

            sample = Sample(
                data_dir=self.data_dir,
                sentrix_id=sentrix_id,
                sentrix_position=sentrix_position,
                **row,
            )

            self.__samples.append(sample)
