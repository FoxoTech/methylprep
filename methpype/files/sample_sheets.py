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

        Method:
            If any row in the file contains these column names, it passes: `{0}`

        Arguments:
            filepath_or_buffer {{file-like}} -- the sample sheet file to parse.

        Returns:
            [boolean] -- Whether the file is a valid sample sheet.
        """.format(REQUIRED_HEADERS)
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
        """Validates and reads a sample sheet file, building a DataFrame from the parsed rows.

        Method:
            It autodetects whether a sample sheet is formatted in Infinium MethylationEPIC style, or without the headers.
            Rows must contain these columns: {0}
            See https://support.illumina.com/downloads/infinium-methylationepic-sample-sheet.html for more information about file formatting.

            Format 1: First row of file contains header data.
            Format 2: header is not the first row. Header begins on the row after [Data] appears in first column.

        Dev notes:
            It loads whole file using pandas.read_csv to better handle whitespace/matching on headers.""".format(REQUIRED_HEADERS)

        LOGGER.info('Parsing sample_sheet')

        if not self.is_sample_sheet(sample_sheet_file):
            columns = ', '.join(REQUIRED_HEADERS)
            raise ValueError(f'Cannot find header with values: {columns}')

        test_sheet = pd.read_csv(
            sample_sheet_file,
            header = None,  # so that is includes row[0] as data -- [this is for looking for the header]
            keep_default_na=False,
            skip_blank_lines=True,
            dtype=str,
        )
        test_sheet = test_sheet.to_dict('records')  # list of dicts
        available_retries = 25
        start_row = None
        if REQUIRED_HEADERS.issubset(set(test_sheet[0].values())):
            start_row = 0
        else:
            for idx,row in enumerate(test_sheet):  # header is not the first row. alt format is that header begins on row after [Data]
                if not available_retries:
                    print(f'DEBUG {cur_line} {line_bits}')
                    raise ValueError('Sample sheet is invalid. Could not find start of data row, assuming there should be a [Data] row to start data, and no more than 25 preceding rows.')
                if '[Data]' in row.values():
                    start_row = idx + 1  # the header begins right after [Data]
                    break
                available_retries -= 1
        if start_row == None:
            raise ValueError("error - did not parse header right")

        # preceding code strips out any non-data rows from sample_sheet_file before loading into dataframe.
        reset_file(sample_sheet_file)
        self.__data_frame = pd.read_csv(
            sample_sheet_file,
            header=start_row,
            keep_default_na=False,
            skip_blank_lines=True,
            dtype=str,
        )

        reset_file(sample_sheet_file)

    def build_samples(self):
        """Builds Sample objects from the processed sample sheet rows.

        Added to Sample as class_method: if the idat file is not in the same folder, (check if exists!) looks recursively for that filename and updates the data_dir for that Sample.
        """

        self.__samples = []

        logging.info('Building samples')

        for _index, row in self.__data_frame.iterrows():
            sentrix_id = row['Sentrix_ID'].strip()
            sentrix_position = row['Sentrix_Position'].strip()

            if not (sentrix_id and sentrix_position):
                continue

            sample = Sample(
                data_dir=self.data_dir,  # this assumes the .idat files are in the same folder with the samplesheet.
                sentrix_id=sentrix_id,
                sentrix_position=sentrix_position,
                **row,
            )

            self.__samples.append(sample)
