# Lib
import logging
from pathlib import Path, PurePath
import pandas as pd
import re
# App
from ..models import Sample
from ..utils import get_file_object, reset_file


__all__ = ['SampleSheet', 'get_sample_sheet',  'get_sample_sheet_s3', 'find_sample_sheet', 'create_sample_sheet']


LOGGER = logging.getLogger(__name__)

REQUIRED_HEADERS = {'Sample_Name', 'Sentrix_ID', 'Sentrix_Position'}
ALT_REQUIRED_HEADERS = {'Sample_Name', 'SentrixBarcode_A', 'SentrixPosition_A'}


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
    LOGGER.debug('Reading sample sheet')

    if not filepath:
        filepath = find_sample_sheet(dir_path)

    data_dir = PurePath(filepath).parent
    return SampleSheet(filepath, data_dir)


def get_sample_sheet_s3(zip_reader):
    """ reads a zipfile and considers all filenames with 'sample_sheet' but will test all csv.
    the zip_reader is an amazon S3ZipReader object capable of reading the zipfile header."""
    ext_matched = [
        file_name
        for file_name in zip_reader.file_names
        if PurePath(file_name).suffix == '.csv'
    ]

    name_matched = [
        file_name
        for file_name in ext_matched
        if 'sample_sheet' in file_name.lower()
        or 'samplesheet' in file_name.lower()
    ]

    candidates = name_matched or ext_matched
    for file_name in candidates:
        sample_sheet_obj = zip_reader.get_file(file_name)
        if SampleSheet.is_sample_sheet(sample_sheet_obj):
            data_dir = PurePath(file_name).parent
            return SampleSheet(sample_sheet_obj, data_dir)
    raise FileNotFoundError('Could not find sample sheet in s3 file.')


def find_sample_sheet(dir_path, return_all=False):
    """Find sample sheet file for Illumina methylation array.

    Notes:
        looks for csv files in {dir_path}.
        If more than one csv file found, returns the one
        that has "sample_sheet" or 'samplesheet' in its name.
        Otherwise, raises error.

    Arguments:
        dir_path {string or path-like} -- Base directory of the sample sheet and associated IDAT files.
        return_all -- if True,
            returns a list of paths to samplesheets, if multiple present, instead of raising an error.

    Raises:
        FileNotFoundError: [description]
        Exception: [description]

    Returns:
        [string] -- Path to sample sheet in base directory
    """
    LOGGER.debug('Searching for sample_sheet in %s', dir_path)

    sample_dir = Path(dir_path)

    if not sample_dir.is_dir():
        raise FileNotFoundError(f'{dir_path} is not a valid directory path')

    csv_files = sample_dir.rglob('*.csv')
    candidates = [
        csv_file for csv_file in csv_files
        if SampleSheet.is_valid_csv(csv_file)
        and SampleSheet.is_sample_sheet(csv_file)
        and 'sample' in str(csv_file).lower()
        and 'sheet' in str(csv_file).lower()
    ]

    num_candidates = len(candidates)

    if num_candidates == 0:
        errors = [
            {'name': csv_file.name,
            'has_headers': SampleSheet.is_valid_csv(csv_file),
            'pandas_can_open': SampleSheet.is_sample_sheet(csv_file)}
            for csv_file in csv_files
        ]
        if errors == []:
            raise FileNotFoundError(f"Could not find sample sheet.")
        else:
            raise FileNotFoundError(f"Could not find sample sheet. (candidate files: {errors})")

    if num_candidates > 1:
        name_matched = [
            file_name
            for file_name in candidates
            if 'sample_sheet' in file_name.stem.lower()
            or 'samplesheet' in file_name.stem.lower()
        ]
        if len(name_matched) == 1:
            pass
        else:
            if return_all:
                return name_matched
            else:
                raise Exception(f"Too many sample sheets in this directory. Move or rename redundant ones. Or specify the path to the one to use with --sample_sheet. (candidate files: {candidates})")

    sample_sheet_file = candidates[0]
    LOGGER.debug('Found sample sheet file: %s', sample_sheet_file)
    return sample_sheet_file


def create_sample_sheet(dir_path, matrix_file=False, output_file='samplesheet.csv',
    sample_type='', sample_sub_type=''):
    """Creates a samplesheet.csv file from the .IDAT files of a GEO series directory

    Arguments:
        dir_path {string or path-like} -- Base directory of the sample sheet and associated IDAT files.
        matrix_file {boolean} -- Whether or not a Series Matrix File should be searched for names. (default: {False})

        ========== | ========= | ==== | =======
        parameter  | required | type | effect
        ========== | =========  ==== | =======
        sample_type | optional | string | label all samples in created sheet as this type (i.e. blood, saliva, tumor cells)
        sample_sub_type |  optional | string | further detail sample type for batch
        controls | optional | list of sample_names | assign all samples in controls list to be "control samples", not treatment samples.
        ========== | ========= | ==== | =======

    Note:
        Because sample_names are only generated from Matrix files, this method won't let you assign controls to samples from CLI.
        Would require all sample names be passed in from CLI as well, a pretty messy endeavor.

    Raises:
        FileNotFoundError: The directory could not be found.
    """

    sample_dir = Path(dir_path)

    if not sample_dir.is_dir():
        raise FileNotFoundError(f'{dir_path} is not a valid directory path')

    idat_files = sample_dir.rglob('*Grn.idat*') #.gz OK

    _dict = {'GSM_ID': [], 'Sample_Name': [], 'Sentrix_ID': [], 'Sentrix_Position': []}

    # additional optional columns
    addl_cols = []
    if sample_type:
        _dict['Sample_Type'] = []
        addl_cols.append('Sample_Type')
    if sample_sub_type:
        _dict['Sample_Sub_Type'] = []
        addl_cols.append('Sample_Sub_Type')

    file_name_error_msg = "This .idat file does not have the right pattern to auto-generate a sample sheet: {0}"
    for idat in idat_files:
        # split string by '/', last element is local file name
        try:
            filename = str(idat).split("/")[-1]
            split_filename = filename.split("_")

            if split_filename[0].startswith('GSM'):
                _dict['GSM_ID'].append(split_filename[0])
                _dict['Sentrix_ID'].append(split_filename[1])
                _dict['Sentrix_Position'].append(split_filename[2])
            elif len(split_filename) == 3:
                _dict['GSM_ID'].append("")
                _dict['Sentrix_ID'].append(split_filename[0])
                _dict['Sentrix_Position'].append(split_filename[1])
            else:
                raise ValueError(file_name_error_msg.format(idat))
        except:
            raise ValueError(file_name_error_msg.format(idat))

        if sample_type:
            _dict['Sample_Type'].append(sample_type)
        if sample_sub_type:
            _dict['Sample_Sub_Type'].append(sample_sub_type)

    if matrix_file:
        _dict['Sample_Name'] = sample_names_from_matrix(dir_path, _dict['GSM_ID'])
    else:
        # generate sample names
        for i in range (1, len(_dict['GSM_ID']) + 1):
            _dict['Sample_Name'].append("Sample_" + str(i))

    df = pd.DataFrame(data=_dict)
    df.to_csv(path_or_buf=(PurePath(dir_path, output_file)),index=False)

    LOGGER.info(f"[!] Created sample sheet: {dir_path}/samplesheet.csv with {len(_dict['GSM_ID'])} GSM_IDs")


def sample_names_from_matrix(dir_path, ordered_GSMs=None):
    """Extracts sample names from a GEO Series Matrix File and returns them in the order of the inputted GSM_IDs

    Arguments:
        dir_path {string or path-like} -- Base directory of the sample sheet and associated IDAT files.
        ordered_GSMs {list of strings} -- List of ordered GSM_IDs

    Raises:
        FileNotFoundError: The Series Matrix file could not be found

    Returns:
        [list of strings] -- Ordered Sample Names
    """

    sample_dir = Path(dir_path)
    matrix_files = list(sample_dir.glob('*matrix.txt'))
    if len(matrix_files) == 0:
        raise FileNotFoundError('No Series Matrix file found')

    f = open(matrix_files[0], "r") # loads the first matching one
    line = f.readline()

    sample_geo_accession = ''
    sample_title = ''
    while line:
        if "!Sample_title" in line:
            sample_title = line
            # print(line)
        if "!Sample_geo_accession" in line:
            sample_geo_accession = line
        if "!series_matrix_table_begin" in line:
            break
        line = f.readline()

    # in the matrix file, two consecutive lines contain quoted strings, separated by spaces with all the sample names and GSM IDs, respectively.
    unordered_Sample_Names = (re.findall(r'"(.*?)"', sample_title))
    unordered_GSMs = (re.findall(r'"(.*?)"', sample_geo_accession))
    GSM_to_name = dict(zip(unordered_GSMs, unordered_Sample_Names))
    if ordered_GSMs:
        ordered_Sample_Names = [GSM_to_name.get(GSM,'') for GSM in ordered_GSMs]
        return ordered_Sample_Names
    else:
        return unordered_Sample_Names


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
        self.fields = {}
        self.renamed_fields = {}

        self.data_dir = data_dir
        self.headers = []
        self.alt_headers = None

        with get_file_object(filepath_or_buffer) as sample_sheet_file:
            self.read(sample_sheet_file)

    @staticmethod
    def is_sample_sheet(filepath_or_buffer):
        """Checks if the provided file-like object is a valid sample sheet.

        Method:
            If any row in the file contains these column names, it passes: `{0}`
            Alternatively, if all of these column names are present instead, it also passes, and processing will expect these: `{1}`

        Arguments:
            filepath_or_buffer {{file-like}} -- the sample sheet file to parse.

        Returns:
            [boolean] -- Whether the file is a valid sample sheet.
        """.format(REQUIRED_HEADERS, ALT_REQUIRED_HEADERS)
        data_frame = pd.read_csv(filepath_or_buffer, header=None, nrows=25)

        reset_file(filepath_or_buffer)

        for _, row in data_frame.iterrows():
            if REQUIRED_HEADERS.issubset(row.values):
                return True
            elif ALT_REQUIRED_HEADERS.issubset(row.values):
                return True

        return False

    @staticmethod
    def is_valid_csv(filepath_or_buffer):
        try:
            data_frame = pd.read_csv(filepath_or_buffer, header=None, nrows=25)
            return True
        except Exception:
            return False

    def get_samples(self):
        """Retrieves Sample objects from the processed sample sheet rows,
        building them if necessary."""
        if not self.__samples:
            self.build_samples()
        return self.__samples

    def get_sample(self, sample_name):
        """ scans all samples for one matching sample_name, if provided.
        If no sample_name, then it returns all samples."""
        # this isn't automatically done, but needed here to work.
        null = self.get_samples()

        candidates = [
            sample
            for sample in self.__samples
            if sample.name == sample_name
        ]
        # or    sample.GSM_ID == sample_name or
        # sample.Sample_Name == sample_name

        num_candidates = len(candidates)
        if num_candidates != 1:
            raise ValueError(f'Expected sample with name `{sample_name}`. Found {num_candidates}')

        return candidates[0]

    def build_samples(self):
        """Builds Sample objects from the processed sample sheet rows.

        Added to Sample as class_method: if the idat file is not in the same folder, (check if exists!) looks recursively for that filename and updates the data_dir for that Sample.
        """

        self.__samples = []

        #logging.info('Building samples')

        for _index, row in self.__data_frame.iterrows():
            if self.alt_headers:
                sentrix_id = row['SentrixBarcode_A'].strip()
                sentrix_position = row['SentrixPosition_A'].strip()
            else:
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
            if sample.renamed_fields != {}:
                self.renamed_fields.update(sample.renamed_fields)
            self.fields.update(sample.fields)
            self.__samples.append(sample)

    def contains_column(self, column_name):
        """ helper function to determine if sample_sheet contains a specific column, such as GSM_ID.
        SampleSheet must already have __data_frame in it."""
        if column_name in self.__data_frame:
            return True
        return False

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

        LOGGER.debug('Parsing sample_sheet')

        if not self.is_sample_sheet(sample_sheet_file):
            columns = ', '.join(REQUIRED_HEADERS)
            alt_columns = ', '.join(ALT_REQUIRED_HEADERS)
            raise ValueError(f'Cannot find header with values: {columns} or {alt_columns}')

        # first, parse headers and reset
        # this puts all the sample_sheet header rows into SampleSheet.headers list.
        rows_to_scan=100
        cur_line = sample_sheet_file.readline()
        while not cur_line.startswith(b'[Data]'):
            if rows_to_scan == 0:
                if self.headers == {}:
                    LOGGING.info("Finished scanning sample_sheet; did not find header info.")
                break
            raw_line = cur_line.decode()
            if raw_line:
                self.headers.append(raw_line)
            cur_line = sample_sheet_file.readline()
            rows_to_scan -= 1
        reset_file(sample_sheet_file)

        test_sheet = pd.read_csv(
            sample_sheet_file,
            header = None,  # this ensures row[0] included as data -- [this is for looking for the header]
            keep_default_na=False,
            skip_blank_lines=True,
            dtype=str,
        )
        test_sheet = test_sheet.to_dict('records')  # list of dicts
        rows_to_scan = 25 # scan first 25 rows of document
        start_row = None
        for idx,row in enumerate(test_sheet):  # header is not the first row. alt format is that header begins on row after [Data]
            if rows_to_scan == 0:
                LOGGER.info(f'DEBUG {cur_line} {line_bits}')
                raise ValueError('Sample sheet is invalid. Could not find start of data row, assuming there should be a [Data] row to start data, and no more than 25 preceding rows.')
            if '[Data]' in row.values():
                # Format 1 parsing: assume the header begins right after [Data]
                start_row = idx + 1
                break
            if REQUIRED_HEADERS.issubset(row.values()):
                # Format 2 parsing: no [Data] and probably first row is header.
                start_row = idx
                self.alt_headers = False
                break
            if ALT_REQUIRED_HEADERS.issubset(row.values()):
                # Format 2 parsing: no [Data] and probably first row is header.
                start_row = idx
                self.alt_headers = True
                break
            rows_to_scan -= 1
        if start_row == None:
            raise ValueError("error - did not parse header right")

        # preceding code uses `start_row` to strip out any non-data rows from sample_sheet_file before loading into dataframe.
        reset_file(sample_sheet_file)
        self.__data_frame = pd.read_csv(
            sample_sheet_file,
            header=start_row,
            keep_default_na=False,
            skip_blank_lines=True,
            dtype=str,
        )
        reset_file(sample_sheet_file)

        # rename ALT columns to standard columns in the sample_sheet dataframe now.
        if self.alt_headers:
            self.rename_alt_headers()

    def rename_alt_headers(self):
        columns = {'SentrixBarcode_A':'Sentrix_ID','SentrixPosition_A':'Sentrix_Position'}
        self.__data_frame = self.__data_frame.rename(columns=columns)
        LOGGER.info(f"Renamed SampleSheet columns {columns}")

    def build_meta_data(self, samples = None):
        """Takes a list of samples and returns a data_frame that can be saved as a pickle. """
        if samples:
            pass
        elif not samples and hasattr(self, '__samples'):
            samples = getattr(self, '__samples')
        else:
            raise ValueError("Either provide a list of samples or run SampleSheet.get_samples() first.")
        field_classattr_lookup = {
            'Sentrix_ID': 'sentrix_id',
            'Sentrix_Position': 'sentrix_position',
            'Sample_Group': 'group',
            'Sample_Name': 'name',
            'Sample_Plate': 'plate',
            'Pool_ID': 'pool',
            'Sample_Well': 'well',
            'GSM_ID': 'GSM_ID',
            'Sample_Type': 'type',
            'Sub_Type': 'sub_type',
            'Control': 'is_control',
        }
        # sample_sheet.fields is a complete mapping of original and renamed_fields
        cols = list(self.fields.values()) + ['Sample_ID']
        meta_frame = pd.DataFrame(columns=cols)
        # row contains the renamed fields, and pulls in the original data from sample_sheet
        for sample in samples:
            row = {}
            for field in self.fields.keys():
                if self.fields[field] in field_classattr_lookup:
                    row[ self.fields[field] ] = getattr(sample, field_classattr_lookup[self.fields[field]] )
                elif field in self.renamed_fields:
                    row[ self.fields[field] ] = getattr(sample, self.renamed_fields[field])
                else:
                    LOGGER.info(f"extra column: {field} ignored")
            # add the UID that matches m_value/beta value pickles
            #... unless there's a GSM_ID too
            row['Sample_ID'] = f"{row['Sentrix_ID']}_{row['Sentrix_Position']}"
            meta_frame = meta_frame.append(row, ignore_index=True)
        return meta_frame
