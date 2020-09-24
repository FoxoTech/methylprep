# Lib
from enum import IntEnum, unique
import pandas as pd
# App
from ..utils import (
    get_file_object,
    read_and_reset,
    read_byte,
    read_char,
    read_int,
    read_long,
    read_results,
    read_short,
    read_string,
    npread,
)


__all__ = ['IdatDataset']


# Constants
# ----------------------------------------------------------------------------

DEFAULT_IDAT_VERSION = 3
DEFAULT_IDAT_FILE_ID = 'IDAT'


@unique
class IdatHeaderLocation(IntEnum):
    """Unique IntEnum defining constant values for byte offsets of IDAT headers.
    """
    FILE_TYPE = 0
    VERSION = 4
    FIELD_COUNT = 12
    SECTION_OFFSETS = 16


@unique
class IdatSectionCode(IntEnum):
    """Unique IntEnum defining constant values for byte offsets of IDAT headers.
    These values come from the field codes of the Bioconductor illuminaio package.

    MM: refer to https://bioconductor.org/packages/release/bioc/vignettes/illuminaio/inst/doc/EncryptedFormat.pdf
    and https://bioconductor.org/packages/release/bioc/vignettes/illuminaio/inst/doc/illuminaio.pdf
    """
    ILLUMINA_ID = 102
    STD_DEV = 103
    MEAN = 104
    NUM_BEADS = 107
    MID_BLOCK = 200
    RUN_INFO = 300
    RED_GREEN = 400
    MOSTLY_NULL = 401
    BARCODE = 402
    CHIP_TYPE = 403
    MOSTLY_A = 404
    UNKNOWN_1 = 405
    UNKNOWN_2 = 406
    UNKNOWN_3 = 407
    UNKNOWN_4 = 408
    UNKNOWN_5 = 409
    UNKNOWN_6 = 410
    UNKNOWN_7 = 510
    NUM_SNPS_READ = 1000


"""Constant set of 'Unknown' IDAT sections."""
UNKNOWN_IDAT_SECTIONS = (
    IdatSectionCode.MOSTLY_NULL,
    IdatSectionCode.MOSTLY_A,
    IdatSectionCode.UNKNOWN_1,
    IdatSectionCode.UNKNOWN_2,
    IdatSectionCode.UNKNOWN_3,
    IdatSectionCode.UNKNOWN_4,
    IdatSectionCode.UNKNOWN_5,
    IdatSectionCode.UNKNOWN_6,
    IdatSectionCode.UNKNOWN_7,
)


# Object Definitions
# ----------------------------------------------------------------------------

class RunInfo():
    """A dataclass defining IDAT run information.

    Arguments:
        idat_file {file-like} -- An open IDAT file.
    """

    __slots__ = [
        'run_time',
        'block_type',
        'block_pars',
        'block_code',
        'code_version',
    ]

    def __init__(self, idat_file):
        """Initializes the RunInfo data class by reading the provided idat_file."""
        self.run_time = read_string(idat_file)
        self.block_type = read_string(idat_file)
        self.block_pars = read_string(idat_file)
        self.block_code = read_string(idat_file)
        self.code_version = read_string(idat_file)

class IdatDataset():
    """Validates and parses an Illumina IDAT file.

    Arguments:
        filepath_or_buffer {file-like} -- the IDAT file to parse.
        channel {Channel} -- the fluorescent channel (Channel.RED or Channel.GREEN)
        that produced the IDAT dataset.

    Keyword Arguments:
        idat_id {string} -- expected IDAT file identifier (default: {DEFAULT_IDAT_FILE_ID})
        idat_version {integer} -- expected IDAT version (default: {DEFAULT_IDAT_VERSION})

    Raises:
        ValueError: The IDAT file has an incorrect identifier or version specifier.
    """
    def __init__(
        self,
        filepath_or_buffer,
        channel,
        idat_id=DEFAULT_IDAT_FILE_ID,
        idat_version=DEFAULT_IDAT_VERSION,
    ):
        """Initializes the IdatDataset, reads and parses the IDAT file."""
        self.channel = channel
        self.barcode = None
        self.chip_type = None
        self.n_beads = 0
        self.n_snps_read = 0
        self.run_info = []

        with get_file_object(filepath_or_buffer) as idat_file:
            # assert file is indeed IDAT format
            if not self.is_idat_file(idat_file, idat_id):
                raise ValueError('Not an IDAT file. Unsupported file type.')

            # assert correct IDAT file version
            if not self.is_correct_version(idat_file, idat_version):
                raise ValueError('Not a version 3 IDAT file. Unsupported IDAT version.')

            self.probe_means = self.read(idat_file)

    @staticmethod
    @read_and_reset
    def is_idat_file(idat_file, expected):
        """Checks if the provided file has the correct identifier.

        Arguments:
            idat_file {file-like} -- the IDAT file to check.
            expected {string} -- expected IDAT file identifier.

        Returns:
            [boolean] -- If the IDAT file identifier matches the expected value
        """
        idat_file.seek(IdatHeaderLocation.FILE_TYPE.value)
        file_type = read_char(idat_file, len(expected))
        return file_type.lower() == expected.lower()

    @staticmethod
    @read_and_reset
    def is_correct_version(idat_file, expected):
        """Checks if the provided file has the correct version.

        Arguments:
            idat_file {file-like} -- the IDAT file to check.
            expected {integer} -- expected IDAT version.

        Returns:
            [boolean] -- If the IDAT file version matches the expected value
        """
        idat_file.seek(IdatHeaderLocation.VERSION.value)
        idat_version = read_long(idat_file)
        return str(idat_version) == str(expected)

    @staticmethod
    @read_and_reset
    def get_section_offsets(idat_file):
        """Parses the IDAT file header to get the byte position
        for the start of each section.

        Arguments:
            idat_file {file-like} -- the IDAT file to process.

        Returns:
            [dict] -- The byte offset for each file section.
        """
        idat_file.seek(IdatHeaderLocation.FIELD_COUNT.value)
        num_fields = read_int(idat_file)

        idat_file.seek(IdatHeaderLocation.SECTION_OFFSETS.value)

        offsets = {}
        for _idx in range(num_fields):
            key = read_short(idat_file)
            offsets[key] = read_long(idat_file)

        return offsets

    def read(self, idat_file):
        """Reads the IDAT file and parses the appropriate sections. Joins the
        mean probe intensity values with their Illumina probe ID.

        Arguments:
            idat_file {file-like} -- the IDAT file to process.

        Returns:
            DataFrame -- mean probe intensity values indexed by Illumina ID.
        """
        section_offsets = self.get_section_offsets(idat_file)

        def seek_to_section(section_code):
            offset = section_offsets[section_code.value]
            idat_file.seek(offset)

        seek_to_section(IdatSectionCode.BARCODE)
        self.barcode = read_string(idat_file)

        seek_to_section(IdatSectionCode.CHIP_TYPE)
        self.chip_type = read_string(idat_file)

        seek_to_section(IdatSectionCode.NUM_SNPS_READ)
        self.n_snps_read = read_int(idat_file)

        seek_to_section(IdatSectionCode.NUM_BEADS)
        self.n_beads = npread(idat_file, '<u1', self.n_snps_read)

        seek_to_section(IdatSectionCode.ILLUMINA_ID)
        illumina_ids = npread(idat_file, '<i4', self.n_snps_read)

        seek_to_section(IdatSectionCode.MEAN)
        probe_means = npread(idat_file, '<u2', self.n_snps_read)

        probe_records = dict(zip(illumina_ids, probe_means))

        data_frame = pd.DataFrame.from_dict(
            data=probe_records,
            orient='index',
            columns=['mean_value'],
            dtype='float32', # int16 could work, and reduce memory by 1/2, but some raw values were > 32127 -- without prenormalization, you get negative values back, which breaks stuff.
        )

        data_frame.index.name = 'illumina_id'

        return data_frame
