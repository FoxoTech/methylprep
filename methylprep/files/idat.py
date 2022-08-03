# Lib
from enum import IntEnum, unique
import pandas as pd
import numpy as np
import struct
from pprint import pprint
import logging

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

LOGGER = logging.getLogger(__name__)
LOGGER.setLevel( logging.WARNING )

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
    and source: https://github.com/snewhouse/glu-genetics/blob/master/glu/lib/illumina.py
    more on encrypted IDAT format here: https://www.bioconductor.org/packages/release/bioc/vignettes/illuminaio/inst/doc/EncryptedFormat.pdf
    """
    ILLUMINA_ID = 102
    STD_DEV = 103
    MEAN = 104
    NUM_BEADS = 107 # how many replicate measurements for each probe
    MID_BLOCK = 200
    RUN_INFO = 300
    RED_GREEN = 400
    MOSTLY_NULL = 401 # manifest
    BARCODE = 402
    CHIP_TYPE = 403 # format
    MOSTLY_A = 404  # label
    UNKNOWN_1 = 405 # opa
    UNKNOWN_2 = 406 # sampleid
    UNKNOWN_3 = 407 # descr
    UNKNOWN_4 = 408 # plate
    UNKNOWN_5 = 409 # well
    UNKNOWN_6 = 410
    UNKNOWN_7 = 510 # unknown
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
""" DEPRECATED: use IdatDataset(... verbose=True) instead.
class RunInfo():
    '''A dataclass defining IDAT run information.

    Arguments:
        idat_file {file-like} -- An open IDAT file.'''

    __slots__ = [
        'run_time',
        'block_type',
        'block_pars',
        'block_code',
        'code_version',
    ]

    def __init__(self, idat_file):
        # Initializes the RunInfo data class by reading the provided idat_file.
        self.run_time = read_string(idat_file)
        self.block_type = read_string(idat_file)
        self.block_pars = read_string(idat_file)
        self.block_code = read_string(idat_file)
        self.code_version = read_string(idat_file)
"""

class IdatDataset():
    """Validates and parses an Illumina IDAT file.

    Arguments:
        filepath_or_buffer {file-like} -- the IDAT file to parse.
        channel {Channel} -- the fluorescent channel (Channel.RED or Channel.GREEN)
        that produced the IDAT dataset.

    Keyword Arguments:
        idat_id {string} -- expected IDAT file identifier (default: {DEFAULT_IDAT_FILE_ID})
        idat_version {integer} -- expected IDAT version (default: {DEFAULT_IDAT_VERSION})
        bit {string, default 'float32'} -- 'float16' will pre-normalize intensities,
            capping max intensity at 32127. This cuts data size in half, but will reduce
            precision on ~0.01% of probes. [effectively downscaling fluorescence]
    Raises:
        ValueError: The IDAT file has an incorrect identifier or version specifier.
    """
    def __init__(
        self,
        filepath_or_buffer,
        channel,
        idat_id=DEFAULT_IDAT_FILE_ID,
        idat_version=DEFAULT_IDAT_VERSION,
        verbose=False,
        std_dev=False,
        nbeads=False,
        bit='float32',
    ):
        """Initializes the IdatDataset, reads and parses the IDAT file."""
        self.verbose = verbose
        self.channel = channel
        self.barcode = None
        self.chip_type = None
        self.n_beads = 0
        self.n_snps_read = 0
        self.run_info = []
        self.include_std_dev = std_dev
        self.include_n_beads = nbeads
        self.bit = bit

        with get_file_object(filepath_or_buffer) as idat_file:
            # assert file is indeed IDAT format
            if not self.is_idat_file(idat_file, idat_id):
                raise ValueError('Not an IDAT file. Unsupported file type.')

            # assert correct IDAT file version
            if not self.is_correct_version(idat_file, idat_version):
                raise ValueError('Not a version 3 IDAT file. Unsupported IDAT version.')

            self.probe_means = self.read(idat_file)
            if self.overflow_check() is False:
                LOGGER.warning("IDAT: contains negative probe values (uint16 overflow error)")
            if self.verbose:
                self.meta(idat_file)

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
        self.n_beads = npread(idat_file, '<u1', self.n_snps_read) # was <u1

        seek_to_section(IdatSectionCode.ILLUMINA_ID)
        illumina_ids = npread(idat_file, '<i4', self.n_snps_read)

        seek_to_section(IdatSectionCode.MEAN)
        probe_means = npread(idat_file, '<u2', self.n_snps_read) # '<u2' reads data as numpy unsigned-float16

        seek_to_section(IdatSectionCode.RUN_INFO)
        runinfo_entry_count, = struct.unpack('<L', idat_file.read(4))
        for i in range(runinfo_entry_count):
            timestamp    = read_string(idat_file)
            entry_type   = read_string(idat_file)
            parameters   = read_string(idat_file)
            codeblock    = read_string(idat_file)
            code_version = read_string(idat_file)
            self.run_info.append( (timestamp, entry_type, parameters, codeblock, code_version) )

        if self.include_std_dev and self.include_n_beads:
            seek_to_section(IdatSectionCode.STD_DEV)
            std_devs = npread(idat_file, '<u2', self.n_snps_read)
            # print(len(probe_means), len(std_devs), len(self.n_beads))
            data_frame = pd.DataFrame(
                data={'mean_value':probe_means, 'std_dev':std_devs, 'n_beads':self.n_beads},
                index=illumina_ids,
                columns=['mean_value','std_dev','n_beads'],
                dtype=self.bit, # int16 could work, and reduce memory by 1/2, but some raw values were > 32127 -- without prenormalization, you get negative values back, which breaks stuff.
            )

        elif self.include_std_dev:
            seek_to_section(IdatSectionCode.STD_DEV)
            std_devs = npread(idat_file, '<u2', self.n_snps_read)
            # print(len(probe_means), len(std_devs))
            data_frame = pd.DataFrame(
                data={'mean_value':probe_means, 'std_dev':std_devs},
                index=illumina_ids,
                columns=['mean_value','std_dev'],
                dtype=self.bit, # int16 could work, and reduce memory by 1/2, but some raw values were > 32127 -- without prenormalization, you get negative values back, which breaks stuff.
            )
        elif self.include_n_beads:
            data_frame = pd.DataFrame(
                data={'mean_value':probe_means, 'n_beads':self.n_beads},
                index=illumina_ids,
                columns=['mean_value','n_beads'],
                dtype=self.bit, # int16 could work, and reduce memory by 1/2, but some raw values were > 32127 -- without prenormalization, you get negative values back, which breaks stuff.
            )
        else:
            # casting astype(int) fixes pandas bug reading uint16, by coercing to float32-compatible format.
            probe_records = dict(zip(illumina_ids, probe_means.astype(int)))
            data_frame = pd.DataFrame.from_dict(
                data=probe_records,
                orient='index',
                columns=['mean_value'],
                dtype=self.bit, # int16 could work, and reduce memory by 1/2, but some raw values were > 32127 -- without prenormalization, you get negative values back, which breaks stuff.
            )
            data_frame.index.name = 'illumina_id'

        if self.bit == 'float16':
            data_frame = data_frame.clip(upper=32127)
            data_frame = data_frame.astype('int16')

        return data_frame


    def meta(self, idat_file):
        """To enable this, initialize idatDataset with verbose=True"""
        section_offsets = self.get_section_offsets(idat_file)
        def seek_to_section(section_code):
            offset = section_offsets[section_code.value]
            idat_file.seek(offset)

        idat_file.seek(IdatHeaderLocation.VERSION.value)
        idat_version = read_long(idat_file)
        print(f"idat version: {idat_version}")
        print("file includes [illumina_id | mean | std_dev | nbeads] tabular data for all probes.")
        #seek_to_section(IdatSectionCode.ILLUMINA_ID)
        #illumina_ids = npread(idat_file, '<i4', self.n_snps_read)
        #print(f"102 (count of) illumina_ids: {len(illumina_ids)}")
        #seek_to_section(IdatSectionCode.STD_DEV)
        #print(f"103 std_dev: {read_string(idat_file)}")
        #seek_to_section(IdatSectionCode.MEAN)
        #print(f"104 mean: {read_string(idat_file)}")
        seek_to_section(IdatSectionCode.STD_DEV)
        probe_std = npread(idat_file, '<u2', self.n_snps_read)
        print(f"103 std_dev average: {round(np.mean(probe_std),1)}, N={len(probe_std)}")
        print(f"107 num_beads average: {round(np.mean(self.n_beads),1)}, N={len(self.n_beads)}")
        seek_to_section(IdatSectionCode.RED_GREEN) #400
        print(f"400 red_green: {read_string(idat_file)}")
        seek_to_section(IdatSectionCode.MOSTLY_NULL) #401
        print(f"401 manifest: {read_string(idat_file)}")
        print(f"402 barcode: {self.barcode}")
        print(f"403 chip_type: {self.chip_type}")
        print(f"1000 n_snps_read: {self.n_snps_read}")

        for line in self.run_info:
            if line[1] == 'Scan':
                print(f"300 run info: {line}")
                break

    def overflow_check(self):
        if hasattr(self, 'probe_means'):
            if (self.probe_means.values < 0).any():
                # n_affected = self.probe_means[self.probe_means.mean_value < 0].count().values[0]
                return False
        return True # passes, no misread probes
