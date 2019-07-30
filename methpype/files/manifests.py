# Lib
import logging
from pathlib import Path
from urllib.parse import urljoin
import numpy as np
import pandas as pd
# App
from ..models import ArrayType, Channel, ProbeType
from ..utils import (
    download_file,
    get_file_object,
    inner_join_data,
    is_file_like,
    reset_file,
)


__all__ = ['Manifest']


LOGGER = logging.getLogger(__name__)

MANIFEST_DIR_NAME = '.methpype_manifest_files'
MANIFEST_DIR_PATH = f'~/{MANIFEST_DIR_NAME}'
MANIFEST_BUCKET_NAME = 'array-manifest-files'
MANIFEST_REMOTE_PATH = f'https://s3.amazonaws.com/{MANIFEST_BUCKET_NAME}/'

ARRAY_TYPE_MANIFEST_FILENAMES = {
    ArrayType.ILLUMINA_27K: 'hm27.hg19.manifest.csv.gz', #'humanmethylation27_270596_v1-2.csv.gz',
    ArrayType.ILLUMINA_450K: 'HumanMethylation450_15017482_v1-2.CoreColumns.csv.gz',
    ArrayType.ILLUMINA_EPIC: 'MethylationEPIC_v-1-0_B4.CoreColumns.csv.gz',
    ArrayType.ILLUMINA_EPIC_PLUS: 'CombinedManifestEPIC.manifest.CoreColumns.csv.gz',
}

MANIFEST_COLUMNS = (
    'IlmnID',
    'AddressA_ID',
    'AddressB_ID',
    'Infinium_Design_Type',
    'Color_Channel',
    'Genome_Build',
    'CHR',
    'MAPINFO',
    'Strand',
)

CONTROL_COLUMNS = (
    'Address_ID',
    'Control_Type',
    'Color',
    'Extended_Type',
)


class Manifest():
    """Provides an object interface to an Illumina array manifest file.

    Arguments:
        array_type {ArrayType} -- The type of array to process.

    Keyword Arguments:
        filepath_or_buffer {file-like} -- a pre-existing manifest filepath (default: {None})

    Raises:
        ValueError: The sample sheet is not formatted properly or a sample cannot be found.
    """

    __genome_df = None
    __probe_type_subsets = None

    def __init__(self, array_type, filepath_or_buffer=None):
        self.array_type = array_type

        if filepath_or_buffer is None:
            filepath_or_buffer = self.download_default(array_type)

        with get_file_object(filepath_or_buffer) as manifest_file:
            self.__data_frame = self.read_probes(manifest_file)
            self.__control_data_frame = self.read_control_probes(manifest_file)

    @property
    def columns(self):
        return MANIFEST_COLUMNS

    @property
    def data_frame(self):
        return self.__data_frame

    @property
    def control_data_frame(self):
        return self.__control_data_frame

    @staticmethod
    def download_default(array_type):
        """Downloads the appropriate manifest file if one does not already exist.

        Arguments:
            array_type {ArrayType} -- The type of array to process.

        Returns:
            [PurePath] -- Path to the manifest file.
        """
        dir_path = Path(MANIFEST_DIR_PATH).expanduser()
        filename = ARRAY_TYPE_MANIFEST_FILENAMES[array_type]
        filepath = Path(dir_path).joinpath(filename)

        if Path.exists(filepath):
            return filepath

        LOGGER.info('Downloading manifest: %s', filename)
        src_url = urljoin(MANIFEST_REMOTE_PATH, filename)
        download_file(filename, src_url, dir_path)

        return filepath

    @staticmethod
    def seek_to_start(manifest_file):
        """ find the start of the data part of the manifest. first left-most column must be "IlmnID" to be found."""
        reset_file(manifest_file)

        current_pos = manifest_file.tell()
        header_line = manifest_file.readline()

        while not header_line.startswith(b'IlmnID'):
            current_pos = manifest_file.tell()
            header_line = manifest_file.readline()

        if current_pos == 0:
            manifest_file.seek(current_pos)
        else:
            manifest_file.seek(current_pos - 1)

    def read_probes(self, manifest_file):
        LOGGER.info('Reading manifest file: %s', manifest_file)

        self.seek_to_start(manifest_file)

        data_frame = pd.read_csv(
            manifest_file,
            comment='[',
            dtype=self.get_data_types(),
            usecols=self.columns,
            nrows=self.array_type.num_probes,
            index_col='IlmnID',
        )

        def get_probe_type(name, infinium_type):
            probe_type = ProbeType.from_manifest_values(name, infinium_type)
            return probe_type.value

        vectorized_get_type = np.vectorize(get_probe_type)

        data_frame['probe_type'] = vectorized_get_type(
            data_frame.index.values,
            data_frame['Infinium_Design_Type'].values,
        )

        return data_frame

    def read_control_probes(self, manifest_file):
        LOGGER.info('Reading control probes from manifest file: %s', manifest_file)

        self.seek_to_start(manifest_file)

        num_headers = 4

        return pd.read_csv(
            manifest_file,
            comment='[',
            header=None,
            index_col=0,
            names=CONTROL_COLUMNS,
            nrows=self.array_type.num_controls,
            skiprows=self.array_type.num_probes + num_headers,
            usecols=range(len(CONTROL_COLUMNS)),
        )

    def map_to_genome(self, data_frame):
        genome_df = self.get_genome_data()
        merged_df = inner_join_data(data_frame, genome_df)
        return merged_df

    def get_genome_data(self):
        if self.__genome_df is not None:
            return self.__genome_df

        LOGGER.info('Building genome data frame')

        genome_columns = [
            'Genome_Build',
            'CHR',
            'MAPINFO',
            'Strand',
        ]

        self.__genome_df = self.data_frame[genome_columns]
        return self.__genome_df

    def get_data_types(self):
        data_types = {
            key: str
            for key in self.columns
        }
        data_types['AddressA_ID'] = 'float64'
        data_types['AddressB_ID'] = 'float64'
        return data_types

    def get_loci_count(self):
        """Returns the number of unique loci/identifiers in the manifest
        """
        return self.data_frame['Name'].nunique()

    def get_loci_names(self):
        """Returns the list of unique loci/identifiers in the manifest
        """
        return self.data_frame['Name'].unique()

    def get_probe_details(self, probe_type, channel=None):
        if not isinstance(probe_type, ProbeType):
            raise Exception('probe_type is not a valid ProbeType')

        if channel and not isinstance(channel, Channel):
            raise Exception('channel not a valid Channel')

        data_frame = self.data_frame
        probe_type_mask = data_frame['probe_type'].values == probe_type.value

        if not channel:
            return data_frame[probe_type_mask]

        channel_mask = data_frame['Color_Channel'].values == channel.value
        return data_frame[probe_type_mask & channel_mask]
