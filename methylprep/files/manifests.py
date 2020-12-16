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

MANIFEST_DIR_NAME = '.methylprep_manifest_files'
MANIFEST_DIR_PATH = f'~/{MANIFEST_DIR_NAME}'
MANIFEST_DIR_PATH_LAMBDA = f'/tmp/{MANIFEST_DIR_NAME}'
MANIFEST_BUCKET_NAME = 'array-manifest-files'
MANIFEST_REMOTE_PATH = f'https://s3.amazonaws.com/{MANIFEST_BUCKET_NAME}/'

ARRAY_TYPE_MANIFEST_FILENAMES = {
    ArrayType.ILLUMINA_27K: 'hm27.hg19.manifest.csv.gz', #'humanmethylation27_270596_v1-2.csv.gz',
    ArrayType.ILLUMINA_450K: 'HumanMethylation450_15017482_v1-2.CoreColumns.csv.gz',
    ArrayType.ILLUMINA_EPIC: 'MethylationEPIC_v-1-0_B4.CoreColumns.csv.gz',
    ArrayType.ILLUMINA_EPIC_PLUS: 'CombinedManifestEPIC.manifest.CoreColumns.csv.gz',
    ArrayType.ILLUMINA_MOUSE: 'LEGX_C20_manifest_mouse_min.csv.gz',
}
ARRAY_FILENAME = {
    '27k': 'hm27.hg19.manifest.csv.gz',
    '450k': 'HumanMethylation450_15017482_v1-2.CoreColumns.csv.gz',
    'epic': 'MethylationEPIC_v-1-0_B4.CoreColumns.csv.gz',
    'epic+': 'CombinedManifestEPIC.manifest.CoreColumns.csv.gz',
    'mouse': 'LEGX_C20_manifest_mouse_min.csv.gz',
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

MOUSE_MANIFEST_COLUMNS = (
    'IlmnID',
    'AddressA_ID',
    'AddressB_ID',
    'Infinium_Design_Type',
    'Color_Channel',
    'Genome_Build',
    'CHR',
    'MAPINFO',
    'Strand',
    'Probe_Type', # additional, needed to identify mouse-specific probes (mu) | and control probe sub_types
)

CONTROL_COLUMNS = (
    'Address_ID',
    'Control_Type',
    'Color',
    'Extended_Type',
    # control probes don't have 'IlmnID' values set -- these probes are not locii specific
    # these column names don't appear in manifest. they are added when importing the control section of rows
)


class Manifest():
    """Provides an object interface to an Illumina array manifest file.

    Arguments:
        array_type {ArrayType} -- The type of array to process.
        values are styled like ArrayType.ILLUMINA_27K, ArrayType.ILLUMINA_EPIC or ArrayType('epic'), ArrayType('mouse')

    Keyword Arguments:
        filepath_or_buffer {file-like} -- a pre-existing manifest filepath (default: {None})

    Raises:
        ValueError: The sample sheet is not formatted properly or a sample cannot be found.
    """

    __genome_df = None
    __probe_type_subsets = None # apparently not used anywhere in methylprep

    def __init__(self, array_type, filepath_or_buffer=None, on_lambda=False):
        self.array_type = array_type
        self.on_lambda = on_lambda # changes filepath to /tmp for the read-only file system

        if filepath_or_buffer is None:
            filepath_or_buffer = self.download_default(array_type, self.on_lambda)

        with get_file_object(filepath_or_buffer) as manifest_file:
            self.__data_frame = self.read_probes(manifest_file)
            self.__control_data_frame = self.read_control_probes(manifest_file)
            self.__snp_data_frame = self.read_snp_probes(manifest_file)
            if self.array_type == ArrayType.ILLUMINA_MOUSE:
                self.__mouse_data_frame = self.read_mouse_probes(manifest_file)
            else:
                self.__mouse_data_frame = pd.DataFrame()

    @property
    def columns(self):
        if self.array_type == ArrayType.ILLUMINA_MOUSE:
            return MOUSE_MANIFEST_COLUMNS
        else:
            return MANIFEST_COLUMNS

    @property
    def data_frame(self):
        return self.__data_frame

    @property
    def control_data_frame(self):
        return self.__control_data_frame

    @property
    def snp_data_frame(self):
        return self.__snp_data_frame

    @property
    def mouse_data_frame(self):
        return self.__mouse_data_frame

    @staticmethod
    def download_default(array_type, on_lambda=False):
        """Downloads the appropriate manifest file if one does not already exist.

        Arguments:
            array_type {ArrayType} -- The type of array to process.

        Returns:
            [PurePath] -- Path to the manifest file.
        """
        dir_path = Path(MANIFEST_DIR_PATH).expanduser()
        if on_lambda:
            dir_path = Path(MANIFEST_DIR_PATH_LAMBDA).expanduser()
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
        LOGGER.info(f'Reading manifest file: {Path(manifest_file.name).stem}')

        self.seek_to_start(manifest_file)

        #print(f"DEBUG read_probes {self.get_data_types()}")
        #print(f"{self.columns}, {self.array_type.num_probes - 1}")
        data_frame = pd.read_csv(
            manifest_file,
            comment='[',
            dtype=self.get_data_types(),
            usecols=self.columns,
            nrows=self.array_type.num_probes - 1, # -1 because every array.num_probes is one more than the total number of rows (dunno why -- Byerly's work) found in manifest.
            # the -1 applies if the manifest has one extra row between the cg and control probes (a [Controls],,,,,, row)
            index_col='IlmnID',
        )

        # AddressB_ID in manifest includes NaNs and INTs and becomes floats, which breaks. forcing back here.
        #data_frame['AddressB_ID'] = data_frame['AddressB_ID'].astype('Int64') # converts floats to ints; leaves NaNs inplace
        # TURNS out, int or float both work for manifests. NOT the source of the error with mouse.
        #LOGGER.info('AddressB_ID')
        #LOGGER.info(f"{data_frame['AddressB_ID']}")

        def get_probe_type(name, infinium_type):
            """what:
        returns one of (I, II, SnpI, SnpII, Control)

        how:
            .from_manifest_values() returns probe type using either
            the Infinium_Design_Type (I or II) or the name (starts with 'rs' == SnpI)
            and 'Control' is none of the above."""
            probe_type = ProbeType.from_manifest_values(name, infinium_type)
            return probe_type.value

        vectorized_get_type = np.vectorize(get_probe_type)

        data_frame['probe_type'] = vectorized_get_type(
            data_frame.index.values,
            data_frame['Infinium_Design_Type'].values,
        )

        #print((f"""DEBUG read_manifest probe types: Control {data_frame[data_frame['probe_type'].str.contains('Control')].shape} I {data_frame[data_frame['probe_type'].str.contains('I')].shape} """
        #      f"""II {data_frame[data_frame['probe_type'].str.contains('II')].shape} SnpI {data_frame[data_frame['probe_type'].str.contains('SnpI')].shape} """
        #      f"""SnpII {data_frame[data_frame['probe_type'].str.contains('SnpII')].shape}"""))

        return data_frame

    def read_control_probes(self, manifest_file):
        """ Unlike other probes, control probes have no IlmnID because they're not locus-specific.
        they also use arbitrary columns, ignoring the header at start of manifest file. """
        #LOGGER.info(f'Reading control probes: {Path(manifest_file.name).stem}')

        self.seek_to_start(manifest_file)

        #num_headers = 4 -- removed on Jan 21 2020.
        # Steve Byerly added this =4, but the manifests seem to start controls at the point specified, with no extra columns based on .num_probes stored.
        num_headers = 0

        #print(f"CONTROLS nrows {self.array_type.num_controls} SKIP {self.array_type.num_probes + num_headers} usecols {range(len(CONTROL_COLUMNS))}")
        return pd.read_csv(
            manifest_file,
            comment='[',
            header=None,
            index_col=0, # illumina_id, not IlmnID here
            names=CONTROL_COLUMNS, # this gives these columns new names, because they have none. loading stuff at end of CSV after probes end.
            nrows=self.array_type.num_controls,
            skiprows=self.array_type.num_probes + num_headers,
            usecols=range(len(CONTROL_COLUMNS)),
        )

    def read_snp_probes(self, manifest_file):
        """ Unlike cpg and control probes, these rs probes are NOT sequential in all arrays. """
        #LOGGER.info(f'Reading snp probes: {Path(manifest_file.name).stem} --> {snp_df.shape[0]} found')
        self.seek_to_start(manifest_file)
        # since these are not sequential, loading everything and filtering by IlmnID.
        snp_df = pd.read_csv(
            manifest_file,
            low_memory=False)
        snp_df = snp_df[snp_df['IlmnID'].str.match('rs', na=False)]
        return snp_df

    def read_mouse_probes(self, manifest_file):
        """ ILLUMINA_MOUSE contains unique probes whose names begin with 'mu' and 'rp'
        for 'murine' and 'repeat', respectively. This creates a dataframe of these probes,
        which are not processed like normal cg/ch probes. """
        self.seek_to_start(manifest_file)
        mouse_df = pd.read_csv(
            manifest_file,
            low_memory=False) # low_memory=Fase is required because control probes create mixed-types in columns.
        mouse_df = mouse_df[(mouse_df['Probe_Type'] == 'rp') | (mouse_df['IlmnID'].str.startswith('uk', na=False)) | (mouse_df['Probe_Type'] == 'mu')] # 'mu' probes now start with 'cg' instead and have 'mu' in another column
        return mouse_df

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
            key: str for key in self.columns
        }
        data_types['AddressA_ID'] = 'Int64' #'float64' -- dtype found only in pandas 0.24 or greater
        data_types['AddressB_ID'] = 'Int64' #'float64'
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
        """given a probe type (I, II, SnpI, SnpII, Control) and a channel (Channel.RED | Channel.GREEN),
        This will return info needed to map probes to their names (e.g. cg0031313 or rs00542420),
        which are NOT in the idat files."""
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
