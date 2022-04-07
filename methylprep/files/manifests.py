# Lib
import logging
from pathlib import Path, PurePath
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

MANIFEST_DIR_NAME = 'methylprep_manifest_files/'
# Download manifests into library so multiple users don't need their own copy
MANIFEST_DIR_PATH = PurePath(Path(__file__).parent.resolve(), MANIFEST_DIR_NAME)
MANIFEST_DIR_PATH_LAMBDA = f'/tmp/{MANIFEST_DIR_NAME}'
MANIFEST_BUCKET_NAME = 'array-manifest-files'
MANIFEST_REMOTE_PATH = f'https://s3.amazonaws.com/{MANIFEST_BUCKET_NAME}/'

ARRAY_FILENAME = {
    '27k': 'hm27.hg19.manifest.csv.gz',
    #'humanmethylation27_270596_v1-2.csv.gz',
    '450k': 'HumanMethylation450k_15017482_v3.csv.gz',
    #'HumanMethylation450_15017482_v1-2.CoreColumns.csv.gz',
    'epic': 'HumanMethylationEPIC_manifest_v2.csv.gz',
    #'MethylationEPIC_v-1-0_B4.CoreColumns.csv.gz',
    'epic+': 'CombinedManifestEPIC_manifest_CoreColumns_v2.csv.gz',
    #'CombinedManifestEPIC.manifest.CoreColumns.csv.gz',
    'mouse': 'MM285_manifest_v3.csv.gz',
    #'MM285_mm39_manifest_v2.csv.gz',
    ###### BE SURE TO ALSO UPDATE arrays.py ArrayType.num_controls if updating a manifest here. #######
}
ARRAY_TYPE_MANIFEST_FILENAMES = {
    ArrayType.ILLUMINA_27K: ARRAY_FILENAME['27k'],
    ArrayType.ILLUMINA_450K: ARRAY_FILENAME['450k'],
    ArrayType.ILLUMINA_EPIC: ARRAY_FILENAME['epic'],
    ArrayType.ILLUMINA_EPIC_PLUS: ARRAY_FILENAME['epic+'],
    ArrayType.ILLUMINA_MOUSE: ARRAY_FILENAME['mouse'],
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
    'OLD_Genome_Build',
    'OLD_CHR',
    'OLD_MAPINFO',
    'OLD_Strand',
)

MOUSE_MANIFEST_COLUMNS = (
    'IlmnID',
    'AddressA_ID',
    'AddressB_ID',
    'Infinium_Design_Type',
    'Color_Channel',
    'design', # replaces Probe_Type in v1.4.6+ with tons of design meta data. only 'Random' and 'Multi' matter in code.
    #'Probe_Type', # pre v1.4.6, needed to identify mouse-specific probes (mu) | and control probe sub_types
    'Genome_Build',
    'CHR',
    'MAPINFO',
    'Strand',
    'OLD_Genome_Build',
    'OLD_CHR',
    'OLD_MAPINFO',
    'OLD_Strand',
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

    def __init__(self, array_name: str, filepath_or_buffer=None, on_lambda=False, verbose=True):
        self.array_name = array_name
        array_str_to_class = dict(zip(list(ARRAY_FILENAME.keys()), list(ARRAY_TYPE_MANIFEST_FILENAMES.keys())))
        # Map array_name as provided as a string to self.array_type which will be of type enum 'ArrayType'
        if array_name in list(array_str_to_class.keys()):
            self.array_type = array_str_to_class[array_name]
        else:
            self.array_type = ArrayType.CUSTOM
        self.on_lambda = on_lambda # changes filepath to /tmp for the read-only file system
        self.verbose = verbose
        # If you don't pass in a filepath, Methylprep uses one of its own manifests which may need downloading
        if filepath_or_buffer is None:
            filepath_or_buffer = self.download_default(self.array_type, self.on_lambda)
            LOGGER.info(f"Using Methylprep manifest {self.array_name} @ {filepath_or_buffer}")
            # This approach to loading a manifest only makes sense if array type is configured in arrays.py
            # which is only true for Methylprep manifests
            with get_file_object(filepath_or_buffer) as manifest_file:
                self.__data_frame = self.read_probes(manifest_file) # DF uses IlmnID as index
                self.__control_data_frame = self.read_control_probes(manifest_file) # DF uses Address_ID as index
                self.__snp_data_frame = self.read_snp_probes(manifest_file) # DF uses neither IlmnID nor Address as index?
                if self.array_type == ArrayType.ILLUMINA_MOUSE:
                    self.__mouse_data_frame = self.read_mouse_probes(manifest_file)
                else:
                    self.__mouse_data_frame = pd.DataFrame()
        else:
            LOGGER.info(f"Using CUSTOM manifest {self.array_name} @ {filepath_or_buffer}")
            self._load_custom_manifest(filepath_or_buffer)

        LOGGER.info(f"""Manifest has {str(len(self.__data_frame))} probes (where {str(len(self.__snp_data_frame))} are snps) and {str(len(self.__control_data_frame))} control probes""")

    def _load_custom_manifest(self, custom_manifest_path):
        optional = ['OLD_CHR', 'OLD_Strand', 'OLD_Genome_Build', 'OLD_MAPINFO']
        minimum_cols = [col for col in self.columns if col not in optional]

        with get_file_object(custom_manifest_path) as manifest_file:
            file_df = pd.read_csv(manifest_file, comment='[', low_memory=False)
            missing_cols = [col for col in minimum_cols if col not in (file_df.columns.to_list()+[file_df.index.name])]
            if missing_cols != []:
                raise ValueError(f"Custom manifest {custom_manifest_path} missing {','.join(missing_cols)} columns")
            # Controls in manifests follow column order of 'Address_ID', 'Control_Type', 'Color', 'Extended_Type',
            # After row starting with [Controls],,,, 
            # Non control probes follow column order of 'IlmnID', 'Name', 'AddressA_ID', 'AddressB_ID'...
            # This means that after loading the file, any row where AddressB_ID is a string is a control probe
            control_probes = file_df[file_df["AddressB_ID"].str.contains('[^0-9.]', regex=True, na=False)]
            control_probes = control_probes[["IlmnID", "Name", "AddressA_ID", "AddressB_ID"]]
            self.__control_data_frame = control_probes.rename(
                columns={
                    "IlmnID": "Address_ID",
                    "Name": "Control_Type",
                    "AddressA_ID": "Color",
                    "AddressB_ID": "Extended_Type"
                }
            )
            # Snps have IlmnID with leading rs
            self.__snp_data_frame = file_df[file_df["IlmnID"].str.match('rs', na=False)].astype(
                {'AddressA_ID':'float64', 'AddressB_ID':'float64'})
            # At this point self.__data_frame should be any probe that is not a control
            self.__data_frame = file_df.drop(self.__control_data_frame.index, axis=0).astype(
                {'AddressA_ID':'float64', 'AddressB_ID':'float64'})
            self.__data_frame = self.__data_frame.set_index('IlmnID')
            # Methylprep wants index of control df to be Address_ID
            self.__control_data_frame = self.__control_data_frame.astype({'Address_ID': 'int64'})
            self.__control_data_frame = self.__control_data_frame.set_index('Address_ID')
            # Do transformation on Infinim_Design_Type column 
            vectorized_get_type = np.vectorize(self.get_probe_type)
            self.__data_frame['probe_type'] = vectorized_get_type(
                self.__data_frame.index.values,
                self.__data_frame['Infinium_Design_Type'].values,
            )
            # IDK it seems like methylsuite wants this variable avaialable?
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
    def get_probe_type(name, infinium_type):
        """returns one of (I, II, SnpI, SnpII, Control)

        .from_manifest_values() returns probe type using either the Infinium_Design_Type (I or II) or the name
        (starts with 'rs' == SnpI) and 'Control' is none of the above."""
        probe_type = ProbeType.from_manifest_values(name, infinium_type)
        return probe_type.value

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

        LOGGER.info(f"Downloading manifest: {Path(filename).stem}")
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
            if not header_line: #EOF
                raise EOFError("The first (left-most) column in your manifest must contain 'IlmnID'. This defines the header row.")
            header_line = manifest_file.readline()

        if current_pos == 0:
            manifest_file.seek(current_pos)
        else:
            manifest_file.seek(current_pos - 1)

    def read_probes(self, manifest_file):
        if self.verbose:
            LOGGER.info(f'Reading manifest file: {Path(manifest_file.name).stem}')

        try:
            data_frame = pd.read_csv(
                manifest_file,
                comment='[',
                dtype=self.get_data_types(),
                usecols=self.columns,
                nrows=self.array_type.num_probes,
                # the -1 applies if the manifest has one extra row between the cg and control probes (a [Controls],,,,,, row) --- fixed in v1.5.6
                index_col='IlmnID',
            )
        except ValueError:
            optional = ['OLD_CHR', 'OLD_Strand', 'OLD_Genome_Build', 'OLD_MAPINFO']
            use_columns = [col for col in self.columns if col not in optional]
            data_frame = pd.read_csv(
                manifest_file,
                comment='[',
                dtype=self.get_data_types(),
                usecols=use_columns,
                nrows=self.array_type.num_probes,
                # the -1 applies if the manifest has one extra row between the cg and control probes (a [Controls],,,,,, row) --- fixed in v1.5.6
                index_col='IlmnID',
            )
            LOGGER.info(f"Some optional genome mapping columns were not found in {Path(manifest_file.name).stem}")
        # AddressB_ID in manifest includes NaNs and INTs and becomes floats, which breaks. forcing back here.
        #data_frame['AddressB_ID'] = data_frame['AddressB_ID'].astype('Int64') # converts floats to ints; leaves NaNs inplace
        # TURNS out, int or float both work for manifests. NOT the source of the error with mouse.
        vectorized_get_type = np.vectorize(self.get_probe_type)
        data_frame['probe_type'] = vectorized_get_type(
            data_frame.index.values,
            data_frame['Infinium_Design_Type'].values,
        )
        return data_frame

    def read_control_probes(self, manifest_file):
        """ Unlike other probes, control probes have no IlmnID because they're not locus-specific.
        they also use arbitrary columns, ignoring the header at start of manifest file. """
        #LOGGER.info(f'Reading control probes: {Path(manifest_file.name).stem}')

        self.seek_to_start(manifest_file)

        return pd.read_csv(
            manifest_file,
            comment='[',
            header=None,
            index_col=0, # illumina_id, not IlmnID here
            names=CONTROL_COLUMNS, # this gives these columns new names, because they have none. loading stuff at end of CSV after probes end.
            nrows=self.array_type.num_controls,
            skiprows=self.array_type.num_probes+1, #without the +1, it includes the last cpg probe in controls and breaks stuff.
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
        # 'O' type columns won't match in SigSet, so forcing float64 here. Also, float32 won't cover all probe IDs; must be float64.
        snp_df = snp_df[snp_df['IlmnID'].str.match('rs', na=False)].astype({'AddressA_ID':'float64', 'AddressB_ID':'float64'})
        return snp_df

    def read_mouse_probes(self, manifest_file):
        """ ILLUMINA_MOUSE contains unique probes whose names begin with 'mu' and 'rp'
        for 'murine' and 'repeat', respectively. This creates a dataframe of these probes,
        which are not processed like normal cg/ch probes. """
        self.seek_to_start(manifest_file)
        mouse_df = pd.read_csv(
            manifest_file,
            low_memory=False) # low_memory=Fase is required because control probes create mixed-types in columns.
        #--- pre v1.4.6: mouse_df = mouse_df[(mouse_df['Probe_Type'] == 'rp') | (mouse_df['IlmnID'].str.startswith('uk', na=False)) | (mouse_df['Probe_Type'] == 'mu')]
        #--- pre v1.4.6: 'mu' probes start with 'cg' instead and have 'mu' in Probe_Type column
        mouse_df = mouse_df[(mouse_df['design'] == 'Multi') | (mouse_df['design'] == 'Random')]
        return mouse_df

    def get_data_types(self):
        data_types = {
            key: str for key in self.columns
        }
        data_types['AddressA_ID'] = 'Int64' #'float64' -- dtype found only in pandas 0.24 or greater
        data_types['AddressB_ID'] = 'Int64' #'float64'
        return data_types

    def get_probe_details(self, probe_type, channel=None):
        """used by infer_channel_switch. Given a probe type (I, II, SnpI, SnpII, Control) and a channel (Channel.RED | Channel.GREEN),
        this will return info needed to map probes to their names (e.g. cg0031313 or rs00542420), which are NOT in the idat files."""
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