# Lib
import logging
from pathlib import PurePath, Path
from urllib.parse import urlparse, urlunparse

LOGGER = logging.getLogger(__name__)
REQUIRED = ['Sentrix_ID', 'Sentrix_Position', 'SentrixBarcode_A', 'SentrixPosition_A', 'Control',
    'Sample_Group', 'Sample_Name', 'Sample_Plate', 'Pool_ID', 'Sample_Well', 'GSM_ID',
    'Sample_Type', 'Sub_Type']

class Sample():
    """Object representing a row in a SampleSheet file

Arguments:
    data_dir {string or path-like} -- Base directory of the sample sheet and associated IDAT files.
    sentrix_id {string} -- The slide number of the processed array.
    sentrix_position {string} -- The position on the processed slide.

Keyword Arguments:
    addl_fields {} -- Additional metadata describing the sample.
    including experiment subject meta data:

        name (sample name, unique id)
        Sample_Type
        Control
        GSM_ID (same as sample name if using GEO public data)

    array meta data:

        group
        plate
        pool
        well
    """

    def __init__(self, data_dir, sentrix_id, sentrix_position, **addl_fields):
        self.data_dir = data_dir
        self.sentrix_id = sentrix_id
        self.sentrix_position = sentrix_position
        self.renamed_fields = {}

        # any OTHER sample_sheet columns are passed in exactly as they appear, if possible, and if column names exist.
        # these will pass into the meta_data pkl created, and any renamed fields must be noted in a lookup.
        for field in addl_fields:
            if field not in REQUIRED:
                new_field_name = field.replace(' ','_')
                if len(field) == 0:
                    continue
                if field[0].isdigit():
                    new_field_name = field[1:]
                if not field.isalnum(): # letters or numbers, or caps. no spaces or unicode
                    import re
                    new_field_name = re.sub(r'\W+', '', new_field_name)
                setattr(self, new_field_name, addl_fields[field])
                self.renamed_fields[field] = new_field_name
        self.group = addl_fields.get('Sample_Group')
        self.name = addl_fields.get('Sample_Name')
        self.plate = addl_fields.get('Sample_Plate')
        self.pool = addl_fields.get('Pool_ID')
        self.well = addl_fields.get('Sample_Well')
        self.GSM_ID = addl_fields.get('GSM_ID') # for GEO published sample compatability
        self.type = addl_fields.get('Sample_Type','Unknown') # from GEO MINiML meta data
        self.sub_type = addl_fields.get('Sub_Type') # from GEO
        self.is_control = True if addl_fields.get('Control') in (1,'1',True, 'True', 'true', 'TRUE') else False
        self.fields = {}
        self.fields.update(self.renamed_fields)
        self.fields.update({
            'Sentrix_ID': 'Sentrix_ID',
            'Sentrix_Position': 'Sentrix_Position', # these will be standardized here, regardless of sample_sheet variation names
            'Sample_Group': 'Sample_Group',
            'Sample_Name': 'Sample_Name',
            'Sample_Plate': 'Sample_Plate',
            'Sample_Type': 'Sample_Type',
            'Sub_Type': 'Sub_Type',
            'Sample_Well': 'Sample_Well',
            'Pool_ID': 'Pool_ID',
            'GSM_ID': 'GSM_ID',
            'Control': 'Control',
        })

    def __str__(self):
        return f'{self.sentrix_id}_{self.sentrix_position}'

    @property
    def base_filename(self):
        return f'{self.sentrix_id}_{self.sentrix_position}'

    @property
    def alternate_base_filename(self):
        """ GEO data sets using this file name convention."""
        if hasattr(self,'GSM_ID') and self.GSM_ID != None:
            return f'{self.GSM_ID}_{self.sentrix_id}_{self.sentrix_position}'
        else:
            return f'{self.sentrix_id}_{self.sentrix_position}'

    def get_filepath(self, extension, suffix=None, verify=True):
        """builds the filepath based on custom file extensions and suffixes during processing.

        Params (verify):
            tests whether file exists, either in data_dir or somewhere in recursive search path of data_dir.
        Export:
            uses this later to fetch the place where a file ought to be created -- but doesn't exist yet, so use verify=False.

        Notes:
            _suffix -- used to create the `<file>_processed` files.
        """
        _suffix = ''
        if suffix is not None:
            _suffix = f'_{suffix}'

        filename = f'{self.base_filename}{_suffix}.{extension}'
        alt_filename = f'{self.alternate_base_filename}{_suffix}.{extension}'
        path = PurePath(self.data_dir, str(self.sentrix_id), filename)
        if verify:
            # confirm this sample IDAT file exists, and update its filepath if different.
            # if filename fails, it will check alt_filename too.
            path = self._build_and_verify_path(filename, alt_filename, allow_compressed=True)
        return str(path)

    def get_file_s3(self, zip_reader, extension, suffix=None):
        """ replaces get_filepath, but for `s3` context. Since these files
        are compressed within a single zipfile in the bucket, they don't
        resolve to PurePaths."""
        _suffix = ''
        if suffix is not None:
            _suffix = f'_{suffix}'

        filename_to_match = f'{self.base_filename}{_suffix}.{extension}'
        for zip_filename in zip_reader.file_names:
            if not zip_filename.endswith('.idat'):
                continue
            if filename_to_match in zip_filename:
                # this is packed within the zipfile still, but zip_reader can fetch it.
                LOGGER.info(zip_reader.get_file_info(zip_filename))
                return zip_reader.get_file(zip_filename, match_partial=False)

    def _build_and_verify_path(self, filename, alt_filename=None, allow_compressed=False):
        """
        Added to Sample as class_method:
            if the matching filename for idat file is not in the same folder.
            check if exists:
            then look recursively for that filename and update the data_dir for that Sample.
            return the complete filepath.

        alt_filename:
            because public data sets on GEO have samplesheets with a different pattern, if the primary file pattern
            fails to match, it will try the alt_filename pattern before returning a FileNotFoundError.

        replaces _build_path
        """
        same_dir_path = PurePath(self.data_dir, str(self.sentrix_id), filename)
        if Path(same_dir_path).is_file():
            # this idat file is in the same folder, no more searching needed.
            return same_dir_path

        if allow_compressed and Path(same_dir_path.with_suffix('.gz')).is_file():
            return same_dir_path

        # otherwise, do a recursive search for this file and return the first path found.
        #file_pattern = f'{self.data_dir}/**/{filename}'
        #file_matches = glob(file_pattern, recursive=True)
        file_matches = list(Path(self.data_dir).rglob(filename))
        if (not file_matches) and allow_compressed:
            file_matches = list(Path(self.data_dir).rglob(filename + '.gz'))
        if file_matches == []:
            if alt_filename != None and alt_filename != filename:
                # Note: both patterns will be identical if GSM_ID missing from sample sheet.
                alt_file_matches = list(Path(self.data_dir).rglob(alt_filename))
                if (not alt_file_matches) and allow_compressed:
                    alt_file_matches = list(Path(self.data_dir).rglob(alt_filename + '.gz'))
                if len(alt_file_matches) > 1:
                    LOGGER.warning(f'Multiple ({len(alt_file_matches)}) files matched {alt_filename} -- saved path to first one: {alt_file_matches[0]}')
                if len(alt_file_matches) > 0:
                    return alt_file_matches[0]
            raise FileNotFoundError(f'No files in {self.data_dir} (or sub-folders) match this sample id: {filename} OR {alt_filename}')
        elif len(file_matches) > 1:
            LOGGER.warning(f'Multiple ({len(file_matches)}) files matched {alt_filename} -- saved path to first one: {file_matches[0]}')
        return file_matches[0]

    def get_export_filepath(self, extension='csv'):
        """ Called by run_pipeline to find the folder/filename to export data as CSV, but CSV file doesn't exist yet."""
        return self.get_filepath(extension, 'processed', verify=False)
