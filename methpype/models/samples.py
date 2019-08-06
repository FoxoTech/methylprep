# Lib
import logging
from pathlib import PurePath, Path
from urllib.parse import urlparse, urlunparse
from glob import glob

LOGGER = logging.getLogger(__name__)

class Sample():
    """Object representing a row in a SampleSheet file

    Arguments:
        data_dir {string or path-like} -- Base directory of the sample sheet and associated IDAT files.
        sentrix_id {string} -- The slide number of the processed array.
        sentrix_position {string} -- The position on the processed slide.

    Keyword Arguments:
        addl_fields {} -- Additional metadata describing the sample.
    """

    def __init__(self, data_dir, sentrix_id, sentrix_position, **addl_fields):
        self.data_dir = data_dir
        self.sentrix_id = sentrix_id
        self.sentrix_position = sentrix_position

        self.group = addl_fields.get('Sample_Group')
        self.name = addl_fields.get('Sample_Name')
        self.plate = addl_fields.get('Sample_Plate')
        self.pool = addl_fields.get('Pool_ID')
        self.well = addl_fields.get('Sample_Well')
        self.GSM_ID = addl_fields.get('GSM_ID') # for GEO published sample compatability

    def __str__(self):
        return f'{self.sentrix_id}_{self.sentrix_position}'

    @property
    def base_filename(self):
        return f'{self.sentrix_id}_{self.sentrix_position}'

    @property
    def alternate_base_filename(self):
        if hasattr(self,'GSM_ID'):
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
            path = self._build_and_verify_path(filename, alt_filename)
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
                print(zip_reader.get_file_info(zip_filename))
                return zip_reader.get_file(zip_filename, match_partial=False)

    def _build_and_verify_path(self, filename, alt_filename=None):
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

        # otherwise, do a recursive search for this file and return the first path found.
        file_pattern = f'{self.data_dir}/**/{filename}'
        file_matches = glob(file_pattern, recursive=True)
        if file_matches == []:
            if alt_filename != None and alt_filename != filename:
                # Note: both patterns will be identical if GSM_ID missing from sample sheet.
                alt_file_pattern = f'{self.data_dir}/**/{alt_filename}'
                alt_file_matches = glob(alt_file_pattern, recursive=True)
                if len(alt_file_matches) > 1:
                    LOGGER.warning(f'Multiple ({len(alt_file_matches)}) files matched {alt_file_pattern} -- saved path to first one: {alt_file_matches[0]}')
                if len(alt_file_matches) > 0:
                    return alt_file_matches[0]
            raise FileNotFoundError(f'No files in {self.data_dir} (or sub-folders) match this sample id: {filename}')
        elif len(file_matches) > 1:
            LOGGER.warning(f'Multiple ({len(file_matches)}) files matched {file_pattern} -- saved path to first one: {file_matches[0]}')
        return file_matches[0]

    def get_export_filepath(self):
        """ Called by run_pipeline to find the folder/filename to export data as CSV, but CSV file doesn't exist yet."""
        return self.get_filepath('csv', 'processed', verify=False)
