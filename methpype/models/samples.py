# Lib
from pathlib import PurePath


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

    def __str__(self):
        return f'{self.sentrix_id}_{self.sentrix_position}'

    @property
    def base_filename(self):
        return f'{self.sentrix_id}_{self.sentrix_position}'

    def get_filepath(self, extension, suffix=None):
        _suffix = ''
        if suffix is not None:
            _suffix = f'_{suffix}'

        filename = f'{self.base_filename}{_suffix}.{extension}'
        path = self._build_path(filename)
        return str(path)

    def get_export_filepath(self):
        return self.get_filepath('csv', 'processed')

    def _build_path(self, filename):
        return PurePath(self.data_dir, str(self.sentrix_id), filename)
