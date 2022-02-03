# Lib
from pathlib import Path
import pytest

# App
from methylprep.files import IdatDataset
from methylprep.models import Channel

class TestIdatModel(object):
    test_data_dir = 'docs/example_data/GSE69852'
    test_idat_file = str(Path(test_data_dir, '9247377093_R02C01_Grn.idat'))

    def test_init_throws_with_no_args(self):
        with pytest.raises(TypeError):
            IdatDataset()

    def test_init_filepath_str(self):
        idat = IdatDataset(self.test_idat_file, channel=Channel.GREEN)

    def test_init_fails_if_passing_in_open_file_object(self):
        with pytest.raises(UnicodeDecodeError):
            with open(self.test_idat_file,'r') as f:
                idat = IdatDataset(f, Channel.GREEN)

    def test_idat_full_meta(self):
        idat = IdatDataset(
            self.test_idat_file,
            channel=Channel.GREEN,
            verbose=True, # <--- calls IdatDataset.meta()
            std_dev=True,
            nbeads=True,
            bit='float16')
