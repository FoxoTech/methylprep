# Lib
import pytest

# App
from methylprep.files import IdatDataset


class TestIdatModel(object):
    def test_init_throws_with_no_args(self):
        with pytest.raises(TypeError):
            IdatDataset()
