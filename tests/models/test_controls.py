# LIb
import pytest
# App
from methylprep.models.controls import ControlType


class TestControlType():

    def test_normalization_green(self):
        assert ControlType.normalization_green() == ('NORM_C', 'NORM_G')

    def test_normalization_red(self):
        assert ControlType.normalization_red() == ('NORM_A', 'NORM_T')
