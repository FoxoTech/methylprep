# LIb
import pytest
# App
from methpype.models import ArrayType

''' NOT IN USE ANYWHERE
class TestArray():
    def test_inits_correctly(self):
        arr = Array('test name', ArrayType.ILLUMINA_450K)
        assert arr.name == 'test name'
        assert arr.array_type == ArrayType.ILLUMINA_450K
'''

class TestArrayType():
    def test_enum_has_correct_values(self):
        assert ArrayType.CUSTOM.value == 'custom'
        assert ArrayType.ILLUMINA_27K.value == '27k'
        assert ArrayType.ILLUMINA_450K.value == '450k'
        assert ArrayType.ILLUMINA_EPIC.value == 'epic'
        assert ArrayType.ILLUMINA_EPIC_PLUS.value == 'epic+'


class TestArrayTypeFromProbeCount():
    def test_no_probe_counts_throws(self):
        with pytest.raises(TypeError):
            ArrayType.from_probe_count()

    def test_out_of_range_probe_counts_throws(self, capsys):
        with pytest.raises(ValueError):
            results = ArrayType.from_probe_count(100)
            assert results is None

            captured = capsys.readouterr()
            assert captured.out == 'Unknown array type\n'

    def test_27k_probe_counts_lower_bound(self):
        array_type = ArrayType.from_probe_count(54000)
        assert array_type == ArrayType.ILLUMINA_27K

    def test_27k_probe_counts_upper_bound(self):
        array_type = ArrayType.from_probe_count(56000)
        assert array_type == ArrayType.ILLUMINA_27K

    def test_450k_probe_counts_lower_bound(self):
        array_type = ArrayType.from_probe_count(622000)
        assert array_type == ArrayType.ILLUMINA_450K

    def test_450k_probe_counts_upper_bound(self):
        array_type = ArrayType.from_probe_count(623000)
        assert array_type == ArrayType.ILLUMINA_450K

    def test_epic_probe_counts_lower_bound(self):
        array_type = ArrayType.from_probe_count(1050000)
        assert array_type == ArrayType.ILLUMINA_EPIC

    def test_epic_probe_counts_upper_bound(self):
        array_type = ArrayType.from_probe_count(1053000)
        assert array_type == ArrayType.ILLUMINA_EPIC

    def test_epic_plus(self):
        array_type = ArrayType.from_probe_count(1055583)
        assert array_type == ArrayType.ILLUMINA_EPIC_PLUS
