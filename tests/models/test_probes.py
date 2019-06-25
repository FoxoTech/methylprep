# LIb
import pytest
# App
from methpype.models import Probe, ProbeType


class TestProbe():
    def test_inits_correctly(self):
        probe = Probe('test_address', 'test_illumina', ProbeType.SNP_ONE)
        assert probe.address == 'test_address'
        assert probe.illumina_id == 'test_illumina'
        assert probe.probe_type is ProbeType.SNP_ONE


class TestProbeType():
    def test_enum_has_correct_values(self):
        assert ProbeType.ONE.value == 'I'
        assert ProbeType.TWO.value == 'II'
        assert ProbeType.CONTROL.value == 'Control'
        assert ProbeType.SNP_ONE.value == 'SnpI'
        assert ProbeType.SNP_TWO.value == 'SnpII'

    def test_str_returns_value(self):
        assert str(ProbeType.ONE) == 'I'
        assert str(ProbeType.TWO) == 'II'
        assert str(ProbeType.CONTROL) == 'Control'
        assert str(ProbeType.SNP_ONE) == 'SnpI'
        assert str(ProbeType.SNP_TWO) == 'SnpII'


class TestProbeTypeFromManifestValues():
    non_snp_name = 'cg1234'
    snp_name = 'rs1234'

    def test_no_args_throws(self):
        with pytest.raises(TypeError):
            ProbeType.from_manifest_values()

        with pytest.raises(TypeError):
            ProbeType.from_manifest_values(self.non_snp_name)

    def test_unknown_non_snp_returns_control(self):
        results = ProbeType.from_manifest_values(self.non_snp_name, 'random')
        assert results is ProbeType.CONTROL

    def test_unknown_is_snp_returns_control(self):
        results = ProbeType.from_manifest_values(self.snp_name, 'random')
        assert results is ProbeType.CONTROL

    def test_type1_non_snp_returns_type1(self):
        results = ProbeType.from_manifest_values(self.non_snp_name, 'I')
        assert results is ProbeType.ONE

    def test_type1_is_snp_returns_type1snp(self):
        results = ProbeType.from_manifest_values(self.snp_name, 'I')
        assert results is ProbeType.SNP_ONE

    def test_type2_non_snp_returns_type2(self):
        results = ProbeType.from_manifest_values(self.non_snp_name, 'II')
        assert results is ProbeType.TWO

    def test_type2_is_snp_returns_type2snp(self):
        results = ProbeType.from_manifest_values(self.snp_name, 'II')
        assert results is ProbeType.SNP_TWO
