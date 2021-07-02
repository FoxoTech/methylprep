# App
from methylprep.models import Channel, Sample, ArrayType, SigSet # MethylationDataset, RawDataset
from methylprep.files import SampleSheet, Manifest, IdatDataset
from pathlib import Path

class TestSigSet():

    @staticmethod
    def test_idat_dataset_control_snp():
        data_dir = 'docs/example_data/GSE69852'
        sample = Sample(data_dir, '9247377093', 'R02C01')
        red_idat = IdatDataset(Path(data_dir, '9247377093_R02C01_Red.idat'), Channel.RED)
        green_idat = IdatDataset(Path(data_dir,'9247377093_R02C01_Grn.idat'), Channel.GREEN)
        manifest = Manifest(ArrayType.ILLUMINA_450K) # this could be REALLY slow, if docker/test env has to download it.
        sigset = SigSet(sample, green_idat, red_idat, manifest)
        if sigset.ctrl_green.loc[27630314].mean_value != 299.0:
            raise AssertionError("fcontrol green mean intensity mismatch")
        if sigset.methylated.loc['rs1019916'].Meth != 7469.0:
            raise AssertionError("SNP meth failed")
        if sigset.methylated.shape[0] != 485577:
            raise AssertionError("meth failed: unexpected number of SNP probes found")
        if sigset.unmethylated.shape[0] != 485577:
            raise AssertionError("unmeth failed: unexpected number of SNP probes found")
        if sigset.II.shape[0] != 350076:
            raise AssertionError("SigSet.II failed: unexpected number of SNP probes found")
        if sigset.IG.shape[0] != 46298:
            raise AssertionError("SigSet.IG failed: unexpected number of SNP probes found")
        if sigset.IR.shape[0] != 89203:
            raise AssertionError("SigSet.IR failed: unexpected number of SNP probes found")
        if sigset.snp_methylated.shape[0] != 65:
            raise AssertionError("SNP meth failed: unexpected number of SNP probes found")
        if sigset.snp_unmethylated.shape[0] != 65:
            raise AssertionError("SNP unmeth failed: unexpected number of SNP probes found")
        if sigset.unmethylated.loc['rs1019916'].Unmeth != 5291.0:
            raise AssertionError("SNP unmeth failed")
        if sigset.snp_unmethylated.loc['rs1019916'].Unmeth != 5291.0:
            raise AssertionError("SNP unmeth failed")

    @staticmethod
    def test_sigset_mouse():
        data_dir = 'docs/example_data/mouse'
        mouse_file = '204879580038_R06C02'
        sample = Sample(data_dir, '204879580038', 'R06C02')
        red_idat = IdatDataset(Path(data_dir, f'{mouse_file}_Red.idat'), Channel.RED)
        green_idat = IdatDataset(Path(data_dir,f'{mouse_file}_Grn.idat'), Channel.GREEN)
        manifest = Manifest(ArrayType.ILLUMINA_MOUSE) # this could be REALLY slow, if docker/test env has to download it.
        sigset = SigSet(sample, green_idat, red_idat, manifest)
        #raw_dataset = RawDataset(sample, green_idat, red_idat) --- class still works, but deprecated; not used in pipeline anywhere.
        #test_raw_dataset = RawDataset.from_sample(sample)
        #fg_controls = raw_dataset.get_fg_controls(manifest, Channel.GREEN)
        #if fg_controls.loc[27630314].mean_value != 299.0:
        if sigset.ctrl_green.loc[27630314].mean_value != 206.0:
            raise AssertionError("control green failed")
        #meth_data = MethylationDataset.snp_methylated(raw_dataset, manifest)
        if sigset.methylated.loc['rs108256820_TC21']['Meth'] != 1889.0:
            raise AssertionError("methylated failed")
        if sigset.snp_methylated.loc['rs108256820_TC21']['Meth'] != 1889.0:
            raise AssertionError("SNP meth failed")
        if sigset.snp_methylated.shape[0] != 1485: # sesame has 1486
            raise AssertionError(f"SNP meth failed: expected 1485 SNP probes, found {sigset.snp_methylated.shape[0]}")
        # added 615 probes into MM285_mm39_v2
        if sigset.methylated.shape[0] != 293194: #292580, 292582 before dropping duplicates in v1
            raise AssertionError(f"meth failed: expected 293194 probes, found {sigset.methylated.shape[0]}")
        if sigset.unmethylated.shape[0] != 293194: #292580, 292583 before dropping duplicates in v1
            raise AssertionError(f"meth failed: expected 293194 probes, found {sigset.unmethylated.shape[0]}")
        if sigset.II.shape[0] != 228271: #227968:
            raise AssertionError(f"SigSet II: expected 228271 probes, found {sigset.II.shape[0]}")
        if sigset.IG.shape[0] != 17697: #17627:
            raise AssertionError(f"SigSet IG: expected 17697 probes, found {sigset.IG.shape[0]}")
        if sigset.IR.shape[0] != 47231: #46990:
            raise AssertionError(f"SigSet IR: expected 47231 probes, found {sigset.IR.shape[0]}")
