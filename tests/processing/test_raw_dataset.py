from io import StringIO
# App
from methylprep.models import Channel, Sample, ArrayType, MethylationDataset
from methylprep.processing import raw_dataset, RawDataset
from methylprep.files import SampleSheet, Manifest, IdatDataset
from pathlib import Path

# TODO:
# get_raw_datasets(sample_sheet, sample_name=None, from_s3=None)
# tset sample_name as string, as list, as None

class TestRawDataset():

    @staticmethod
    def test_get_raw_dataset_list():
        sample_sheet = SampleSheet('tests/processing/dummy_sheet.csv', '.')
        # too slow, and done elsewhere already as part of test run_pipeline
        #datasets = raw_dataset.get_raw_datasets(sample_sheet, sample_name=['Sample_1','Sample_2'], from_s3=None)
        #print(datasets)

    @staticmethod
    def test_methylation_dataset_control_snp():
        data_dir = 'docs/example_data/GSE69852'
        sample = Sample(data_dir, '9247377093', 'R02C01')
        red_idat = IdatDataset(Path(data_dir, '9247377093_R02C01_Red.idat'), Channel.RED)
        green_idat = IdatDataset(Path(data_dir,'9247377093_R02C01_Grn.idat'), Channel.GREEN)
        raw_dataset = RawDataset(sample, green_idat, red_idat)
        manifest = Manifest(ArrayType.ILLUMINA_450K) # this could be REALLY slow, if docker/test env has to download it.
        fg_controls = raw_dataset.get_fg_controls(manifest, Channel.GREEN)
        if fg_controls.loc[27630314].mean_value != 299.0:
            raise AssertionError("fg_controls failed")
        meth_data = MethylationDataset.snp_methylated(raw_dataset, manifest)
        if meth_data.data_frame.loc['rs1019916'].mean_value != 7469.0:
            raise AssertionError("SNP meth failed")
        if meth_data.data_frame.shape[0] != 65:
            raise AssertionError("SNP meth failed: unexpected number of SNP probes found")
        unmeth_data = MethylationDataset.snp_unmethylated(raw_dataset, manifest)
        if unmeth_data.data_frame.loc['rs1019916'].mean_value != 5291.0:
            raise AssertionError("SNP unmeth failed")
        if unmeth_data.data_frame.shape[0] != 65:
            raise AssertionError("SNP unmeth failed: unexpected number of SNP probes found")
