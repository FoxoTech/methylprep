from io import StringIO
# App
from methylprep.processing import raw_dataset
from methylprep.files import SampleSheet


# TODO:
# get_raw_datasets(sample_sheet, sample_name=None, from_s3=None)
# tset sample_name as string, as list, as None

class TestRawDataset():

    @staticmethod
    def test_get_raw_dataset_list():
        sample_sheet = SampleSheet('tests/processing/dummy_sheet.csv', '.')
        #datasets = raw_dataset.get_raw_datasets(sample_sheet, sample_name=['Sample_1','Sample_2'], from_s3=None)
        #print(datasets)
