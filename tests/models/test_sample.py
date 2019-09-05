# App
from methylprep.models import Channel, Sample


class TestSample():

    def _make_test_sample(self):
        return Sample(
            data_dir='//test_data_dir',
            sentrix_id='test_sentrix_id',
            sentrix_position='test_sentrix_position',
            Sample_Name='test_name',
        )

    def test_inits_correctly(self):
        sample = self._make_test_sample()

        assert sample.data_dir == '//test_data_dir'
        assert sample.name == 'test_name'
        assert sample.sentrix_id == 'test_sentrix_id'
        assert sample.sentrix_position == 'test_sentrix_position'

    def test_str_returns_sentrix_id_position(self):
        sample = self._make_test_sample()

        assert str(sample) == 'test_sentrix_id_test_sentrix_position'

    def test_base_filename_returns_sentrix_id_position(self):
        sample = self._make_test_sample()

        assert sample.base_filename == 'test_sentrix_id_test_sentrix_position'

    def test_get_filepath_no_suffix_returns_path_no_suffix(self):
        sample = self._make_test_sample()
        # use verify=False for when the file doesn't actually exist.
        result = sample.get_filepath('txt', verify=False)
        assert str(result) == '//test_data_dir/test_sentrix_id/test_sentrix_id_test_sentrix_position.txt'

    def test_get_filepath_with_suffix_returns_path_with_suffix(self):
        sample = self._make_test_sample()
        # use verify=False for when the file doesn't actually exist.
        result = sample.get_filepath('txt', 'suffix', verify=False)
        assert str(result) == '//test_data_dir/test_sentrix_id/test_sentrix_id_test_sentrix_position_suffix.txt'

    @staticmethod
    def _make_real_test_sample():
        return Sample(
            data_dir='docs/example_data/GSE100850',
            sentrix_id='200526210010',
            sentrix_position='R01C01',
            Sample_Name='200526210010_R01C01_Grn.idat',
        )

    def test_get_filepath_and_file_must_exist(self):
        realsample = self._make_real_test_sample()
        # use verify=True for when the file must exist.
        result = realsample.get_filepath('idat', 'Grn', verify=True)
        assert str(result) == 'docs/example_data/GSE100850/200526210010_R01C01_Grn.idat'
