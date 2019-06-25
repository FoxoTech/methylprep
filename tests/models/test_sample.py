# App
from methpype.models import Channel, Sample


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

        result = sample.get_filepath('txt')
        assert str(result) == '//test_data_dir/test_sentrix_id/test_sentrix_id_test_sentrix_position.txt'

    def test_get_filepath_with_suffix_returns_path_with_suffix(self):
        sample = self._make_test_sample()

        result = sample.get_filepath('txt', 'suffix')
        assert str(result) == '//test_data_dir/test_sentrix_id/test_sentrix_id_test_sentrix_position_suffix.txt'
