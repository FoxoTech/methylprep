# Lib
import pytest
# App
from methylprep.files.sample_sheets import find_sample_sheet


class TestGetSampleSheet():

    def test_raises_with_missing_dir(self, tmp_path):
        with pytest.raises(FileNotFoundError):
            missing_path = tmp_path.joinpath('fake_path')
            find_sample_sheet(missing_path)

    def test_raises_with_valid_path_not_dir(self, tmp_path):
        with pytest.raises(FileNotFoundError):
            missing_path = tmp_path.joinpath('fake_path.txt')
            missing_path.touch()
            find_sample_sheet(missing_path)

    def test_raises_if_no_csv_files(self, mocker, tmp_path):
        with pytest.raises(FileNotFoundError):
            sample_sheet_path = tmp_path.joinpath('sample_sheet.txt')
            sample_sheet_path.touch()

            mock_is_sample_sheet = mocker.patch('methylprep.files.SampleSheet.is_sample_sheet')
            mock_is_sample_sheet.return_value = True

            find_sample_sheet(tmp_path)

    def test_raises_if_no_valid_csv_files(self, mocker, tmp_path):
        with pytest.raises(FileNotFoundError):
            sample_sheet_path = tmp_path.joinpath('sample_sheet.csv')
            sample_sheet_path.touch()

            mock_is_sample_sheet = mocker.patch('methylprep.files.SampleSheet.is_sample_sheet')
            mock_is_sample_sheet.return_value = False

            find_sample_sheet(tmp_path)

    def test_raises_if_too_many_valid_files(self, mocker, tmp_path):
        with pytest.raises(Exception):
            path_1 = tmp_path.joinpath('sample_sheet1.csv')
            path_1.touch()

            path_2 = tmp_path.joinpath('sample_sheet2.csv')
            path_2.touch()

            mock_is_sample_sheet = mocker.patch('methylprep.files.SampleSheet.is_sample_sheet')
            mock_is_sample_sheet.return_value = True

            find_sample_sheet(tmp_path)

    def test_returns_samplesheet_with_validpath(self, mocker, tmp_path):
        sample_sheet_path = tmp_path.joinpath('sample_sheet.csv')
        sample_sheet_path.touch()

        mock_is_sample_sheet = mocker.patch('methylprep.files.SampleSheet.is_sample_sheet')
        mock_is_sample_sheet.return_value = True

        result = find_sample_sheet(tmp_path)

        assert result == sample_sheet_path
