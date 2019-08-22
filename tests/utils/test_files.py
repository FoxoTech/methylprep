# LIb
from io import StringIO, BytesIO
from pathlib import Path
from unittest.mock import patch
import pytest
# App
from methylprep.utils.files import (
    download_file,
    ensure_directory_exists,
    read_and_reset,
)


@pytest.fixture(scope='function')
def file_obj():
    # setup code
    infile = StringIO()

    def make_infile(write_bytes):
        infile.write(write_bytes)
        infile.seek(0)
        return infile

    yield make_infile

    # teardown code
    infile.close()


class TestDecoratorReadAndReset():
    def test_returns_file_to_same_position(self, file_obj):
        infile = file_obj('test text')
        expected_position = 5
        infile.seek(expected_position)

        assert infile.tell() == expected_position

        def reader(passed_infile):
            return passed_infile.read()

        wrapped = read_and_reset(reader)
        result = wrapped(infile)

        assert result == 'text'
        assert infile.tell() == expected_position


class TestEnsureDirectoryExists():
    def test_noops_if_path_exists(self, mocker):
        mock_exists = mocker.patch('pathlib.Path.exists')
        mock_mkdir = mocker.patch('pathlib.Path.mkdir')
        mock_exists.return_value = True

        ensure_directory_exists('test_path')

        assert mock_mkdir.called is False

    def test_mkdir_if_path_not_exists(self, mocker):
        mock_exists = mocker.patch('pathlib.Path.exists')
        mock_mkdir = mocker.patch('pathlib.Path.mkdir')
        mock_exists.return_value = False

        ensure_directory_exists('test_path')

        assert mock_mkdir.call_count == 1
        assert mock_mkdir.called_with('test_path', parents=True)


@patch('methylprep.utils.files.urlopen')
@patch('methylprep.utils.files.shutil')
class TestDownloadFile():
    mock_filename = 'mock filename.txt'
    mock_src_url = 'src_url'
    tmpdir = None

    @pytest.fixture(autouse=True)
    def transact(self, tmpdir):
        self.tmpdir = Path(tmpdir)

    def test_raises_if_file_exists_overwrite_false(self, *_args):
        with pytest.raises(FileExistsError):
            self.tmpdir.joinpath(self.mock_filename).touch()
            download_file(self.mock_filename, self.mock_src_url, self.tmpdir, overwrite=False)

    def test_doesnt_raise_if_file_exists_overwrite_true(self, *_args):
        self.tmpdir.joinpath(self.mock_filename).touch()
        download_file(self.mock_filename, self.mock_src_url, self.tmpdir, overwrite=True)

    def test_makes_dir_if_dest_dir_not_exists(self, *_args):
        new_dir = self.tmpdir.joinpath('not_exists')
        assert new_dir.exists() is False

        download_file(self.mock_filename, self.mock_src_url, new_dir)
        assert new_dir.exists() is True

    def test_skips_mkdir_if_dest_dir_exists(self, *_args):
        new_dir = self.tmpdir.joinpath('not_exists')
        new_dir.mkdir()
        assert new_dir.exists() is True

        download_file(self.mock_filename, self.mock_src_url, new_dir)
        assert new_dir.exists() is True

    def test_saves_content_to_dest_path(self, mock_shutil, mock_urlopen, *_args):
        mock_shutil.stop()
        expected_content = b'response content'
        response = BytesIO(expected_content)
        mock_urlopen.return_value = response
        expected_filepath = self.tmpdir.joinpath(self.mock_filename)
        assert expected_filepath.exists() is False

        download_file(self.mock_filename, self.mock_src_url, self.tmpdir)
        assert expected_filepath.exists() is True
        assert mock_shutil.copyfileobj.call_count == 1
