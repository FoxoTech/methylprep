# LIb
from io import BytesIO
from unittest.mock import patch
import pytest
# App
from methpype.utils.parsing import (
    bytes_to_int,
    read_byte,
    read_char,
    read_int,
    read_long,
    read_results,
    read_short,
    read_string,
)


@pytest.fixture(scope='function')
def file_obj():
    # setup code
    infile = BytesIO()

    def make_infile(write_bytes):
        infile.write(write_bytes)
        infile.seek(0)
        return infile

    yield make_infile

    # teardown code
    infile.close()


class TestConvertToInt():
    @patch('builtins.int')
    def test_bytes_to_int_call_builtin_int(self, mock_int):
        byte_val = bytes([17])
        bytes_to_int(byte_val)
        assert mock_int.from_bytes.call_count == 1

    @patch('builtins.int')
    def test_bytes_to_int_uses_default_values(self, mock_int):
        byte_val = bytes([17])
        bytes_to_int(byte_val)
        assert mock_int.from_bytes.call_count == 1
        assert mock_int.from_bytes.called_with(byte_val, signed=False)

    @patch('builtins.int')
    def test_bytes_to_int_uses_provided_values(self, mock_int):
        byte_val = bytes([17])
        bytes_to_int(byte_val, signed=True)
        assert mock_int.from_bytes.call_count == 1
        assert mock_int.from_bytes.called_with(byte_val, signed=True)


class TestReadBinary():
    def test_read_results(self, file_obj):
        infile = file_obj(b'\xce\x00\xce\x00\x90\xbf\xa1\x00')
        result = read_results(infile, read_int, 2)
        assert result == [13500622, 10600336]


class TestReadByte():
    def test_read_byte(self, file_obj):
        infile = file_obj(b'\xff')
        result = read_byte(infile)
        assert result == 255


class TestReadShort():
    def test_read_short(self, file_obj):
        infile = file_obj(b'\xce\x00')
        result = read_short(infile)
        assert result == 206


class TestReadInt():
    def test_read_int(self, file_obj):
        infile = file_obj(b'\xce\x00\xce\x00')
        result = read_int(infile)
        assert result == 13500622


class TestReadLong():
    def test_read_long(self, file_obj):
        infile = file_obj(b'\xce\x00\xce\x00')
        result = read_long(infile)
        assert result == 13500622


class TestReadChar():
    def test_read_char(self, file_obj):
        infile = file_obj(b'test file extra')
        result = read_char(infile, 9)
        assert result == 'test file'


class TestReadString():
    def test_read_string(self, file_obj):
        test_text = b'test text for read_string'
        num_chars = len(test_text)
        encoded_num_chars = bytes([num_chars])
        test_input = encoded_num_chars + test_text
        infile = file_obj(test_input)
        result = read_string(infile)
        assert result == 'test text for read_string'
