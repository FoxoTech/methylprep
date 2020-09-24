# Lib
import logging
import numpy as np


LOGGER = logging.getLogger(__name__)


__all__ = [
    'read_byte',
    'read_char',
    'read_int',
    'read_long',
    'read_results',
    'read_short',
    'read_string',
    'npread',
]


def bytes_to_int(input_bytes, signed=False):
    """Returns the integer represented by the given array of bytes.
    Pre-sets the byteorder to be little-endian.

    Arguments:
        input_bytes -- Holds the array of bytes to convert.  The argument must either
            support the buffer protocol or be an iterable object producing bytes.
            Bytes and bytearray are examples of built-in objects that support the
            buffer protocol.

    Keyword Arguments:
        signed {bool} -- Indicates whether two's complement is used to represent the integer. (default: {False})

    Returns:
        [integer] -- Integer value converted from the supplied bytes.
    """
    return int.from_bytes(input_bytes, byteorder='little', signed=signed)


def read_results(infile, parser, num_elements, allow_early_end=False):
    """DEPRECATED function, replaced by npread() which runs faster.
    Parses a binary file multiple times, allowing for control if the
    file ends prematurely.

    Arguments:
        infile {file-like} -- The binary file to read the select number of bytes.
        parser {[type]} -- Parsing function to apply iteratively over the infile.
        num_bytes {integer} -- The number of elements to parse.

    Keyword Arguments:
        allow_early_end {bool} -- Whether it is ok to reach the end of the file early. (default: {False})

    Raises:
        EOFError: If the end of the file is reached before the number of elements have
            been processed.

    Returns:
        [list(any)] -- A list of the parsed values.
    """
    results = []

    while len(results) < num_elements:
        try:
            results.append(parser(infile))
            continue
        except EOFError:
            if allow_early_end:
                break

            raise EOFError('End of file reached before number of results parsed')

    return results


def read_byte(infile):
    """Converts a single byte to an integer value.

    Arguments:
        infile {file-like} -- The binary file to read the select number of bytes.

    Returns:
        [integer] -- Unsigned integer value converted from the supplied bytes.
    """
    return bytes_to_int(infile.read(1), signed=False)


def read_short(infile):
    """Converts a two-byte element to an integer value.

    Arguments:
        infile {file-like} -- The binary file to read the select number of bytes.

    Returns:
        [integer] -- Unsigned integer value converted from the supplied bytes.
    """
    return bytes_to_int(infile.read(2), signed=False)


def read_int(infile):
    """Converts a four-byte element to an integer value.

    Arguments:
        infile {file-like} -- The binary file to read the select number of bytes.

    Returns:
        [integer] -- Signed integer value converted from the supplied bytes.
    """
    return bytes_to_int(infile.read(4), signed=True)


def read_long(infile):
    """Converts an eight-byte element to an integer value.

    Arguments:
        infile {file-like} -- The binary file to read the select number of bytes.

    Returns:
        [integer] -- Signed integer value converted from the supplied bytes.
    """
    return bytes_to_int(infile.read(8), signed=True)


def read_char(infile, num_bytes):
    """Converts an array of bytes to a string.

    Arguments:
        infile {file-like} -- The binary file to read the select number of bytes.
        num_bytes {integer} -- The number of bytes to read and parse.

    Returns:
        [string] -- UTF-8 decoded string value.
    """
    return infile.read(num_bytes).decode('utf-8')


def read_string(infile):
    """Converts an array of bytes to a string.

    Arguments:
        infile {file-like} -- The binary file to read the select number of bytes.

    Returns:
        [string] -- UTF-8 decoded string value.
    """
    num_bytes = read_byte(infile)
    num_chars = num_bytes % 128
    shift = 0

    while num_bytes // 128 == 1:
        num_bytes = read_byte(infile)
        shift += 7
        offset = (num_bytes % 128) * (2 ** shift)
        num_chars += offset

    return read_char(infile, num_chars)

def npread(file_like, dtype, n):
    """Parses a binary file multiple times, allowing for control if the
    file ends prematurely. This replaces read_results() and runs faster.
    And it provides support for reading gzipped idat files without decompressing.

    Arguments:
        infile {file-like} -- The binary file to read the select number of bytes.
        dtype {data type} -- used within idat files
        n {number of snps read} -- see files/idat.py for how this function is applied.

    Raises:
        EOFError: If the end of the file is reached before the number of elements have
            been processed.

    Returns:
        [list(any)] -- A list of the parsed values.
    """
    dtype=np.dtype(dtype)
    # np.readfile is not able to read from gzopene-d file
    a=file_like.read(dtype.itemsize*n)
    if len(a) != dtype.itemsize*n:
        raise EOFError('End of file reached before number of results parsed')
    r=np.frombuffer(a, dtype, n)
    if r.size != n:
        raise EOFError('End of file reached before number of results parsed')
    return r
