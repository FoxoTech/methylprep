# Lib
import gzip
import logging
from pathlib import Path, PurePath
import shutil
from urllib.request import urlopen
from urllib.error import URLError
import ssl


__all__ = [
    'download_file',
    'ensure_directory_exists',
    'get_file_object',
    'is_file_like',
    'read_and_reset',
    'reset_file',
]


LOGGER = logging.getLogger(__name__)


def read_and_reset(inner):
    """Decorator that resets a file-like object back to the original
    position after the function has been called."""

    def wrapper(infile, *args, **kwargs):
        current_position = infile.tell()
        rval = inner(infile, *args, **kwargs)
        infile.seek(current_position)
        return rval

    return wrapper


def make_path_like(path_like):
    """Attempts to convert a string to a Path instance."""

    if isinstance(path_like, Path):
        return path_like

    try:
        return Path(path_like)
    except TypeError:
        raise TypeError(f'could not convert to Path: {path_like}')


def require_path(inner):
    """Decorator that ensure the argument provided to the inner function
    is a Path instance."""

    def wrapped(orig_path, *args, **kwargs):
        path_like = make_path_like(orig_path)
        return inner(path_like, *args, **kwargs)

    return wrapped


@require_path
def ensure_directory_exists(path_like):
    """Ensures the ancestor directories of the provided path
    exist, making them if they do not."""
    if path_like.exists():
        return

    parent_dir = path_like
    if path_like.suffix:
        parent_dir = path_like.parent

    parent_dir.mkdir(parents=True, exist_ok=True)


def download_file(filename, src_url, dest_dir, overwrite=False):
    """download_file now defaults to non-SSL if SSL fails, with warning to user.
    MacOS doesn't have ceritifi installed by default."""
    dir_path = make_path_like(dest_dir)
    dest_path = dir_path.joinpath(filename)

    if not dest_path.exists():
        ensure_directory_exists(dest_dir)
    elif not overwrite:
        # check if file already exists, and return if it is there.
        LOGGER.info(f'File exists: {dest_path} Set overwrite=True to overwrite the file.')
        #raise FileExistsError(f'File exists: {dest_path}') # -- raising an error here terminates lambda, except that pipeline_s3 catches it.
        return
    try:
        with urlopen(src_url) as response:
            with open(dest_path, 'wb') as out_file:
                shutil.copyfileobj(response, out_file)
    except URLError as e:
        LOGGER.error(e)
        LOGGER.info("If you got [SSL: CERTIFICATE_VERIFY_FAILED] error and you're using MacOS, go to folder /Applications/Python 3.X and run 'Install Certificates.command' to fix this. It cannot download from https.")
        # <urlopen error [SSL: CERTIFICATE_VERIFY_FAILED] certificate verify failed: unable to get local issuer certificate (_ssl.c:1056)>
        try:
            LOGGER.info("retrying without SSL")
            context = ssl._create_unverified_context()
            with urlopen(src_url, context=context) as response:
                with open(dest_path, 'wb') as out_file:
                    shutil.copyfileobj(response, out_file)
        except URLError as e:
            raise URLError(e)

def is_file_like(obj):
    """Check if the object is a file-like object.
    For objects to be considered file-like, they must be an iterator AND have either a
    `read` and/or `write` method as an attribute.
    Note: file-like objects must be iterable, but iterable objects need not be file-like.

    Arguments:
        obj {any} --The object to check.

    Returns:
        [boolean] -- [description]

    Examples:
    --------
    >>> buffer(StringIO("data"))
    >>> is_file_like(buffer)
    True
    >>> is_file_like([1, 2, 3])
    False
    """

    if not (hasattr(obj, 'read') or hasattr(obj, 'write')):
        return False

    if not hasattr(obj, '__iter__'):
        return False

    return True


def get_file_object(filepath_or_buffer):
    """Returns a file-like object based on the provided input.
    If the input argument is a string, it will attempt to open the file
    in 'rb' mode.
    """
    if is_file_like(filepath_or_buffer):
        return filepath_or_buffer

    if PurePath(filepath_or_buffer).suffix == '.gz':
        return gzip.open(filepath_or_buffer, 'rb')

    return open(filepath_or_buffer, 'rb')


def reset_file(filepath_or_buffer):
    """Attempts to return the open file to the beginning if it is seekable."""
    if not hasattr(filepath_or_buffer, 'seek'):
        return

    filepath_or_buffer.seek(0)
