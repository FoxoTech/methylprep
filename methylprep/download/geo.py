# Lib
import logging
from ftplib import FTP
import ftplib
import socket
from pathlib import Path, PurePath
import os
from tarfile import ReadError
import re
import zipfile
import gzip
import tarfile
import zlib
import shutil
from bs4 import BeautifulSoup
import dateutil # python-dateutil, non-built-in
import datetime
import pickle
from urllib.request import urlopen
import pandas as pd
from tqdm import tqdm
from collections import Counter
import sys

# unique to find_betas_any_source()...
try:
    import methylcheck # not required generally, but needed for one of these functions; handled within.
    # geo.py uses .load, .read_geo, and .read_geo_processed
except ImportError:
    pass

try:
    from pandas.errors import ParserError
except ImportError:
    try:
        from pandas.io.parsers import ParserError
    except Exception:
        class ParserError(ValueError):
            """
            Exception that is raised by an error encountered in parsing file contents.

            This is a generic error raised for errors encountered when functions like
            `read_csv` or `read_html` are parsing contents of a file.

            See Also
            --------
            read_csv : Read CSV (comma-separated) file into a DataFrame.
            read_html : Read HTML table into a DataFrame."""
            pass
import io
import json
import time
import tempfile
import requests
import zipfile
import random
import subprocess

# app
from .miniml import sample_sheet_from_miniml, sample_sheet_from_idats, convert_miniml
# cannot relative-import here because process_data uses geo.py.
#from .process_data import confirm_dataset_contains_idats, get_attachment_info, run_series
import methylprep.download.process_data

__all__ = [
    'geo_download',
    'pipeline_find_betas_any_source',
    'geo_metadata',
    'download_geo_processed',
    'search',
]


#logging.basicConfig(level=logging.DEBUG) # always verbose
LOGGER = logging.getLogger(__name__)
LOGGER.setLevel( logging.INFO )

def geo_download(geo_id, series_path, geo_platforms, clean=True, decompress=True, meta_only=False):
    """Downloads the IDATs and metadata for a GEO series

    Arguments:
        geo_id [required]
            the GEO Accension for the desired series (e.g. GSE134293)
        series_path [required]
            the directory to download the data to
        geo_platforms [required]
            the list of supported GEO platforms
        clean
            whether or not to delete files once they are no longer need (True by default)

    Note about GEO IDs:
        You can use the NIH online search to find data sets, then click "Send to:" at the button of a results page,
        and export a list of unique IDs as text file. These IDs are not GEO_IDs used here. First, remove the first
        three digits from the number, so Series ID: 200134293 is GEO accension ID: 134293, then include the GSE part,
        like "GSE134293" in your CLI parameters.

    meta_only: set to True if you only want it to download the meta data, not the RAW data file.

    This function returns True or False, depending on whether the downloaded data is correct."""
    success = True
    series_dir = Path(series_path)
    raw_filename = f"{geo_id}_RAW.tar"
    miniml_packaged_file = f"{geo_id}_family.xml"

    if not os.path.exists(series_path):
        raise FileNotFoundError(f'{geo_id} directory not found.')

    for platform in geo_platforms:
        if not Path(f"{series_path}/{platform}").exists():
            Path(f"{series_path}/{platform}").mkdir(parents=True, exist_ok=True)

    ftp = FTP('ftp.ncbi.nlm.nih.gov', timeout=120) # 2 mins
    ftp.login()
    ftp.cwd(f"geo/series/{geo_id[:-3]}nnn/{geo_id}")

    # some miniml files are stored in multiple parts; only need the last part
    try:
        miniml_files = list(ftp.mlsd('miniml'))
        miniml_files = {k:int(v['size']) for k,v in miniml_files if int(v['size']) > 0}
        # if there are multiple matching files, then the LAST file will always have the xml part we need.
        # the rest are table-tbl files we don't need.
        if len(miniml_files) > 1:
            pattern = f"{geo_id}_family\.xml(.*)(\.tgz)?"
            miniml_matches = {miniml_filename: re.match(pattern, miniml_filename).groups()[0] for miniml_filename in miniml_files}
            for k,v in miniml_matches.items():
                if re.search('(\d)', v):
                    miniml_matches[k] = int(re.search('(\d)', v).groups()[0])
                else:
                    miniml_matches[k] = 0
            # sort so highest filename is first, then return first in list
            miniml_filename, miniml_part = sorted(miniml_matches.items(), key=lambda i:i[1], reverse=True)[0]
            filesize = miniml_files[miniml_filename]
        else:
            miniml_filename = f"{geo_id}_family.xml.tgz"
            ftp.voidcmd('TYPE I') # reset to binary
            filesize = ftp.size(f"miniml/{miniml_filename}")
    except Exception as e:
        LOGGER.error(f"miniml multipart check ERROR: {e}")
        try: # fallback will always get a file, but file may not contain the XML we need
            miniml_filename = f"{geo_id}_family.xml.tgz"
            filesize = ftp.size(f"miniml/{miniml_filename}") # -- gives 550 error because CWD puts it in ASCII mode.
        except Exception as e:
            LOGGER.error(f"ftp.size ERROR: {e}")
            ftp.voidcmd('TYPE I') # from https://stackoverflow.com/a/22093848/536538
            filesize = ftp.size(f"miniml/{miniml_filename}") # -- gives 550 error because CWD puts it in ASCII mode.

    def is_valid_tgz(filepath=f"{series_path}/{miniml_filename}"):
        if Path(filepath).exists():
            try:
                min_tar = tarfile.open(filepath)
                for file in min_tar.getnames():
                    pass
                return True
            except gzip.BadGzipFile:
                LOGGER.info(f"Found a family series file, but appears to be corrupted, so re-downloading.")
                Path(filepath).unlink()
                return False
            except zlib.error:
                LOGGER.info(f"Found a family series file, but appears to be corrupted (invalid stored block lengths), so re-downloading.")
                Path(filepath).unlink()
                return False
        else:
            return False # file does not exist yet, so download

    if not Path(f"{series_path}/{miniml_filename}").exists():
        # test opening it; may be corrupted
        attempt = 0
        while is_valid_tgz(f"{series_path}/{miniml_filename}") is False:
            LOGGER.info(f"Downloading {miniml_filename}; expecting {int(filesize/1048576)} MB") # ({current} of {len(miniml_files)})
            miniml_file = open(f"{series_path}/{miniml_filename}", 'wb')
            try:
                #for filename,filestats in ftp.mlsd(path="miniml", facts=["size"]):
                #    if filename == miniml_filename:
                #        filesize = filestats['size']
                #        break
                with tqdm(unit = 'b', unit_scale = True, unit_divisor=1024, leave = False, miniters = 1, desc = geo_id, total = filesize) as tqdm_instance:
                    def tqdm_callback(data):
                        tqdm_instance.update(len(data))
                        miniml_file.write(data)

                        # check to see if download has completed
                        if (miniml_file.tell() == filesize):
                            # close down the progress bar
                            # prevents an extra ftp call which generates a timed out exception
                            tqdm_instance.close()

                    ftp.retrbinary(f"RETR miniml/{miniml_filename}", tqdm_callback, 8192)
            except Exception as e:
                LOGGER.error(e)
                LOGGER.warning('Failed to create a progress bar, but retrying. It should be downloading...')
                ftp.retrbinary(f"RETR miniml/{miniml_filename}", miniml_file.write, 8192)
            miniml_file.close()
            attempt += 1
            if attempt > 1:
                LOGGER.info(f"[Attempt #{attempt}]")
            if attempt == 5:
                LOGGER.error(f"Could not download after {attempt} attempts. Aborting.")
                break

        min_tar = tarfile.open(f"{series_path}/{miniml_filename}")
        file_found = False
        for file in min_tar.getnames():
            if file == miniml_packaged_file:
                min_tar.extract(file, path=series_path)
                file_found = True
        if not file_found:
            LOGGER.warning(f"Did not find {miniml_packaged_file} within {series_path}/{miniml_filename}")
        min_tar.close()
        if clean:
            Path(f"{series_path}/{miniml_filename}").unlink()
    ftp.quit()
    if meta_only is True:
        return success

    if list(series_dir.glob('*.idat.gz')) == [] and list(series_dir.glob('**/*.idat')) == []:
        if not Path(f"{series_path}/{raw_filename}").exists():
            ftp = FTP('ftp.ncbi.nlm.nih.gov',
                      timeout=59)  # see issue https://bugs.python.org/issue30956 (must be <60s because of a bug)
            ftp.login()
            ftp.cwd(f"geo/series/{geo_id[:-3]}nnn/{geo_id}")
            raw_file = open(f"{series_path}/{raw_filename}", 'wb')
            filesize = ftp.size(f"suppl/{raw_filename}")

            import requests
            nih_root = 'https://www.ncbi.nlm.nih.gov/geo/download/'
            url = f'{nih_root}?acc={geo_id}&format=file'
            with requests.get(url, stream=True) as ss:
                ss.raise_for_status()
                filesize = int(ss.headers['Content-Length'])
                with tqdm(unit = 'b', unit_scale = True, unit_divisor=1024, leave=False, miniters=1, desc=raw_file.name, total=filesize) as tqdm_instance:
                    for chunk in ss.raw.stream(8192, decode_content=False):
                        if chunk:
                            raw_file.write(chunk)
                            tqdm_instance.update(len(chunk))
            """
            try:
                with tqdm(unit = 'b', unit_scale = True, unit_divisor=1024, leave = False, miniters = 1, desc = geo_id, total = filesize) as tqdm_instance:
                    def tqdm_callback(data):
                        # update progress bar and write out data
                        tqdm_instance.update(len(data))
                        raw_file.write(data)

                        # check to see if download has completed
                        if (raw_file.tell() == filesize):
                            # close down the progress bar
                            # prevents an extra ftp call which generates a timed out exception
                            tqdm_instance.close()

                    # make ftp data request
                    ftp.retrbinary(f"RETR suppl/{raw_filename}", tqdm_callback, 8192)

                # quit the current ftp session
                ftp.quit()

            except ftplib.all_errors as e:
                LOGGER.warning(f"FTP error: {e}")
            except Exception as e:
                LOGGER.warning(f"FTP exception: {e}")
            """
            # check the downloaded file size
            rawsize = raw_file.tell()
            if rawsize != filesize:
                LOGGER.info(f'geo_download sizes: FTP:{filesize} Raw:{rawsize} Diff:{filesize-rawsize}')
                LOGGER.info(f"geo_download: File {raw_filename} download failed.")
                return False
            raw_file.close()

        LOGGER.info(f"Unpacking {raw_filename}")
        try:
            tar = tarfile.open(f"{series_path}/{raw_filename}")
            # let user know if this lack idats
            if not any([(True if '.idat' in member.name else False) for member in list(tar.getmembers())]):
                file_endings = Counter([tuple(PurePath(member.name).suffixes) for member in list(tar.getmembers())])
                file_endings = [(k,v) for k,v in file_endings.most_common() if v > 1]
                LOGGER.warning(f'No idat files found in {raw_filename}. {len(list(tar.getmembers()))} files found: {file_endings}.')
                success = False
            for member in tar.getmembers():
                if re.match('.*.idat.gz', member.name):
                    tar.extract(member, path=series_path)
        except ReadError as e:
            raise ReadError(f"There appears to be an incomplete download of {geo_id}. Please delete those files and run this again.")
            success = False
        tar.close()
        if clean:
            os.remove(f"{series_path}/{raw_filename}")

    if decompress:
        for gz in series_dir.glob("*.idat.gz"):
            gz_string = str(gz)
            with gzip.open(gz_string, 'rb') as f_in:
                with open(gz_string[:-3], 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
            if clean:
                gz.unlink()
        LOGGER.info(f"Downloaded and unpacked {geo_id} idats")
    elif not decompress:
        LOGGER.info(f"Downloaded {geo_id} idats without decompressing")
    return success


def geo_metadata(geo_id, series_path, geo_platforms, path):
    """Reads the metadata for the given series (MINiML file) and creates a metadata dictionary and sample sheet for each platform

    Arguments:
        geo_id [required]
            the GEO Accension for the desired series
        series_path [required]
            the directory containing the series data
        geo_platforms [required]
            the list of supported GEO platforms
        path
            the path to the directory containing dictionaries for each platform

    Returns:
        A list of platforms that the series contains samples of"""
    pipeline_kwargs = {}
    miniml_filename = f"{geo_id}_family.xml"
    with open(f"{series_path}/{miniml_filename}", 'r') as fp:
        soup = BeautifulSoup(fp, "xml")

    meta_dicts = {}
    samples_dict = {}

    for platform in geo_platforms:
        meta_dicts[platform] = {}
        samples_dict[platform] = {}

    samples = soup.MINiML.find_all("Sample")
    for sample in samples:
        platform = sample.find('Platform-Ref')['ref']
        accession = sample.Accession.text
        title = sample.Title.text
        attributes_dir = {}
        attributes_dir['source'] = sample.Source.text
        attributes_dir['platform'] = platform
        attributes_dir['title'] = title
        for char in sample.find_all('Characteristics'):
            attributes_dir[char['tag']] = char.text.strip()
        if sample.find('Description'):
            attributes_dir['description'] = sample.find('Description').text.strip()

        # only some MINiML files have this.
        try:
            split_idat = sample.find('Supplementary-Data').text.split("/")[-1].split("_")
            attributes_dir['Sample_ID'] = f"{split_idat[1]}_{split_idat[2]}" # matches beta_value column names
            attributes_dir['Sentrix_ID'] = f"{split_idat[1]}"
            attributes_dir['Sentrix_Position'] = f"{split_idat[2]}"
        except:
            LOGGER.info( "MINiML file does not provide (Sentrix_ID_R00C00)" )

        if platform in geo_platforms:
            for idat in sample.find_all('Supplementary-Data'):
                if idat['type'] == 'IDAT':
                    file_name = (idat.text.split("/")[-1]).strip()[:-3]
                    if (not(Path(f"{series_path}/{file_name}").is_file())) and Path(f"{series_path}/{file_name}.gz").is_file():
                        file_name = file_name+".gz"
                    try:
                        shutil.move(f"{series_path}/{file_name}", f"{series_path}/{platform}/{file_name}")
                    except FileNotFoundError:
                        # this doesn't throw an error if file is already in the right folder
                        if not Path(f"{series_path}/{platform}/{file_name}").is_file():
                            raise FileNotFoundError (f"Could not move file {series_path}/{file_name} after downloading.")

            meta_dicts[platform][accession] = attributes_dir
            samples_dict[platform][accession] = title
        else:
            # this ought to keep other idat files from being included in the package.
            LOGGER.warning(f'Sample: {title[:40]} has unrecognized platform: {platform}; not moving data file')
    LOGGER.info(f"Found {len(attributes_dir)} tags for {len(samples)} samples: {attributes_dir}")
    fp.close()

    seen_platforms = []

    for platform in geo_platforms:
        if meta_dicts[platform]:
            meta_dict_filename = f"{geo_id}_{platform}_dict.pkl"
            pickle.dump(meta_dicts[platform], open(f"{series_path}/{meta_dict_filename}", 'wb'))
            if not os.path.exists(f"{path}/{platform}_dictionaries/{geo_id}_dict.pkl"):
                shutil.move(f"{series_path}/{meta_dict_filename}", f"{path}/{platform}_dictionaries/{geo_id}_dict.pkl")
            try:
                sample_sheet_from_miniml(geo_id, series_path, platform, samples_dict[platform], meta_dicts[platform], save_df=True)
            except ValueError:
                # in this case, the samplesheet meta data wasn't consistent, so using idat filenames instead
                try:
                    sample_sheet_from_idats(geo_id, series_path, platform, samples_dict[platform])
                except ValueError as e:
                    LOGGER.info(f"{e}; will try run_pipeline's create samplesheet method.")
                    pipeline_kwargs['make_sample_sheet'] = True
            if platform not in seen_platforms:
                seen_platforms.append(platform)

    return seen_platforms, pipeline_kwargs


def pipeline_find_betas_any_source(**kwargs):
    """beta_bake: Sets up a script to run methylprep that saves directly to path or S3.
The slowest part of processing GEO datasets is downloading, so this handles that.

STEPS
    - uses `methylprep alert -k <keywords>` to curate a list of GEO IDs worth grabbing.
        note that version 1 will only process idats.
        also runs methylcheck.load on processed files, if installed.
    - downloads a zipfile, uncompresses it,
    - creates a samplesheet,
    - moves it into foxo-test-pipeline-raw for processing.
    - You get back a zipfile with all the output data.

required kwargs:
    - project_name: string, like GSE123456, to specify one GEO data set to download.
        to initialize, specify one GEO id as an input when starting the function.
        - beforehand, you can use `methylprep alert` to verify the data exists.
        - OR you can pass in a string of GEO_ID separated by commas without any spaces and it will split them.
optional kwargs:
    - function: 'geo' (optional, ignored; used to specify this pipeline to run from command line)
    - data_dir:
        - default is current working directory ('.') if omitted
        - use to specify where all files will be downloaded, processed, and finally stored, unless `--cleanup=False`.
        - if using AWS S3 settings below, this will be ignored.
    - verbose: False, default is minimal logging messages.
    - save_source: if True, it will retain .idat and/or -tbl-1.txt files used to generate beta_values dataframe pkl files.
    - compress: if True, it will package everything together in a {geo_id}.zip file, or use gzip if files are too big for zip.
        - default is False
    - clean: If True, removes files from folder, except the compressed output zip file. (Requires compress to the True too)

It will use local disk by default, but if you want it to run in AWS batch + efs provide these:
    - efs (AWS elastic file system name, for lambda or AWS batch processing)
    - clean: default True. If False, does not explicitly remove the tempfolder files at end, or move files into data_dir output filepath/folder.
        - if you need to keep folders in working/efs folder instead of moving them to the data_dir.
        - use cleanup=False when embedding this in an AWS/batch/S3 context,
        then use the `working tempfolder` path and filenames returned to copy these files into S3.

returns:
    - if a single GEO_ID, returns a dict with "error", "filenames", and "tempdir" keys.
    - if mulitple GEO_IDs, returns a dict with "error", "geo_ids" (nested dict), and "tempdir" keys. Uses same tempdir for everything, so clean should be set to True.
    - "error" will be None if it worked okay.
    - "filenames" will be a list of filenames that were created as outputs (type=string)
    - "tempdir" will be the python tempfile tempory-directory object. Passing this out prevents
        garbage collector from removing it when the function ends, so you can retrive these files and
        run tempdir.cleanup() manually. Otherwise, python will remove the tempdir for you when python closes,
        so copy whatever you want out of it first. This makes it possible to use this function with AWS EFS (elastic file systems)
        as part of a lambda or aws-batch function where disk space is more limited.

NOTE: v1.3.0 does NOT support multiple GEO IDs yet.
    """
    if 'methylcheck' not in sys.modules:
        assert ImportError("You cannot run this without `methylcheck` installed.")
    if not kwargs.get('project_name'):
        return {"filenames": [], "tempdir": None, "error": "`project_name` is required (to specify one GEO_ID or a list of them as a comma separated string without spaces)"}
    import zipfile # FOR SOME REASON, importing zipfile at top of file doesn't work in this function :( -- prob because I imported this function without loading the whole file. Or I reassigned var 'zipfile' by accident somewhere.
    BATCH_SIZE=33 # for processing idats
    geo_ids = kwargs['project_name'].split(',') # always a list.
    if not kwargs.get('compress'):
        kwargs['compress'] = False
    if not kwargs.get('move'):
        # CLI always sets move to True
        kwargs['move'] = True # DEBUG for EFS where it might be better to NOT move out of workdir at end
    if kwargs.get('verbose') == False:
        LOGGER.setLevel( logging.WARNING )
    LOGGER.info(f"DEBUG: find_betas: Your command line inputs: {kwargs}")

    #1: set working folder, and make sure /mnt/efs exists if used.
    if kwargs.get('efs') and kwargs['efs'] is not None:
        EFS = kwargs['efs']
        mounted = False
        try:
            #efs_files = list(Path(EFS).glob('*'))
            mounted = Path(EFS).exists()
            LOGGER.info(f"{EFS} exists: {mounted}")
        except Exception as e:
            LOGGER.error(f"{EFS} error: {e}")
        if not mounted and os.name != 'nt':
            result = subprocess.run(["df", "-h"], stdout=subprocess.PIPE)
            LOGGER.warning(f"EFS mount [df -h]: {result.stdout.decode('utf-8')}")
            raise FileNotFoundError("This batch function has no {EFS} mounted. Cannot proceed.")
        working = tempfile.TemporaryDirectory(dir=EFS)
    elif not kwargs.get('efs') and kwargs.get('data_dir') != None:
        if not Path(kwargs.get('data_dir')).exists():
            Path(kwargs.get('data_dir')).mkdir(parents=True, exist_ok=True)
        working = tempfile.TemporaryDirectory(dir=kwargs['data_dir'])
        EFS = working.name
    elif not kwargs.get('efs') and kwargs.get('data_dir') == None:
        working = tempfile.TemporaryDirectory(dir='.')
        EFS = working.name
    if kwargs.get('data_dir') == None:
        kwargs['data_dir'] = '.' # CLI seems to pass in None and the get() doesn't get it.

    LOGGER.info(f"Starting batch GEO pipeline processor for {geo_ids}.")
    #1: run EACH geo series through downloader/processor
    each_geo_result = {}
    for geo_id in geo_ids:
        # this download processed data, or idats, or -tbl-1.txt files inside the _family.xml.tgz files.
        # also downloads meta_data
        # download_geo_processed() downloads nothing if idats exist.
        result = download_geo_processed(geo_id, working, verbose=kwargs.get('verbose', False), use_headers=True)
        LOGGER.info(f"DEBUG: download_geo_processed result: {result}")
        extracted_files = [_file.name for _file in list(Path(working.name).rglob('*')) if not _file.is_dir()]
        zipfile_paths = [_file for _file in list(Path(working.name).rglob('*')) if not _file.is_dir()]
        # LOGGER.info(f"DEBUG: zipfile_paths 465: {zipfile_paths}")

        #2: CHECK: Are all files in working.dir or a sub-folder? Must return path to each file off of workdir.
        zipfile_names = []
        for zipfile_path in zipfile_paths:
            #LOGGER.info(f"DEBUG: checking {zipfile_path.name}")
            working_name_last_part = Path(working.name).parts[-1]
            _file = Path(working.name, zipfile_path.name)
            if _file.exists() and _file.is_file(): # no sub-path used
                zipfile_names.append(zipfile_path.name)
                #LOGGER.info(f"DEBUG: found {zipfile_path.name}")
            elif Path(zipfile_path).exists() and Path(zipfile_path).is_file():
                for part in reversed(Path(zipfile_path).parts):
                    if path == working_name_last_part:
                        new_zipfile_name = f"{working_name_last_part}/{zipfile_path.name}"
                    zipfile_names.append(new_zipfile_name)
                LOGGER.info(f"Found file {zipfile_path.name} was in some sub-folder, so returning as {new_zipfile_name}.")
            else:
                #raise FileNotFoundError(f"ERROR: Downloaded a file but cannot parse the path to it: {zipfile_path}")
                return {"filenames": zipfile_names, "tempdir": working, "error": f"Downloaded a file but cannot parse the path to it: {zipfile_path}"}
        # LOGGER.info(f"DEBUG zipfile_names 485: {zipfile_names}")

        #3: memory check / debug
        if os.name != 'nt' and kwargs.get('verbose') == True:
            total_disk_space = subprocess.check_output(['du','-sh', EFS]).split()[0].decode('utf-8')
            LOGGER.info(f"Tempfolder {EFS} contains {len(list([str(k) for k in Path(EFS).rglob('*')]))} files, {total_disk_space} total.")
            LOGGER.info(f"DEBUG: {EFS} all files: {list([str(k) for k in Path(working.name).rglob('*')])}")

        #4: LOGGING OUT -- result dict tells this function what it found.
        if result['found_idats'] == False and result['processed_files'] == True and result['tbl_txt_files'] == False:
            LOGGER.info(f"Found {len(extracted_files)} files for {geo_id}.")
        if result['found_idats'] == False and result['processed_files'] == False and result['tbl_txt_files'] == False:
            LOGGER.warning(f"No downloadable methylation data found for {geo_id}.")
            continue
        if result['tbl_txt_files'] == True:
            LOGGER.info(f"Found -tbl-1.txt files with meta data for {geo_id}.")
        if result['tbl_txt_files'] == True:
            LOGGER.info(f"Found processed csv files for {geo_id}.")


        if result['found_idats'] == True:
            LOGGER.info(f"Found {len(extracted_files)} IDATs for {geo_id}.")
            try:
                # dict_only: should download IDATs without processing them.
                download_success = methylprep.download.process_data.run_series(geo_id,
                    working.name,
                    dict_only=True,
                    batch_size=BATCH_SIZE,
                    clean=True,
                    abort_if_no_idats=True)
            except Exception as e:
                LOGGER.error(f"ERROR run_series: {e}")
                import errno
                if hasattr(e,'errno') and e.errno == errno.ENOSPC:
                    import traceback;print('DEBUG run_series No Space:', traceback.format_exc())
                if hasattr(e,'errno') and e.errno == errno.ENOENT:
                    import traceback;print('DEBUG run_series File not found:', traceback.format_exc())
                return {"error":e, "filenames": zipfile_names, "tempdir": working}
            if download_success is False:
                return {"error":"IDATs detected but failed to download IDAT files", "filenames": zipfile_names, "tempdir": working}
            LOGGER.info(f"Ran series {geo_id}")
            extracted_files = [_file.name for _file in list(Path(working.name).rglob('*')) if not _file.is_dir()]
            zipfile_paths = [_file for _file in list(Path(working.name).rglob('*')) if not _file.is_dir()]
            #LOGGER.info(f"DEBUG 521: UPDATED the extracted files {zipfile_paths} || {extracted_files}")

        #5: compress and package files, and move out of working.
        if kwargs.get('compress') is True:
            # check if any one file is too big, and switch to gzip if so.
            use_gzip = False if result['found_idats'] else True # True yields separate files
            for k in Path(working.name).rglob('*'):
                if k.stat().st_size >= zipfile.ZIP64_LIMIT:
                    # this next line assumes only one GSE ID in function.
                    LOGGER.info(f"Switching to gzip because {str(k)} is greater than 2GB. This probably breaks if processing IDATS. Found_idats = {result['found_idats']}")
                    use_gzip = True
                    break

            # save all process files and move
            # IDATs will also be gzipped separately if they're too big.
            debug_zip_paths = []
            if use_gzip:
                for k in zipfile_paths:
                    if '.gz' in k.suffixes:
                        LOGGER.info(f"DEBUG: skip gzipping {k}")
                        continue
                    LOGGER.info(f"DEBUG: gzipping {k}")
                    # gzip each file in-place, then upload them. These are big pkl files.
                    # this will also catch the GSExxxxx-tbl-1.txt pickled beta_values.pkl dataframe, if it exists.
                    with open(k, 'rb') as file_in:
                        gzip_name = Path(working.name, f"{str(k.name)}.gz")
                        LOGGER.info(f"DEBUG: gzip outfile: {gzip_name}")
                        with gzip.open(gzip_name, 'wb') as file_out:
                            shutil.copyfileobj(file_in, file_out)
                            # --- REMOVE the prev file here ???
                            zipfile_names.remove(k.name)
                            #file_size = file_out.seek(0, io.SEEK_END)
                            #LOGGER.info(f"File: {gzip_name.name} -- {round(file_size/1000000)} MB")
                            if os.name != 'nt':
                                presult = subprocess.run(["du", "-chs", EFS], stdout=subprocess.PIPE)
                                presult = presult.stdout.decode('utf-8').split('\t')[0]
                                LOGGER.info(f"{gzip_name.name} -- {presult} total")
                    debug_zip_paths.append(gzip_name.name)
                    zipfile_names.append(gzip_name.name)
            else:
                zipfile_name = f"{geo_id}.zip"
                with zipfile.ZipFile(Path(working.name, zipfile_name),
                    mode='w',
                    compression=zipfile.ZIP_DEFLATED,
                    allowZip64=True,
                    compresslevel=9) as zip:

                    #for k in Path(working.name).rglob('*'):
                    for k in zipfile_paths:
                        LOGGER.info(f"DEBUG: zipping {k}")
                        if k.is_dir():
                            continue
                        if k.name == zipfile_name:
                            continue # there is an empty file created in the same folder I'm zipping up, so need to skip this guy.
                        zip.write(str(k), k.name) # 2nd arg arcname will drop folder structure in zipfile (the /mnt/efs/tmpfolder stuff)
                        zipinfo = zip.getinfo(k.name)
                        LOGGER.info(f"{k.name} ({round(zipinfo.file_size/1000)} --> {round(zipinfo.compress_size/1000)} KB)")
                    LOGGER.info(f"In ZipFile {Path(working.name, zipfile_name)}: {zip.namelist()}")
                zipfile_names.append(zipfile_name)

        if kwargs.get('move') == True:
            if result['found_idats'] == True:
                files = Path(working.name).rglob('*')
            else:
                files = zipfile_paths
            for _file in files:
                if '.DS_Store' in _file.name:
                    continue
                # print(_file)
                try:
                    shutil.move(str(_file), Path(kwargs['data_dir']))
                except Exception as e:
                    LOGGER.warning(f"couldn't move file: {e}")

        #7: delete/move tempfile/efs stuff
        #efs_files = list(Path(working.name).glob('*'))
        # NOTE: if running a bunch of GEO_IDs, and clean is False, everything will be in the same temp workdir, so make sure to clean/move each one.
        # skip this step if you need to keep folders in working/efs folder instead of moving them to the data_dir.
        if kwargs.get('clean') is True and kwargs.get('move') is False: # and kwargs.get('compress') is True:
            # moving important files out of working folder, then returning these as a list.
            # this would break in lambda/aws-batch because efs is huge but the host drive is not.
            #zipfiles = [_zipfile for _zipfile in list(Path(working.name).glob('*')) if _zipfile in zipfile_names]
            for zipfile_name in zipfile_names:
                if kwargs.get('verbose') == True:
                    LOGGER.info(f"Copying {zipfile_name} to {kwargs.get('data_dir','.')}")
                shutil.copy(Path(working.name, zipfile_name), kwargs.get('data_dir','.'))
            efs_exists = Path(working.name).exists()
            if len(geo_ids) > 1 and geo_id != geo_ids[-1]:
                # remove files, but don't let folder go away with cleanup(), unless this is the last one.
                for _file in Path(working.name).rglob('*'):
                    if _file.is_file():
                        _file.unlink()
                LOGGER.info(f"Removed temp files; left working dir intact.")
            else:
                working.cleanup()
            final_files = list(Path(kwargs.get('data_dir','.')).glob('*'))
            if kwargs.get('verbose') == True:
                LOGGER.info(f"Removing temp_files: (exists: {efs_exists} --> exists: {Path(working.name).exists()})")
                LOGGER.info(f"Files saved: {final_files} vs zipfile_names: {zipfile_names}")

        if kwargs.get('move') is True:
            working = kwargs['data_dir'] # return the correct folder location if moving, after clean step

        result['filenames'] = zipfile_names
        each_geo_result[geo_id] = result

    if len(geo_ids) > 1:
        LOGGER.warning("Returning a different, nested DICT structure because more than one GEO_ID was processed in this function.")
        return {"error": None, "geo_ids": each_geo_result, "tempdir": working}
    return {
        "error":None, "filenames": zipfile_names, "tempdir": working,
        "found_idats": result["found_idats"],
        "processed_files": result["processed_files"],
        "tbl_txt_files": result["tbl_txt_files"],
        }


def download_geo_processed(geo_id, working, verbose=False, use_headers=False):
    """Uses methylprep/methylcheck to get processed beta values.
    use_headers: if True, it will use the series_matrix headers to create a samplesheet instead of MiNiML file, which is faster,
    but lacks Sentrix_ID and Sentrix_Position data. But this is only needed for methylprep process. In this case, the
    meta data and betas are keyed using GSM_IDs instead."""
    filename_keywords = ['matrix', 'processed', 'signals', 'intensities', 'normalized', 'intensity', 'raw_data', 'mean', 'average', 'beta']
    filename_exclude_keywords = ['RNA','Illumina']
    if verbose:
        LOGGER.info(f"Searching GEO NIH for {geo_id}")
    result_df = search(geo_id, filepath=working.name, verbose=verbose)
    """ note: each df row columns are:
    ROW = {'title': 'GSE:  'url': 'samples':
        'date': 'platform': 'idats':
        file_name_1, file_size_1, file_link_1 ... }
    """
    # confirm no idats, then
    # find the links and download if match pattern.
    found_idats = False
    found_headers = False
    downloaded_files = False
    tbl_txt_files = False
    for idx, row in result_df.iterrows():
        # first, skip rows if IDATs exist for the series.
        if verbose:
            LOGGER.info(f"{idx}: {dict(row)}")
        if row.get('idats') == '1':
            if verbose:
                LOGGER.info(f"IDATs exist for {geo_id}")
            found_idats = True
            return {'found_idats': found_idats, 'processed_files': downloaded_files, 'tbl_txt_files': tbl_txt_files}

        # second, see if the series_matrix file exists, and use that.
        series_server = "https://ftp.ncbi.nlm.nih.gov"
        series_path = f"geo/series/{geo_id[:-3]}nnn/{geo_id}/matrix"
        file_name = f"{geo_id}_series_matrix.txt.gz"
        series_matrix_path = f"{series_server}/{series_path}/{file_name}"
        series_matrix_path2 = f"{series_server}/geo/series/{geo_id[:-3]}nnn/matrix/{geo_id}/{file_name}"
        # some series paths include the GPL platform code too, if multiple platforms are present
        series_matrix_path3 = f"{series_server}/{series_path}/{geo_id}-GPL21145_series_matrix.txt.gz"
        series_matrix_path4 = f"{series_server}/{series_path}/{geo_id}-GPL13534_series_matrix.txt.gz"
        try:
            saved_file = file_name.replace(' ','_')
            response = requests.head(series_matrix_path)
            response2 = requests.head(series_matrix_path2)
            response3 = requests.head(series_matrix_path3)
            response4 = requests.head(series_matrix_path4)
            if response.status_code == 200:
                url = series_matrix_path
            elif response2.status_code == 200:
                url = series_matrix_path2
            elif response3.status_code == 200:
                url = series_matrix_path3
            elif response4.status_code4 == 200:
                url = series_matrix_path4
            else:
                url = None
            if url != None:
                if verbose:
                    LOGGER.info(f"Downloading {series_matrix_path}")
                with requests.get(url, stream=True) as r:
                    total_size_in_bytes= int(r.headers.get('content-length', 0))
                    block_size = 1024 #1 Kibibyte
                    progress_bar = tqdm(desc=saved_file, total=total_size_in_bytes, unit='iB', unit_scale=True)
                    with open(Path(working.name, saved_file), 'wb') as f:
                        #shutil.copyfileobj(r.raw, f) --- needed when saving to EFS
                        for data in r.iter_content(block_size):
                            progress_bar.update(len(data))
                            f.write(data)
                    progress_bar.close()
                if total_size_in_bytes != 0 and progress_bar.n != total_size_in_bytes:
                    LOGGER.error("ERROR, something went wrong")
                else:
                    # TEST the downloaded file. possibly avoid downloading miniml file because this has the meta data already.
                    saved_file_path = Path(working.name,saved_file)
                    if url == None or saved_file_path.exists() == False:
                        LOGGER.info("no series_matrix file downloaded")
                    elif '.gz' in saved_file_path.suffixes:
                        LOGGER.info("Un-gzipping...")
                        unzipped_file = str(saved_file_path)[:-3]
                        with gzip.open(saved_file_path, 'rb') as f_in:
                            with open(unzipped_file, 'wb') as f_out:
                                shutil.copyfileobj(f_in, f_out)
                        # delete gzip
                        if Path(unzipped_file).exists():
                            saved_file_path.unlink()
                        if Path(unzipped_file).stat().st_size < 100000:
                            LOGGER.warning(f"Series Matrix file size ({Path(unzipped_file).stat().st_size/1000}K) is too small to contain beta_values.")
                        import methylcheck
                        data = methylcheck.read_geo_processed.read_series_matrix(unzipped_file, include_headers_df=True)
                        if verbose:
                            LOGGER.info(f"{geo_id} data: {data['df'].shape} headers: {data['headers_df'].shape}")
                        downloaded_files = True
                        if isinstance(data.get('headers_df'),pd.DataFrame):
                            found_headers = True
                        try:
                            samplesheet = samplesheet_from_series_matrix(data['headers_df'])
                            samplesheet.to_csv(Path(Path(unzipped_file).parent, f"{geo_id}_samplesheet.csv"), index=False)
                        except Exception as e:
                            LOGGER.error(f"Could not create samplesheet from series_matrix headers: {e}")
                        if data.get('series_dict'):
                            with open(Path(Path(unzipped_file).parent, f"{geo_id}_series_summary.json"), 'w', encoding='utf8') as f:
                                json.dump(data['series_dict'],f)
                        if isinstance(data.get('df'), pd.DataFrame): # betas
                            if len(data['df']) == 0:
                                LOGGER.error(f"beta values DataFrame is empty: {data['df'].shape}")
                                downloaded_files = False # the processed data file (series matrix) was empty, so setting to False
                            else:
                                data['df'].to_pickle(Path(Path(unzipped_file).parent, f"{geo_id}_beta_values.pkl"))
                        # if not compressing later, move this file out of tempfolder
                        if compress != True:
                            current = Path(unzipped_file).parent
                            parent = Path(unzipped_file).parent.parent
                            # shutil.move(unzipped_file, str(parent)) --- the .txt file is no longer needed. everything is repackaged.
                            shutil.move(str(Path(current, f"{geo_id}_samplesheet.csv")), str(parent))
                            shutil.move(str(Path(current, f"{geo_id}_series_summary.json")), str(parent))
                            shutil.move(str(Path(current, f"{geo_id}_beta_values.pkl")), str(parent))
                        continue
        except Exception as e:
            LOGGER.info(f"Series_matrix download failed: {e}, trying other saved files")
            import traceback
            LOGGER.info(traceback.format_exc())
            pass

        # third, follow the other file links to processed data from the search DF.
        for i in range(3):
            if row.get(f'file_name_{i}') and row.get(f'file_size_{i}') and row.get(f'file_link_{i}'):
                file_name = row.get(f'file_name_{i}')
                file_size = row.get(f'file_size_{i}')
                file_link = row.get(f'file_link_{i}')
                if (any([keyword.lower() in file_name.lower() for keyword in filename_keywords]) and
                    all([keyword.lower() not in file_name.lower() for keyword in filename_exclude_keywords])
                    ):
                    if verbose:
                        LOGGER.info(f"Matched {file_name}; size: {file_size}")
                    if 'ftp://' in file_link:
                        raw_filename = file_link
                        series_path = working.name
                        saved_file = file_name.replace(' ','_')
                        ftp = FTP('ftp.ncbi.nlm.nih.gov', timeout=120) # 2 mins
                        ftp.login()
                        ftp.cwd(f"geo/series/{geo_id[:-3]}nnn/{geo_id}")
                        raw_file = open(f"{series_path}/{saved_file}", 'wb')
                        filesize = ftp.size(f"suppl/{file_name}")
                        LOGGER.info(f"FTPing {saved_file} -- {round(filesize/1000000)} MB") #  from geo/series/{geo_id[:-3]}nnn/{geo_id}
                        try:
                            #ftp.retrbinary(f"RETR suppl/{raw_filename}", raw_file.write)
                            ftp.retrbinary(f"RETR suppl/{file_name}", raw_file.write)
                        except Exception as e:
                            LOGGER.error("error: {e}, trying {file_link} instead of geo/series/{geo_id[:-3]}nnn/{geo_id}/suppl/{raw_filename}")
                            ftp.retrbinary(f"RETR {file_link}", raw_file.write)
                        raw_file.close()
                        if verbose:
                            LOGGER.info(f"Downloaded {raw_filename}")
                        downloaded_files = True

                    elif 'https://' in file_link:
                        saved_file = file_name.replace(' ','_')
                        if verbose:
                            LOGGER.info(f"Downloading {saved_file} from {file_link}")
                        with requests.get(file_link, stream=True) as r:
                            with open(Path(working.name, saved_file), 'wb') as f:
                                shutil.copyfileobj(r.raw, f)
                        if verbose:
                            LOGGER.info(f"Downloaded {saved_file}")
                        downloaded_files = True
                    else:
                        LOGGER.error(f"Unrecognized protocol in {file_link}")

                    this = Path(working.name, saved_file)
                    try:
                        if verbose:
                            LOGGER.info(f"Trying read_geo() on {this}")
                        import methylcheck
                        beta_df = methylcheck.read_geo(this, verbose=verbose)
                        is_df = isinstance(beta_df, pd.DataFrame)
                        if is_df and verbose:
                            LOGGER.info(f"Result is a dataframe.")
                        else:
                            if verbose:
                                LOGGER.info(f"Result is NOT a dataframe.")
                    except Exception as e:
                        import traceback
                        LOGGER.info(f"ERROR: {e}")
                        LOGGER.info(traceback.format_exc())
                        return {'found_idats': found_idats, 'processed_files': downloaded_files, 'tbl_txt_files': tbl_txt_files}

                    if is_df:
                        # save to disk, then load again. don't overwriting pre-existing file of same name
                        df_file = Path(working.name, f"{geo_id}_beta_values.pkl")
                        if df_file.exists(): # could use a while loop here...
                            if verbose:
                                LOGGER.info(f"{df_file} already exists. Trying an alternate name.")
                            alt_name = Path(saved_file).stem.replace(geo_id,'')
                            df_file = Path(working.name, f"{geo_id}_beta_values_from_{alt_name}.pkl")
                            if verbose:
                                LOGGER.info(f"Alt name: {df_file}")
                            if df_file.exists(): # guaranteed to work, but less informative.
                                random_id = ''.join([random.choice('abcdefghjiklmnopqrstuvwxyz1234567890') for i in range(16)])
                                df_file = Path(working.name, f"{geo_id}_beta_values_{random_id}.pkl")
                        beta_df.to_pickle(df_file)
                        if verbose:
                            LOGGER.info(f"Saved {df_file.name}")
                        try:
                            if verbose:
                                LOGGER.info(f"reopening {df_file.name}")
                            if 'methylcheck' in sys.modules:
                                beta_df = methylcheck.load(df_file, verbose=verbose, silent=(not verbose))
                            else:
                                beta_df = pd.read_pickle(df_file)
                            if isinstance(beta_df, pd.DataFrame) and verbose:
                                LOGGER.info(f"df shape: {beta_df.shape}")
                        except Exception as e:
                            import traceback
                            LOGGER.info(f"ERROR: {e}")
                            #LOGGER.info(traceback.format_exc())
                else:
                    if verbose:
                        LOGGER.info(f"Skipped {file_name}")

    if found_headers == True and use_headers == True:
        pass
    elif downloaded_files == True and use_headers == False:
        if verbose:
            LOGGER.info(f"Getting MINiML meta data")
        ftp = FTP('ftp.ncbi.nlm.nih.gov', timeout=120) # 2 mins
        ftp.login()
        ftp.cwd(f"geo/series/{geo_id[:-3]}nnn/{geo_id}/miniml")
        # look for additional parts of xml files
        file_parts = []
        for filename,filestats in ftp.mlsd(facts=["size"]):
            if geo_id in filename:
                #LOGGER.info(f"DEBUG {filename} -- {round(int(filestats['size'])/1000000)} MB")
                file_parts.append(filename)
        if len(file_parts) > 1:
            LOGGER.info(f"The {geo_id}_family.xml miniml file for this series contains multiple ({len(file_parts)}) files.")

        try:
            # methylprep.download.convert_miniml
            local_files = convert_miniml(geo_id, data_dir=working.name, merge=False, download_it=True, extract_controls=False, require_keyword=None, sync_idats=False, verbose=verbose)
            # HERE -- copy all of these (.xml files) into the S3 output folder
            gsm_files = []
            for local_file in local_files:
                if '.tgz' in Path(local_file).suffixes:
                    if verbose == True:
                        LOGGER.info(f"Unpacking: {local_file}")
                    shutil.unpack_archive(local_file)
                    all_files = list(Path(working.name).rglob('*'))
                    gsm_files = list(Path(working.name).rglob('GSM*.txt'))
                    non_gsm_files = [file for file in all_files if file not in gsm_files]
                    if verbose == True:
                        LOGGER.info(f"{len(all_files)} local_files, {len(gsm_files)}, non-GSM: {len(non_gsm_files)} | GSMs: {gsm_files}")
                else:
                    if verbose == True:
                        LOGGER.info(f"local_file skipped: {local_file}")
            if len(gsm_files) > 0:
                tbl_txt_files = True
                LOGGER.info(f"DEBUG: gsm_files detected: {gsm_files} -- to_df()")
                beta_df = betas_from_tbl_txt_files(gsm_files) # overwrites each time, but last cycle should be a complete list.
                beta_file = f"{geo_id}_beta_values-tbl-1.pkl"
                beta_df.to_pickle(beta_file)
                if verbose:
                    LOGGER.info(f"{beta_file} written, exists: {Path(working.name,beta_file).exists()}")
                    # uploading to s3 will happen in outside function, because it is a .pkl file in working.name dir.

            if Path(working.name, f"{geo_id}_family.xml").exists():
                # check again, as unpack_archive should create this from multi-part _family files
                # this should create the meta_data.pkl file that gets auto-saved in pipeline.
                local_files = convert_miniml(geo_id,
                    data_dir=working.name,
                    merge=False,
                    download_it=False,
                    extract_controls=False,
                    require_keyword=None,
                    sync_idats=False)

        except Exception as e:
            LOGGER.error(f"convert_miniml error: {e}")
            import traceback
            LOGGER.error(traceback.format_exc())
    return {'found_idats': found_idats, 'processed_files': downloaded_files, 'tbl_txt_files': tbl_txt_files}


def betas_from_tbl_txt_files(file_list, remove_after=True):
    """input: list of file paths to be converted into one beta DF and saved, returning saved file path."""
    samples = [] # list of dfs to merge, with cols being sample names and probes as index
    for file in file_list:
        FILE = Path(file)
        if FILE.suffix == '.txt' and FILE.name.startswith('GSM') and '-tbl' in FILE.name:
            sample_id = FILE.name.split('-')[0]
            sample = pd.read_csv(FILE, sep='\t', header=0, names=['IlmnID', sample_id])
            sample = sample.set_index('IlmnID')
            samples.append(sample)
            #LOGGER.info(f"{len(samples)}: {sample_id} -> {len(sample)}")
    df_ok = False
    try:
        df = pd.concat(samples, axis=1, sort=False)
        df_ok = True
    except Exception as e:
        LOGGER.error("Could not concat these samples into a dataframe. The probe names in rows don't line up.")
        LOGGER.error(e)
    if df_ok == False:
        return
    if remove_after == True:
        for file in file_list: # deleting source files
            FILE = Path(file)
            if FILE.suffix == '.txt' and FILE.name.startswith('GSM') and '-tbl' in FILE.name:
                FILE.unlink()
    LOGGER.info(f"Found {df.shape[1]} samples with {df.shape[0]} probes from {len(file_list)} GSM-txt files.")
    return df


def search(keyword, filepath='.', verbose=True):
    """ CLI/cron function to check for new datasets.
    set up as a weekly cron.
    uses a local storage file to compare with old datasets in <pattern>_meta.csv.
    saves the dates of each dataset from GEO; calculates any new ones as new rows. updates csv.

options:
    pass in -k keyword
    verbose (True|False) --- reports to page; saves csv too

returns:
    saves a CSV to disk and returns a dataframe of results"""
    # FIRST: get any previous search results
    filename = Path(filepath, f'geo_alert {keyword}.csv')
    if filename.exists():
        prev_data = pd.read_csv(filename)
        if verbose:
            LOGGER.info(f'Previous search: {len(prev_data)} results')
    else:
        prev_data = None

    ROOT = 'http://www.ncbi.nlm.nih.gov/'
    # NOTE: query_page was limited to 20 results, so had to use summary table page instead
    #query_page = 'http://www.ncbi.nlm.nih.gov/gds/?term=GPL13534'
    summary_page = 'http://www.ncbi.nlm.nih.gov/geo/browse/?view=series&display=500&zsort=date&search='  #methylation'
    if keyword:
        for word in keyword.split():
            #query_page += f'+{word}'
            summary_page += f'+{word}'
    #query_page += '+AND+%22gse%22%5BEntry+Type%5D' # limits to datasets
    if verbose:
        LOGGER.info(summary_page)
    summary_html = urlopen(summary_page).read()
    #query_html = urlopen(query_page).read()
    try:
        soup = BeautifulSoup(summary_html, 'html.parser')
        total_series = soup.find(id='count')
        if verbose:
            LOGGER.info(total_series.text) # e.g. 1 series
        table = soup.find(id="geo_data")
        data = []
        for row in tqdm(table.find_all('tr'), desc='Checking for idats', disable=(not verbose)):
            gse = row.find('a') # always first row
            gse = gse.text if gse else 'None'
            if gse == 'None':
                continue # header
            title = row.find('td', {'class': 'title'})
            title = title.text.strip() if title else 'No title'
            samples = [i.text for i in row.find_all('a') if 'samples' in i.get('href','')]
            samples = samples[0] if len(samples) > 0 else 0
            try:
                samples = int(samples)
            except:
                if samples == 'Samples':
                    continue
            date = row.find('td', {'class': 'date'})
            if date:
                date = dateutil.parser.parse(date.text)
                date = date.strftime('%Y-%m-%d') # YEAR-month-day
            else:
                date = ''
            ROW = {
                'title': title,
                'GSE': gse,
                'url': ROOT + row.find('a').get('href', ''),
                'samples': samples,
                'date': date,
                'platform': None,
                'idats': None,
            }
            geo_acc_page = f"http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={gse}"
            html = urlopen(geo_acc_page).read()
            PLATFORMS = ['GPL21145', 'GPL13534', 'GPL23976', 'GPL8490', 'GPL16304', 'GPL18809'] # GPL16304, GPL28271 (Horvath)
            for platform in PLATFORMS:
                if platform in str(html):
                    ROW['platform'] = platform
                    break
            usable = methylprep.download.process_data.confirm_dataset_contains_idats(gse) # True|False
            info = methylprep.download.process_data.get_attachment_info(gse) # list of dicts for each file with name, link, size
            # prefill, so csv is square; always 3 file spots saved
            for i in range(3):
                ROW[f'file_name_{i}'] = ''
                ROW[f'file_size_{i}'] = ''
                ROW[f'file_link_{i}'] = ''
            for i, filedata in enumerate(info):
                ROW[f'file_name_{i}'] = filedata['name']
                ROW[f'file_size_{i}'] = filedata['size']
                ROW[f'file_link_{i}'] = filedata['link']
                if i == 2:
                    break
            ROW['idats'] = '1' if usable else '0'
            #soup = BeautifulSoup(query_html, 'html.parser')
            data.append(ROW)
    except Exception as e:
        LOGGER.error(e)
    if type(prev_data) != type(None):
        new_results = len(data) - len(prev_data)
        if verbose:
            LOGGER.info(f'{new_results} new results found.')
        # do stuff here with only the new results.
        # need a hook API here (i.e. slackbot) - call another CLI function
    # overwrite results either way.
    df = pd.DataFrame(data)
    df.to_csv(filename)
    if verbose:
        LOGGER.info(f"{filename} written")
    return df


def samplesheet_from_series_matrix(df):
    """input: header_df from methylcheck.read_series_matrix.
    This approach matches meta-data with sample betas without needing sentrix_ids. Key is GSMxxxx.
    This won't support methylprep process function, but fine for all post-process functions needing phenotype data.
    This parses multiple Characteristics columns into separate colummns in dataframe."""
    missing = ['Sentrix_ID', 'Sentrix_Position']
    columns = {
    '!Sample_geo_accession':'Sample_ID',
    '!Sample_source_name_ch1': 'source',
    '!Sample_platform_id':'platform',
    '!Sample_description':'Description', # may have more than one row in df
    '!Sample_characteristics_ch1':'Characteristics', # may have more than one row in df
    }
    #test: samplesheet = samplesheet_from_series_matrix(data['headers_df'])
    samplesheet_rows = []
    for gsm in df.columns:
        data = df[gsm]
        new = {}
        missing_data = Counter()
        for column,field in columns.items():
            if column == '!Sample_characteristics_ch1':
                continue # parse each row as a separate column, with : as key-value separator.
            if column not in data:
                missing_data[column] += 1
            else:
                if isinstance(data.loc[column], pd.Series):
                    # merge values
                    new[field] = " | ".join(list(data.loc[column].values))
                else:
                    new[field] = data.loc[column]
        if missing_data:
            LOGGER.warning(f"missing from meta data: {missing_data.most_common()}")
        column = '!Sample_characteristics_ch1'
        OVERWRITE_WARNINGS = Counter()
        if isinstance(data.loc[column], pd.Series): # only detects multi-line data here
            # parse each key-value as a separate column; only works if the labels in every sample are exactly the same.
            for item in data.loc[column]:
                if isinstance(item, str) and ':' in item:
                    try:
                        key,value = item.split(':')
                    except:
                        LOGGER.warning("could not split key-value pair because extra : present")
                        continue
                    if key.strip() in new: # possible that same key appears twice, or matches some other meta data. Data loss.
                        OVERWRITE_WARNINGS[key.strip()] += 1
                    else:
                        new[key.strip()] = value.strip()
                elif item == '':
                    continue
                else:
                    LOGGER.warning(f"Characteristic '{item}' not understood")
        if len(OVERWRITE_WARNINGS) > 0:
            LOGGER.warning(f"These Sample Characteristics had the same labels and were lost: {OVERWRITE_WARNINGS.most_common()}")
        new['GSM_ID'] = new.get('Sample_ID')
        for column,value in data.items():
            if column in columns.keys():
                continue
            field = column.replace('!Sample_','')
            # some fields are multiple rows, so represented by multiple index rows in this headers_df
            if field not in new and isinstance(data.loc[column], pd.Series):
                value = " | ".join(list(data.loc[column].values))
            if field not in new:
                new[field] = value
        new.update({field:None for field in missing})
        samplesheet_rows.append(new)
    # finally, make sure Sample_ID (index column) does not appear more than once. Some GEOs have this problem.
    # Avoid creating files with two Sample_ID columns, or any other duplicate columns:
    sample_sheet_meta_data = pd.DataFrame(data=samplesheet_rows)
    if any(sample_sheet_meta_data.columns.duplicated()):
        sample_sheet_meta_data = sample_sheet_meta_data.loc[:, ~sample_sheet_meta_data.columns.duplicated()]
    return sample_sheet_meta_data
