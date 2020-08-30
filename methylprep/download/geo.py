# Lib
import logging
from ftplib import FTP
import socket
from pathlib import Path, PurePath
import os
import tarfile
from tarfile import ReadError
import re
import gzip
import shutil
from bs4 import BeautifulSoup
import pickle
import pandas as pd
from tqdm import tqdm
from collections import Counter
# app
from .miniml import sample_sheet_from_miniml, sample_sheet_from_idats

#logging.basicConfig(level=logging.DEBUG) # always verbose
LOGGER = logging.getLogger(__name__)
LOGGER.setLevel( logging.INFO )

def geo_download(geo_id, series_path, geo_platforms, clean=True):
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

    This function returns True or False, depending on whether the downloaded data is correct."""
    success = True
    series_dir = Path(series_path)
    raw_filename = f"{geo_id}_RAW.tar"
    miniml_filename = f"{geo_id}_family.xml"

    if not os.path.exists(series_path):
        raise FileNotFoundError(f'{geo_id} directory not found.')

    for platform in geo_platforms:
        if not Path(f"{series_path}/{platform}").exists():
            Path(f"{series_path}/{platform}").mkdir()

    ftp = FTP('ftp.ncbi.nlm.nih.gov', timeout=120) # 2 mins
    ftp.login()
    ftp.cwd(f"geo/series/{geo_id[:-3]}nnn/{geo_id}")

    if not Path(f"{series_path}/{miniml_filename}").exists():
        if not Path(f"{series_path}/{miniml_filename}.tgz").exists():
            LOGGER.info(f"Downloading {miniml_filename}")
            miniml_file = open(f"{series_path}/{miniml_filename}.tgz", 'wb')
            try:
                filesize = ftp.size(f"miniml/{miniml_filename}.tgz")
                with tqdm(unit = 'b', unit_scale = True, leave = False, miniters = 1, desc = geo_id, total = filesize) as tqdm_instance:
                    def tqdm_callback(data):
                        tqdm_instance.update(len(data))
                        miniml_file.write(data)
                    ftp.retrbinary(f"RETR miniml/{miniml_filename}.tgz", tqdm_callback)
            except Exception as e:
                print(e)
                LOGGER.info('tqdm: Failed to create a progress bar, but it is downloading...')
                ftp.retrbinary(f"RETR miniml/{miniml_filename}.tgz", miniml_file.write)
            miniml_file.close()
            LOGGER.info(f"Downloaded {miniml_filename}")
        ftp.close()
        LOGGER.info(f"Unpacking {miniml_filename}")
        min_tar = tarfile.open(f"{series_path}/{miniml_filename}.tgz")
        for file in min_tar.getnames():
            if file == miniml_filename:
                min_tar.extract(file, path=series_path)
        min_tar.close()
        if clean:
            Path(f"{series_path}/{miniml_filename}.tgz").unlink()

    ftp = FTP('ftp.ncbi.nlm.nih.gov', timeout=59) # see issue https://bugs.python.org/issue30956 (must be <60s because of a bug)
    ftp.login()
    ftp.cwd(f"geo/series/{geo_id[:-3]}nnn/{geo_id}")

    if not list(series_dir.glob('**/*.idat')):
        if not list(series_dir.glob('*.idat.gz')):
            if not Path(f"{series_path}/{raw_filename}").exists():
                raw_file = open(f"{series_path}/{raw_filename}", 'wb')
                filesize = ftp.size(f"suppl/{raw_filename}")
                try:
                    try:
                        with tqdm(unit = 'b', unit_scale = True, leave = False, miniters = 1, desc = geo_id, total = filesize) as tqdm_instance:
                            def tqdm_callback(data):
                                tqdm_instance.update(len(data))
                                raw_file.write(data)
                            ftp.retrbinary(f"RETR suppl/{raw_filename}", tqdm_callback)
                    except Exception as e:
                        LOGGER.info('tqdm: Failed to create a progress bar, but it is downloading...')
                        ftp.retrbinary(f"RETR suppl/{raw_filename}", raw_file.write)
                except socket.timeout as e:
                    LOGGER.warning(f"FTP timeout error.")
                    # seems to happen AFTER download is done, so just ignoring it.
                LOGGER.info(f"Closing file {raw_filename}")
                raw_file.close()
                LOGGER.info(f"Downloaded {raw_filename}")
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
        LOGGER.info(f"Decompressing {geo_id} IDAT files")
        for gz in series_dir.glob("*.idat.gz"):
            gz_string = str(gz)
            with gzip.open(gz_string, 'rb') as f_in:
                with open(gz_string[:-3], 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
            if clean:
                os.remove(gz_string)

    LOGGER.info(f"Downloaded and unpacked {geo_id}")
    ftp.quit()
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
                    try:
                        shutil.move(f"{series_path}/{file_name}", f"{series_path}/{platform}/{file_name}")
                    except FileNotFoundError:
                        # this doesn't throw an error if file is already in the right folder
                        if not Path(f"{series_path}/{platform}/{file_name}").is_file():
                            raise FileNotFoundError ("Could not move file after downloading.")

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
