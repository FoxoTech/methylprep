# Lib
import logging
from ftplib import FTP
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

LOGGER = logging.getLogger(__name__)

def geo_download(geo_id, series_path, geo_platforms, clean=True):
    """Downloads the IDATs and metadata for a GEO series

    Arguments:
        geo_id [required]
            the GEO Accension for the desired series
        series_path [required]
            the directory to download the data to
        geo_platforms [required]
            the list of supported GEO platforms
        clean
            whether or not to delete files once they are no longer need (True by default)"""
    series_dir = Path(series_path)
    raw_filename = f"{geo_id}_RAW.tar"
    miniml_filename = f"{geo_id}_family.xml"

    if not os.path.exists(series_path):
        raise FileNotFoundError(f'{geo_id} directory not found.')

    for platform in geo_platforms:
        if not os.path.exists(f"{series_path}/{platform}"):
            os.mkdir(f"{series_path}/{platform}")

    ftp = FTP('ftp.ncbi.nlm.nih.gov')
    ftp.login()
    ftp.cwd(f"geo/series/{geo_id[:-3]}nnn/{geo_id}")

    #LOGGER.info(f"Downloading {geo_id}")
    if not list(series_dir.glob('**/*.idat')):
        if not list(series_dir.glob('*.idat.gz')):
            if not os.path.exists(f"{series_path}/{raw_filename}"):
                #LOGGER.info(f"Downloading {raw_filename}")
                raw_file = open(f"{series_path}/{raw_filename}", 'wb')
                filesize = ftp.size(f"suppl/{raw_filename}")
                try:
                    with tqdm(unit = 'b', unit_scale = True, leave = False, miniters = 1, desc = geo_id, total = filesize) as tqdm_instance:
                        def tqdm_callback(data):
                            tqdm_instance.update(len(data))
                            raw_file.write(data)
                        ftp.retrbinary(f"RETR suppl/{raw_filename}", tqdm_callback)
                except Exception as e:
                    LOGGER.info('tqdm: Failed to create a progress bar, but it is downloading...')
                    ftp.retrbinary(f"RETR suppl/{raw_filename}", raw_file.write)
                raw_file.close()
                LOGGER.info(f"Downloaded {raw_filename}")
            LOGGER.info(f"Unpacking {raw_filename}")
            try:
                tar = tarfile.open(f"{series_path}/{raw_filename}")
                for member in tar.getmembers():
                    if re.match('.*.idat.gz', member.name):
                        tar.extract(member, path=series_path)
            except ReadError as e:
                raise ReadError(f"There appears to be an incomplete download of {geo_id}. Please delete those files and run this again.")
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

    if not os.path.exists(f"{series_path}/{miniml_filename}"):
        if not os.path.exists(f"{series_path}/{miniml_filename}.tgz"):
            LOGGER.info(f"Downloading {miniml_filename}")
            miniml_file = open(f"{series_path}/{miniml_filename}.tgz", 'wb')
            ftp.retrbinary(f"RETR miniml/{miniml_filename}.tgz", miniml_file.write)
            miniml_file.close()
            LOGGER.info(f"Downloaded {miniml_filename}")
        LOGGER.info(f"Unpacking {miniml_filename}")
        min_tar = tarfile.open(f"{series_path}/{miniml_filename}.tgz")
        for file in min_tar.getnames():
            if file == miniml_filename:
                min_tar.extract(file, path=series_path)
        if clean:
            os.remove(f"{series_path}/{miniml_filename}.tgz")

    LOGGER.info(f"Downloaded and unpacked {geo_id}")

    ftp.quit()

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
        split_idat = sample.find('Supplementary-Data').text.split("/")[-1].split("_")
        attributes_dir['methylprep_name'] = f"{split_idat[1]}_{split_idat[2]}"

        if platform in geo_platforms:
            for idat in sample.find_all('Supplementary-Data'):
                if idat['type'] == 'IDAT':
                    file_name = (idat.text.split("/")[-1]).strip()[:-3]
                    shutil.move(f"{series_path}/{file_name}", f"{series_path}/{platform}/{file_name}")

            meta_dicts[platform][accession] = attributes_dir
            samples_dict[platform][accession] = title
        else:
            raise ValueError(f'Sample: {title} has unrecognized platform: {platform}')
    fp.close()

    seen_platforms = []

    for platform in geo_platforms:
        if meta_dicts[platform]:
            meta_dict_filename = f"{geo_id}_{platform}_dict.pkl"
            pickle.dump(meta_dicts[platform], open(f"{series_path}/{meta_dict_filename}", 'wb'))
            if not os.path.exists(f"{path}/{platform}_dictionaries/{geo_id}_dict.pkl"):
                shutil.copyfile(f"{series_path}/{meta_dict_filename}", f"{path}/{platform}_dictionaries/{geo_id}_dict.pkl")
            sample_sheet_from_min(geo_id, series_path, platform, samples_dict[platform])
            if platform not in seen_platforms:
                seen_platforms.append(platform)

    return seen_platforms

def sample_sheet_from_min(geo_id, series_path, platform, platform_samples_dict):
    """Creates a sample_sheet for all samples of a particular platform for a given series

    Arguments:
        geo_id [required]
            the GEO Accension for the desired series
        series_path [required]
            the directory containing the series data
        platform [required]
            the platform to generate a sample sheet for
        platform_samples_dict
            the dictionary of samples for the given platform"""
    series_dir = Path(f"{series_path}/{platform}")
    idat_files = series_dir.glob('*Grn.idat')
    _dict = {'GSM_ID': [], 'Sample_Name': [], 'Sentrix_ID': [], 'Sentrix_Position': []}
    for idat in idat_files:
        filename = str(idat).split("/")[-1]
        if re.match('(GSM[0-9]+_[0-9a-zA-Z]+_R0[0-9].0[0-9].Grn.idat)', filename):
            split_filename = filename.split("_")
            _dict['GSM_ID'].append(split_filename[0])
            _dict['Sentrix_ID'].append(split_filename[1])
            _dict['Sentrix_Position'].append(split_filename[2])
        else:
            raise ValueError(f"{filename} has an unexpected naming format")

    # confusing logic here, names_dir is unordered, but retrieving by key in correct order alleviates
    for key in _dict['GSM_ID']:
        _dict['Sample_Name'].append(platform_samples_dict[key])

    df = pd.DataFrame(data=_dict)
    df.to_csv(path_or_buf=(PurePath(f"{series_path}/{platform}", 'samplesheet.csv')),index=False)
