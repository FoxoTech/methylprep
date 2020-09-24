# Lib
import logging
from ftplib import FTP
from pathlib import Path, PurePath
import os
import zipfile
import csv
import shutil
import pickle
import re
import pandas as pd
from tqdm import tqdm


LOGGER = logging.getLogger(__name__)

def ae_download(ae_id, series_path, ae_platforms, clean=True):
    """Downloads the IDATs and metadata for an ArrayExpress series

    Arguments:
        ae_id [required]
            the ArrayExpress Accension for the desired series
        series_path [required]
            the directory to download the data to
        ae_platforms [required]
            the list of supported ArrayExpress platforms
        clean
            whether or not to delete files once they are no longer need (True by default)"""
    success = True
    series_dir = Path(series_path)
    sdrf_filename = f"{ae_id}.sdrf.txt"

    if not os.path.exists(series_path):
        raise FileNotFoundError(f'{ae_id} directory not found')

    for platform in ae_platforms:
        if not os.path.exists(f"{series_path}/{platform}"):
            os.mkdir(f"{series_path}/{platform}")

    ftp = FTP('ftp.ebi.ac.uk', timeout=59)
    ftp.login()
    ae_split = ae_id.split("-")[1] # MTAB, GEOD, etc
    ftp.cwd(f"/pub/databases/arrayexpress/data/experiment/{ae_split}/{ae_id}")

    #LOGGER.info(f"Downloading {ae_id}")
    if not list(series_dir.glob('**/*.idat')):
        for file in ftp.nlst():
            if file.split(".")[1] == 'raw':
                if not os.path.exists(f"{series_path}/{file}"):
                    #LOGGER.info(f"Downloading {file} from ArrayExpress")
                    raw_file = open(f"{series_path}/{file}", 'wb')
                    filesize = ftp.size(file)
                    try:
                        with tqdm(unit = 'b', unit_scale = True, leave = False, miniters = 1, desc = ae_id, total = filesize) as tqdm_instance:
                            def tqdm_callback(data):
                                tqdm_instance.update(len(data))
                                raw_file.write(data)
                            ftp.retrbinary(f"RETR {file}", tqdm_callback)
                    except Exception as e:
                        LOGGER.error(e)
                        LOGGER.info('tqdm: Failed to create a progress bar, but it is downloading...')
                        ftp.retrbinary(f"RETR {file}", raw_file.write)
                    raw_file.close()
                    LOGGER.info(f"Downloaded {file}")
                LOGGER.info(f"Unpacking {file}")
                zips = series_dir.glob('*.raw.*.zip')
                for zip in zips:
                    zip_obj = zipfile.ZipFile(str(zip))
                    zip_obj.extractall(path=series_path)
                    if clean:
                        LOGGER.info(f"Removing {zip}")
                        Path(zip).unlink()

    if not os.path.exists(f"{series_path}/{sdrf_filename}"):
        LOGGER.info(f"Downloading {sdrf_filename}")
        sdrf_file = open(f"{series_path}/{sdrf_filename}", 'wb')
        ftp.retrbinary(f"RETR {sdrf_filename}", sdrf_file.write)
        sdrf_file.close()

    LOGGER.info(f"Downloaded and unpacked {ae_id}")
    ftp.quit()
    return success


def ae_metadata(ae_id, series_path, ae_platforms, path):
    """Reads the metadata for the given series (SDRF file) and creates a metadata dictionary and sample sheet for each platform

    Arguments:
        ae_id [required]
            the ArrayExpress Accension for the desired series
        series_path [required]
            the directory containing the series data
        ae_platforms [required]
            the list of supported ArrayExpress platforms
        path
            the path to the directory containing dictionaries for each platform

    Returns:
        A list of platforms that the series contains samples of"""
    pipeline_kwargs = {}
    with open(f"{series_path}/{ae_id}.sdrf.txt") as fp:
        reader = csv.reader(fp, delimiter="\t")
        d = list(reader)

    meta_dicts = {}

    for platform in ae_platforms:
        meta_dicts[platform] = {}

    num_cols = len(d[0])
    num_rows = len(d)

    for row_num in range(1, num_rows):
        if (row_num % 2) == 0:
            sample_dict = {}
            for i in range (0, num_cols):
                sample_dict[d[0][i]] = d[row_num][i]
            sample_name = f"{ae_id}-sample-{str(int(row_num / 2))}"

            platform = sample_dict['Array Design REF']
            if platform in ae_platforms:
                red_idat = sample_dict['Array Data File']
                green_idat = red_idat[:-8] + "Grn.idat"

                shutil.move(f"{series_path}/{red_idat}", f"{series_path}/{platform}/{red_idat}")
                shutil.move(f"{series_path}/{green_idat}", f"{series_path}/{platform}/{green_idat}")

                meta_dicts[platform][sample_name] = sample_dict
            else:
                raise ValueError(f'Sample: {sample_name} has unrecognized platform: {platform}')
    fp.close()

    seen_platforms = []

    for platform in ae_platforms:
        if meta_dicts[platform]:
            meta_dict_filename = f"{ae_id}_{platform}_dict.pkl"
            pickle.dump(meta_dicts[platform], open(f"{series_path}/{meta_dict_filename}", 'wb'))
            if not os.path.exists(f"{path}/{platform}_dictionaries/{meta_dict_filename}"):
                shutil.copyfile(f"{series_path}/{meta_dict_filename}", f"{path}/{platform}_dictionaries/{ae_id}_dict.pkl")
            sample_sheet_from_sdrf(ae_id, series_path, platform, meta_dicts[platform])
            if platform not in seen_platforms:
                seen_platforms.append(platform)

    return seen_platforms, pipeline_kwargs


def sample_sheet_from_sdrf(ae_id, series_path, platform, platform_samples_dict):
    """Creates a sample_sheet for all samples of a particular platform for a given series

    Arguments:
        ae_id [required]
            the ArrayExpress Accension for the desired series
        series_path [required]
            the directory containing the series data
        platform [required]
            the platform to generate a sample sheet for
        platform_samples_dict
            the dictionary of samples for the given platform"""
    _dict = {'Sample_Name': [], 'Sentrix_ID': [], 'Sentrix_Position': []}
    for sample_name in platform_samples_dict:
        filename = platform_samples_dict[sample_name]['Array Data File']
        if re.match('[0-9a-zA-Z]+_R0[0-9].0[0-9].Red.idat', filename):
                split_assay_name = filename.split("_")

                _dict['Sample_Name'].append(sample_name)
                _dict['Sentrix_ID'].append(split_assay_name[0])
                _dict['Sentrix_Position'].append(split_assay_name[1])
        else:
            raise ValueError(f"{sample_name} Array Data File column has unexpected format")

    df = pd.DataFrame(data=_dict)
    df.to_csv(path_or_buf=(PurePath(f"{series_path}/{platform}", 'samplesheet.csv')),index=False)
