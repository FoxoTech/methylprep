# Lib
import os
import logging
from pathlib import Path, PurePath
from urllib.request import urlopen
import shutil
import pandas as pd
from bs4 import BeautifulSoup
import re
# App
from .geo import (
    geo_download,
    geo_metadata
)
from .array_express import(
    ae_download,
    ae_metadata
)
from methylprep import run_pipeline


LOGGER = logging.getLogger(__name__)

GEO_PLATFORMS = ['GPL21145', 'GPL13534', 'GPL8490'] # GPL23976 -- is an 850k genotyping array that is bundled with some datasets
AE_PLATFORMS = ['A-MEXP-2255', 'A-GEOD-21145', 'A-GEOD-13534']
PLATFORMS = GEO_PLATFORMS + AE_PLATFORMS
BATCH_SIZE = 100


def run_series(id, path, dict_only=False, batch_size=BATCH_SIZE, clean=True, abort_if_no_idats=True,
               decompress=True):
    """Downloads the IDATs and metadata for a series then generates one metadata dictionary and one beta value matrix for each platform in the series

    Arguments:
        id [required]
            the series ID (can be a GEO or ArrayExpress ID)
        path [required]
            the path to the directory to download the data to. It is assumed a dictionaries and beta values
            directory has been created for each platform (and will create one for each if not)
        dict_only
            if True, downloads idat files and meta data and creates data dictionaries for each platform, but does not process them further.
        batch_size
            the batch_size to use when processing samples (number of samples run at a time).
            By default is set to the constant 100.
        clean
            if True, removes intermediate processing files"""
    if not Path(f"{str(path)}/{PLATFORMS[0]}_dictionaries").exists():
        initialize(str(path))

    path = str(path)
    series_path = f"{path}" #"/{id}"

    if not os.path.exists(series_path):
        LOGGER.info(f"Creating directory for {id}")
        os.mkdir(series_path)

    download_success = False
    if id[:3] == 'GSE':
        series_type = 'GEO'
        if confirm_dataset_contains_idats(id) == False:
            LOGGER.error(f"[!] Geo data set {id} probably does NOT contain usable raw data (in .idat format). Not downloading.") # Press CTRL-C to cancel the download.")
        if abort_if_no_idats and confirm_dataset_contains_idats(id) == False:
            download_success = False
        else:
            download_success = geo_download(id, series_path, GEO_PLATFORMS, clean=clean, decompress=decompress)
    elif id[:7] == 'E-MTAB-':
        series_type = 'AE'
        download_success = ae_download(id, series_path, AE_PLATFORMS, clean=clean)
    else:
        raise ValueError(f"[ERROR] Series type not recognized. (The ID should begin with GSE or E-MTAB-)")

    if download_success == True:
        dicts = list(Path(series_path).rglob(f'{id}_dict.pkl'))
        if not dicts:
            if series_type == 'GEO':
                seen_platforms, pipeline_kwargs = geo_metadata(id, series_path, GEO_PLATFORMS, str(path))
            elif series_type == 'AE':
                seen_platforms, pipeline_kwargs = ae_metadata(id, series_path, AE_PLATFORMS, str(path))
        else:
            pipeline_kwargs = {} # ambigious whether {'make_sample_sheet':True} is needed here
            seen_platforms = []
            for d in dicts:
                for platform_name in PLATFORMS:
                    if platform_name in str(d): # case sensitive, and Path().match fails
                        seen_platforms.append(platform_name)  #str(d).split("/")[-1].split("_")[1])

    cleanup(str(path))
    if not dict_only and download_success:
        if pipeline_kwargs.get('make_sample_sheet') == True:
            pipeline_kwargs['meta_data_frame'] = True # otherwise, don't make a second meta_data pkl
        process_series(id, str(path), seen_platforms, batch_size, **pipeline_kwargs)
    if download_success == False:
        LOGGER.error("Series failed to download successfully.")
    return download_success


def process_series(id, path, seen_platforms, batch_size, **kwargs):
    """Processes the samples for each platform for the specified series, saving one pickled dataframe for each platform

    Arguments:
        id [required]
            the Accension for the desired series
        series_path [required]
            the directory containing the series data
        seen_platforms [required]
            the platforms the series has samples of
        batch_size
            the number of samples to process at a time"""
    for platform in seen_platforms:
        if not Path(path, platform, f"{id}_beta_values.pkl").exists():
            data_dir = f"{path}/{platform}"
            LOGGER.info(f"Processing {id} -- {platform} samples")
            LOGGER.info(kwargs)
            run_pipeline(data_dir,
                betas=True,
                save_control=True,
                poobah=True,
                export_poobah=True,
                quality_mask=True,
                batch_size=batch_size,
                make_sample_sheet=kwargs.get('make_sample_sheet',False),
                meta_data_frame=kwargs.get('meta_data_frame',False)
                ) #make_sample_sheet handled within miniml.py logic


def run_series_list(list_file, path, dict_only=False, batch_size=BATCH_SIZE, **kwargs):
    """Downloads the IDATs and metadata for a list of series, creating metadata dictionaries and dataframes of sample beta_values

    Arguments:
        list_file [required]
            the name of the file containing a list of GEO_IDS and/or Array Express IDs to download and process.
            This file must be located in the directory data is downloaded to.
            Each line of the file should contain the name of one data series ID.
        path [required]
            the path to the directory to download the data to. It is assumed a dictionaries and beta values
            directory has been created for each platform (and will create one for each if not)
        dict_only
            if True, only downloads data and creates dictionaries for each platform
        batch_size
            the batch_size to use when processing samples (number of samples run at a time).
            By default is set to the constant 100."""
    path = str(path)

    #if not os.path.exists(f"{path}/{PLATFORMS[0]}_beta_values"):
    #    initialize(str(path))

    try:
        fp = open(f"{path}/{str(list_file)}", 'r')
    except FileNotFoundError:
        LOGGER.error("""Specify your list of GEO series IDs to download using a text file in the folder where data should be saved. Put one ID on each line.""")
        return
    for series_id in fp:
        series_id = series_id.strip()
        series_path = Path(path, series_id) # run_series and geo_download get confused if idats already present, so this avoids that confusion
        try:
            series_path.mkdir(parents=True, exist_ok=True)
            LOGGER.info(f"Running {series_id}")
            run_series(series_id, series_path, dict_only=dict_only, batch_size=batch_size, **kwargs)
        except (ValueError, FileNotFoundError) as e:
            LOGGER.info(f"Error with {series_id}: {e}")
            with open("problem_series.txt", "a+") as fp:
                fp.write(f"{series_id} ({e})\n")
            fp.close()


def initialize(path):
    """Creates one directory for dictionaries and one directory for beta_values per platform

    Arguments:
        path [required]
            the path to the directory to create the platform directories"""
    if not Path(path).is_dir():
        #LOGGER.debug(f"Created {path} directory.")
        Path(path).mkdir(parents=True, exist_ok=True)
    for platform in PLATFORMS:
        #if not os.path.exists(f"{path}/{platform}_beta_values"):
        #    #LOGGER.debug(f"Created {platform} beta_values directory")
        #    os.mkdir(f"{path}/{platform}_beta_values")
        if not os.path.exists(f"{path}/{platform}_dictionaries"):
            #LOGGER.debug(f"Created {platform} dictionaries directory")
            os.mkdir(f"{path}/{platform}_dictionaries")

def confirm_dataset_contains_idats(geo_id):
    """ quickly scans the GEO accession viewer page for this dataset. if IDATs are mentioned, the file probably contains idats.
    Also - ensures that the geoxxx_RAW.ZIP file is large enough to contain data and not just manifest files."""
    geo_acc_page = f"http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={geo_id}"
    html = urlopen(geo_acc_page).read()
    idat = True if 'TAR (of IDAT)' in str(html) else False
    bigzip = False
    try:
        soup = BeautifulSoup(html, 'html.parser')
        table = [i for i in soup.find_all('table') if 'Supplementary file' in i.text]
        filesizes = [i for i in table[0].find_all('td') if 'Mb' in i.text or 'Gb' in i.text]
        # MB or GB?
        GB = True if len([i for i in table[0].find_all('td') if 'Gb' in i.text]) > 0 else False
        filesizes = [int(re.search(r'(\d+).*',i.text).group(1)) for i in filesizes if re.search(r'(\d+)',i.text)]
        if filesizes != []:
            if not GB and max(filesizes) > 195:
                bigzip = True
            elif GB and max(filesizes) > 0:
                bigzip = True
    except:
        pass
    if idat and bigzip:
        return True
    else:
        return False


def get_attachment_info(geo_id):
    """for a given GEO page, get file names, sizes, types, links"""
    geo_acc_page = f"http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={geo_id}"
    html = urlopen(geo_acc_page).read()
    info = [] # list of file data dicts
    try:
        soup = BeautifulSoup(html, 'html.parser')
        table = [i for i in soup.find_all('table') if 'Supplementary file' in i.text]
        if len(table) > 1:
            table = [table[-1]] # there are no meta tags for this table, but it should be the last one. having problems with nested tables getting stored.
            # alt method: findChildren( recursive=False )
        filesizes = [i for i in table[0].find_all('td') if 'Mb' in i.text or 'Gb' in i.text]
        # MB or GB?
        GB = True if len([i for i in table[0].find_all('td') if 'Gb' in i.text]) > 0 else False
        filesizes = [int(re.search(r'(\d+).*',i.text).group(1)) for i in filesizes if re.search(r'(\d+)',i.text)]
        if filesizes != []:
            for i,row in enumerate(table[0].find_all('tr')):
                if row.find_all('td')[0].text == 'Supplementary file':
                    continue
                filename = row.find_all('td')[0].text if len(row.find_all('td')) > 0 else ''
                filesize = row.find_all('td')[1].text if len(row.find_all('td')) > 1 else ''
                if len(row.find_all('td')) >= 2 and len(row.find_all('td')[2].find_all('a')) > 0:
                    filelink = row.find_all('td')[2].find_all('a')[0].get('href')
                    if filelink.startswith('/geo/'):
                        filelink = 'https://www.ncbi.nlm.nih.gov' + filelink
                else:
                    filelink = ''
                info.append({
                    'name': filename,
                    'size': filesize,
                    'link': filelink,
                })
    except Exception as e:
        LOGGER.error(f"Error parsing file data: {e}")
    return info


def cleanup(path):
    """removes unused/empty directories
    Arguments:
        path [required]
            the root path to check recursively"""
    if not Path(path).is_dir():
        raise ValueError(f"{path} doesn't exist")
    # _dictionaries are not needed after meta_data created.
    for platform in PLATFORMS:
        if Path(f"{path}/{platform}_dictionaries").is_dir():
            for file in Path(f"{path}/{platform}_dictionaries").rglob('*_dict.pkl'):
                file.unlink()
    #folders = [f"{path}/{platform}_beta_values" for platform in PLATFORMS]
    folders = [f"{path}/{platform}_dictionaries" for platform in PLATFORMS]
    folders.extend([f"{path}/{platform}" for platform in PLATFORMS]) # if no data, remove it.
    for folder in folders:
        if not Path(folder).is_dir():
            continue
        non_empty_dirs = {str(p.parent) for p in Path(folder).rglob('*') if p.is_file()}
        if non_empty_dirs == set():
            Path(folder).rmdir()
