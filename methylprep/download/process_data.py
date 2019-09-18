# Lib
import os
import logging
from pathlib import Path, PurePath
import shutil
import pandas as pd
# App
from .geo import (
    geo_download,
    geo_metadata
)
from .array_express import(
    ae_download,
    ae_metadata
)
import methylprep


LOGGER = logging.getLogger(__name__)

GEO_PLATFORMS = ['GPL21145', 'GPL13534']
AE_PLATFORMS = ['A-MEXP-2255', 'A-GEOD-21145']
PLATFORMS = GEO_PLATFORMS + AE_PLATFORMS
BATCH_SIZE = 100


def run_series(id, path, dict_only=False, batch_size=BATCH_SIZE, clean=True):
    """Downloads the IDATs and metadata for a series then generates one metadata dictionary and one beta value matrix for each platform in the series

    Arguments:
        id [required]
            the series ID (can be a GEO or ArrayExpress ID)
        path [required]
            the path to the directory to download the data to. It is assumed a dictionaries and beta values
            directory has been created for each platform (and will create one for each if not)
        dict_only
            if True, only downloads data and creates dictionaries for each platform
        batch_size
            the batch_size to use when processing samples (number of samples run at a time).
            By default is set to the constant 100."""
    if not os.path.exists(f"{str(path)}/{PLATFORMS[0]}_beta_values"):
        initialize(str(path))

    path = str(path)
    series_path = f"{path}/{id}"
    series_dir = Path(series_path)

    if not os.path.exists(series_path):
        LOGGER.info(f"Creating directory for {id}")
        os.mkdir(series_path)

    if id[:3] == 'GSE':
        series_type = 'GEO'
        geo_download(id, series_path, GEO_PLATFORMS, clean=clean)
    elif id[:7] == 'E-MTAB-':
        series_type = 'AE'
        ae_download(id, series_path, AE_PLATFORMS, clean=clean)
    else:
        raise ValueError(f"[ERROR] Series type not recognized. (The ID should begin with GSE or E-MTAB-)")

    dicts = list(series_dir.glob('*_dict.pkl'))
    if not dicts:
        if series_type == 'GEO':
            seen_platforms = geo_metadata(id, series_path, GEO_PLATFORMS, str(path))
        elif series_type == 'AE':
            seen_platforms = ae_metadata(id, series_path, AE_PLATFORMS, str(path))
    else:
        seen_platforms = []
        for d in list(dicts):
            seen_platforms.append(str(d).split("/")[-1].split("_")[1])

    if not dict_only:
        process_series(id, str(path), seen_platforms, batch_size)


def process_series(id, path, seen_platforms, batch_size):
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
        if not os.path.exists(f"{path}/{platform}_beta_values/{id}_beta_values.pkl"):
            LOGGER.info(f"Processing {id} {platform} samples")
            methylprep.run_pipeline(f"{path}/{id}/{platform}", betas=True, batch_size=batch_size)

            dir = Path('./')
            dfs = []
            betas = dir.glob('beta_values_*.pkl')
            betas_list = list(betas)

            for i in range(1,len(betas_list) + 1):
                df = pd.read_pickle(f"beta_values_{str(i)}.pkl")
                dfs.append(df)
            if len(dfs) > 1:
                joined_df = pd.concat(dfs, axis=1)
            else:
                joined_df = dfs[0]

            joined_df.to_pickle(f"{path}/{platform}_beta_values/{id}_beta_values.pkl")
            for beta in betas_list:
                os.remove(beta)
            LOGGER.info(f"Consolidated {id} {platform} samples, saved to {id}_beta_values.pkl")


def run_series_list(list_file, path, dict_only=False, batch_size=BATCH_SIZE):
    """Downloads the IDATs and metadata for a list of series, creating metadata dictionaries and dataframes of sample beta_values

    Arguments:
        list_file [required]
            the name of the file containing the series to download and process.
            This file must be located in the directory data is downloaded to (path).
            Each line of the file contains the name of one series.
        path [required]
            the path to the directory to download the data to. It is assumed a dictionaries and beta values
            directory has been created for each platform (and will create one for each if not)
        dict_only
            if True, only downloads data and creates dictionaries for each platform
        batch_size
            the batch_size to use when processing samples (number of samples run at a time).
            By default is set to the constant 100."""
    path = str(path)

    if not os.path.exists(f"{path}/{PLATFORMS[0]}_beta_values"):
        initialize(str(path))

    fp = open(f"{path}/{str(list_file)}", 'r')
    for series_id in fp:
        try:
            LOGGER.info(f"Running {series_id.strip()}")
            run_series(series_id.strip(), path, dict_only=dict_only, batch_size=batch_size)
        except (ValueError, FileNotFoundError) as e:
            LOGGER.info(f"Error with {series_id.strip()}: {e}")
            with open("problem_series.txt", "a") as fp:
                fp.write(f"{series_id.strip()} ({e})\n")
            fp.close()


def initialize(path):
    """Creates one directory for dictionaries and one directory for beta_values per platform

    Arguments:
        path [required]
            the path to the directory to create the platform directories"""

    if not os.path.exists(path):
        os.mkdir(f"{path}")
        LOGGER.info(f"Created the {path} directory.")

    for platform in PLATFORMS:
        if not os.path.exists(f"{path}/{platform}_beta_values"):
            LOGGER.info(f"Creating {platform} beta_values directory")
            os.mkdir(f"{path}/{platform}_beta_values")

        if not os.path.exists(f"{path}/{platform}_dictionaries"):
            LOGGER.info(f"Creating {platform} dictionaries directory")
            os.mkdir(f"{path}/{platform}_dictionaries")

        if not os.path.exists(f"{path}/problem_series.txt"):
            open('problem_series.txt', 'a').close()
