import re
from pathlib import Path
from bs4 import BeautifulSoup
import pandas as pd
import tarfile
from tarfile import ReadError
from tqdm import tqdm
from ftplib import FTP
import logging
# app
from .samplesheet_sync_idat import remove_idats_not_in_samplesheet

LOGGER = logging.getLogger(__name__)
LOGGER.setLevel( logging.INFO )

geo_platforms=['GPL13534', 'GPL21145']

### FIX: make extract_controls FALSE by default, cli to turn on in meta_data ###

def convert_miniml(geo_id, data_dir='.', merge=True, download_it=True, extract_controls=False, require_keyword=None, sync_idats=False):
    """This scans the datadir for an xml file with the geo_id in it.
    Then it parses it and saves the useful stuff to a dataframe called "sample_sheet_meta_data.pkl".
    DOES NOT REQUIRE idats.

Arguments:
    merge:
        If merge==True and there is a file with 'samplesheet' in the folder, and that sheet has GSM_IDs,
        merge that data into this samplesheet. Useful for when you have idats and want one combined samplesheet for the dataset.

    download_it:
        if miniml file not in data_dir path, it will download it from web.

    extract_controls [experimental]:
        if you only want to retain samples from the whole set that have certain keywords,
        such as "control" or "blood", this experimental flag will rewrite the samplesheet with only the parts you want,
        then feed that into run_pipeline with named samples.
    require_keyword [experimental]:
        another way to eliminate samples from samplesheets, before passing into the processor.
        if specified, this keyword needs to appear somewhere in the values of a samplesheet.
    sync_idats:
        If flagged, this will search `data_dir` for idats and remove any of those that are not found in the filtered samplesheet.
        Requires you to run the download function first to get all idats, before you run this `meta_data` function.

    Attempts to also read idat filenames, if they exist, but won't fail if they don't.
    """
    file_pattern = f'{geo_id}*.xml'
    files = list(Path(data_dir).rglob(file_pattern))
    if files:
        file = Path(files[0])
    else:
        if download_it == True:
            file = download_miniml(geo_id, data_dir)
        else:
            print('did not find any miniml XML files to convert.')
            return
    with open(file, 'r') as fp:
        soup = BeautifulSoup(fp, "xml")
    meta_dicts = {'GPL13534':{}, 'GPL21145':{}}
    samples_dict = {'GPL13534':{}, 'GPL21145':{}}
    samples = soup.MINiML.find_all("Sample")
    missing_methylprep_names = 0
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
            attributes_dir['Sample_ID'] = f"{split_idat[1]}_{split_idat[2]}" # {accession / GSM} is not used in methylprep data columns
            attributes_dir['Sentrix_ID'] = f"{split_idat[1]}"
            attributes_dir['Sentrix_Position'] = f"{split_idat[2]}"
        except:
            missing_methylprep_names += 1
        if platform in geo_platforms:
            meta_dicts[platform][accession] = attributes_dir
            samples_dict[platform][accession] = title
    if missing_methylprep_names > 0:
        LOGGER.info( f"MINiML file does not provide `methylprep_name` (sentrix_id_R00C00) for {missing_methylprep_names}/{len(samples)} samples." )

    for platform in geo_platforms:
        if meta_dicts[platform]:
            #Path(f"{data_dir}/{platform}").mkdir(parents=True, exist_ok=True)
            meta_dicts[platform] = meta_from_idat_filenames(data_dir, meta_dicts[platform])
            if merge == True:
                # find samplesheet.csv with idat ID/Position and merge in here.
                meta_dicts[platform] = merge_sample_sheets(data_dir, meta_dicts[platform])
            sample_sheet_from_miniml(geo_id, data_dir, platform, samples_dict[platform], meta_dicts[platform],
                save_df=True, extract_controls=extract_controls, require_keyword=require_keyword)
            if sync_idats:
                paths = list(Path(data_dir).rglob(f'{geo_id}_{platform}_samplesheet.csv'))
                if len(paths) > 0:
                    remove_idats_not_in_samplesheet(paths[0], data_dir)
                else:
                    LOGGER.warning(f"Could not locate file ({f'{geo_id}_{platform}_samplesheet.csv'}) in path of {data_dir}.")
    cleanup(data_dir)

def download_miniml(geo_id, series_path):
    """Downloads the MINIML metadata for a GEO series

    Arguments:
        geo_id [required]
            the GEO Accension for the desired series (e.g. GSE134293)
        series_path [required]
            the directory to download the data to

    Note about GEO IDs:
        You can use the NIH online search to find data sets, then click "Send to:" at the button of a results page,
        and export a list of unique IDs as text file. These IDs are not GEO_IDs used here. First, remove the first
        three digits from the number, so Series ID: 200134293 is GEO accension ID: 134293, then include the GSE part,
        like "GSE134293" in your CLI parameters. """
    series_dir = Path(series_path)
    miniml_filename = f"{geo_id}_family.xml"

    if not Path(series_path).exists():
        Path(series_path).mkdir()

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
        LOGGER.info(f"Unpacking {miniml_filename}")
        if not Path(f"{series_path}/{miniml_filename}.tgz").exists():
            print( f"{series_path}/{miniml_filename}.tgz missing" )
        min_tar = tarfile.open(f"{series_path}/{miniml_filename}.tgz")
        for file in min_tar.getnames():
            if file == miniml_filename:
                min_tar.extract(file, path=series_path)
            if Path(f"{series_path}/{miniml_filename}.tgz").exists():
                Path(f"{series_path}/{miniml_filename}.tgz").unlink()
    LOGGER.info(f"Downloaded and unpacked {geo_id}")
    ftp.quit()
    return f"{series_path}/{miniml_filename}"


def merge_sample_sheets(data_dir, meta_dict):
    # potentially grab Sentrix_ ID and Position from samplesheet of idats.
    files = list(Path(data_dir).rglob(f'*samplesheet*.csv'))
    if files:
        for file in files:
            try:
                samplesheet = pd.read_csv(file)
            except:
                LOGGER.warning('could not open {file}')
                continue
            if 'GSM_ID' in samplesheet:
                if set(samplesheet['GSM_ID']) & set(meta_dict) != set():
                    # meta_dict is a dict of dicts with GSM_IDS as keys
                    # add lists of Sentrix_ID and Sentrix_Position to meta_dict.
                    gsm_sid = {row.GSM_ID:row.Sentrix_ID for idx,row in samplesheet.iterrows() if row.get('Sentrix_ID')}
                    gsm_spos = {row.GSM_ID:row.Sentrix_Position for idx,row in samplesheet.iterrows() if row.get('Sentrix_Position')}
                    updated = 0
                    not_found = 0
                    for GSM_ID in samplesheet['GSM_ID']:
                        if GSM_ID in meta_dict and GSM_ID in gsm_sid and GSM_ID in gsm_spos:
                            if not meta_dict[GSM_ID].get('Sentrix_ID'):
                                meta_dict[GSM_ID]['Sentrix_ID'] = gsm_sid[GSM_ID]
                                updated += 1
                            if not meta_dict[GSM_ID].get('Sentrix_Position'):
                                meta_dict[GSM_ID]['Sentrix_Position'] = gsm_spos[GSM_ID]
                        else:
                            not_found += 1
                    if updated > 0:
                        LOGGER.info(f'{file}: meta_dict updated with {updated} Sentrix IDs/Positions')
                    if not_found > 0:
                        LOGGER.info(f'{file}: {not_found} GSM_IDs in samplesheet where not found in meta_dict')
    return meta_dict


def meta_from_idat_filenames(data_dir, meta_dict):
    gsm_id_pos = {}
    for idx, idat in enumerate(Path(data_dir).rglob('*Grn.idat')):
        if len(idat.stem.split("_")) == 4: # and 'Grn' in idat.stem:
            gsm_id = idat.stem.split("_")[0]
            sentrix_id = idat.stem.split("_")[1]
            sentrix_position = idat.stem.split("_")[2]
            gsm_id_pos[gsm_id] = {'Sentrix_ID': sentrix_id, 'Sentrix_Position': sentrix_position}
        elif len(idat.stem.split("_")) == 3: # and 'Grn' in idat.stem:
            gsm_id = idx
            sentrix_id = idat.stem.split("_")[0]
            sentrix_position = idat.stem.split("_")[1]
            gsm_id_pos[gsm_id] = {'Sentrix_ID': sentrix_id, 'Sentrix_Position': sentrix_position}

    if len(gsm_id_pos) == len(meta_dict):
        for gsm_id, data in gsm_id_pos.items():
            if gsm_id in meta_dict:
                meta_dict[gsm_id]['Sentrix_ID'] = data['Sentrix_ID']
                meta_dict[gsm_id]['Sentrix_Position'] = data['Sentrix_Position']
        LOGGER.info(f"found {len(gsm_id_pos)} idat files; updated meta_data.")
    return meta_dict


def sample_sheet_from_miniml(geo_id, series_path, platform, samp_dict, meta_dict, save_df=False, extract_controls=False, require_keyword=None):
    """Creates a sample_sheet for all samples of a particular platform for a given series
    Arguments:
        geo_id [required]
            the GEO Accension for the desired series
        series_path [required]
            the directory containing the series data
        platform [required]
            the platform to generate a sample sheet for
        samp_dict
            the dictionary of {sample GSM_ID: title} for the given platform
        meta_dict
            the dictionary of {sample GSM_ID: meta_dict} for the given platform

    Notes:
        cleans up redundant columns, including Sample_Name == title and platform
            """
    # id and position are NOT in the miniml file, must match with GSMID later
    _dict = {'GSM_ID': [], 'Sample_Name': [], 'Sentrix_ID': [], 'Sentrix_Position': []}
    for GSM_ID, sample_title in samp_dict.items():
        _dict['GSM_ID'].append(GSM_ID)
        _dict['Sample_Name'].append(sample_title)
        for key, value in meta_dict[GSM_ID].items():
            if key not in _dict:
                _dict[key] = []
            _dict[key].append(value)
    # arrays must be all same length
    out = _dict.copy()
    for column in _dict.keys():
        # Sentrix_ID and Sentrix_Position may be missing.
        if len(_dict[column]) == 0:
            LOGGER.info(f"dropped {column} (empty)")
            out.pop(column,None)
        if set(_dict[column]) & set(geo_platforms) != set():
            LOGGER.info(f"dropped `{column}` ({set(_dict[column]) & set(geo_platforms)})")
            out.pop(column,None) # don't need platform, since saved in folder this way.
        if column in ('title','source'):
            # drop these columns first, if redundant
            for other_column in _dict.keys():
                if set(_dict[column]) == set(_dict[other_column]) and column != other_column:
                    LOGGER.info(f"{column} == {other_column}; dropping {column}")
                    out.pop(column,None)

    try:
        df = pd.DataFrame(data=out)
    except ValueError as e: # arrays must all be same length
        from collections import Counter
        LOGGER.info(f"ValueError - array lengths vary in sample meta data: {[(key, len(val)) for key,val in out.items()]}")
        ## would be HARD to salvage it by filling in blanks for missing rows ##
        raise ValueError(f"{e}; this happens when a samplesheet is missing descriptors for one or more samples.")
    # filter: only retain control samples
    if extract_controls:
        keywords = ['control','ctrl']
        filtered = df.copy()
        for idx, row in df.iterrows():
            # some column needs to match one of the keywords to be retained
            tested = []
            for keyword in keywords:
                for val in row.values:
                    if re.search(keyword, val, re.I):
                        tested.append(val)
            if tested == []:
                filtered.drop(idx, inplace=True)
        LOGGER.info(f"Filtering controls: Retained {len(filtered.index)}/{len(df.index)} samples matching {keywords}.")
        del df
        df = filtered
    if require_keyword:
        filtered = df.copy()
        for idx, row in df.iterrows():
            # some column value needs to match this keyword for sample to be retained
            tested = []
            for val in row.values:
                if re.search(require_keyword, val, re.I):
                    tested.append(val)
            if tested == []:
                filtered.drop(idx, inplace=True)
        LOGGER.info(f"Filtering keyword: Retained {len(filtered.index)}/{len(df.index)} samples matching `{require_keyword}``.")
        del df
        df = filtered

    #LOGGER.info(df.head())
    LOGGER.info(f"Final samplesheet contains {df.shape[0]} rows and {df.shape[1]} columns")
    if len(df.columns) < 30:
        LOGGER.info(f"{list(df.columns)}")
    if Path(series_path, platform).exists():
        df.to_csv(path_or_buf=(Path(series_path, platform, f'{geo_id}_{platform}_samplesheet.csv')),index=False)
    else:
        df.to_csv(path_or_buf=(Path(f"{series_path}", f'{geo_id}_{platform}_samplesheet.csv')),index=False)
    if save_df:
        if Path(series_path, platform).exists():
            df.to_pickle(Path(series_path, platform, f'{geo_id}_{platform}_meta_data.pkl'))
        else:
            df.to_pickle(Path(f"{series_path}", f'{geo_id}_{platform}_meta_data.pkl'))


def sample_sheet_from_idats(geo_id, series_path, platform, platform_samples_dict, save_df=False):
    """This is a simpler "fallback" parser of miniml.sample_sheet_from_miniml().

    Creates a sample_sheet for all samples of a particular platform for a given series

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
    df.to_csv(path_or_buf=(Path(f"{series_path}/{platform}", 'samplesheet.csv')),index=False)
    if save_df:
        df.to_pickle(Path(f"{series_path}/{platform}", 'sample_sheet_meta_data.pkl'))


def cleanup(path):
    """removes unused/empty directories
    Arguments:
        path [required]
            the root path to check recursively"""
    if not Path(path).is_dir():
        raise ValueError(f"{path} doesn't exist")
    folders = [p for p in list(Path(path).glob('*')) if p.is_dir()]
    for folder in folders:
        filled_folders = {str(p.parent) for p in Path(folder).rglob('*') if p.is_file()}
        if filled_folders == set():
            Path(folder).rmdir()
