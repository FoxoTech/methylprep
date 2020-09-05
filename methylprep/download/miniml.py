import re
from pathlib import Path
import shutil
from bs4 import BeautifulSoup
import pandas as pd
import tarfile
from tarfile import ReadError
from tqdm import tqdm
from ftplib import FTP
import logging
import time
# app
from .samplesheet_sync_idat import remove_idats_not_in_samplesheet
from methylprep import run_pipeline
#from methylprep.download.process_data import run_series --- breaks, dunno why!
import methylprep.download.process_data

LOGGER = logging.getLogger(__name__)
LOGGER.setLevel( logging.INFO )

geo_platforms=['GPL13534', 'GPL21145']

### FIX: make extract_controls FALSE by default, cli to turn on in meta_data ###

__all__ = [
    'convert_miniml',
    'download_miniml',
    'merge_sample_sheets',
    'meta_from_idat_filenames',
    'sample_sheet_from_miniml',
    'sample_sheet_from_idats',
    'cleanup',
    'build_composite_dataset',
]


def convert_miniml(geo_id, data_dir='.', merge=True, download_it=True, extract_controls=False, require_keyword=None, sync_idats=False, remove_tgz=False, verbose=False):
    """This scans the datadir for an xml file with the geo_id in it.
    Then it parses it and saves the useful stuff to a dataframe called "sample_sheet_meta_data.pkl".
    DOES NOT REQUIRE idats.

CLI version:
    python -m meta_data -i GSExxxxx -d <my_folder>

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
    local_files = []
    files = list(Path(data_dir).rglob(file_pattern))
    if files:
        file = Path(files[0])
        result = {'local_files': [file]}
    else:
        if download_it == True:
            # option to have it retain the .xml.tgz files, which -- if multipart -- contain each smaple betas
            result = download_miniml(geo_id, data_dir, remove_tgz=remove_tgz, verbose=verbose) ### methylprep v1.3.0 will change value return from string to list. Only used here within methylprep.
            file = result['meta_data']
            # also returns a list of local_files found, if there were samples tgz'ed in with meta_data
            local_files = result['meta_data']
            if verbose:
                LOGGER.info(f"Attempting to read xml: {file}")
        else:
            LOGGER.info('Did not find any miniml XML files to convert.')
            return local_files
    with open(file, 'r') as fp:
        soup = BeautifulSoup(fp, "xml")
    meta_dicts = {'GPL13534':{}, 'GPL21145':{}}
    samples_dict = {'GPL13534':{}, 'GPL21145':{}}
    samples = soup.MINiML.find_all("Sample")
    missing_methylprep_names = 0
    warned_once = False
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
        elif warned_once is False:
            LOGGER.warning(f"{platform} is not among the supported platforms {geo_platforms}")
            warned_once = True
    if missing_methylprep_names > 0:
        LOGGER.info( f"MINiML file does not provide (Sentrix_id_R00C00) for {missing_methylprep_names}/{len(samples)} samples." )

    for platform in geo_platforms:
        if verbose:
            LOGGER.info(f"Platform {platform} --> {len(meta_dicts[platform])} samples")
        if meta_dicts[platform]:
            #Path(f"{data_dir}/{platform}").mkdir(parents=True, exist_ok=True)
            meta_dicts[platform] = meta_from_idat_filenames(data_dir, meta_dicts[platform])
            if merge == True:
                # find samplesheet.csv with idat ID/Position and merge in here.
                meta_dicts[platform] = merge_sample_sheets(data_dir, meta_dicts[platform])
            sample_sheet_from_miniml(geo_id, data_dir, platform, samples_dict[platform], meta_dicts[platform],
                save_df=True, extract_controls=extract_controls, require_keyword=require_keyword, verbose=verbose)
            if sync_idats:
                paths = list(Path(data_dir).rglob(f'{geo_id}_{platform}_samplesheet.csv'))
                if len(paths) > 0:
                    remove_idats_not_in_samplesheet(paths[0], data_dir)
                else:
                    LOGGER.warning(f"Could not locate file ({f'{geo_id}_{platform}_samplesheet.csv'}) in path of {data_dir}.")
    cleanup(data_dir)
    return result['local_files']


def download_miniml(geo_id, series_path, remove_tgz=False, verbose=False):
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
        like "GSE134293" in your CLI parameters.
    Note on alt filenames:
        GSE50660_family.xml was actually three files:
            GSE50660_family.xml.part3.tgz
            GSE50660_family.xml.part2.tgz
            GSE50660_family.xml.tgz
        so v1.3.0 adds support for multi-part miniml files.
            (these often include a lot of additional processed idats data)

note: returns a LIST instead of a single filename now. update convert_miniml() to use that.
        """
    series_dir = Path(series_path)
    expected_miniml_filename = f"{geo_id}_family.xml"

    if not Path(series_path).exists():
        Path(series_path).mkdir()

    ftp = FTP('ftp.ncbi.nlm.nih.gov', timeout=120) # 2 mins
    ftp.login()
    ftp.cwd(f"geo/series/{geo_id[:-3]}nnn/{geo_id}")

    # look for additional parts of xml files, so long as each part contains the geo_id.
    file_parts = []
    for filename,filestats in ftp.mlsd(path=f"miniml", facts=["size"]):
        if geo_id in filename:
            file_parts.append(filename)

    local_files = []
    for miniml_filename in file_parts:
        LOGGER.info(f"Downloading {miniml_filename}")
        miniml_ftp_path = f"miniml/{miniml_filename}"
        miniml_local_path = f"{series_path}/{miniml_filename}"
        miniml_file = open(miniml_local_path, 'wb')
        try:
            # Alternatively, as a hacky workaround, you might try download a file -- any file -- before sending your SIZE command. With request.UseBinary = true for that request, it should cause your client to send the "TYPE I" command to the FTP server. (And it won't matter if that download request fails; the TYPE command will still have been sent.)
            filesize = ftp.size(miniml_ftp_path)
            with tqdm(unit = 'b', unit_scale = True, leave = False, miniters = 1, desc = geo_id, total = filesize) as tqdm_instance:
                def tqdm_callback(data):
                    tqdm_instance.update(len(data))
                    miniml_file.write(data)
                ftp.retrbinary(f"RETR {miniml_ftp_path}", tqdm_callback)
        except Exception as e:
            if not str(e).startswith('550'):
                LOGGER.info(e)
            if verbose:
                LOGGER.info('Failed to create a progress bar, but it is downloading...')
            ftp.retrbinary(f"RETR {miniml_ftp_path}", miniml_file.write)
        miniml_file.close()
        if not Path(miniml_local_path).exists():
            LOGGER.error(f"{miniml_local_path} missing")
        min_tar = tarfile.open(miniml_local_path)
        for file in min_tar.getnames(): # this only extracts the one .xml file we need.
            #LOGGER.info(f"DEBUG miniml_filename {expected_miniml_filename} == tar_file {file} | {file == expected_miniml_filename}")
            if file == expected_miniml_filename:
                min_tar.extract(file, path=series_path)
        if Path(miniml_local_path).exists() and remove_tgz:
            Path(miniml_local_path).unlink()
        if verbose:
            LOGGER.info(f"Downloaded and unpacked {miniml_filename}")

    for extracted_file in Path(series_path).glob('*'):
        if verbose:
            LOGGER.info(f"Extracted {extracted_file}")
        if expected_miniml_filename == extracted_file.name:
            miniml_local_path = f"{series_path}/{extracted_file.name}"
            if verbose:
                LOGGER.info(f"*Identified* meta_data: {extracted_file}")
        local_files.append(str(extracted_file))
    ftp.quit()
    return {'meta_data': miniml_local_path, 'local_files': local_files}


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


def sample_sheet_from_miniml(geo_id, series_path, platform, samp_dict, meta_dict, save_df=False, extract_controls=False, require_keyword=None, verbose=True):
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
            LOGGER.debug(f"dropped {column} (empty)")
            out.pop(column,None)
        if set(_dict[column]) & set(geo_platforms) != set():
            LOGGER.debug(f"dropped `{column}` ({set(_dict[column]) & set(geo_platforms)})")
            out.pop(column,None) # don't need platform, since saved in folder this way.
        if column in ('title','source'):
            # drop these columns first, if redundant
            for other_column in _dict.keys():
                if set(_dict[column]) == set(_dict[other_column]) and column != other_column:
                    LOGGER.debug(f"{column} == {other_column}; dropping {column}")
                    out.pop(column,None)
    try:
        df = pd.DataFrame(data=out)
    except ValueError as e: # arrays must all be same length
        from collections import Counter
        LOGGER.info(f"ValueError - array lengths vary in sample meta data: {[(key, len(val)) for key,val in out.items()]}")
        LOGGER.info(f"Trying another method to save partial sample data into csv/pkl file.")
        ## alt method to salvage it by filling in blanks for missing rows -- but doesn't seem to capture Sentrix_ID / Sentrix_Position ##
        column_counts = {'GSM_ID': [], 'Sample_Name': []} # column_name : [GSM_IDs included]
        out = {} # list of dicts as sample rows, keyed to GSM_IDs
        for GSM_ID, sample_title in samp_dict.items():
            out[GSM_ID] = {'GSM_ID': GSM_ID, 'Sample_Name': sample_title}
            column_counts['GSM_ID'].append(GSM_ID)
            column_counts['Sample_Name'].append(GSM_ID)
            for key, value in meta_dict[GSM_ID].items():
                out[GSM_ID][key] = value
                if key not in column_counts:
                    column_counts[key] = [GSM_ID]
                else:
                    column_counts[key].append(GSM_ID)
        # check if any fields are missing GSM_IDs, and fill in with blanks
        for col_name, gsm_id_list in column_counts.items():
            if len(out) != len(gsm_id_list):
                for GSM_ID in out.keys():
                    if GSM_ID not in gsm_id_list:
                        out[GSM_ID][col_name] = ''
        try:
            df = pd.DataFrame(data=list(out.values()))
        except:
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
    if len(df.columns) < 30 and verbose:
        LOGGER.info(f"{list(df.columns)}")
    if Path(series_path, platform).exists():
        df.to_csv(path_or_buf=(Path(series_path, platform, f'{geo_id}_{platform}_samplesheet.csv')),index=False)
    else:
        df.to_csv(path_or_buf=(Path(series_path, f'{geo_id}_{platform}_samplesheet.csv')),index=False)
    if save_df:
        if Path(series_path, platform).exists():
            df.to_pickle(Path(series_path, platform, f'{geo_id}_{platform}_meta_data.pkl'))
        else:
            df.to_pickle(Path(series_path, f'{geo_id}_{platform}_meta_data.pkl'))


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


def build_composite_dataset(geo_id_list, data_dir, merge=True, download_it=True, extract_controls=False, require_keyword=None, sync_idats=True, betas=False, m_value=False, export=False):
    """A wrapper function for convert_miniml() to download a list of GEO datasets
    and process only those samples that meet criteria. Specifically - grab the "control" or "normal" samples
    from a bunch of experiments for one tissue type (e.g. "blood"), process them, and put all the resulting
    beta_values and/or m_values pkl files in one place, so that you can run `methylize.load_both()` to
    create a combined reference dataset for QC, analysis, or meta-analysis.

    Arguments:
        geo_id_list (required):
            A list of GEO "GSEnnn" ids. From command line, pass these in as separate values
        data_dir:
            folder to save data
        merge (True):
            If merge==True and there is a file with 'samplesheet' in the folder, and that sheet has GSM_IDs,
            merge that data into this samplesheet. Useful for when you have idats and want one combined samplesheet for the dataset.

        download_it (True):
            if miniml file not in data_dir path, it will download it from web.

        extract_controls (True)):
            if you only want to retain samples from the whole set that have certain keywords,
            such as "control" or "blood", this experimental flag will rewrite the samplesheet with only the parts you want,
            then feed that into run_pipeline with named samples.
        require_keyword (None):
            another way to eliminate samples from samplesheets, before passing into the processor.
            if specified, the "keyword" string passed in must appear somewhere in the values of a samplesheet
            for sample to be downloaded, processed, retained.
        sync_idats:
            If flagged, this will search `data_dir` for idats and remove any of those that are not found in the filtered samplesheet.
            Requires you to run the download function first to get all idats, before you run this `meta_data` function.
        betas:
            process beta_values
        m_value:
            process m_values

        - Attempts to also read idat filenames, if they exist, but won't fail if they don't.
        - removes unneeded files as it goes, but leaves the xml MINiML file and folder there as a marker if a geo dataset fails to download. So it won't try again on resume.
    """
    def remove_unused_files(geo_id, geo_folder):
        if list(Path(data_dir, geo_id).rglob('*.idat')) == []:
            for file in Path(data_dir, geo_id).glob("*_samplesheet.csv"):
                file.unlink()
            for file in Path(data_dir, geo_id).glob("*_meta_data.pkl"):
                file.unlink()
            # the XML and folder is used to mark failed downloads on resume, so it skips them next time
            #for file in Path(data_dir, geo_id).glob("*_family.xml"):
            #    file.unlink()
            #try:
            #    Path(data_dir, geo_id).rmdir()
            #except Exception as e:
            #    LOGGER.error(f"Path {data_dir/geo_id} is not empty. Could not remove.")
            return True
        return False
    start_time = time.process_time()
    # note: parser uses VERBOSE setting to show/suppress INFO and DEBUG level messages. WARNING/ERROR msgs always appear.
    # get the ids from file
    try:
        with open(Path(data_dir,geo_id_list), 'r') as fp:
            geo_ids = [series_id.strip() for series_id in fp]
    except FileNotFoundError:
        LOGGER.error("""File not found: Specify your list of GEO series IDs to download using a text file in the folder where data should be saved. Put one ID on each line. """)
        return
    except ValueError as e:
        LOGGER.error(f"Error with {fp.name}: {e}")
        fp.close()
        return
    fp.close()

    geo_folders = {}
    for geo_id in geo_ids:
        # exclude failed folders: if the folder exists and only contains the miniml family.xml file, skip it.
        if Path(data_dir,geo_id).exists() and Path(data_dir,geo_id, f'{geo_id}_family.xml').exists():
            if len(list(Path(data_dir, geo_id).rglob('*'))) == 1:
                LOGGER.info(f"Skipping {geo_id}; appears to be a prior run that didn't match filters, or was missing data.")
                continue
        # exclude geo series whose HTML pages don't say TAR (of idat).
        if methylprep.download.process_data.confirm_dataset_contains_idats(geo_id) == False:
            LOGGER.error(f"[!] Geo data set {geo_id} probably does NOT contain usable raw data (in .idat format). Not downloading.")
            continue

        geo_folder = Path(data_dir,geo_id)
        geo_folders[geo_id] = geo_folder

        # download meta_data
        LOGGER.info(f"Running {geo_id}")
        try:
            convert_miniml(
                geo_id,
                data_dir=geo_folder,
                merge=merge,
                download_it=download_it,
                extract_controls=extract_controls,
                require_keyword=require_keyword,
                sync_idats=False) #no idat files exist yet.
        except Exception as e:
            LOGGER.error(f'Processing meta_data failed: {e}')
            continue
        ## if the samplesheet is empty, stop.
        abort = False
        for platform in geo_platforms:
            if Path(data_dir, geo_id, f"{geo_id}_{platform}_samplesheet.csv").is_file():
                samplesheet = pd.read_csv(Path(data_dir, geo_id, f"{geo_id}_{platform}_samplesheet.csv"))
                if len(samplesheet.index) == 0:
                    LOGGER.warning(f"Aborting {geo_id}: No samples match filters (control:{extract_controls}, keyword: {require_keyword})")
                    # TODO: cleanup this first.
                    abort = True
                    remove_unused_files(geo_id, geo_folder)
                    geo_folders.pop(geo_id)
                    continue
        if abort:
            continue

        # check if no idats yet
        if list(Path(geo_folder).rglob('*.idat')) != []:
            pass
        else:
            try:
                # download idats
                methylprep.download.process_data.run_series(
                    geo_id,
                    geo_folder,
                    dict_only=True,
                    batch_size=200,
                    clean=True, # do later in this function
                    )
            except Exception as e:
                LOGGER.error(f'Downloading IDATs failed: {e}')
                # THIS DOES NOT CLEAN THE FOLDERS IF ERROR
                continue

        if sync_idats:
            for platform in geo_platforms:
                paths = list(Path(data_dir, geo_id).rglob(f'{geo_id}_{platform}_samplesheet.csv'))
                if len(paths) > 0:
                    remove_idats_not_in_samplesheet(paths[0], Path(data_dir, geo_id))
                #else:
                #    LOGGER.error(f"Could not locate file ({geo_id}_{platform}_samplesheet.csv}) in path of {data_dir}/{geo_id}")

        # there are probably two versions of samplesheets. if they're the same, remove one.
        # if they differ, keep the one created by miniml (where filtering happens)
        # this should only proceed this far if both samplesheets contain samples.
        for platform in geo_platforms:
            samp_A = Path(data_dir, geo_id, f"{geo_id}_{platform}_samplesheet.csv")
            samp_B = Path(data_dir, geo_id, platform, f"{geo_id}_{platform}_samplesheet.csv")
            meta_B = Path(data_dir, geo_id, platform, f"{geo_id}_{platform}_meta_data.pkl")
            if (samp_A.is_file() and samp_B.is_file()):
                samplesheet = pd.read_csv(samp_A)
                basic_samplesheet = pd.read_csv(samp_B)
                #if samplesheet == basic_samplesheet: OR in every case, keep the miniml-filtered version.
                #shutil.move(samp_A, samp_B)
                samp_B.unlink()
                if meta_B.is_file():
                    meta_B.unlink()
        # next, remove all files if there are no idats in the folder.
        if remove_unused_files(geo_id, geo_folder):
            geo_folders.pop(geo_id)


    if export == False and betas == False and m_value == False:
        LOGGER.info("Not processing data, because no output types are specified.")
    else:
        for geo_id, geo_folder in geo_folders.items():
            try:
                run_pipeline(geo_folder, #maybe pass in sample_sheet_filepath if it gets confused here.
                    betas=betas,
                    m_value=m_value,
                    batch_size=200,
                    export=export,
                    meta_data_frame=False
                    )
            except Exception as e:
                LOGGER.warning(f'Processing IDATs failed for {geo_id}: {e}')

    LOGGER.info('[!] Consolidating data files [!]')
    # consoldate data into one folder and remove rest.
    file_patterns = [
        'beta_values_*.pkl' if betas else None,
        'm_values_*.pkl' if m_value else None
        ]
    for pattern in file_patterns:
        if pattern == None:
            continue
        datas = []
        samples = 0
        for file in Path(data_dir).rglob(pattern):
            if file.parts[-2] in geo_ids:
                data = pd.read_pickle(file)
                # data should have probes in rows, samples in columns.
                if data.shape[1] > data.shape[0]:
                    data = data.transpose()
                datas.append(data)
                samples += data.shape[1]
                print(f"-- {len(datas)} {file.parts[-2]} {data.shape}")
        if datas:
            big_data = pd.concat(datas, axis=1, ignore_index=False, sort=True, copy=False)
            del datas
            LOGGER.info(f"[!] saved {samples} samples to disk from {pattern}.")
            if 'beta' in pattern:
                big_data.to_pickle(Path(data_dir,'beta_values.pkl'))
            if 'm_value' in pattern:
                big_data.to_pickle(Path(data_dir,'m_values.pkl'))
    end_time = time.process_time()
    LOGGER.info(f"[*] Process time: {round((start_time - end_time)/60.0,1)} min")
