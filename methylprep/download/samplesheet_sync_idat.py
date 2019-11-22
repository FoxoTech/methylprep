from pathlib import Path
import logging
import pandas as pd
# app

LOGGER = logging.getLogger(__name__)
LOGGER.setLevel( logging.INFO )

def remove_idats_not_in_samplesheet(samplesheet_filepath, sample_path):
    """ compares all idats against a whitelist of those in samplesheet, removing the rest from sample_path.
    idats can end with .idat.gz or .idat """
    samples = pd.read_csv(samplesheet_filepath)
    all_idats = list(Path(sample_path).rglob('*.idat')) + list(Path(sample_path).rglob('*.idat.gz'))
    all_idats_names = [i.name for i in all_idats]
    # these are VALID idats to retain
    save_list = []
    try:
        idat_fileparts = [f"{row['GSM_ID']}_{row['Sentrix_ID']}_{row['Sentrix_Position']}" for (idx,row) in samples.iterrows()]
    except KeyError as e:
        LOGGER.error(f"Samplesheet is missing {e}.")
        return 
    for file in idat_fileparts:
        files = [f"{file}_Grn.idat", f"{file}_Grn.idat.gz", f"{file}_Red.idat", f"{file}_Red.idat.gz"]
        for idat in files:
            if idat in all_idats_names:
                save_list.append(idat)
        #files = [f"{file}_Grn.idat", f"{file}_Grn.idat.gz", f"{file}_Red.idat", f"{file}_Red.idat.gz"]
        #if Path(idat).exists():
    remove_list = [idat for idat in all_idats if idat.name not in save_list]
    #LOGGER.info(f"removing {len(remove_list)} idats out of a total of {len(all_idats)} found,")
    worked = 'OK' if len(samples.index) == len(save_list)/2 else 'ERROR'
    if worked != 'OK':
        return
    removed = 0
    for idat in all_idats:
        if idat.name in save_list:
            continue
        if Path(idat).exists():
            Path(idat).unlink()
            #print('-',idat)
            removed += 1
    #LOGGER.info(f'removed {removed} idat files not in samplesheet. ready to process remaining ones.')
    LOGGER.info(f"retaining {len(save_list)} files for {len(samples.index)} samples ({worked}). (Dropped {len(remove_list)} idats)")

# remove_idats_not_in_samplesheet('GSE89278/GSE89278_GPL13534_samplesheet.csv','GSE89278')
