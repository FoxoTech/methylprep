import numpy as np
import pandas as pd
from pathlib import Path
from importlib import resources # py3.7+
pkg_namespace = 'methylprep.models'

def sketchy_probes_warning(filepath):
    """ not implemented anywhere yet """
    # used to warn user if running a Catalog_ID that contains sketchy illumina probes

    with resources.path(pkg_namespace, 'illumina_sketchy_probes_996.npy') as probe_filepath:
        sketchy_probes = np.load(probe_filepath)
    # "If the first 8 numbers of Sentrix_ID (i.e. xxxxxxxx0001) are greater or equal to 20422033,
    # then the BeadChip originates from production batches using the new manufacturing process."
    new_manufacturing_cutoff_id = 20422033

    if Path(filepath).name.startswith('GSM'):
        catalog_id = int(Path(filepath).name.split('_')[1][:8])
    else:
        catalog_id = int(Path(filepath).name[:8])
    if catalog_id  >= new_manufacturing_cutoff_id:
        return Path(filepath).name
    else:
        return None

with resources.path(pkg_namespace, 'qualityMask450.txt.gz') as probe_filepath:
    qualityMask450 = pd.read_csv(probe_filepath)['x']
with resources.path(pkg_namespace, 'qualityMaskEPIC.txt.gz') as probe_filepath:
    qualityMaskEPIC = pd.read_csv(probe_filepath)['x']
