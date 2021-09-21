import numpy as np
import pandas as pd
from pathlib import Path
try:
    from importlib import resources # py >= 3.7
except ImportError: # py < 3.7
    import pkg_resources
    loader = pkg_resources.resource_filename
pkg_namespace = 'methylprep.models'

try:
    if pkg_resources:
        # resource_filename context manager does not support "with"
        probe_filepath = loader(pkg_namespace, 'qualityMask450.txt.gz')
        qualityMask450 = pd.read_csv(probe_filepath)['x']
        probe_filepath = loader(pkg_namespace, 'qualityMaskEPIC.txt.gz')
        qualityMaskEPIC = pd.read_csv(probe_filepath)['x']
        probe_filepath = loader(pkg_namespace, 'qualityMaskEPICPLUS.txt.gz')
        qualityMaskEPICPLUS = pd.read_csv(probe_filepath)['x']
        probe_filepath = loader(pkg_namespace, 'qualityMaskmouse.txt.gz')
        qualityMaskmouse = pd.read_csv(probe_filepath)['x']
except:
    with resources.path(pkg_namespace, 'qualityMask450.txt.gz') as probe_filepath:
        qualityMask450 = pd.read_csv(probe_filepath)['x']
    with resources.path(pkg_namespace, 'qualityMaskEPIC.txt.gz') as probe_filepath:
        qualityMaskEPIC = pd.read_csv(probe_filepath)['x']
    with resources.path(pkg_namespace, 'qualityMaskEPICPLUS.txt.gz') as probe_filepath:
        qualityMaskEPICPLUS = pd.read_csv(probe_filepath)['x']
    with resources.path(pkg_namespace, 'qualityMaskmouse.txt.gz') as probe_filepath:
        qualityMaskmouse = pd.read_csv(probe_filepath)['x']


"""
def sketchy_probes_warning(filepath):
    ''' not implemented anywhere yet '''
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
"""

# mouse had 293199 probes in mask from wanding: sesameData::sesameDataGet("MM285.address")$ordering
#qualityMaskEPICPLUS = pd.concat([qualityMask450, qualityMaskEPIC], axis=0)
#qualityMaskEPICPLUS = qualityMaskEPICPLUS.drop_duplicates()
# need a lookup here to rename probes to match EPICPLUS
# here, I combined the other masks and renamed probes to match epic+. 3889 probes in epic_plus don't match either of these.
# qualityMaskEPICPLUS_r.to_csv('/Users/mmaxmeister/methylprep/methylprep/models/qualityMaskEPICPLUS.txt', index=False)
# qualityMaskEPICPLUS = qualityMaskEPICPLUS.apply( lambda row: row.split('_')[0] )
# in preprocess.py, I'm simply passing through all 3889 probes whose EPIC names don't match EPIC+ cg.... names.
