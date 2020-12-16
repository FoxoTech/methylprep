# Lib
from enum import Enum, unique

import logging
LOGGER = logging.getLogger(__name__)

@unique
class ArrayType(Enum):
    """This class stores meta data about array types, such as numbers of probes of each type, and how to guess the array from probes in idat files."""

    CUSTOM = 'custom'
    ILLUMINA_27K = '27k'
    ILLUMINA_450K = '450k'
    ILLUMINA_EPIC = 'epic'
    ILLUMINA_EPIC_PLUS = 'epic+'
    ILLUMINA_MOUSE = 'mouse'

    def __str__(self):
        return self.value

    @classmethod
    def from_probe_count(cls, probe_count):
        """Determines array type using number of probes counted in raw idat file. Returns array string."""
        if probe_count == 1055583 or probe_count == 868578:
            return cls.ILLUMINA_EPIC_PLUS

        if 622000 <= probe_count <= 623000:
            return cls.ILLUMINA_450K

        if 1050000 <= probe_count <= 1053000:
            return cls.ILLUMINA_EPIC

        if 54000 <= probe_count <= 56000:
            return cls.ILLUMINA_27K

        if 315000 <= probe_count <= 362000: #V1 actual count from idat: 315639
            return cls.ILLUMINA_MOUSE
            #B1 V1 274390 actual probes == rows in manifest
            #B3 V2 299344 actual probes == rows in manifest file; 361821 count from idat

        if 56000 <= probe_count <= 1100000:
            LOGGER.warning(f'Probe count ({probe_count}) falls outside of normal range. Setting to newest array type: EPIC')
            return cls.ILLUMINA_EPIC

        raise ValueError(f'Unknown array type: ({probe_count} probes detected)')

    @property
    def num_probes(self):
        """ used to load normal cg+ch probes from start of manifest until this point. """
        probe_counts = {
            ArrayType.ILLUMINA_27K: 27578,
            ArrayType.ILLUMINA_450K: 485578,
            ArrayType.ILLUMINA_EPIC: 865919,
            ArrayType.ILLUMINA_EPIC_PLUS: 868699, # was 868699 until Jan 21, 2020. corrected.
            # if EPIC+ is not set to 868699, noob fails downstream. but there are only 868698 probes by my count.
            #ArrayType.ILLUMINA_MOUSE: 268833, #274390 rows in manifest RND1 on 2020-03-25.
            # this includes all types. so ch+cg types == 268832
            # test: added +1 because mouse controls were short by one. and this fixed it. prob
            # need to test them ALL and add +1 header in manifest.py code.
            #B1 V2 = 273757, # test: list where all control probes start
            #ArrayType.ILLUMINA_MOUSE: 298710, #B3 V2 row number for first control probe (after [Controls],,,,,,header)
            ArrayType.ILLUMINA_MOUSE: 297414,
            #287054 #287054 is first control row; no header row
            #297415 # row number for first control probe (after header [Controls],,,, ) with row count starting at zero.
        }
        return probe_counts.get(self)

    @property
    def num_controls(self):
        probe_counts = {
            ArrayType.ILLUMINA_27K: 144, # the manifest does not contain control probe data (illumina's site included)
            ArrayType.ILLUMINA_450K: 850,
            ArrayType.ILLUMINA_EPIC: 635,
            ArrayType.ILLUMINA_EPIC_PLUS: 635,
            ArrayType.ILLUMINA_MOUSE: 1966,
        }
        return probe_counts.get(self)

    @property
    def num_snps(self):
        probe_counts = {
            ArrayType.ILLUMINA_27K: 0,
            ArrayType.ILLUMINA_450K: 65,
            ArrayType.ILLUMINA_EPIC: 59,
            ArrayType.ILLUMINA_EPIC_PLUS: 120,
            ArrayType.ILLUMINA_MOUSE: 1353, #in v2: 536, #was at end of file, now before control (testing)
        }
        return probe_counts.get(self)

""" # doesn't appear to be used anywhere. -- and pipeline Array model conflicts with it.
class Array():
    __slots__ = ['name', 'array_type']

    def __init__(self, name, array_type):
        self.name = name
        self.array_type = array_type
"""
