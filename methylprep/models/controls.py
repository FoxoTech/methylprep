# Lib
from enum import Enum, unique


@unique
class ControlType(Enum):
    STAINING = 'STAINING' # NOT FOUND IN MOUSE MANIFEST
    EXTENSION = 'EXTENSION'
    HYBRIDIZATION = 'HYBRIDIZATION'
    TARGET_REMOVAL = 'TARGET REMOVAL'
    BISULFITE_CONVERSION_I = 'BISULFITE CONVERSION I'
    BISULFITE_CONVERSION_II = 'BISULFITE CONVERSION II'
    SPECIFICITY_I = 'SPECIFICITY I'
    SPECIFICITY_II = 'SPECIFICITY II'
    NON_POLYMORPHIC = 'NON-POLYMORPHIC' # changed from 'NON - POLYMORPHIC' on 2020-12-01 to match manifests
    NEGATIVE = 'NEGATIVE'
    RESTORATION = 'RESTORATION'
    NORM_A = 'NORM_A'
    NORM_C = 'NORM_C'
    NORM_G = 'NORM_G'
    NORM_T = 'NORM_T'

    @classmethod
    def normalization_green(cls):
        return (
            cls.NORM_C.value,
            cls.NORM_G.value,
        )

    @classmethod
    def normalization_red(cls):
        return (
            cls.NORM_A.value,
            cls.NORM_T.value,
        )


class ControlProbe():
    """ NOT USED ANYWHERE """
    __slots__ = [
        'address',
        'control_type',
        'color',
        'extended_type',
    ]

    def __init__(self, address, control_type, color, extended_type):
        self.address = address
        self.control_type = control_type
        self.color = color
        self.extended_type = extended_type
