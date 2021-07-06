# Lib
from enum import Enum, unique


@unique
class Channel(Enum):
    """ idat probes measure either a red or green fluorescence.
    This specifies which to return within idat.py: red_idat or green_idat."""
    RED = 'Red'
    GREEN = 'Grn'

    def __str__(self):
        return self.value

    @property
    def is_green(self):
        return self == self.GREEN

    @property
    def is_red(self):
        return self == self.RED


@unique
class ProbeAddress(Enum):
    """AddressA_ID and AddressB_ID are columns in the manifest csv that contain internal Illumina probe identifiers.

    Type II probes use AddressA_ID; Type I uses both, because there are two probes, two colors.

    probe intensities in .idat files are keyed to one of these ids, but processed data is always keyed to the
    IlmnID probe "names" -- so this is used in converting between IDs. It is used to define probe sets below
    in this probes.py"""
    A = 'A'
    B = 'B'

    @property
    def header_name(self):
        if self == self.A:
            return 'AddressA_ID'
        return 'AddressB_ID'


@unique
class ProbeType(Enum):
    """ probes can either be type I or type II for CpG or Snp sequences. Control probes are used for background
    correction in different fluorescence ranges and staining efficiency.
    Type I probes record EITHER a red or a green value.
    Type II probes record both values together.
    NOOB uses the red fluorescence on a green probe and vice versa to calculate background fluorescence."""
    ONE = 'I'
    TWO = 'II'
    SNP_ONE = 'SnpI'
    SNP_TWO = 'SnpII'
    CONTROL = 'Control'
    # I was separating out mouse probes EARLY, here, but found they need to be processed like all other probes, THEN removed in post-processing stage.
    #MOUSE_ONE = 'MouseI'
    #MOUSE_TWO = 'MouseII'

    def __str__(self):
        return self.value

    @staticmethod
    def from_manifest_values(name, infinium_type):
        """ this function determines which of four kinds of probe goes with this name, using either
        the Infinium_Design_Type (I or II) or the name (starts with 'rs')
        and decides 'Control' is non of the above."""
        is_control = any([name.startswith('rs'),
                        name.startswith('ctl'),
                        name.startswith('neg'),
                        name.startswith('BSC'),
                        name.startswith('NON'),
                        ])
        is_snp = name.startswith('rs')

        if is_control and is_snp:
            if infinium_type == 'I':
                return ProbeType.SNP_ONE
            elif infinium_type == 'II':
                return ProbeType.SNP_TWO
            else:
                return ProbeType.CONTROL
        elif is_control:
            return ProbeType.CONTROL

        elif infinium_type == 'I':
            return ProbeType.ONE

        elif infinium_type == 'II':
            return ProbeType.TWO

        elif infinium_type in ('IR','IG'): # mouse only -- these are type I probes but Bret's files label them this way
            return ProbeType.ONE

        return ProbeType.CONTROL


class Probe():
    """ this doesn't appear to be instantiated anywhere in methylprep """
    __slots__ = [
        'address',
        'illumina_id',
        'probe_type',
    ]

    def __init__(self, address, illumina_id, probe_type):
        self.address = address
        self.illumina_id = illumina_id
        self.probe_type = probe_type


class ProbeSubset():
    """ used below in probes.py to define sub-sets of probes:
    foreground-(red|green|all), or (un)methylated probes
    """
    __slots__ = [
        'data_channel',
        'probe_address',
        'probe_channel',
        'probe_type',
    ]

    def __init__(self, data_channel, probe_address, probe_channel, probe_type):
        self.data_channel = data_channel
        self.probe_address = probe_address
        self.probe_channel = probe_channel
        self.probe_type = probe_type

    def __str__(self):
        return f'{self.probe_type}-{self.probe_channel}'

    @property
    def is_green(self):
        return self.data_channel.is_green

    @property
    def is_red(self):
        return self.data_channel.is_red

    @property
    def column_name(self):
        return self.probe_address.header_name

    #def get_probe_details(self, manifest):
    #    return manifest.get_probe_details(
    #        probe_type=self.probe_type,
    #        channel=self.probe_channel,
    #    )




METHYLATED_SNP_PROBES = (
    ProbeSubset(
        data_channel=Channel.GREEN,
        probe_address=ProbeAddress.A,
        probe_channel=None,
        probe_type=ProbeType.SNP_TWO,
    ),
    ProbeSubset(
        data_channel=Channel.GREEN,
        probe_address=ProbeAddress.B,
        probe_channel=Channel.GREEN,
        probe_type=ProbeType.SNP_ONE,
    ),
    ProbeSubset(
        data_channel=Channel.RED,
        probe_address=ProbeAddress.B,
        probe_channel=Channel.RED,
        probe_type=ProbeType.SNP_ONE,
    ),
)

UNMETHYLATED_SNP_PROBES = (
    ProbeSubset(
        data_channel=Channel.RED,
        probe_address=ProbeAddress.A,
        probe_channel=None,
        probe_type=ProbeType.SNP_TWO,
    ),
    ProbeSubset(
        data_channel=Channel.GREEN,
        probe_address=ProbeAddress.A,
        probe_channel=Channel.GREEN,
        probe_type=ProbeType.SNP_ONE,
    ),
    ProbeSubset(
        data_channel=Channel.RED,
        probe_address=ProbeAddress.A,
        probe_channel=Channel.RED,
        probe_type=ProbeType.SNP_ONE,
    ),
)

FG_GREEN_PROBE_SUBSETS = (
    ProbeSubset(
        data_channel=Channel.GREEN,
        probe_address=ProbeAddress.A,
        probe_channel=None,
        probe_type=ProbeType.TWO,
    ),
    ProbeSubset(
        data_channel=Channel.GREEN,
        probe_address=ProbeAddress.A,
        probe_channel=Channel.GREEN,
        probe_type=ProbeType.ONE,
    ),
    ProbeSubset(
        data_channel=Channel.GREEN,
        probe_address=ProbeAddress.B,
        probe_channel=Channel.GREEN,
        probe_type=ProbeType.ONE,
    ),
)

FG_RED_PROBE_SUBSETS = (
    ProbeSubset(
        data_channel=Channel.RED,
        probe_address=ProbeAddress.A,
        probe_channel=None,
        probe_type=ProbeType.TWO,
    ),
    ProbeSubset(
        data_channel=Channel.RED,
        probe_address=ProbeAddress.A,
        probe_channel=Channel.RED,
        probe_type=ProbeType.ONE,
    ),
    ProbeSubset(
        data_channel=Channel.RED,
        probe_address=ProbeAddress.B,
        probe_channel=Channel.RED,
        probe_type=ProbeType.ONE,
    ),
)

FG_PROBE_SUBSETS = {
    Channel.GREEN: FG_GREEN_PROBE_SUBSETS,
    Channel.RED: FG_RED_PROBE_SUBSETS,
}

METHYLATED_PROBE_SUBSETS = ( # == IG[AddressB] + II[green] + IR[AddressB]
    ProbeSubset(
        data_channel=Channel.GREEN,
        probe_address=ProbeAddress.A,
        probe_channel=None,
        probe_type=ProbeType.TWO,
    ),
    ProbeSubset(
        data_channel=Channel.GREEN,
        probe_address=ProbeAddress.B,
        probe_channel=Channel.GREEN,
        probe_type=ProbeType.ONE,
    ),
    ProbeSubset(
        data_channel=Channel.RED,
        probe_address=ProbeAddress.B,
        probe_channel=Channel.RED,
        probe_type=ProbeType.ONE,
    ),
)

UNMETHYLATED_PROBE_SUBSETS = ( # == IG [AddressA] + II[red] + IR[AddressA]
    ProbeSubset(
        data_channel=Channel.RED,
        probe_address=ProbeAddress.A,
        probe_channel=None,
        probe_type=ProbeType.TWO,
    ),
    ProbeSubset(
        data_channel=Channel.GREEN,
        probe_address=ProbeAddress.A,
        probe_channel=Channel.GREEN,
        probe_type=ProbeType.ONE,
    ),
    ProbeSubset(
        data_channel=Channel.RED,
        probe_address=ProbeAddress.A,
        probe_channel=Channel.RED,
        probe_type=ProbeType.ONE,
    ),
)

''' DISABLED because these faster to process with normals then remove during postprocessing.
METHYLATED_MOUSE_PROBES = (
    ProbeSubset(
        data_channel=Channel.GREEN,
        probe_address=ProbeAddress.A,
        probe_channel=None,
        probe_type=ProbeType.MOUSE_TWO,
    ),
    ProbeSubset(
        data_channel=Channel.GREEN,
        probe_address=ProbeAddress.B,
        probe_channel=Channel.GREEN,
        probe_type=ProbeType.MOUSE_ONE,
    ),
    ProbeSubset(
        data_channel=Channel.RED,
        probe_address=ProbeAddress.B,
        probe_channel=Channel.RED,
        probe_type=ProbeType.MOUSE_ONE,
    ),
)

UNMETHYLATED_MOUSE_PROBES = (
    ProbeSubset(
        data_channel=Channel.RED,
        probe_address=ProbeAddress.A,
        probe_channel=None,
        probe_type=ProbeType.MOUSE_TWO,
    ),
    ProbeSubset(
        data_channel=Channel.GREEN,
        probe_address=ProbeAddress.A,
        probe_channel=Channel.GREEN,
        probe_type=ProbeType.MOUSE_ONE,
    ),
    ProbeSubset(
        data_channel=Channel.RED,
        probe_address=ProbeAddress.A,
        probe_channel=Channel.RED,
        probe_type=ProbeType.MOUSE_ONE,
    ),
)


# SNP-PREP-NOTES from JAN 2020
# notebook equivalent of SNP_PROBES below.
#manifest = data_containers[0].manifest.data_frame
#snpProbesI = manifest[manifest['probe_type']=='SnpI']
#snpProbesI_Grn = snpProbesI[snpProbesI['Color_Channel']=='Grn']
#snpProbesI_Red = snpProbesI[snpProbesI['Color_Channel']=='Red']
#snpProbesII = manifest[manifest['probe_type']=='SnpII']

SNP_PROBES = (
    # SNP_II_PROBES
    ProbeSubset(
        data_channel=Channel.GREEN,
        probe_address=ProbeAddress.A,
        probe_channel=None,
        probe_type=ProbeType.SNP_TWO,
    ),
    # SNP_I_GREEN_PROBES
    ProbeSubset(
        data_channel=Channel.GREEN,
        probe_address=ProbeAddress.A,
        probe_channel=Channel.GREEN,
        probe_type=ProbeType.SNP_ONE,
    ),
    # SNP_I_RED_PROBES
    ProbeSubset(
        data_channel=Channel.RED,
        probe_address=ProbeAddress.B,
        probe_channel=Channel.RED,
        probe_type=ProbeType.SNP_ONE,
    )
)
'''
