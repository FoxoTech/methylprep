# Lib
from enum import Enum, unique


@unique
class Channel(Enum):
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
    A = 'A'
    B = 'B'

    @property
    def header_name(self):
        if self == self.A:
            return 'AddressA_ID'
        return 'AddressB_ID'


@unique
class ProbeType(Enum):
    ONE = 'I'
    TWO = 'II'
    SNP_ONE = 'SnpI'
    SNP_TWO = 'SnpII'
    CONTROL = 'Control'

    def __str__(self):
        return self.value

    @staticmethod
    def from_manifest_values(name, infinium_type):
        is_snp = name.startswith('rs')

        if infinium_type == 'I':
            return ProbeType.SNP_ONE if is_snp else ProbeType.ONE

        if infinium_type == 'II':
            return ProbeType.SNP_TWO if is_snp else ProbeType.TWO

        return ProbeType.CONTROL


class Probe():
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

    def get_probe_details(self, manifest):
        return manifest.get_probe_details(
            probe_type=self.probe_type,
            channel=self.probe_channel,
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

METHYLATED_PROBE_SUBSETS = (
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

UNMETHYLATED_PROBE_SUBSETS = (
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
