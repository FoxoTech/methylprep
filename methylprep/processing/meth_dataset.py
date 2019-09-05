# Lib
import logging
import pandas as pd
# App
from ..models import METHYLATED_PROBE_SUBSETS, UNMETHYLATED_PROBE_SUBSETS


__all__ = ['MethylationDataset']


LOGGER = logging.getLogger(__name__)


class MethylationDataset():
    """Wrapper for a collection of methylated or unmethylated probes and their mean intensity values,
    providing common functionality for the subset of probes.

    Arguments:
        raw_dataset {RawDataset} -- A sample's RawDataset for a single well on the processed array.
        manifest {Manifest} -- The Manifest for the correlated RawDataset's array type.
        probe_subsets {list(ProbeSubset)} -- Collection of ProbeSubsets that correspond tot he probe type
            (methylated or unmethylated).
    """
    __bg_corrected = False
    __preprocessed = False

    def __init__(self, raw_dataset, manifest, probe_subsets):
        LOGGER.info('Preprocessing methylation dataset: %s', raw_dataset.sample)

        self.probe_subsets = probe_subsets
        self.raw_dataset = raw_dataset

        self.data_frames = {
            probe_subset: self._get_subset_means(manifest, probe_subset)
            for probe_subset in probe_subsets
        }

        self.data_frame = self.build_data_frame()

    @classmethod
    def methylated(cls, raw_dataset, manifest):
        return cls(raw_dataset, manifest, METHYLATED_PROBE_SUBSETS)

    @classmethod
    def unmethylated(cls, raw_dataset, manifest):
        return cls(raw_dataset, manifest, UNMETHYLATED_PROBE_SUBSETS)

    def build_data_frame(self):
        return pd.concat(self.data_frames.values())

    def _get_subset_means(self, manifest, probe_subset):
        channel_means = self.raw_dataset.get_channel_means(probe_subset.data_channel)
        channel_means = channel_means.assign(Channel=probe_subset.data_channel.value)

        probe_details = probe_subset.get_probe_details(manifest)

        return probe_details.merge(
            channel_means,
            how='inner',
            left_on=probe_subset.probe_address.header_name,
            right_index=True,
            suffixes=(False, False),
        )

    def set_bg_corrected(self, green_corrected, red_corrected):
        for probe_subset in self.data_frames:
            if probe_subset.is_red:
                corrected_values = red_corrected
            elif probe_subset.is_green:
                corrected_values = green_corrected
            else:
                raise ValueError('No data_channel for probe_subset')

            self._set_subset_bg_corrected(probe_subset, corrected_values)

        self.data_frame = self.build_data_frame()
        self.__bg_corrected = True

    def _set_subset_bg_corrected(self, probe_subset, corrected_values):
        original = self.data_frames[probe_subset]
        column = probe_subset.column_name

        filtered_corrected = corrected_values.loc[original[column]]

        updated = original.merge(
            filtered_corrected[['bg_corrected']],
            how='inner',
            left_on=column,
            right_index=True,
            suffixes=(False, False),
        )

        self.data_frames[probe_subset] = updated

    def set_noob(self, red_factor):
        for probe_subset, data_frame in self.data_frames.items():
            if probe_subset.is_red:
                data_frame = data_frame.assign(noob=data_frame['bg_corrected'] * red_factor)
            else:
                data_frame = data_frame.assign(noob=data_frame['bg_corrected'])

            self.data_frames[probe_subset] = data_frame

        self.data_frame = self.build_data_frame()
        self.__preprocessed = True
