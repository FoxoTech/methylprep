# Lib
import logging
import pandas as pd
# App
from ..models import (
    METHYLATED_PROBE_SUBSETS,
    UNMETHYLATED_PROBE_SUBSETS,
    METHYLATED_SNP_PROBES,
    UNMETHYLATED_SNP_PROBES,
    #METHYLATED_MOUSE_PROBES,
    #UNMETHYLATED_MOUSE_PROBES
)


__all__ = ['MethylationDataset']


LOGGER = logging.getLogger(__name__)


class MethylationDataset():
    """Wrapper for a collection of methylated or unmethylated probes and their mean intensity values,
    providing common functionality for the subset of probes.

    Arguments:
        raw_dataset {RawDataset} -- A sample's RawDataset for a single well on the processed array.
        manifest {Manifest} -- The Manifest for the correlated RawDataset's array type.
        probe_subsets {list(ProbeSubset)} -- Collection of ProbeSubsets that correspond to the probe type
        (methylated or unmethylated).

    note: self.methylated.data_frame 'bg_corrected' and 'noob' values will be same under preprocess_sesame_noob,
    but different under minfi/legacy pre-v1.4.0 results. And this 'noob' will not match SampleDataContainer.dataframe
    because dye-bias correction happens later in processing.
    """
    __bg_corrected = False
    __preprocessed = False # AKA NOOB CORRECTED
    __dye_bias_corrected = False

    def __init__(self, raw_dataset, manifest, probe_subsets):
        #LOGGER.info('Preprocessing methylation dataset: %s', raw_dataset.sample)

        self.probe_subsets = probe_subsets
        self.raw_dataset = raw_dataset # __init__ uses red_idat and green_idat IdatDatasets

        self.data_frames = {
            probe_subset: self._get_subset_means(manifest, probe_subset)
            for probe_subset in probe_subsets
        }

        self.data_frame = self.build_data_frame()

    @classmethod
    def methylated(cls, raw_dataset, manifest):
        """ convenience method that feeds in a pre-defined list of methylated CpG locii probes """
        return cls(raw_dataset, manifest, METHYLATED_PROBE_SUBSETS)

    @classmethod
    def unmethylated(cls, raw_dataset, manifest):
        """ convenience method that feeds in a pre-defined list of UNmethylated CpG locii probes """
        return cls(raw_dataset, manifest, UNMETHYLATED_PROBE_SUBSETS)

    @classmethod
    def snp_methylated(cls, raw_dataset, manifest):
        """ convenience method that feeds in a pre-defined list of methylated Snp locii probes """
        return cls(raw_dataset, manifest, METHYLATED_SNP_PROBES)

    @classmethod
    def snp_unmethylated(cls, raw_dataset, manifest):
        """ convenience method that feeds in a pre-defined list of UNmethylated Snp locii probes """
        return cls(raw_dataset, manifest, UNMETHYLATED_SNP_PROBES)

    def build_data_frame(self):
        return pd.concat(self.data_frames.values())

    def _get_subset_means(self, manifest, probe_subset):
        """ nearly the same as raw_data.get_subset_means, but this index is IlmnID and raw_data index is illumina_id."""
        channel_means_df = self.raw_dataset.get_channel_means(probe_subset.data_channel)
        channel_means_df = channel_means_df.assign(Channel=probe_subset.data_channel.value)

        #probe_details = probe_subset.get_probe_details(manifest)
        probe_details = manifest.get_probe_details(probe_subset.probe_type, probe_subset.probe_channel)

        # check here for probes that are missing data in manifest, and drop them if they are (better to be imperfect with warnings)
        if probe_details[probe_subset.probe_address.header_name].isna().sum() > 0:
            print('These probes are probably incorrect in your manifest; processing cannot continue.')
            print( probe_details.loc[ probe_details[probe_subset.probe_address.header_name].isna() ].index )
            pre_shape = probe_details.shape
            probe_details = probe_details.drop( probe_details[ probe_details[probe_subset.probe_address.header_name].isna() ].index )
            print(f"{pre_shape[0] - probe_details.shape[0]} removed; {probe_details[probe_subset.probe_address.header_name].isna().sum()} nan remaining; but downstream steps will not work.")
            # this still won't fix it, because OOB also does some filtering.

        # check here for probes that are missing data in manifest, and drop them if they are (better to be imperfect with warnings)
        if probe_details[probe_subset.probe_address.header_name].isna().sum() > 0:
            print('These probes are probably incorrect in your manifest; processing cannot continue.')
            print( probe_details.loc[ probe_details[probe_subset.probe_address.header_name].isna() ].index )
            pre_shape = probe_details.shape
            probe_details = probe_details.drop( probe_details[ probe_details[probe_subset.probe_address.header_name].isna() ].index )
            print(f"{pre_shape[0] - probe_details.shape[0]} removed; {probe_details[probe_subset.probe_address.header_name].isna().sum()} nan remaining; but downstream steps will not work.")
            # this still won't fix it, because OOB also does some filtering.

        return probe_details.merge(
            channel_means_df,
            how='inner',
            left_on=probe_subset.probe_address.header_name, # AddressA_ID or AddressB_ID
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
        column = probe_subset.column_name # AddressA_ID or AddressB_ID

        filtered_corrected = corrected_values.loc[original[column]] # adds the 'bg_corrected' column

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

            #data_frame = data_frame.drop('bg_corrected', axis='columns') # no longer needed
            self.data_frames[probe_subset] = data_frame

        self.data_frame = self.build_data_frame()
        self.__preprocessed = True
