# Lib
import logging
import pandas as pd
# App
from ..models import (
    FG_PROBE_SUBSETS,
    ArrayType,
    Channel,
    ProbeType,
)
from ..files import IdatDataset
from ..utils import inner_join_data
from ..utils.progress_bar import * # checks environment and imports tqdm appropriately.
from collections import Counter

__all__ = ['RawDataset', 'RawMetaDataset', 'get_raw_datasets', 'get_raw_meta_datasets', 'get_array_type']


LOGGER = logging.getLogger(__name__)


def get_raw_datasets(sample_sheet, sample_name=None, from_s3=None, meta_only=False):
    """Generates a collection of RawDataset instances for the samples in a sample sheet.

    Arguments:
        sample_sheet {SampleSheet} -- The SampleSheet from which the data originates.

    Keyword Arguments:
        sample_name {string} -- Optional: one sample to process from the sample_sheet. (default: {None})
        from_s3 {zip_reader} -- pass in a S3ZipReader object to extract idat files from a zipfile hosted on s3.
        meta_only {True/False} -- doesn't read idat files, only parses the meta data about them.
        (RawMetaDataset is same as RawDataset but has no idat probe values stored in object, because not needed in pipeline)

    Raises:
        ValueError: If the number of probes between raw datasets differ.

    Returns:
        [RawDataset] -- A RawDataset instance.

    """

    LOGGER.debug('Generating raw datasets from sample sheet')

    if not sample_name:
        samples = sample_sheet.get_samples()
    elif type(sample_name) is list:
        samples = [
            sample_sheet.get_sample(sample)
            for sample in sample_name
        ]
    else:
        samples = [sample_sheet.get_sample(sample_name)]
        LOGGER.info("Found sample in SampleSheet: {0}".format(sample_name))

    if from_s3 and meta_only:
        parser = RawMetaDataset
        raw_datasets = [parser(sample) for sample in samples]
    elif from_s3 and not meta_only:
        parser = RawDataset.from_sample_s3
        zip_reader = from_s3
        raw_datasets = tqdm([parser(zip_reader, sample) for sample in samples], total=len(samples), desc='Getting raw datasets')
    elif not from_s3 and not meta_only:
        parser = RawDataset.from_sample
        raw_datasets = tqdm([parser(sample) for sample in samples], total=len(samples), desc='Getting raw datasets')

    if not meta_only:
        # ensure all idat files have same number of probes
        probe_counts = {
            dataset.n_snps_read
            for dataset in raw_datasets
        }

        if len(probe_counts) != 1:
            # also explain which samples have which probes -- for splitting samples up
            probe_sample_counts = Counter([dataset.n_snps_read for dataset in raw_datasets])
            samples_by_probe_count = {probe_count:[] for probe_count in list(probe_counts)}
            for dataset in raw_datasets:
                sample_name = f"{dataset.sample.sentrix_id}_{dataset.sample.sentrix_id}"
                samples_by_probe_count[dataset.n_snps_read].append(sample_name)
            LOGGER.error(f'Samples grouped by probe count: {probe_sample_counts.most_common()}')
            LOGGER.error(f'{samples_by_probe_count}')
            raise ValueError(f'IDATs with varying number of probes: {probe_counts}')

    return raw_datasets


def get_array_type(raw_datasets):
    """ provide a list of raw_datasets and it will return the array type by counting probes """
    array_types = {dataset.array_type for dataset in raw_datasets}
    if len(array_types) == 0:
        raise ValueError('could not identify array type from IDATs')
    elif len(array_types) != 1:
        raise ValueError('IDATs with varying array types')
    array_type = array_types.pop()
    return array_type


class RawDataset():
    """Wrapper for a sample and its pair of raw IdatDataset values.

    Arguments:
        sample {Sample} -- A Sample parsed from the sample sheet.
        green_idat {IdatDataset} -- The sample's GREEN channel IdatDataset.
        red_idat {IdatDataset} -- The sample's RED channel IdatDataset.

    Raises:
        ValueError: If the IDAT file pair have differing number of probes.
        TypeError: If an invalid Channel is provided when parsing an IDAT file.
    """

    def __init__(self, sample, green_idat, red_idat):
        snps_read = {green_idat.n_snps_read, red_idat.n_snps_read}

        if len(snps_read) > 1:
            raise ValueError('IDAT files have a varying number of probes (compared Grn to Red channel)')

        self.n_snps_read = snps_read.pop()
        self.green_idat = green_idat
        self.red_idat = red_idat
        self.sample = sample
        self.array_type = ArrayType.from_probe_count(self.n_snps_read)

    @classmethod
    def from_sample(cls, sample):
        green_filepath = sample.get_filepath('idat', Channel.GREEN)
        green_idat = IdatDataset(green_filepath, channel=Channel.GREEN)

        red_filepath = sample.get_filepath('idat', Channel.RED)
        red_idat = IdatDataset(red_filepath, channel=Channel.RED)
        return cls(sample, green_idat, red_idat)

    @classmethod
    def from_sample_s3(cls, zip_reader, sample):
        green_filepath = sample.get_file_s3(zip_reader, 'idat', suffix=Channel.GREEN)
        green_idat = IdatDataset(green_filepath, channel=Channel.GREEN)

        red_filepath = sample.get_file_s3(zip_reader, 'idat', suffix=Channel.RED)
        red_idat = IdatDataset(red_filepath, channel=Channel.RED)
        return cls(sample, green_idat, red_idat)

    def get_channel_means(self, channel):
        if not isinstance(channel, Channel):
            raise TypeError('channel is not a valid Channel')
        if channel is Channel.GREEN:
            return self.green_idat.probe_means
        return self.red_idat.probe_means

    def get_fg_controls(self, manifest, channel):
        #LOGGER.info('Preprocessing %s foreground controls dataset: %s', channel, self.sample)
        control_probes = manifest.control_data_frame
        channel_means = self.get_channel_means(channel).astype('float16')
        return inner_join_data(control_probes, channel_means)

    def get_oob_controls(self, manifest):
        """ Out-of-bound controls are the mean intensity values for the
        channel in the opposite channel's probes """
        oob_green = self.filter_oob_probes(Channel.RED, manifest, self.green_idat)
        oob_red = self.filter_oob_probes(Channel.GREEN, manifest, self.red_idat)

        oob_green['Channel'] = Channel.GREEN.value
        oob_red['Channel'] = Channel.RED.value

        return {
            Channel.GREEN: oob_green,
            Channel.RED: oob_red,
        }

    def filter_oob_probes(self, channel, manifest, idat_dataset):
        """ this is the step where it appears that illumina_id (internal probe numbers)
        are matched to the AddressA_ID / B_IDs from manifest,
        which allows for 'cgXXXXXXX' probe names to be used later. """
        probe_details = manifest.get_probe_details(
            probe_type=ProbeType.ONE,
            channel=channel,
        )
        # 2020-03-25: probe_details was returning an empty DataFrame with mouse,
        # because two new probe types existed (IR, IG) -- note that new types results
        # in this null issue and a huber ZeroDivisionError ultimately in CLI.

        probe_details = probe_details[['AddressA_ID', 'AddressB_ID']]

        probe_means = idat_dataset.probe_means

        set_a = probe_details.merge(
            probe_means,
            how='inner',
            left_on='AddressA_ID',
            right_index=True,
            suffixes=(False, False),
        )

        set_b = probe_details.merge(
            probe_means,
            how='inner',
            left_on='AddressB_ID',
            right_index=True,
            suffixes=(False, False),
        )

        oob_probes = set_a.append(set_b) # will contain duplicates for probes that have both red and grn channels (II)
        return oob_probes

    def get_fg_values(self, manifest, channel):
        """ appears to only be used in NOOB function """
        #LOGGER.info('Preprocessing %s foreground datasets: %s', channel, self.sample)

        probe_subsets = FG_PROBE_SUBSETS[channel]

        channel_foregrounds = [
            self.get_subset_means(probe_subset, manifest)
            for probe_subset in probe_subsets
        ]

        # debug - trying to locate the SNP signal
        #for probe_subset in probe_subsets:
        #    print(probe_subset.probe_address, probe_subset.probe_type, probe_subset.data_channel, probe_subset.probe_channel)
        #    # this has both ProbeAddress.A and IlmnID -- check for rs.
        #    test = pd.concat(channel_foregrounds)
        #   print('get_fg_values', test.shape, test.index.duplicated())
        #    print([rs for rs in test['IlmnID'] if 'rs' in rs])

        return pd.concat(channel_foregrounds)

    def get_subset_means(self, probe_subset, manifest):
        """apparently, not called anywhere """
        channel_means_df = self.get_channel_means(probe_subset.data_channel)
        probe_details = probe_subset.get_probe_details(manifest)
        column_name = probe_subset.column_name

        merge_df = probe_details[[column_name, 'probe_type']]
        merge_df = merge_df.reset_index()
        merge_df = merge_df.set_index(column_name)

        return inner_join_data(channel_means_df, merge_df)


class RawMetaDataset():
    """Wrapper for a sample and meta data, without its pair of raw IdatDataset values.

    Arguments:
        sample {Sample} -- A Sample parsed from the sample sheet.

    each Sample contains (at a minimum):
        data_dir=self.data_dir
        sentrix_id=sentrix_id
        sentrix_position=sentrix_position
    """

    def __init__(self, sample):
        self.sample = sample
