# Lib
import logging
import pandas as pd
# App
from ..models import (
    ArrayType,
    Channel,
    ProbeType,
)
from ..models.probes import FG_PROBE_SUBSETS
from ..files import IdatDataset
from ..utils import inner_join_data
# from ..utils.progress_bar import * # checks environment and imports tqdm appropriately.
from collections import Counter

__all__ = ['RawDataset', 'MethylationDataset', 'get_raw_datasets', 'get_array_type']


LOGGER = logging.getLogger(__name__)


# moved from pipeline.py as no longer needed.
def get_manifest(idat_datasets, array_type=None, manifest_filepath=None):
    """Return a Manifest, given a list of idat_datasets (from idats).

    Arguments:
        raw_datasets {list(RawDataset)} -- Collection of RawDataset instances that
            require a manifest file for the related array_type.

    Keyword Arguments:
        array_type {ArrayType} -- The type of array to process. If not provided, it
            will be inferred from the number of probes in the IDAT file. (default: {None})
        manifest_filepath {path-like} -- Path to the manifest file. If not provided,
            it will be inferred from the array_type and downloaded if necessary (default: {None})

    Returns:
        [Manifest] -- A Manifest instance.
    """

    """ provide a list of raw_datasets and it will return the array type by counting probes """
    if array_type is None:
        array_type = get_array_type(idat_datasets)
    return Manifest(array_type, manifest_filepath)


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
        raw_datasets = [parser(zip_reader, sample) for sample in samples]
    elif not from_s3 and not meta_only:
        parser = RawDataset.from_sample
        raw_datasets = [parser(sample) for sample in samples]

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
                sample_name = f"{dataset.sample.sentrix_id}_{dataset.sample.sentrix_position}"
                samples_by_probe_count[dataset.n_snps_read].append(sample_name)
            LOGGER.error(f'Samples grouped by probe count: {probe_sample_counts.most_common()}')
            LOGGER.error(f'{samples_by_probe_count}')
            raise ValueError(f'IDATs with varying number of probes: {probe_counts}')

    return raw_datasets

# moved to meth_dataset
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
    """Wrapper for a sample and its pair of raw IdatDataset values. Not used in methylprep since v1.5+

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
        channel_means = self.get_channel_means(channel).astype('float32') #.astype('float16')
        return inner_join_data(control_probes, channel_means)

    def get_oob_controls(self, manifest, include_rs=True):
        """ Out-of-bound controls are the mean intensity values for the
        channel in the opposite channel's probes (IG oob and IR oob)"""
        oobG = self.filter_oob_probes(Channel.RED, manifest, self.green_idat, include_rs=include_rs) # manIR + green values
        oobR = self.filter_oob_probes(Channel.GREEN, manifest, self.red_idat, include_rs=include_rs) # manIG + red values

        oobG['Channel'] = Channel.GREEN.value
        oobR['Channel'] = Channel.RED.value

        return {
            Channel.GREEN: oobG,
            Channel.RED: oobR,
        }


    def get_infer_channel_probes(self, manifest, debug=False):
        """ like filter_oob_probes, but returns two dataframes for green and red channels with meth and unmeth columns
        effectively criss-crosses the red-oob channels and appends to green, and appends green-oob to red
        returns a dict with 'green' and 'red' channel probes

        MOVED TO infer_channel_switch.py as a standalone function."""
        probe_details_IR = manifest.get_probe_details(
            probe_type=ProbeType.ONE,
            channel=Channel.RED,
        )
        probe_details_IG = manifest.get_probe_details(
            probe_type=ProbeType.ONE,
            channel=Channel.GREEN,
        )
        # need: illumina_id in index, (green)'meth', (red)'unmeth'
        idat_meth_unmeth = (
            self.green_idat.probe_means
            .rename(columns={'mean_value':'meth'})
            .sort_index()
            .merge(
            self.red_idat.probe_means.rename(columns={'mean_value':'unmeth'}).sort_index(),
            sort=False,
            left_index=True,
            right_index=True)
            )

        # OOB PROBE values are IR(unmeth) and IG(meth); I'll replace IR(meth) and IG(unmeth) below
        # RED channel; uses AddressA_ID for oob IR(unmeth)
        oobR = probe_details_IG.merge(
            idat_meth_unmeth,
            how='inner',
            left_on='AddressA_ID',
            right_index=True,
            suffixes=(False, False),
        )
        oobR = oobR[['meth','unmeth']].sort_index()
        # GREEN channel; AddressB_ID for oob IG(meth)
        oobG = probe_details_IR.merge(
            idat_meth_unmeth,
            how='inner',
            left_on='AddressB_ID',
            right_index=True,
            suffixes=(False, False),
        )
        oobG = oobG[['meth','unmeth']].sort_index()

        # IN BAND probes should be IR(meth) and IG(unmeth)
        # NOTE: below uses same idat DF as before, but the probe_details from manifest are swapped.
        red_in_band = probe_details_IR.merge(
            idat_meth_unmeth,
            how='inner',
            left_on='AddressA_ID',
            right_index=True,
            suffixes=(False, False),
        )
        red_in_band = red_in_band[['meth','unmeth']].sort_index()
        green_in_band = probe_details_IG.merge(
            idat_meth_unmeth,
            how='inner',
            left_on='AddressB_ID',
            right_index=True,
            suffixes=(False, False),
        )
        green_in_band = green_in_band[['meth','unmeth']].sort_index()

        ## HACK: I can't read/get idats to match sesame exactly, so moving columns around to match
        # - swap oob-green[unmeth] with red[meth]
        # - swap oob-red[meth] with green[unmeth]
        oobG_unmeth = oobG['unmeth'].copy()
        oobR_meth = oobR['meth'].copy()
        oobG['unmeth'] = red_in_band['meth'].copy()
        oobR['meth'] = green_in_band['unmeth'].copy()
        red_in_band['meth'] = oobG_unmeth
        green_in_band['unmeth'] = oobR_meth

        # next, add the green-in-band to oobG and red-in-band to oobR
        oobG_IG = oobG.append(green_in_band).sort_index()
        oobR_IR = oobR.append(red_in_band).sort_index()

        # channel swap requires a way to update idats with illumina_ids
        lookupIR = probe_details_IR.merge(
            idat_meth_unmeth,
            how='inner',
            left_on='AddressA_ID',
            right_index=True,
            suffixes=(False, False),
        )[['AddressA_ID','AddressB_ID']]
        lookupIG = probe_details_IG.merge(
            idat_meth_unmeth,
            how='inner',
            left_on='AddressB_ID',
            right_index=True,
            suffixes=(False, False),
        )[['AddressA_ID','AddressB_ID']]
        lookup = lookupIG.append(lookupIR).sort_index()

        if debug:
            return {'green': oobG_IG, 'red': oobR_IR, 'oobG': oobG, 'oobR':oobR, 'IG': green_in_band, 'IR': red_in_band, 'lookup': lookup}
        return {'green': oobG_IG, 'red': oobR_IR, 'IR': red_in_band, 'IG': green_in_band, 'lookup': lookup}


    def filter_oob_probes(self, channel, manifest, idat_dataset, include_rs=True):
        """ adds both channels, and rs probes """
        # channel should be methylprep.models.Channel.RED or methylprep.models.Channel.GREEN
        probe_means = idat_dataset.probe_means # index matches AddressA_ID or AddressB_ID, depending on RED/GREEN channel

        probes = manifest.get_probe_details(
            probe_type=ProbeType.ONE, # returns IR or IG cgxxxx probes only
            channel=channel,
        )[['AddressA_ID', 'AddressB_ID']]
        if include_rs:
            snp_probes = manifest.get_probe_details(
                probe_type=ProbeType.SNP_ONE,
                channel=channel,
            )[['AddressA_ID', 'AddressB_ID']]
            probes = pd.concat([probes, snp_probes])

        if channel == Channel.RED:
            oobG = probes.merge(
                probe_means, # green channel X AddresB (meth) channel
                how='inner',
                left_on='AddressB_ID',
                right_index=True,
                suffixes=(False, False),
            ).rename(columns={'mean_value': 'meth'})
            oobG = oobG.merge(
                probe_means, # green channel X AddresA (unmeth) channel
                how='inner',
                left_on='AddressA_ID',
                right_index=True,
                suffixes=(False, False),
            ).rename(columns={'mean_value': 'unmeth'}).sort_values('IlmnID')
            return oobG.drop(['AddressA_ID', 'AddressB_ID'], axis=1)

        if channel == Channel.GREEN:
            oobR = probes.merge(
                probe_means, # red channel X AddressB for (meth)
                how='inner',
                left_on='AddressB_ID',
                right_index=True,
                suffixes=(False, False),
            ).rename(columns={'mean_value': 'meth'}).sort_values('IlmnID')
            oobR = oobR.merge(
                probe_means, # red channel X AddressA for (unmeth)
                how='inner',
                left_on='AddressA_ID',
                right_index=True,
                suffixes=(False, False),
            ).rename(columns={'mean_value': 'unmeth'})
            return oobR.drop(['AddressA_ID', 'AddressB_ID'], axis=1)


    def _old_filter_oob_probes(self, channel, manifest, idat_dataset):
        """ this is the step where it appears that illumina_id (internal probe numbers)
        are matched to the AddressA_ID / B_IDs from manifest,
        which allows for 'cgXXXXXXX' probe names to be used later.

        infer_channel_switch() requires both meth and unmeth in the probes returned."""
        probe_details = manifest.get_probe_details(
            probe_type=ProbeType.ONE,
            channel=channel,
        )
        # 2020-03-25: probe_details was returning an empty DataFrame with mouse,
        # because two new probe types existed (IR, IG) -- note that new types results
        # in this null issue and a huber ZeroDivisionError ultimately in CLI.
        # 2021-02-16: IG and IR probes have both AddressA_ID and AddressB_ID, but only one should be used.
        # so filtering here based on channel.
        probe_details = probe_details[['AddressA_ID', 'AddressB_ID']]
        probe_means = idat_dataset.probe_means # index matches AddressA_ID or AddressB_ID, depending on RED/GREEN channel
        if channel == Channel.RED:
            set_a = probe_details.merge(
                probe_means, # green channel for oob
                how='inner',
                left_on='AddressA_ID',
                right_index=True,
                suffixes=(False, False),
            )
            return set_a.drop(['AddressA_ID', 'AddressB_ID'], axis=1)

        if channel == Channel.GREEN:
            set_b = probe_details.merge(
                probe_means, # red channel for oob
                how='inner',
                left_on='AddressB_ID',
                right_index=True,
                suffixes=(False, False),
            )
            return set_b.drop(['AddressA_ID', 'AddressB_ID'], axis=1)
        # below: would contain duplicate probes in index, both the red-in-channel and green-oob, so commented it out.
        #oob_probes = set_a.append(set_b)
        #return oob_probes

    def get_fg_values(self, manifest, channel, index_by='illumina_id'):
        """ appears to only be used in bg_correct part of NOOB function; index is illumina_id, not IlmnID.
        also called in preprocess.preprocess_noob()."""
        #LOGGER.info('Preprocessing %s foreground datasets: %s', channel, self.sample)

        probe_subsets = FG_PROBE_SUBSETS[channel]

        channel_foregrounds = [
            self.get_subset_means(probe_subset, manifest, index_by=index_by)
            for probe_subset in probe_subsets
        ]

        # debug - trying to locate the duplicate mouse probes
        #for idx,probe_subset in enumerate(probe_subsets):
        #    print(probe_subset.probe_address, probe_subset.probe_type, probe_subset.data_channel, probe_subset.probe_channel)
            # this has both ProbeAddress.A and IlmnID -- check for rs.
            #test = pd.concat(channel_foregrounds)
        #    test = channel_foregrounds[idx]
        #    print(test.shape, test.index.duplicated().sum())
            #print([rs for rs in test['IlmnID'] if 'rs' in rs])
            #if test.index.duplicated().sum() > 0:
        #        import pdb;pdb.set_trace()
        #        print(test[test.index.duplicated()])

        return pd.concat(channel_foregrounds)

    def get_subset_means(self, probe_subset, manifest, index_by='illumina_id'):
        """ called by get_fg_values for each of 6 probe subsets """
        channel_means_df = self.get_channel_means(probe_subset.data_channel)
        probe_details = manifest.get_probe_details(probe_subset.probe_type, probe_subset.probe_channel)
        column_name = probe_subset.column_name
        # instead of dropping duplicates, use index_by='IlmnID' with mouse
        if probe_details[column_name].duplicated().sum() > 0:
            LOGGER.warning(f'filtered {probe_details[column_name].duplicated().sum()} duplicate probes ({column_name})')
            probe_details = probe_details[ ~probe_details[column_name].duplicated() ]

        probe_details = probe_details.reset_index()
        merge_df = probe_details[[column_name, 'probe_type', 'IlmnID']]
        merge_df = merge_df.reset_index()
        merge_df = merge_df.set_index(column_name)

        merged_data =  inner_join_data(channel_means_df, merge_df)
        if index_by == 'IlmnID': # this is only used by preprocess_noob()
            merged_data.index.name = column_name
            merged_data = merged_data.reset_index()
            merged_data = merged_data.set_index('IlmnID')
            merged_data = merged_data.drop('index', axis=1)
            merged_data = merged_data.rename(columns={'AddressA_ID':'illumina_id', 'AddressB_ID': 'illumina_id'})
        else:
            merged_data.index.name = column_name
            merged_data = merged_data.drop(['IlmnID','index'], axis=1)
        return merged_data


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
            #probe_subset: self.raw_dataset.get_subset_means(probe_subset, manifest, index_by='IlmnID')
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
        """ nearly the same as raw_data.get_subset_means, but this index is IlmnID and raw_data index is illumina_id.
        this is called for each probe_subset using the @classmethods above."""
        channel_means_df = self.raw_dataset.get_channel_means(probe_subset.data_channel)
        channel_means_df = channel_means_df.assign(Channel=probe_subset.data_channel.value)
        # [illumina_id (index) | probe_mean_value | Grn or Red]

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

        print(f"_set_subset_bg_corrected {str(probe_subset)} {column} {probe_subset.data_channel} --> {filtered_corrected.shape}")

        updated = original.merge(
            filtered_corrected[['bg_corrected']],
            how='inner',
            left_on=column,
            right_index=True,
            suffixes=(False, False),
        )

        print(f"_set_subset_bg_corrected {original.shape} --> {updated.shape}")

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


class RawMetaDataset():
    """moved to sigset.py
    Wrapper for a sample and meta data, without its pair of raw IdatDataset values.

    Arguments:
        sample {Sample} -- A Sample parsed from the sample sheet.

    each Sample contains (at a minimum):
        data_dir=self.data_dir
        sentrix_id=sentrix_id
        sentrix_position=sentrix_position
    """

    def __init__(self, sample):
        self.sample = sample
