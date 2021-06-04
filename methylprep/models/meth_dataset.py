# Lib
import logging
import pandas as pd
import numpy as np
# App
from ..models import (
    METHYLATED_PROBE_SUBSETS,
    UNMETHYLATED_PROBE_SUBSETS,
    METHYLATED_SNP_PROBES,
    UNMETHYLATED_SNP_PROBES,
    FG_PROBE_SUBSETS,
    ArrayType,
    Channel,
    ProbeType,
)
from ..files import IdatDataset
from ..utils.progress_bar import * # checks environment and imports tqdm appropriately.
from collections import Counter


__all__ = ['MethylationDataset', 'SigSet', 'parse_sample_sheet_into_idat_datasets', 'RawMetaDataset']


LOGGER = logging.getLogger(__name__)



class RawMetaDataset():
    """Wrapper for a sample and meta data, without its pair of raw IdatDataset values."""
    def __init__(self, sample):
        self.sample = sample

def parse_sample_sheet_into_idat_datasets(sample_sheet, sample_name=None, from_s3=None, meta_only=False):
    """Generates a collection of IdatDatasets from samples in a sample sheet.

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
        [RawDatasets] -- A list of idat data pairs, each a dict like {'green_idat': green_idat, 'red_idat': red_idat}
    """
    # now idat_datasets is not a class, but just a list of dicts, with each dict being a pair of red_idat and green_idat Objects.

    LOGGER.debug('Reading IDATs from sample sheet')

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
        idat_datasets = [parser(sample) for sample in samples]
    elif from_s3 and not meta_only:
        #parser = RawDataset.from_sample_s3
        zip_reader = from_s3
        def parser(zip_reader, sample):
            green_filepath = sample.get_filepath('idat', Channel.GREEN)
            green_idat = IdatDataset(green_filepath, channel=Channel.GREEN)
            red_filepath = sample.get_filepath('idat', Channel.RED)
            red_idat = IdatDataset(red_filepath, channel=Channel.RED)
            return {'green_idat': green_idat, 'red_idat': red_idat, 'sample': sample}
        idat_datasets = tqdm([parser(zip_reader, sample) for sample in samples], total=len(samples), desc='Reading IDATs')
    elif not from_s3 and not meta_only:
        #parser = RawDataset.from_sample
        def parser(sample):
            green_filepath = sample.get_filepath('idat', Channel.GREEN)
            green_idat = IdatDataset(green_filepath, channel=Channel.GREEN)
            red_filepath = sample.get_filepath('idat', Channel.RED)
            red_idat = IdatDataset(red_filepath, channel=Channel.RED)
            return {'green_idat': green_idat, 'red_idat': red_idat, 'sample': sample}
        idat_datasets = tqdm([parser(sample) for sample in samples], total=len(samples), desc='Reading IDATs')

    if not meta_only:
        # ensure all idat files have same number of probes
        batch_probe_counts = set()
        counts_per_sample = Counter()
        for dataset in idat_datasets:
            snps_read = {dataset['green_idat'].n_snps_read, dataset['red_idat'].n_snps_read}
            if len(snps_read) > 1:
                raise ValueError('IDAT files have a varying number of probes (compared Grn to Red channel)')
            n_snps_read = snps_read.pop()
            batch_probe_counts.add(n_snps_read)
            counts_per_sample[n_snps_read] += 1
        if len(batch_probe_counts) != 1:
            LOGGER.error(f'Samples grouped by probe count: {counts_per_sample.most_common()}')
            raise ValueError(f'IDATs with varying number of probes: {probe_counts}')
    return idat_datasets


class SigSet():
    """
    I’m gonna try to create a fresh methylprep “SigSet” to replace our methylationDataset and RawDataset objects, which are redundant, and even have redundant functions within them. Part of why I have been frustrated/confused by our code.
    Central to the SeSAMe platform is the SigSet data structure, an S4 class with slots containing signals for six different classes of probes:
    [x] II - Type-II probes;
    [x] IR - Type-I Red channel probes;
    [x] IG - Type-I Grn channel probes;
    [x] oobG - Out-of-band Grn channel probes (matching Type-I Red channel probes in number);
    [x] oobR - Out-of-band Red channel probes (matching Type-I Grn channel probes in number);
    [-] ctl - control probes.
    [x] methylated, unmethylated, snp_methylated, snp_unmethylated
    [x] fg_green, fg_red (opposite of oobG and oobR)

    - just tidying up how we access this stuff, and trying to stick to IlmnID everywhere because the illumina_id within IDAT files is no longer unique as a ref.
    - I checked again, and no other array breaks these rules. But sounds like Bret won’t stick to this pattern going forward with designs. So I suspect other software will break with new arrays, unless they rewrite for this too.

    - this combines every layer of objects between IdatDatasets and SampleDataContainers.
    - this avoids looping through probe subsets, instead referring to a lookup-dataframe of how these relate.
    - avoids probes.py
        probe_type is a derived label, not in manifest (I, II, SnpI, SnpII, control)
    """
    __bg_corrected = False
    __preprocessed = False # AKA NOOB CORRECTED
    __dye_bias_corrected = False

    # we found the SECRET DECODER RING!
    #name | data_channel | probe_address | probe_channel | probe_type | fg_green | fg_red | meth | unmeth | snp_meth | snp_unmeth
    # 'data_channel' refers to either the green_idat red_idat values; 'Color_Channel' + probe_address is used by type-I probes to refer to one of two measurements within that probe.
    df_columns=['name','data_channel','probe_address','Color_Channel','probe_type',
        'fg_green', 'fg_red', 'meth', 'unmeth', 'snp_meth', 'snp_unmeth', 'foreground', 'snp', 'Infinium_Design_Type']
    data=[ # in-band probes
        ['II-None-Green-Meth', 'GREEN', 'AddressA_ID', None, 'II', 1, 0, 1, 0, 0, 0, 1, 0, 'II'],
        ['IG-A-Unmeth', 'GREEN', 'AddressA_ID', 'Grn', 'I', 1, 0, 0, 1, 0, 0, 1, 0, 'I'],
        ['IG-B-Meth', 'GREEN', 'AddressB_ID', 'Grn', 'I', 1, 0, 1, 0, 0, 0, 1, 0, 'I'],
        ['II-None-Red-Unmeth', 'RED', 'AddressA_ID', None, 'II', 0, 1, 0, 1, 0, 0, 1, 0, 'II'],
        ['IR-A-Unmeth', 'RED', 'AddressA_ID', 'Red', 'I', 1, 0, 0, 1, 0, 0, 1, 0, 'I'],
        ['IR-B-Meth', 'RED', 'AddressB_ID', 'Red', 'I', 1, 0, 1, 0, 0, 0, 1, 0, 'I'],
        # in-band SNPS
        ['SnpII-None-Green-Meth', 'GREEN', 'AddressA_ID', None, 'SnpII', 1, 0, 0, 0, 1, 0, 1, 1, 'II'],
        ['SnpIG-A-Unmeth', 'GREEN', 'AddressA_ID', 'Grn', 'SnpI', 1, 0, 0, 0, 0, 1, 1, 1, 'I'],
        ['SnpIG-B-Meth', 'GREEN', 'AddressB_ID', 'Grn', 'SnpI', 1, 0, 0, 0, 1, 0, 1, 1, 'I'],
        ['SnpII-None-Red-Unmeth', 'RED', 'AddressA_ID', None, 'SnpII', 0, 1, 0, 0, 0, 1, 1, 1, 'II'],
        ['SnpIR-A-Unmeth', 'RED', 'AddressA_ID', 'Red', 'SnpI', 1, 0, 0, 0, 0, 1, 1, 1, 'I'],
        ['SnpIR-B-Meth', 'RED', 'AddressB_ID', 'Red', 'SnpI', 1, 0, 1, 0, 0, 0, 1, 1, 'I'],
        # oob subsets -- these should include Snp-oob probes, but they need to be addressed separately.
        ['oobG-Unmeth', 'GREEN', 'AddressA_ID', 'Red', 'I', 0, 0, 0, 1, 0, 0, 0, 0, 'I'],
        ['oobG-Meth',   'GREEN', 'AddressB_ID', 'Red', 'I', 0, 0, 1, 0, 0, 0, 0, 0, 'I'],
        ['oobR-Unmeth', 'RED', 'AddressA_ID', 'Grn', 'I', 0, 0, 0, 1, 0, 0, 0, 0, 'I'],
        ['oobR-Meth',   'RED', 'AddressB_ID', 'Grn', 'I', 0, 0, 1, 0, 0, 0, 0, 0, 'I'],
        ['SnpIG-oobG-Unmeth', 'GREEN', 'AddressA_ID', 'Red', 'SnpI', 0, 0, 0, 0, 0, 1, 0, 1, 'I'],
        ['SnpIG-oobG-Meth',   'GREEN', 'AddressB_ID', 'Red', 'SnpI', 0, 0, 0, 0, 1, 0, 0, 1, 'I'],
        ['SnpIR-oobR-Unmeth', 'RED', 'AddressA_ID', 'Grn', 'SnpI', 0, 0, 0, 0, 0, 1, 0, 1, 'I'],
        ['SnpIR-oobR-Meth',   'RED', 'AddressB_ID', 'Grn', 'SnpI', 0, 0, 0, 0, 1, 0, 0, 1, 'I']]

    idat_decoder = pd.DataFrame(data=data, columns=df_columns).set_index('name')

    # how to decode the idat into logical chunks for processing; sesame II/IG/IR includes snps by default; mprep excludes snps from meth/unmeth/fg_red/fg_green
    subsets = {
        'II': ['II-None-Green-Meth', 'II-None-Red-Unmeth', 'SnpII-None-Green-Meth', 'SnpII-None-Red-Unmeth'],
        'IG': ['IG-A-Unmeth', 'IG-B-Meth', 'SnpIG-A-Unmeth', 'SnpIG-B-Meth'],
        'IR': ['IR-A-Unmeth', 'IR-B-Meth', 'SnpIR-A-Unmeth', 'SnpIR-B-Meth'],
        'oobG': ['oobG-Unmeth','oobG-Meth', 'SnpIG-oobG-Unmeth', 'SnpIG-oobG-Meth'],
        'oobR': ['oobR-Unmeth','oobR-Meth', 'SnpIR-oobR-Unmeth', 'SnpIR-oobR-Meth'],
        'methylated': ['II-None-Green-Meth', 'IG-B-Meth', 'IR-B-Meth'],
        'unmethylated': ['II-None-Red-Unmeth', 'IG-A-Unmeth', 'IR-A-Unmeth'],
        'snp_methylated': ['SnpII-None-Green-Meth', 'SnpIG-B-Meth', 'SnpIR-B-Meth'],
        'snp_unmethylated': ['SnpII-None-Red-Unmeth', 'SnpIG-A-Unmeth', 'SnpIR-A-Unmeth'],
        'fg_green': ['II-None-Green-Meth', 'IG-A-Unmeth', 'IG-B-Meth'],
        'fg_red': ['II-None-Red-Unmeth', 'IR-A-Unmeth', 'IR-B-Meth'],
    }
    # after __init__, SigSet will have class variables for each of the keys in subsets above.

    def __init__(self, sample, green_idat, red_idat, manifest, debug=False):
        """ green_idat has .probe_means and .meta as main functions
        and for extra info, use extra kwargs:
        red= m.files.IdatDataset('9247377093_R02C01_Red.idat', m.models.Channel.RED, verbose=True, std_dev=True, nbeads=True)
        """
        snps_read = {green_idat.n_snps_read, red_idat.n_snps_read}
        if len(snps_read) > 1:
            raise ValueError('IDAT files have a varying number of probes (comparing Grn to Red channel)')
        if not (str(green_idat.channel) == 'Grn' and str(red_idat.channel) == 'Red'):
            raise ValueError("The IDAT files you supplied seem to be reversed. Check the order of your inputs to SigSet")
        self.n_snps_read = snps_read.pop()
        # DEBUG
        #self.array_type = ArrayType.from_probe_count(self.n_snps_read)

        # these next two should be unnecessary, because nothing should be reading idats downstream; use self.data_channel instead
        #self.green_idat = green_idat
        #self.red_idat = red_idat
        self.data_channel = {'GREEN': green_idat.probe_means, 'RED': red_idat.probe_means} # indexed to illumina_ids
        # illumina_ids are all II means, plus a stacked list of type-I-AddressA and type-I-AddressB means
        self.sample = sample
        self.man = manifest.data_frame # relevant columns are 'probe_type', AddressA_ID, AddressB_ID, index, Color_Channel
        self.man = self.man[ ~self.man.index.str.startswith('rs') ] # snp_man covers these
        self.snp_man = manifest.snp_data_frame.set_index('IlmnID')
        self.ctl_man = manifest.control_data_frame
        self.address_code = {'AddressA_ID':'A', 'AddressB_ID':'B', 'A':'AddressA_ID', 'B':'AddressB_ID'}

        """
        ## SigSet EPIC
        ##  - @IG probes: 49989 - 332 4145 70 7094 599 2958 ...
        ##  - @IR probes: 92294 - 183 8040 1949 6152 833 89 ...
        ##  - @II probes: 724612 - 6543 1596 3133 1011 3035 2837 ...
        ##  - @oobG probes: 92294 - 138 277 107 218 232 80 ...
        ##  - @oobR probes: 49989 - 1013 150 81 910 448 183 ...
        ##  - @ctl probes: 635 ...
        ##  - @pval: 866895 - 0.005141179 0.04914081 0.002757492 ...
        """

        for subset, decoder_parts in self.subsets.items():
            data_frames = {}
            for part in decoder_parts:
                i = self.idat_decoder.loc[part]
                ref = self.snp_man if i['snp'] == 1 else self.man
                # can't merge on NAType, so filling in -1s. No probe_means illumina_ids will match -1
                # using -1 instead of NaN throughout solves a lot of problems!
                if ref['AddressA_ID'].isna().sum() > 0:
                    ref['AddressA_ID'] = ref['AddressA_ID'].fillna(-1)
                if ref['AddressB_ID'].isna().sum() > 0:
                    print('filling')
                    #ref.loc[ (ref['AddressB_ID'].isnull()), 'AddressB_ID'] = -1
                    ref['AddressB_ID'].fillna(-1, inplace=True)
                # and pandas won't compare NaN to NaN... so need this extra color_channel filter
                color_channel = ref['Color_Channel'].isna() if i['Color_Channel'] is None else ref['Color_Channel'] == i['Color_Channel']
                probe_ids = ref[ (ref['Infinium_Design_Type'] == i['Infinium_Design_Type']) & (color_channel) ][i['probe_address']]
                probe_means = self.data_channel[i['data_channel']]
                probe_means = probe_means[ probe_means.index.isin(probe_ids) ]
                if debug:
                    print(f"{subset} -- {part}: {probe_ids.shape}, {probe_means.shape}")
                    if len(probe_ids) == 0:
                        print('no probes matched')
                        import pdb;pdb.set_trace()
                # merge and establish IlmnIDs from illumina_ids here
                probe_subset_data = ref[['AddressA_ID', 'AddressB_ID']].merge(probe_means,
                    left_on=i['probe_address'],
                    right_index=True)
                probe_subset_data['used'] = self.address_code[i['probe_address']]
                probe_subset_data = probe_subset_data.rename(
                    columns={'mean_value': 'Meth' if 'Meth' in part else 'Unmeth'})
                #data_frames['Unmeth' if 'Meth' in part else 'Meth'] = None
                data_frames[part] = probe_subset_data

            try:
                # here, put the meth and unmeth parts into separate columns as we combine
                meth_parts = [frame for frame in data_frames.values() if 'Meth' in frame.columns]
                unmeth_parts = [frame for frame in data_frames.values() if 'Unmeth' in frame.columns]
                if unmeth_parts == []:
                    data_frame = pd.concat(meth_parts)
                    data_frame['Unmeth'] = None
                elif meth_parts == []:
                    data_frame = pd.concat(unmeth_parts)
                    data_frame['Meth'] = None
                else:
                    data_frame = pd.concat(meth_parts)
                    data_frame = data_frame.merge(pd.concat(unmeth_parts)[['Unmeth']], left_index=True, right_index=True)
                if debug:
                    print(subset, data_frame.shape)
                setattr(self, subset, data_frame)
            except Exception as e:
                print(e)
                import pdb;pdb.set_trace()

    # originally was `set_bg_corrected` from MethylationDataset | called by NOOB
    def update_probe_means(self, noob_green, noob_red, red_factor=None):
        """ pass in two dataframes (green and red) with IlmnIDs in index and a 'bg_corrected' column in each.

        because __init__ has created each subset as a dataframe with IlmnID in index, this matches to index.
        and uses decoder to parse whether 'Meth' or 'Unmeth' values get updated.

        upstream: container.sigset.update_probe_means(noob_green, noob_red)
        """

        for probe_subset, decoder_parts in self.subsets.items():
            if self.debug: print(f'--- probe_subset {probe_subset} ---')
            df = getattr(self, probe_subset) # contains multiple parts in one dataframe, so need to update with logic each time
            df = df.assign(bg_corrected_Meth=None)
            df = df.assign(bg_corrected_Unmeth=None)
            # need to assign new values to red and green channels separately, and decode into meth/unmeth column
            for part in decoder_parts:
                #**** NOOB is matching using tango A/B IDs upstream, must change. ****#
                # assuming probes are passed in indexed to IlmnIDs, we can merge/update whatever matches
                i = self.idat_decoder.loc[part]
                if i['meth'] == 1:
                    column = 'Meth'
                    bg_column = 'bg_corrected_Meth'
                if i['unmeth'] == 1:
                    column = 'Unmeth'
                    bg_column = 'bg_corrected_Unmeth'
                if i['snp_meth'] == 1:
                    column = 'Meth'
                    bg_column = 'bg_corrected_Meth'
                if i['snp_unmeth'] == 1:
                    column = 'Unmeth'
                    bg_column = 'bg_corrected_Unmeth'

                if i['data_channel'] == 'GREEN':
                    updated = green_corrected.copy()
                elif i['data_channel'] == 'RED':
                    updated = red_corrected.copy()
                    if red_factor is not None:
                        # NOTE: this changes None to NaN and dtype is Object
                        updated = updated['bg_corrected'] * red_factor

                # what IlmnIDs are in this subset?
                ref = self.snp_man if i['snp'] == 1 else self.man
                # and pandas won't compare NaN to NaN... so need this extra color_channel filter
                color_channel = ref['Color_Channel'].isna() if i['Color_Channel'] is None else ref['Color_Channel'] == i['Color_Channel']
                IlmnIDs = ref[ (ref['Infinium_Design_Type'] == i['Infinium_Design_Type']) & (color_channel) ].index

                # use update; cannot merge if the IlmnIDs overlap, to avoid NaNs and bg_corrected_x/_y cols
                updated = updated.rename(columns={'bg_corrected':bg_column})
                updated = updated.loc[ updated.index.isin( IlmnIDs )]
                debug_pre = df[bg_column].isna().sum()
                df.update(updated) # matches on index
                debug_post = df[bg_column].isna().sum()
                num_updated = debug_pre - debug_post
                if self.debug: print(f"{part} {df.shape[0]}, {num_updated} updated")
                if num_updated > 0:
                    print(df.notna().sum())
            # overwrite!
            setattr(self, probe_subset, df)
        self.__preprocessed = True # applied by set_noob
        self.__bg_corrected = True

    # from MethylationDataset | called by self.set_bg_corrected for each probe subset | now part of set_bg_corrected
    def _set_subset_bg_corrected(self, probe_subset, corrected_values):
        pass
    # from MethylationDataset | now part of update_probe_means
    def set_noob(self, red_factor):
        """ same method as update_probe_means, but simply applies a linear correction to all RED channel values """
        #update_probe_means(self, noob_green, noob_red, red_factor)
        raise KeyError("set_noob replaced by update_probe_means in v1.5+")


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
