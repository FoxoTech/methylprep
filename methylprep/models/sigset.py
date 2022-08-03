# Lib
import logging
import pandas as pd
import numpy as np
# App
from ..models import (
    ArrayType,
    Channel,
    ProbeType,
)
from ..files import IdatDataset
from ..utils.progress_bar import * # checks environment and imports tqdm appropriately.
from collections import Counter


__all__ = ['SigSet', 'parse_sample_sheet_into_idat_datasets', 'RawMetaDataset']


LOGGER = logging.getLogger(__name__)

def get_array_type(idat_dataset_pairs):
    """ provide a list of idat_dataset_pairs and it will return the array type, confirming probe counts match in batch. """
    array_types = {dataset['array_type'] for dataset in idat_dataset_pairs}
    if len(array_types) == 0:
        raise ValueError('could not identify array type from IDATs')
    elif len(array_types) != 1:
        raise ValueError('IDATs with varying array types')
    array_type = array_types.pop()
    return array_type

class RawMetaDataset():
    """Wrapper for a sample and meta data, without its pair of raw IdatDataset values."""
    def __init__(self, sample):
        self.sample = sample

def parse_sample_sheet_into_idat_datasets(sample_sheet, sample_name=None, from_s3=None, meta_only=False, bit='float32'):
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

    #LOGGER.info(f'Reading {len(samples)} IDATs from sample sheet')
    if from_s3 and meta_only:
        parser = RawMetaDataset
        idat_datasets = [parser(sample) for sample in samples]
    elif from_s3 and not meta_only:
        #parser = RawDataset.from_sample_s3
        zip_reader = from_s3
        def parser(zip_reader, sample):
            green_filepath = sample.get_filepath('idat', Channel.GREEN)
            green_idat = IdatDataset(green_filepath, channel=Channel.GREEN, bit=bit)
            red_filepath = sample.get_filepath('idat', Channel.RED)
            red_idat = IdatDataset(red_filepath, channel=Channel.RED, bit=bit)
            return {'green_idat': green_idat, 'red_idat': red_idat, 'sample': sample}
        idat_datasets = []
        for sample in tqdm(samples, total=len(samples), desc='Reading IDATs'):
            idat_datasets.append(parser(zip_reader, sample))
    elif not from_s3 and not meta_only:
        #parser = RawDataset.from_sample
        def parser(sample):
            green_filepath = sample.get_filepath('idat', Channel.GREEN)
            green_idat = IdatDataset(green_filepath, channel=Channel.GREEN, bit=bit)
            red_filepath = sample.get_filepath('idat', Channel.RED)
            red_idat = IdatDataset(red_filepath, channel=Channel.RED, bit=bit)
            return {'green_idat': green_idat, 'red_idat': red_idat, 'sample': sample}
        idat_datasets = []
        for sample in tqdm(samples, total=len(samples), desc='Reading IDATs'):
            idat_datasets.append(parser(sample))

    if not meta_only:
        idat_datasets = list(idat_datasets) # tqdm objects are not subscriptable, not like a real list
        # ensure all idat files have same number of probes
        batch_probe_counts = set()
        counts_per_sample = Counter()
        for idx,dataset in enumerate(idat_datasets):
            snps_read = {dataset['green_idat'].n_snps_read, dataset['red_idat'].n_snps_read}
            if len(snps_read) > 1:
                raise ValueError('IDAT files have a varying number of probes (compared Grn to Red channel)')
            n_snps_read = snps_read.pop()
            batch_probe_counts.add(n_snps_read)
            counts_per_sample[n_snps_read] += 1
            idat_datasets[idx]['array_type'] = ArrayType.from_probe_count(n_snps_read)
        if len(batch_probe_counts) != 1:
            array_types = Counter([dataset['array_type'] for dataset in idat_datasets])
            LOGGER.warning(f"These IDATs have varying numbers of probes: {counts_per_sample.most_common()} for these array types: {array_types.most_common()}")
            LOGGER.warning(f"(Processing will drop any probes that are not found across all samples for a given array type.)")
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
    [x] ctrl_green, ctrl_red - control probes.
    [x] methylated, unmethylated, snp_methylated, snp_unmethylated
    [x] fg_green, fg_red (opposite of oobG and oobR) AKA ibG, ibR for in-band probes.

    - just tidying up how we access this stuff, and trying to stick to IlmnID everywhere because the illumina_id within IDAT files is no longer unique as a ref.
    - I checked again, and no other array breaks these rules. But sounds like Bret won’t stick to this pattern going forward with designs. So I suspect other software will break with new arrays, unless they rewrite for this too.

    - this combines every layer of objects between IdatDatasets and SampleDataContainers.
    - this avoids looping through probe subsets, instead referring to a lookup-dataframe of how these relate.
    - avoids probes.py
        probe_type is a derived label, not in manifest (I, II, SnpI, SnpII, control)
    """
    __bg_corrected = False
    __minfi_noob = False # linear_dye applied
    __dye_bias_corrected = False
    __preprocessed = False

    # we found the SECRET DECODER RING!
    #name | data_channel | probe_address | probe_channel | probe_type | fg_green | fg_red | meth | unmeth | snp_meth | snp_unmeth
    # 'data_channel' refers to either the green_idat red_idat values; 'Color_Channel' + probe_address is used by type-I probes to refer to one of two measurements within that probe.
    df_columns=['name','data_channel','probe_address','Color_Channel','probe_type',
        'fg_green', 'fg_red', 'meth', 'unmeth', 'snp_meth', 'snp_unmeth', 'foreground', 'snp', 'Infinium_Design_Type']
    data=[ # in-band probes
        ['II-None-Green-Meth', 'GREEN', 'AddressA_ID', None, 'II', 1, 0, 1, 0, 0, 0, 1, 0, 'II'],
        ['IG-A-Unmeth', 'GREEN', 'AddressA_ID', 'Grn', 'I', 1, 0, 0, 1, 0, 0, 1, 0, 'I'],
        ['IG-B-Meth', 'GREEN', 'AddressB_ID', 'Grn', 'I', 1, 0, 1, 0, 0, 0, 1, 0, 'I'],
        ['II-None-Red-Unmeth', 'RED',   'AddressA_ID', None, 'II', 0, 1, 0, 1, 0, 0, 1, 0, 'II'],
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

    # how to decode the idat into logical chunks for processing; sesame II/IG/IR includes snps by default; mprep separates snps from meth/unmeth/fg_red/fg_green
    subsets = {
        'II': ['II-None-Green-Meth', 'II-None-Red-Unmeth', 'SnpII-None-Green-Meth', 'SnpII-None-Red-Unmeth'],
        'IG': ['IG-A-Unmeth', 'IG-B-Meth', 'SnpIG-A-Unmeth', 'SnpIG-B-Meth'],
        'IR': ['IR-A-Unmeth', 'IR-B-Meth', 'SnpIR-A-Unmeth', 'SnpIR-B-Meth'],
        'oobG': ['oobG-Unmeth','oobG-Meth', 'SnpIG-oobG-Unmeth', 'SnpIG-oobG-Meth'],
        'oobR': ['oobR-Unmeth','oobR-Meth', 'SnpIR-oobR-Unmeth', 'SnpIR-oobR-Meth'],
        'methylated': ['II-None-Green-Meth', 'IG-B-Meth', 'IR-B-Meth', 'SnpII-None-Green-Meth', 'SnpIG-B-Meth', 'SnpIR-B-Meth'],
        'unmethylated': ['II-None-Red-Unmeth', 'IG-A-Unmeth', 'IR-A-Unmeth', 'SnpII-None-Red-Unmeth', 'SnpIG-A-Unmeth', 'SnpIR-A-Unmeth'],
        'snp_methylated': ['SnpII-None-Green-Meth', 'SnpIG-B-Meth', 'SnpIR-B-Meth'],
        'snp_unmethylated': ['SnpII-None-Red-Unmeth', 'SnpIG-A-Unmeth', 'SnpIR-A-Unmeth'],
        #'fg_green': ['II-None-Green-Meth', 'IG-A-Unmeth', 'IG-B-Meth'],
        #'fg_red':   ['II-None-Red-Unmeth', 'IR-A-Unmeth', 'IR-B-Meth'],
        # ibG is fg_green plus the SNPs
        'ibG': ['II-None-Green-Meth', 'IG-A-Unmeth', 'IG-B-Meth', 'SnpIG-A-Unmeth', 'SnpIG-B-Meth', 'SnpII-None-Green-Meth'],
        'ibR': ['II-None-Red-Unmeth', 'IR-A-Unmeth', 'IR-B-Meth', 'SnpIR-A-Unmeth', 'SnpIR-B-Meth', 'SnpII-None-Red-Unmeth'],
        # 'ctrl_green' and 'ctrl_red' are defined attributes below, because they use a separate manifest and addressing system.
    }
    # after __init__, SigSet will have class variables for each of the keys in subsets above.

    def __init__(self, sample, green_idat, red_idat, manifest, debug=False):
        """ green_idat has .probe_means and .meta as main functions
        and for extra info, use extra kwargs:
        red= m.files.IdatDataset('9247377093_R02C01_Red.idat', m.models.Channel.RED, verbose=True, std_dev=True, nbeads=True)
        """
        self.debug = debug
        snps_read = {green_idat.n_snps_read, red_idat.n_snps_read}
        if len(snps_read) > 1:
            raise ValueError('IDAT files have a varying number of probes (comparing Grn to Red channel)')
        if (str(green_idat.channel) != 'Grn' or str(red_idat.channel) != 'Red'):
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
        self.ctrl_green = self.ctl_man.merge(
            green_idat.probe_means.astype('float32'),
            how='inner', left_index=True, right_index=True)
        self.ctrl_red = self.ctl_man.merge(
            red_idat.probe_means.astype('float32'),
            how='inner', left_index=True, right_index=True)
        self.array_type = manifest.array_type
        if self.array_type == ArrayType.ILLUMINA_MOUSE:
            self.mouse_probes_mask = ( (self.man['design'] == 'Multi')  | (self.man['design'] == 'Random') )
        else:
            self.mouse_probes_mask = None
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

        SigSet 450k
        @II 350076 ................... methylated   485512
        @IG 46298 ... oobR 46298 ..... unmethylated 485512
        @IR 89203 ... oobG 89203 ..... snp_methylated   65
        .............................. snp_unmethylated 65
        fg_green 396325 |vs| ibG 396374 (incl 40 + 9 SNPs)  --(flattened)--> 442672
        fg_red   439223 |vs| ibR 439279 (incl 40 + 16 SNPs) --(flattened)--> 528482
        """

        if debug: print('DEBUG comparing [manifest probe_IDs vs idat probe_means]')

        for subset, decoder_parts in self.subsets.items():
            data_frames = {}
            for part in decoder_parts:
                i = self.idat_decoder.loc[part]
                ref = self.snp_man.copy() if i['snp'] == 1 else self.man.copy()
                # can't merge on NAType, so filling in -1s. No probe_means illumina_ids will match -1
                # using -1 instead of NaN throughout solves a lot of problems!
                if ref['AddressA_ID'].isna().sum() > 0:
                    ref.loc[:, 'AddressA_ID'].fillna(-1, inplace=True)
                if ref['AddressB_ID'].isna().sum() > 0:
                    #ref.loc[ (ref['AddressB_ID'].isnull()), 'AddressB_ID'] = -1
                    ref.loc[:, 'AddressB_ID'].fillna(-1, inplace=True)
                # and pandas won't compare NaN to NaN... so need this extra color_channel filter
                color_channel = ref['Color_Channel'].isna() if i['Color_Channel'] is None else (ref['Color_Channel'] == i['Color_Channel'])
                probe_ids = ref[ (ref['Infinium_Design_Type'] == i['Infinium_Design_Type']) & (color_channel) ][i['probe_address']]
                probe_means = self.data_channel[i['data_channel']] # starts with all 361821 mouse probes here, keyed to illumina_ids
                probe_means = probe_means.reset_index() # index is Nth row; illumina_id now a column that can be redundant
                probe_means = probe_means[ probe_means.illumina_id.isin(probe_ids) ]
                if len(probe_ids) == 0:
                    LOGGER.error(f"SigSet.init(): no probes matched for {subset}:{part}")
                #************ DEBUG ***********#
                if debug:
                    #print(f"DEBUG duplicated probe_ids (from manifest): {len( probe_ids[probe_ids.duplicated(keep=False)] )}")
                    duped = len( probe_ids[probe_ids.duplicated(keep=False)] )
                    dupe_msg = f"-- {duped} multiprobes" if duped != 0 else ''
                    means_msg = probe_means.shape if probe_means.shape[0] != probe_ids.shape[0] else 'OK'
                    print(f"DEBUG {subset} -- {part}: {probe_ids.shape} -- {means_msg} {dupe_msg}")
                # 2021-11-29: confirmed that all 361821 mouse means in IDAT DO get read. 4622 of these are control probes, but
                # methylprep only uses 633 of them (matching 635 EPIC probes for QC).
                # 919 of these probes are duplicates having the same illumina_id but different IlmnIDs (TC11, TC12, TC13 etc..) that dont merge right.
                #if subset == 'methylated' and part == 'IR-B-Meth':

                #************ DEBUG ***********#
                # merge and establish IlmnIDs from illumina_ids here
                ref_probe_names = ref[['AddressA_ID', 'AddressB_ID']].reset_index()
                probe_subset_data = ref_probe_names.merge(probe_means, how='inner',
                    left_on=i['probe_address'],
                    right_on='illumina_id')
                probe_subset_data['used'] = self.address_code[i['probe_address']]
                mean_col_name = 'Meth' if 'Meth' in part else 'Unmeth'
                probe_subset_data = probe_subset_data.rename(columns={'mean_value': mean_col_name})
                probe_subset_data = probe_subset_data.drop(['illumina_id'], axis='columns')
                # looks like duplicated probes DO have different mean_values, so it works. nothing is lost.
                #if (probe_subset_data.IlmnID.duplicated().sum() > 0
                #    or probe_subset_data.IlmnID.isna().sum() > 0
                #    or probe_subset_data[mean_col_name].isna().sum() > 0
                #    or probe_subset_data['AddressA_ID'].duplicated().sum() > 0
                #    or ('II' not in part and probe_subset_data['AddressB_ID'].duplicated().sum() > 0)):
                #    print('ERROR')
                #    import pdb;pdb.set_trace()
                probe_subset_data = probe_subset_data.set_index('IlmnID')
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
                    # need to keep NaNs in Meth / Unmeth when merging, so 'outer'
                    test_data_frame = data_frame.merge(pd.concat(unmeth_parts)[['Unmeth']], left_index=True, right_index=True, how='inner')
                    data_frame = data_frame.merge(pd.concat(unmeth_parts)[['Unmeth']], left_index=True, right_index=True, how='outer')
                    #if len(data_frame) != len(test_data_frame):
                    #    print(f"---- DEBUG {subset} {data_frame.shape} vs {test_data_frame.shape}")
                    # -- this explained by having NaNs in either Meth/Unmeth channel
                if debug:
                    print(subset, len(data_frame))
                setattr(self, subset, data_frame)
            except Exception as e:
                raise Exception(f"SigSet: {e}")

        self.starting_probe_counts = {subset: getattr(self, subset).shape[0] for subset in self.subsets.keys()} # DEBUGGING
        self.detect_and_drop_duplicates()
        if debug: self.check_for_probe_loss()

    # originally was `set_bg_corrected` from MethylationDataset | called by NOOB
    def update_probe_means(self, noob_green, noob_red, red_factor=None):
        """ pass in two dataframes (green and red) with IlmnIDs in index and a 'bg_corrected' column in each.

        because __init__ has created each subset as a dataframe with IlmnID in index, this matches to index.
        and uses decoder to parse whether 'Meth' or 'Unmeth' values get updated.

        upstream: container.sigset.update_probe_means(noob_green, noob_red)

        replaces 'bg_corrected' column with 'noob_Meth' or 'noob_Unmeth' column.

        does NOT update ctrl_red or ctrl_green; these are updated within the NOOB function because structually different.
        """

        for probe_subset, decoder_parts in self.subsets.items():
            if self.debug: print(f'--- probe_subset {probe_subset} ---')
            df = getattr(self, probe_subset) # contains multiple parts in one dataframe, so need to update with logic each time
            df = df.assign(noob_Meth=None)
            df['noob_Meth'] = df['noob_Meth'].astype('float32')
            df = df.assign(noob_Unmeth=None)
            df['noob_Unmeth'] = df['noob_Unmeth'].astype('float32')
            # need to assign new values to red and green channels separately, and decode into meth/unmeth column
            for part in decoder_parts:
                try:
                    #**** NOOB is matching using tango A/B IDs upstream, must change. ****#
                    # assuming probes are passed in indexed to IlmnIDs, we can merge/update whatever matches
                    i = self.idat_decoder.loc[part]
                    if i['meth'] == 1:
                        column = 'Meth'
                        bg_column = 'noob_Meth'
                    if i['unmeth'] == 1:
                        column = 'Unmeth'
                        bg_column = 'noob_Unmeth'
                    if i['snp_meth'] == 1:
                        column = 'Meth'
                        bg_column = 'noob_Meth'
                    if i['snp_unmeth'] == 1:
                        column = 'Unmeth'
                        bg_column = 'noob_Unmeth'

                    if i['data_channel'] == 'GREEN':
                        updated = noob_green.copy()
                    elif i['data_channel'] == 'RED':
                        updated = noob_red.copy()
                        if red_factor is not None:
                            # NOTE: this changes None to NaN and dtype is Object
                            updated['bg_corrected'] = (updated['bg_corrected'] * red_factor).round(0)

                    # what IlmnIDs are in this subset?
                    ref = self.snp_man if i['snp'] == 1 else self.man
                    # and pandas won't compare NaN to NaN... so need this extra color_channel filter
                    color_channel = ref['Color_Channel'].isna() if i['Color_Channel'] is None else ref['Color_Channel'] == i['Color_Channel']
                    IlmnIDs = ref[ (ref['Infinium_Design_Type'] == i['Infinium_Design_Type']) & (color_channel) ].index

                    # NOT SURE ABOUT THIS HACK. IT WORKS, but why?
                    if probe_subset in ('oobG','oobR'):
                        # swap the data_channels
                        updated = noob_red.copy() if probe_subset == 'oobG' else noob_green.copy()
                        # ... and optionally, grab the IlmnIDs from df and match them instead. Wasn't necessary. <-- this is prone to mismatching across probe_subsets.
                        #IlmnIDs = set(list(df.sort_index().index)) & set(list(updated['IlmnID']))

                    # use update; cannot merge if the IlmnIDs overlap, to avoid NaNs and bg_corrected_x/_y cols
                    updated = updated.rename(columns={'bg_corrected':bg_column})
                    updated = updated.loc[ updated['IlmnID'].isin(IlmnIDs)]
                    # update() matches on IlmnID index, but noob_green/red have multiple meth/unmeth values per IlmnID so have to split it apart.
                    if updated['IlmnID'].duplicated().sum() != 0:
                        #raise AssertionError(f"{probe_subset} -- {part} contains duplicate IlmnIDs, after filtering")
                        # --- there are duplicate IlmnIDs here, so can't index to it until I split the M and U derived parts out
                        if 'Meth' in part:
                            used = 'M'
                        elif 'Unmeth' in part:
                            used = 'U'
                        else:
                            raise ValueError(f"Not sure how to read {part}; is it Meth or Unmeth?")
                        updated = updated.loc[ updated['used'] == used ]
                        if updated['IlmnID'].duplicated().sum() != 0:
                            raise AssertionError(f"{probe_subset} -- {part} -- {used} contains duplicate IlmnIDs, after filtering twice")
                    updated = updated.set_index('IlmnID')
                    debug_pre = df[bg_column].isna().sum()
                    df.update(updated[[bg_column]])
                    debug_post = df[bg_column].isna().sum()
                    num_updated = debug_pre - debug_post
                    if self.debug: print(f"{part} {df.shape[0]} (+{num_updated})")
                    # overwrite!
                    setattr(self, probe_subset, df)
                except Exception as e:
                    print(f"**** Sigset: error in update_probe_means: {e} ***")
            if self.starting_probe_counts.get(probe_subset) != getattr(self, probe_subset).shape[0]:
                if self.debug: LOGGER.warning(f"Update probes: {probe_subset} count changed from {self.starting_probe_counts.get(probe_subset)} to {getattr(self, probe_subset).shape[0]}")

        self.__preprocessed = True # applied by set_noob
        self.__bg_corrected = True
        self.__minfi_noob = False
        self.__linear_dye = True if red_factor is not None else False

    """
    # from raw_dataset; may no longer be needed, but kept for testing against new approach 2021
    def get_oob_controls(self, green_idat, red_idat, manifest, include_rs=True):
        ''' Out-of-bound controls are the mean intensity values for the
        channel in the opposite channel's probes (IG oob and IR oob)

.. todo::
    TEST -- does this give same output as SigSet.oobG and oobR?
        '''
        param_sets = [
            {'channel': Channel.RED, 'idat': green_idat, 'manifest': manifest},
            {'channel': Channel.GREEN, 'idat': red_idat, 'manifest': manifest},
        ]

        for channel_params in param_sets:
            channel = channel_params['channel']
            idat_dataset = channel_params['idat']
            manifest = channel_params['manifest']
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
                ).rename(columns={'mean_value': 'Meth'})
                oobG = oobG.merge(
                    probe_means, # green channel X AddresA (unmeth) channel
                    how='inner',
                    left_on='AddressA_ID',
                    right_index=True,
                    suffixes=(False, False),
                ).rename(columns={'mean_value': 'Unmeth'}).sort_values('IlmnID')
                oobG.drop(['AddressA_ID', 'AddressB_ID'], axis=1)

            if channel == Channel.GREEN:
                oobR = probes.merge(
                    probe_means, # red channel X AddressB for (meth)
                    how='inner',
                    left_on='AddressB_ID',
                    right_index=True,
                    suffixes=(False, False),
                ).rename(columns={'mean_value': 'Meth'}).sort_values('IlmnID')
                oobR = oobR.merge(
                    probe_means, # red channel X AddressA for (unmeth)
                    how='inner',
                    left_on='AddressA_ID',
                    right_index=True,
                    suffixes=(False, False),
                ).rename(columns={'mean_value': 'Unmeth'})
                oobR.drop(['AddressA_ID', 'AddressB_ID'], axis=1)
        return (oobG.sort_index(), oobR.sort_index())
    """

    # from raw_dataset
    def filter_oob_probes(self, channel, manifest, idat_dataset, include_rs=True):
        raise KeyError("filter_oob_probes replaced by (is part of) SigSet.get_oob_controls in v1.5+")

    # from MethylationDataset | called by self.set_bg_corrected for each probe subset | now part of set_bg_corrected
    # from MethylationDataset | now part of update_probe_means
    def set_noob(self, red_factor):
        """ same method as update_probe_means, but simply applies a linear correction to all RED channel values """
        #update_probe_means(self, noob_green, noob_red, red_factor)
        raise KeyError("set_noob replaced by update_probe_means in v1.5+")

    def detect_and_drop_duplicates(self):
        """ as of v1.5.0, mouse manifest includes a few probes that cause duplicate values, and breaks processing.
        So this removes them. About 5 probes in all.

        Note: This runs during SigSet__init__,
        and might fail if any of these probes are affected by inter_type_I_probe_switch(),
        which theoretically should never happen in mouse. But infer-probes affects the idat probe_means directly,
        and runs before SigSet is created in SampleDataContainer, to avoid double-reading confusion.
        """
        probe_count = 0
        # (1) look for dupes within a subset; mouse.methylated has 2 to drop
        for subset in self.subsets:
            this = getattr(self, subset)
            if this.index.duplicated().sum() > 0:
                pre = this.index.duplicated().sum()
                this = this.loc[ ~this.index.duplicated() ]
                setattr(self, subset, this)
                this = getattr(self, subset)
                probe_count += pre
                if self.debug:
                    LOGGER.info(f"Dropped duplicate probes from SigSet.{subset}: {pre} --> {this.index.duplicated().sum()}")

        # (2) look between paired subsets; the index probe names should match exactly.
        # but if idat probe_means is missing for one or the other (AddressA_ID / AddressB_ID error?)
        # drop these. mouse has 3 to drop.
        matched_sets = [
            ('methylated','unmethylated')
        ]
        # either remove the mismatched ones, or add in missing values to other datasets (assume min fluor of 1.0)
        for partA,partB in matched_sets:
            if set(getattr(self, partA).index) - set(getattr(self, partB).index) != set():
                if self.debug:
                    LOGGER.info(f"mismatched probes ({partA} - {partB}): {set(getattr(self, partA).index) - set(getattr(self, partB).index)}")
                this = getattr(self, partA)
                mismatched = list(set(getattr(self, partA).index) - set(getattr(self, partB).index))
                this = this.loc[ ~this.index.isin(mismatched) ]
                setattr(self, partA, this)
            if set(getattr(self, partB).index) - set(getattr(self, partA).index) != set():
                probe_count += len(set(getattr(self, partB).index) - set(getattr(self, partA).index))
                if self.debug:
                    LOGGER.info(f"mismatched probes ({partB} - {partA}): {set(getattr(self, partB).index) - set(getattr(self, partA).index)}")
                this = getattr(self, partB)
                mismatched = list(set(getattr(self, partB).index) - set(getattr(self, partA).index))
                this = this.loc[ ~this.index.isin(mismatched) ]
                setattr(self, partB, this)

    def check_for_probe_loss(self, stage=''):
        """Debugger runs this during processing to see where mouse probes go missing or get duplicated."""
        if stage != '' and self.debug:
            LOGGER.info(f"[{stage}]")
        for subset in self.subsets:
            if getattr(self, subset).index.duplicated().sum() > 0:
                if self.debug:
                    LOGGER.info(f"[ {getattr(self, subset).index.duplicated().sum()} duplicate probes on SigSet.{subset} ]")
            if getattr(self, subset).shape[0] != self.starting_probe_counts[subset]:
                if self.debug:
                    count_lost = self.starting_probe_counts[subset] - getattr(self, subset).shape[0]
                    LOGGER.info(f"[ {count_lost} probes lost from SigSet.{subset} ]")
