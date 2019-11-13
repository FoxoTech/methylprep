# Lib
import logging
import numpy as np
import pandas as pd
from tqdm import tqdm
from collections import Counter
from pathlib import Path
# App
from ..files import Manifest, get_sample_sheet, create_sample_sheet
from ..models import Channel
from ..utils import ensure_directory_exists
from .meth_dataset import MethylationDataset
from .postprocess import (
    calculate_beta_value,
    calculate_copy_number,
    calculate_m_value,
    consolidate_values_for_sheet
)
from .preprocess import preprocess_noob
from .raw_dataset import get_raw_datasets

__all__ = ['SampleDataContainer', 'get_manifest', 'run_pipeline', 'consolidate_values_for_sheet']


LOGGER = logging.getLogger(__name__)


def get_manifest(raw_datasets, array_type=None, manifest_filepath=None):
    """Generates a SampleSheet instance for a given directory of processed data.

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
    if array_type is None:
        array_types = {dataset.array_type for dataset in raw_datasets}
        if len(array_types) == 0:
            raise ValueError('could not identify array type from IDATs')
        elif len(array_types) != 1:
            raise ValueError('IDATs with varying array types')

        array_type = array_types.pop()

    return Manifest(array_type, manifest_filepath)


def run_pipeline(data_dir, array_type=None, export=False, manifest_filepath=None,
                 sample_sheet_filepath=None, sample_name=None,
                 betas=False, m_value=False, make_sample_sheet=False, batch_size=None,
                 save_uncorrected=False, meta_data_frame=True):
    """The main CLI processing pipeline. This does every processing step and returns a data set.

    Arguments:
        data_dir [required]
            path where idat files can be found, and samplesheet csv.
        array_type [default: autodetect]
            27k, 450k, EPIC, EPIC+
            If omitted, this will autodetect it.
        export [default: False]
            if True, exports a CSV of the processed data for each idat file in sample.
        betas
            if True, saves a pickle (beta_values.pkl) of beta values for all samples
        m_value
            if True, saves a pickle (m_values.pkl) of beta values for all samples
        manifest_filepath [optional]
            if you want to provide a custom manifest, provide the path. Otherwise, it will download
            the appropriate one for you.
        sample_sheet_filepath [optional]
            it will autodetect if ommitted.
        sample_name [optional, list]
            if you don't want to process all samples, you can specify individual as a list.
            if sample_names are specified, this will not also do batch sizes (large batches must process all samples)
        make_sample_sheet [optional]
            if True, generates a sample sheet from idat files called 'samplesheet.csv', so that processing will work.
            From CLI pass in "--no_sample_sheet" to trigger sample sheet auto-generation.
        batch_size [optional]
            if set to any integer, samples will be processed and saved in batches no greater than
            the specified batch size. This will yield multiple output files in the format of
            "beta_values_1.pkl ... beta_values_N.pkl".

    Returns:
        By default, if called as a function, a list of SampleDataContainer objects is returned.

        betas
            if True, will return a single data frame of betavalues instead of a list of SampleDataContainer objects.
            Format is a "wide matrix": columns contain probes and rows contain samples.
        m_value
            if True, will return a single data frame of m_factor values instead of a list of SampleDataContainer objects.
            Format is a "wide matrix": columns contain probes and rows contain samples.

        if batch_size is set to more than 200 samples, nothing is returned but all the files are saved. You can recreate the output by loading the files.

    Processing note:
        The sample_sheet parser will ensure every sample has a unique name and assign one (e.g. Sample1) if missing, or append a number (e.g. _1) if not unique.
        This may cause sample_sheets and processed data in dataframes to not match up. Will fix in future version."""
    LOGGER.info('Running pipeline in: %s', data_dir)
    if sample_name:
        LOGGER.info('Sample names: {0}'.format(sample_name))

    if make_sample_sheet:
        create_sample_sheet(data_dir)
    sample_sheet = get_sample_sheet(data_dir, filepath=sample_sheet_filepath)

    samples = sample_sheet.get_samples()
    if sample_sheet.renamed_fields != {}:
        show_fields = []
        for k,v in sample_sheet.renamed_fields.items():
            if v != k:
                show_fields.append(f"{k} --> {v}\n")
            else:
                show_fields.append(f"{k}\n")
        LOGGER.info(f"Found {len(show_fields)} additional fields in sample_sheet: {''.join(show_fields)}")

    batches = []
    batch = []
    sample_id_counter = 1
    if batch_size:
        if type(batch_size) != int or batch_size < 1:
            raise ValueError('batch_size must be an integer greater than 0')
        for sample in samples:
            if sample_name and sample.name not in sample_name:
                continue

            # batch uses Sample_Name, so ensure these exist
            if sample.name in (None,''):
                sample.name = f'Sample_{sample_id_counter}'
                sample_id_counter += 1
            # and are unique.
            if Counter((s.name for s in samples)).get(sample.name) > 1:
                sample.name = f'{sample.name}_{sample_id_counter}'
                sample_id_counter += 1

            if len(batch) < batch_size:
                batch.append(sample.name)
            else:
                batches.append(batch)
                batch = []
                batch.append(sample.name)
        batches.append(batch)
    else:
        for sample in samples:
            if sample_name and sample.name not in sample_name:
                continue

            # batch uses Sample_Name, so ensure these exist
            if sample.name in (None,''):
                sample.name = f'Sample_{sample_id_counter}'
                sample_id_counter += 1
            # and are unique.
            if Counter((s.name for s in samples)).get(sample.name) > 1:
                sample.name = f'{sample.name}_{sample_id_counter}'
                sample_id_counter += 1

            batch.append(sample.name)
        batches.append(batch)

    data_containers = [] # returned when this runs in interpreter
    for batch_num, batch in enumerate(batches, 1):
        raw_datasets = get_raw_datasets(sample_sheet, sample_name=batch)
        manifest = get_manifest(raw_datasets, array_type, manifest_filepath)

        batch_data_containers = []
        export_paths = set() # inform CLI user where to look
        for raw_dataset in tqdm(raw_datasets):
            data_container = SampleDataContainer(
                raw_dataset=raw_dataset,
                manifest=manifest,
                retain_uncorrected_probe_intensities=save_uncorrected,
            )

            data_container.process_all()
            batch_data_containers.append(data_container)

            if export:
                output_path = data_container.sample.get_export_filepath()
                data_container.export(output_path)
                export_paths.add(output_path)

        if betas:
            df = consolidate_values_for_sheet(batch_data_containers, postprocess_func_colname='beta_value')
            if not batch_size:
                pkl_name = 'beta_values.pkl'
            else:
                pkl_name = f'beta_values_{batch_num}.pkl'
            if df.shape[1] > df.shape[0]:
                df = df.transpose() # put probes as columns for faster loading.
            pd.to_pickle(df, Path(data_dir,pkl_name))
            LOGGER.info(f"saved {pkl_name}")
        if m_value:
            df = consolidate_values_for_sheet(batch_data_containers, postprocess_func_colname='m_value')
            if not batch_size:
                pkl_name = 'm_values.pkl'
            else:
                pkl_name = f'm_values_{batch_num}.pkl'
            if df.shape[1] > df.shape[0]:
                df = df.transpose() # put probes as columns for faster loading.
            pd.to_pickle(df, Path(data_dir,pkl_name))
            LOGGER.info(f"saved {pkl_name}")
        if export:
            LOGGER.info(f"[!] Exported results (csv) to: {export_paths}")

        # consolidating data_containers this will break with really large sample sets, so skip here.
        if batch_size and batch_size >= 200:
            continue
        data_containers.extend(batch_data_containers)

    if meta_data_frame == True:
        #sample_sheet.fields is a complete mapping of original and renamed_fields
        cols = list(sample_sheet.fields.values()) + ['Sample_ID']
        meta_frame = pd.DataFrame(columns=cols)
        field_classattr_lookup = {
        'Sentrix_ID': 'sentrix_id',
        'Sentrix_Position': 'sentrix_position',
        'Sample_Group': 'group',
        'Sample_Name': 'name',
        'Sample_Plate': 'plate',
        'Pool_ID': 'pool',
        'Sample_Well': 'well',
        'GSM_ID': 'GSM_ID',
        'Sample_Type': 'type',
        'Sub_Type': 'sub_type',
        'Control': 'is_control',
        }
        # row contains the renamed fields, and pulls in the original data from sample_sheet
        for sample in samples:
            row = {}
            for field in sample_sheet.fields.keys():
                if sample_sheet.fields[field] in field_classattr_lookup:
                    row[ sample_sheet.fields[field] ] = getattr(sample, field_classattr_lookup[sample_sheet.fields[field]] )
                elif field in sample_sheet.renamed_fields:
                    row[ sample_sheet.fields[field] ] = getattr(sample, sample_sheet.renamed_fields[field])
                else:
                    LOGGER.info(f"extra column: {field} ignored")
                #    row[ sample_sheet.fields[field] ] = getattr(sample, field)
            # add the UID that matches m_value/beta value pickles
            #... unless there's a GSM_ID too
            # appears that methylprep m_value and beta files only include ID_Position as column names.
            #if row.get('GSM_ID') != None:
            #    row['Sample_ID'] = f"{row['GSM_ID']}_{row['Sentrix_ID']}_{row['Sentrix_Position']}"
            #else:
            row['Sample_ID'] = f"{row['Sentrix_ID']}_{row['Sentrix_Position']}"
            meta_frame = meta_frame.append(row, ignore_index=True)
        meta_frame_filename = f'sample_sheet_meta_data.pkl'
        meta_frame.to_pickle(Path(data_dir,meta_frame_filename))
        LOGGER.info(f"[!] Exported meta_data to {meta_frame_filename}")

    # batch processing done; consolidate and return data. This uses much more memory, but not called if in batch mode.
    if batch_size and batch_size >= 200:
        print("Because the batch size was >200 samples, files are saved but no data objects are returned.")
        return
    elif betas:
        return consolidate_values_for_sheet(data_containers, postprocess_func_colname='beta_value')
    elif m_value:
        return consolidate_values_for_sheet(data_containers, postprocess_func_colname='m_value')
    else:
        return data_containers


class SampleDataContainer():
    """Wrapper that provides easy access to slices of data for a Sample,
    its RawDataset, and the pre-configured MethylationDataset subsets of probes.

    Arguments:
        raw_dataset {RawDataset} -- A sample's RawDataset for a single well on the processed array.
        manifest {Manifest} -- The Manifest for the correlated RawDataset's array type.
    """

    __data_frame = None

    def __init__(self, raw_dataset, manifest, retain_uncorrected_probe_intensities=False):
        self.manifest = manifest
        self.raw_dataset = raw_dataset
        self.sample = raw_dataset.sample
        self.retain_uncorrected_probe_intensities=retain_uncorrected_probe_intensities

        self.methylated = MethylationDataset.methylated(raw_dataset, manifest)
        self.unmethylated = MethylationDataset.unmethylated(raw_dataset, manifest)
        self.oob_controls = raw_dataset.get_oob_controls(manifest)

    @property
    def fg_green(self):
        return self.raw_dataset.get_fg_values(self.manifest, Channel.GREEN)

    @property
    def fg_red(self):
        return self.raw_dataset.get_fg_values(self.manifest, Channel.RED)

    @property
    def ctrl_green(self):
        return self.raw_dataset.get_fg_controls(self.manifest, Channel.GREEN)

    @property
    def ctrl_red(self):
        return self.raw_dataset.get_fg_controls(self.manifest, Channel.RED)

    @property
    def oob_green(self):
        return self.oob_controls[Channel.GREEN]

    @property
    def oob_red(self):
        return self.oob_controls[Channel.RED]

    def preprocess(self):
        """ combines the methylated and unmethylated columns from the SampleDataContainer. """
        if not self.__data_frame:
            if self.retain_uncorrected_probe_intensities == True:
                uncorrected_meth = self.methylated.data_frame.copy()
                uncorrected_unmeth = self.unmethylated.data_frame.copy()

            preprocess_noob(self) # apply corrections: bg subtract, then noob

            methylated = self.methylated.data_frame[['noob']]
            unmethylated = self.unmethylated.data_frame[['noob']]

            self.__data_frame = methylated.join(
                unmethylated,
                lsuffix='_meth',
                rsuffix='_unmeth',
            )

            if self.retain_uncorrected_probe_intensities == True:
                self.__data_frame['meth'] = uncorrected_meth['mean_value']
                self.__data_frame['unmeth'] = uncorrected_unmeth['mean_value']

            self.__data_frame = self.__data_frame.round(4)

        return self.__data_frame

    def save_uncorrected_probe_intensities(self, input_dataframe):
        vectorized_func = np.vectorize(postprocess_func)

        input_dataframe[header] = vectorized_func(
            input_dataframe['noob_meth'].values,
            input_dataframe['noob_unmeth'].values,
        )
        return self._postprocess(input_dataframe, None, 'meth')

    def process_m_value(self, input_dataframe):
        """Calculate M value from methylation data"""
        return self._postprocess(input_dataframe, calculate_m_value, 'm_value')

    def process_beta_value(self, input_dataframe):
        """Calculate Beta value from methylation data"""
        return self._postprocess(input_dataframe, calculate_beta_value, 'beta_value')

    def process_copy_number(self, input_dataframe):
        """Calculate copy number value from methylation data"""
        return self._postprocess(input_dataframe, calculate_copy_number, 'cm_value')

    def process_all(self):
        """Runs all pre and post-processing calculations for the dataset."""
        data_frame = self.preprocess()
        data_frame = self.process_m_value(data_frame)
        data_frame = self.process_beta_value(data_frame)
        self.__data_frame = data_frame
        return data_frame

    def export(self, output_path):
        ensure_directory_exists(output_path)
        self.__data_frame.to_csv(output_path)

    def _postprocess(self, input_dataframe, postprocess_func, header):
        vectorized_func = np.vectorize(postprocess_func)

        input_dataframe[header] = vectorized_func(
            input_dataframe['noob_meth'].values,
            input_dataframe['noob_unmeth'].values,
        )

        return input_dataframe
