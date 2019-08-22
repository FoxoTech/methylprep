# Lib
import logging
import numpy as np
import pandas as pd
from tqdm import tqdm
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
                 betas=False, m_value=False, make_sample_sheet=False):
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
            if you want to not process all samples, you can specify them as a list.
        make_sample_sheet [optional]
            if True, generates a sample sheet from idat files called 'samplesheet.csv', so that processing will work.
            From CLI pass in "--no_sample_sheet" to trigger sample sheet auto-generation.

    Returns:
        By default, a list of SampleDataContainer objects are returned.

        betas
            if True, will return a single data frame of betavalues instead of a list of SampleDataContainer objects.
            Format is a "wide matrix": columns contain probes and rows contain samples.
        m_factor
            if True, will return a single data frame of m_factor values instead of a list of SampleDataContainer objects.
            Format is a "wide matrix": columns contain probes and rows contain samples."""

    LOGGER.info('Running pipeline in: %s', data_dir)

    if make_sample_sheet:
        create_sample_sheet(data_dir)
    sample_sheet = get_sample_sheet(data_dir, filepath=sample_sheet_filepath)
    raw_datasets = get_raw_datasets(sample_sheet, sample_name=sample_name)
    manifest = get_manifest(raw_datasets, array_type, manifest_filepath)

    data_containers = []
    export_paths = set() # inform CLI user where to look
    for raw_dataset in tqdm(raw_datasets):
        data_container = SampleDataContainer(
            raw_dataset=raw_dataset,
            manifest=manifest,
        )

        data_container.process_all()
        data_containers.append(data_container)

        if export:
            output_path = data_container.sample.get_export_filepath()
            data_container.export(output_path)
            export_paths.add(output_path)

    if betas:
        df = consolidate_values_for_sheet(data_containers, postprocess_func_colname='beta_value')
        pd.to_pickle(df, 'beta_values.pkl')
        LOGGER.info("saved beta_values.pkl")
        return df
    if m_value:
        df = consolidate_values_for_sheet(data_containers, postprocess_func_colname='m_value')
        pd.to_pickle(df,'m_values.pkl')
        LOGGER.info("saved m_values.pkl")
        return df
    if export:
        # not using LOGGER because this should appear regardless of verbose flag.
        # print(f"[!] Exported results (csv) to: {export_paths}")
        # requires --verbose too.
        LOGGER.info(f"[!] Exported results (csv) to: {export_paths}")
    return data_containers


class SampleDataContainer():
    """Wrapper that provides easy access to slices of data for a Sample,
    its RawDataset, and the pre-configured MethylationDataset subsets of probes.

    Arguments:
        raw_dataset {RawDataset} -- A sample's RawDataset for a single well on the processed array.
        manifest {Manifest} -- The Manifest for the correlated RawDataset's array type.
    """

    __data_frame = None

    def __init__(self, raw_dataset, manifest):
        self.manifest = manifest
        self.raw_dataset = raw_dataset
        self.sample = raw_dataset.sample

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
            preprocess_noob(self)

            methylated = self.methylated.data_frame[['noob']]
            unmethylated = self.unmethylated.data_frame[['noob']]

            self.__data_frame = methylated.join(
                unmethylated,
                lsuffix='_meth',
                rsuffix='_unmeth',
            )

        return self.__data_frame

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
