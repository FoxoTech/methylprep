# Lib
import logging
import numpy as np
# App
from ..files import Manifest, get_sample_sheet
from ..models import Channel
from ..utils import ensure_directory_exists
from .meth_dataset import MethylationDataset
from .postprocess import (
    calculate_beta_value,
    calculate_copy_number,
    calculate_m_value,
)
from .preprocess import preprocess_noob
from .raw_dataset import get_raw_datasets


__all__ = ['SampleDataContainer', 'get_manifest', 'run_pipeline']


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

        if len(array_types) != 1:
            raise 'IDATs with varying array types'

        array_type = array_types.pop()

    return Manifest(array_type, manifest_filepath)


def run_pipeline(data_dir, array_type=None, export=False, manifest_filepath=None,
                 sample_sheet_filepath=None, sample_names=None):
    LOGGER.info('Running pipeline in: %s', data_dir)

    sample_sheet = get_sample_sheet(data_dir, filepath=sample_sheet_filepath)
    raw_datasets = get_raw_datasets(sample_sheet, sample_names=sample_names)
    manifest = get_manifest(raw_datasets, array_type, manifest_filepath)

    data_containers = []
    for raw_dataset in raw_datasets:
        data_container = SampleDataContainer(
            raw_dataset=raw_dataset,
            manifest=manifest,
        )

        data_container.process_all()
        data_containers.append(data_container)

        if export:
            output_path = data_container.sample.get_export_filepath()
            data_container.export(output_path)

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
