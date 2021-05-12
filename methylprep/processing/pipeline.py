# Lib
import logging
import numpy as np
import pandas as pd
from ..utils.progress_bar import * # checks environment and imports tqdm appropriately.
from collections import Counter
from pathlib import Path
import pickle
# App
from ..files import Manifest, get_sample_sheet, create_sample_sheet
from ..models import (
    Channel,
    MethylationDataset,
    ArrayType,
    get_raw_datasets,
    get_array_type
)
from .postprocess import (
    calculate_beta_value,
    calculate_m_value,
    calculate_copy_number,
    consolidate_values_for_sheet,
    consolidate_control_snp,
    one_sample_control_snp,
    consolidate_mouse_probes,
    merge_batches,
)
from ..utils import ensure_directory_exists, is_file_like
from .preprocess import preprocess_noob, preprocess_noob_sesame, _apply_sesame_quality_mask
from .p_value_probe_detection import _pval_sesame_preprocess
from .infer_channel_switch import infer_type_I_probes
from .dye_bias import nonlinear_dye_bias_correction


__all__ = ['SampleDataContainer', 'get_manifest', 'run_pipeline', 'consolidate_values_for_sheet', 'make_pipeline']

LOGGER = logging.getLogger(__name__)


def get_manifest(raw_datasets, array_type=None, manifest_filepath=None):
    """Return a Manifest, given a list of raw_datasets (from idats).

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
        array_type = get_array_type(raw_datasets)
    return Manifest(array_type, manifest_filepath)


def run_pipeline(data_dir, array_type=None, export=False, manifest_filepath=None,
                 sample_sheet_filepath=None, sample_name=None,
                 betas=False, m_value=False, make_sample_sheet=False, batch_size=None,
                 save_uncorrected=False, save_control=True, meta_data_frame=True,
                 bit='float32', poobah=False, export_poobah=False,
                 poobah_decimals=3, poobah_sig=0.05, low_memory=True,
                 sesame=True, quality_mask=None, **kwargs):
    """The main CLI processing pipeline. This does every processing step and returns a data set.

    Required Arguments:
        data_dir [required]
            path where idat files can be found, and a samplesheet csv.

    Optional file and sub-sampling inputs:
        manifest_filepath [optional]
            if you want to provide a custom manifest, provide the path. Otherwise, it will download
            the appropriate one for you.
        sample_sheet_filepath [optional]
            it will autodetect if ommitted.
        make_sample_sheet [optional]
            if True, generates a sample sheet from idat files called 'samplesheet.csv', so that processing will work.
            From CLI pass in "--no_sample_sheet" to trigger sample sheet auto-generation.
        sample_name [optional, list]
            if you don't want to process all samples, you can specify individual samples as a list.
            if sample_names are specified, this will not also do batch sizes (large batches must process all samples)

    Optional processing arguments:
        sesame [default: True]
            If True, applies offsets, poobah, noob, infer_channel_switch, nonlinear-dye-bias-correction, and qualityMask to imitate the output of openSesame function.
            If False, outputs will closely match minfi's processing output.
            Prior to version 1.4.0, file processing matched minfi.
        array_type [default: autodetect]
            27k, 450k, EPIC, EPIC+
            If omitted, this will autodetect it.
        batch_size [optional]
            if set to any integer, samples will be processed and saved in batches no greater than
            the specified batch size. This will yield multiple output files in the format of
            "beta_values_1.pkl ... beta_values_N.pkl".
        bit [default: float32]
            You can change the processed output files to one of: {float16, float32, float64}.
            This will make files & memory usage smaller, often with no loss in precision.
            However, using float16 masy cause an overflow error, resulting in "inf" appearing instead of numbers, and numpy/pandas functions do not universally support float16.
        low_memory [default: True]
            If False, pipeline will not remove intermediate objects and data sets during processing.
            This provides access to probe subsets, foreground, and background probe sets in the
            SampleDataContainer object returned when this is run in a notebook (not CLI).
        quality_mask [default: None]
            If False, process will NOT remove sesame's list of unreliable probes.
            If True, removes probes.
            The default None will defer to sesamee, which defaults to true. But if explicitly set, it will override sesame setting.

    Optional export files:
        meta_data_frame [default: True]
            if True, saves a file, "sample_sheet_meta_data.pkl" with samplesheet info.
        export [default: False]
            if True, exports a CSV of the processed data for each idat file in sample.
        save_uncorrected [default: False]
            if True, adds two additional columns to the processed.csv per sample (meth and unmeth),
            representing the raw fluorescence intensities for all probes.
            It does not apply NOOB correction to values in these columns.
        save_control [default: False]
            if True, adds all Control and SnpI type probe values to a separate pickled dataframe,
            with probes in rows and sample_name in the first column.
            These non-CpG probe names are excluded from processed data and must be stored separately.
        poobah [default: False]
            If specified as True, the pipeline will run Sesame's p-value probe detection method (poobah)
            on samples to remove probes that fail the signal/noise ratio on their fluorescence channels.
            These will appear as NaNs in the resulting dataframes (beta_values.pkl or m_values.pkl).
            All probes, regardless of p-value cutoff, will be retained in CSVs, but there will be a 'poobah_pval'
            column in CSV files that methylcheck.load uses to exclude failed probes upon import at a later step.
        poobah_sig [default: 0.05]
            the p-value level of significance, above which, will exclude probes from output (typical range of 0.001 to 0.1)
        poobah_decimals [default: 3]
            The number of decimal places to round p-value column in the processed CSV output files.
        mouse probes
            Mouse-specific will be saved if processing a mouse array.

    Optional final estimators:
        betas
            if True, saves a pickle (beta_values.pkl) of beta values for all samples
        m_value
            if True, saves a pickle (m_values.pkl) of beta values for all samples
        Note on meth/unmeth:
            if either betas or m_value is True, this will also save two additional files:
            'meth_values.pkl' and 'unmeth_values.pkl' with the same dataframe structure,
            representing raw, uncorrected meth probe intensities for all samples. These are useful
            in some methylcheck functions and load/produce results 100X faster than loading from
            processed CSV output.

    Returns:
        By default, if called as a function, a list of SampleDataContainer objects is returned, with the following execptions:

        betas
            if True, will return a single data frame of betavalues instead of a list of SampleDataContainer objects.
            Format is a "wide matrix": columns contain probes and rows contain samples.
        m_value
            if True, will return a single data frame of m_factor values instead of a list of SampleDataContainer objects.
            Format is a "wide matrix": columns contain probes and rows contain samples.
        if batch_size is set to more than ~600 samples, nothing is returned but all the files are saved. You can recreate/merge output files by loading the files using methylcheck.load().

    Processing notes:
        The sample_sheet parser will ensure every sample has a unique name and assign one (e.g. Sample1) if missing, or append a number (e.g. _1) if not unique.
        This may cause sample_sheets and processed data in dataframes to not match up. Will fix in future version.

        pipeline steps:
            1 make sample sheet or read sample sheet into a list of samples' data
            2 split large projects into batches, if necessary, and ensure unique sample names
            3 read idats
            4 select and read manifest
            5 put everything into SampleDataContainer class objects
            6 process everything, using the pipeline steps specified
                idats -> channel_swaps -> poobah -> quality_mask -> noob -> dye_bias
            7 apply the final estimator function (beta, m_value, or copy number) to all data
            8 export all the data into multiple files, as defined by pipeline
        """
    # support for the make_pipeline wrapper function here; a more structured way to pass in args like sklearn.
    # unexposed flags all start with 'do_': (None will retain default settings)
    do_infer_channel_switch = None # defaults to sesame(True)
    do_noob = None # defaults to True
    do_dye_bias = None # defaults to sesame(True)
    do_save_noob = None
    do_mouse = True
    hidden_kwargs = ['pipeline_steps', 'pipeline_exports', 'debug']
    if kwargs != {}:
        for kwarg in kwargs:
            if kwarg not in hidden_kwargs:
                raise ValueError(f"One of your parameters ({kwarg}) was not recognized. Did you misspell it?")
    if sesame == True:
        poobah = True

    if kwargs != {} and 'pipeline_steps' in kwargs:
        pipeline_steps = kwargs.get('pipeline_steps')
        if 'all' in pipeline_steps:
            do_infer_channel_switch = True
            poobah = True
            quality_mask = True
            do_noob = True
            do_dye_bias = True
            sesame = None # prevent this from overriding elsewhere
        else:
            do_infer_channel_switch = True if 'infer_channel_switch' in pipeline_steps else False
            poobah = True if 'poobah' in pipeline_steps else False
            quality_mask = True if 'quality_mask' in pipeline_steps else False
            do_noob = True if 'noob' in pipeline_steps else False
            do_dye_bias = True if 'dye_bias' in pipeline_steps else False
            sesame = None
    if kwargs != {} and 'pipeline_exports' in kwargs:
        pipeline_exports = kwargs.get('pipeline_exports')
        if 'all' in pipeline_exports:
            export = True # csv
            save_uncorrected = True # meth, unmeth
            do_save_noob = True # noob_meth, noob_unmeth
            export_poobah = True # poobah
            meta_data_frame = True # sample_sheet_meta_data
            do_mouse = True # 'mouse' -- only if array_type matches; False will suppress export
            save_control = True # control
        else:
            export = True if 'csv' in pipeline_exports else False
            save_uncorrected = True if ('meth' in pipeline_exports or 'unmeth' in pipeline_exports) else False
            do_save_noob = True if ('noob_meth' in pipeline_exports or 'noob_unmeth' in pipeline_exports) else False
            export_poobah = True if 'poobah' in pipeline_exports else False
            meta_data_frame = True if 'sample_sheet_meta_data' in pipeline_exports else False
            # mouse is determined by the array_type match, but you can suppress creating this file here
            do_mouse = True if 'mouse' in pipeline_exports else False
            save_control = True if 'control' in pipeline_exports else False


    LOGGER.info('Running pipeline in: %s', data_dir)
    if bit not in ('float64','float32','float16'):
        raise ValueError("Input 'bit' must be one of ('float64','float32','float16') or ommitted.")
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
                show_fields.append(f"{k} --> {v}")
            else:
                show_fields.append(f"{k}")
        LOGGER.info(f"Found {len(show_fields)} additional fields in sample_sheet:\n{' | '.join(show_fields)}")

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

    temp_data_pickles = []
    control_snps = {}
    #data_containers = [] # returned when this runs in interpreter, and < 200 samples
    # v1.3.0 memory fix: save each batch_data_containers object to disk as temp, then load and combine at end.
    # 200 samples still uses 4.8GB of memory/disk space (float64)
    missing_probe_errors = {'noob': [], 'raw':[]}

    for batch_num, batch in enumerate(batches, 1):
        raw_datasets = get_raw_datasets(sample_sheet, sample_name=batch)
        manifest = get_manifest(raw_datasets, array_type, manifest_filepath) # this allows each batch to be a different array type; but not implemented yet. common with older GEO sets.

        batch_data_containers = []
        export_paths = set() # inform CLI user where to look
        for raw_dataset in tqdm(raw_datasets, total=len(raw_datasets), desc="Processing samples"):
            data_container = SampleDataContainer(
                raw_dataset=raw_dataset,
                manifest=manifest,
                retain_uncorrected_probe_intensities=save_uncorrected,
                bit=bit,
                switch_probes=(do_infer_channel_switch or sesame), # this applies all sesame-specific options
                quality_mask= (quality_mask or sesame or False), # this applies all sesame-specific options (beta / noob offsets too)
                do_noob=(do_noob if do_noob != None else True), # None becomes True, but make_pipeline can override with False
                pval=poobah, #defaults to False as of v1.4.0
                poobah_decimals=poobah_decimals,
                poobah_sig=poobah_sig,
                correct_dye_bias=(do_dye_bias or sesame or False), # this applies all sesame-specific options
                debug=kwargs.get('debug',False),
                sesame=sesame,
            )
            data_container.process_all()

            if export:
                output_path = data_container.sample.get_export_filepath()
                data_container.export(output_path)
                export_paths.add(output_path)
                # this tidies-up the tqdm by moving errors to end of batch warning.
                if data_container.noob_processing_missing_probe_errors != []:
                    missing_probe_errors['noob'].extend(data_container.noob_processing_missing_probe_errors)
                if data_container.raw_processing_missing_probe_errors != []:
                    missing_probe_errors['raw'].extend(data_container.raw_processing_missing_probe_errors)

            if save_control: # Process and consolidate now. Keep in memory. These files are small.
                sample_id = f"{data_container.sample.sentrix_id}_{data_container.sample.sentrix_position}"
                control_df = one_sample_control_snp(data_container)
                control_snps[sample_id] = control_df

            # now I can drop all the unneeded stuff from each SampleDataContainer (400MB per sample becomes 92MB)
            # these are stored in SampleDataContainer.__data_frame for processing.
            if low_memory is True:
                # use data_frame values instead of class objects
                del data_container.manifest
                del data_container.raw_dataset
                del data_container.methylated
                del data_container.unmethylated
            batch_data_containers.append(data_container)

            #if str(data_container.sample) == '200069280091_R01C01':
            #    print(f"200069280091_R01C01 -- cg00035864 -- meth -- {data_container._SampleDataContainer__data_frame['meth']['cg00035864']}")
            #    print(f"200069280091_R01C01 -- cg00035864 -- unmeth -- {data_container._SampleDataContainer__data_frame['unmeth']['cg00035864']}")

        LOGGER.info('[finished SampleDataContainer processing]')

        if betas:
            df = consolidate_values_for_sheet(batch_data_containers, postprocess_func_colname='beta_value', bit=bit, poobah=poobah)
            if not batch_size:
                pkl_name = 'beta_values.pkl'
            else:
                pkl_name = f'beta_values_{batch_num}.pkl'
            if df.shape[1] > df.shape[0]:
                df = df.transpose() # put probes as columns for faster loading.
            df = df.astype('float32')
            df = df.sort_index().reindex(sorted(df.columns), axis=1)
            pd.to_pickle(df, Path(data_dir,pkl_name))
            LOGGER.info(f"saved {pkl_name}")
        if m_value:
            df = consolidate_values_for_sheet(batch_data_containers, postprocess_func_colname='m_value', bit=bit, poobah=poobah)
            if not batch_size:
                pkl_name = 'm_values.pkl'
            else:
                pkl_name = f'm_values_{batch_num}.pkl'
            if df.shape[1] > df.shape[0]:
                df = df.transpose() # put probes as columns for faster loading.
            df = df.astype('float32')
            df = df.sort_index().reindex(sorted(df.columns), axis=1)
            pd.to_pickle(df, Path(data_dir,pkl_name))
            LOGGER.info(f"saved {pkl_name}")
        if (do_save_noob is not False) or betas or m_value:
            df = consolidate_values_for_sheet(batch_data_containers, postprocess_func_colname='noob_meth', bit=bit)
            if not batch_size:
                pkl_name = 'noob_meth_values.pkl'
            else:
                pkl_name = f'noob_meth_values_{batch_num}.pkl'
            if df.shape[1] > df.shape[0]:
                df = df.transpose() # put probes as columns for faster loading.
            df = df.astype('float32') if df.isna().sum().sum() > 0 else df.astype('uint16')
            df = df.sort_index().reindex(sorted(df.columns), axis=1)
            pd.to_pickle(df, Path(data_dir,pkl_name))
            LOGGER.info(f"saved {pkl_name}")
            # TWO PARTS
            df = consolidate_values_for_sheet(batch_data_containers, postprocess_func_colname='noob_unmeth', bit=bit)
            if not batch_size:
                pkl_name = 'noob_unmeth_values.pkl'
            else:
                pkl_name = f'noob_unmeth_values_{batch_num}.pkl'
            if df.shape[1] > df.shape[0]:
                df = df.transpose() # put probes as columns for faster loading.
            df = df.astype('float32') if df.isna().sum().sum() > 0 else df.astype('uint16')
            df = df.sort_index().reindex(sorted(df.columns), axis=1)
            pd.to_pickle(df, Path(data_dir,pkl_name))
            LOGGER.info(f"saved {pkl_name}")

        #if (betas or m_value) and save_uncorrected:
        if save_uncorrected:
            df = consolidate_values_for_sheet(batch_data_containers, postprocess_func_colname='meth', bit=bit)
            if not batch_size:
                pkl_name = 'meth_values.pkl'
            else:
                pkl_name = f'meth_values_{batch_num}.pkl'
            if df.shape[1] > df.shape[0]:
                df = df.transpose() # put probes as columns for faster loading.
            df = df.astype('float32') if df.isna().sum().sum() > 0 else df.astype('uint16')
            df = df.sort_index().reindex(sorted(df.columns), axis=1)
            pd.to_pickle(df, Path(data_dir,pkl_name))
            LOGGER.info(f"saved {pkl_name}")
            # TWO PARTS
            df = consolidate_values_for_sheet(batch_data_containers, postprocess_func_colname='unmeth', bit=bit)
            if not batch_size:
                pkl_name = 'unmeth_values.pkl'
            else:
                pkl_name = f'unmeth_values_{batch_num}.pkl'
            if df.shape[1] > df.shape[0]:
                df = df.transpose() # put probes as columns for faster loading.
            df = df.astype('float32') if df.isna().sum().sum() > 0 else df.astype('uint16')
            df = df.sort_index().reindex(sorted(df.columns), axis=1)
            pd.to_pickle(df, Path(data_dir,pkl_name))
            LOGGER.info(f"saved {pkl_name}")

        if manifest.array_type == ArrayType.ILLUMINA_MOUSE and do_mouse:
            # save mouse specific probes
            if not batch_size:
                mouse_probe_filename = f'mouse_probes.pkl'
            else:
                mouse_probe_filename = f'mouse_probes_{batch_num}.pkl'
            consolidate_mouse_probes(batch_data_containers, Path(data_dir, mouse_probe_filename))
            LOGGER.info(f"saved {mouse_probe_filename}")

        if export:
            export_path_parents = list(set([str(Path(e).parent) for e in export_paths]))
            LOGGER.info(f"[!] Exported results (csv) to: {export_path_parents}")

        if export_poobah:
            if all(['poobah_pval' in e._SampleDataContainer__data_frame.columns for e in batch_data_containers]):
                # this option will save a pickled dataframe of the pvalues for all samples, with sample_ids in the column headings and probe names in index.
                # this sets poobah to false in kwargs, otherwise some pvalues would be NaN I think.
                df = consolidate_values_for_sheet(batch_data_containers, postprocess_func_colname='poobah_pval', bit=bit, poobah=False, poobah_sig=poobah_sig)
                if not batch_size:
                    pkl_name = 'poobah_values.pkl'
                else:
                    pkl_name = f'poobah_values_{batch_num}.pkl'
                if df.shape[1] > df.shape[0]:
                    df = df.transpose() # put probes as columns for faster loading.
                df = df.sort_index().reindex(sorted(df.columns), axis=1)
                pd.to_pickle(df, Path(data_dir,pkl_name))
                LOGGER.info(f"saved {pkl_name}")

        # v1.3.0 fixing mem probs: pickling each batch_data_containers object then reloading it later.
        # consolidating data_containers this will break with really large sample sets, so skip here.
        #if batch_size and batch_size >= 200:
        #    continue
        #data_containers.extend(batch_data_containers)

        pkl_name = f"_temp_data_{batch_num}.pkl"
        with open(Path(data_dir,pkl_name), 'wb') as temp_data:
            pickle.dump(batch_data_containers, temp_data)
            temp_data_pickles.append(pkl_name)

    del batch_data_containers

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
        LOGGER.info(f"saved {meta_frame_filename}")

    # FIXED in v1.3.0
    # moved consolidate_control_snp() from this spot to earlier in pipeline, because it uses
    # raw_dataset and this gets removed before pickling _temp files. Here I pickle.dump the SNPS.
    if save_control:
        control_filename = f'control_probes.pkl'
        with open(Path(data_dir, control_filename), 'wb') as control_file:
            pickle.dump(control_snps, control_file)
        LOGGER.info(f"saved {control_filename}")

    # summarize any processing errors
    if missing_probe_errors['noob'] != []:
        avg_missing_per_sample = int(round(sum([item[1] for item in missing_probe_errors['noob']])/len(missing_probe_errors['noob'])))
        samples_affected = len(set([item[0] for item in missing_probe_errors['noob']]))
        LOGGER.warning(f"{samples_affected} samples were missing (or had infinite values) NOOB meth/unmeth probe values (average {avg_missing_per_sample} per sample)")
    if missing_probe_errors['raw'] != []:
        avg_missing_per_sample = int(round(sum([item[1] for item in missing_probe_errors['raw']])/len(missing_probe_errors['raw'])))
        samples_affected = len(set([item[0] for item in missing_probe_errors['raw']]))
        LOGGER.warning(f"{samples_affected} samples were missing (or had infinite values) RAW meth/unmeth probe values (average {avg_missing_per_sample} per sample)")

    # batch processing done; consolidate and return data. This uses much more memory, but not called if in batch mode.
    if batch_size and batch_size >= 200:
        LOGGER.warning("Because the batch size was >=200 samples, files are saved but no data objects are returned.")
        del batch_data_containers
        for temp_data in temp_data_pickles:
            temp_file = Path(data_dir, temp_data)
            temp_file.unlink(missing_ok=True) # delete it
        return

    # consolidate batches and delete parts, if possible
    for file_type in ['beta_values', 'm_values', 'meth_values', 'unmeth_values',
        'noob_meth_values', 'noob_unmeth_values', 'mouse_probes', 'poobah_values']: # control_probes.pkl not included yet
        test_parts = list([str(temp_file) for temp_file in Path(data_dir).rglob(f'{file_type}*.pkl')])
        num_batches = len(test_parts)
        # ensures that only the file_types that appear to be selected get merged.
        #print(f"DEBUG num_batches {num_batches}, batch_size {batch_size}, file_type {file_type}")
        if batch_size and num_batches >= 1: #--- if the batch size was larger than the number of total samples, this will still drop the _1
            merge_batches(num_batches, data_dir, file_type)

    # reload all the big stuff -- after everything important is done.
    # attempts to consolidate all the batch_files below, if they'll fit in memory.
    data_containers = []
    for temp_data in temp_data_pickles:
        temp_file = Path(data_dir, temp_data)
        if temp_file.exists(): #possibly user deletes file while processing, since these are big
            with open(temp_file,'rb') as _file:
                batch_data_containers = pickle.load(_file)
                data_containers.extend(batch_data_containers)
                del batch_data_containers
            temp_file.unlink() # delete it after loading.

    if betas:
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
        bit (default: float64) -- option to store data as float16 or float32 to save space.
        pval (default: False) -- whether to apply p-value-detection algorithm to remove
            unreliable probes (based on signal/noise ratio of fluoresence)
            uses the sesame method (pOOBah) based on out of band background levels

    Jan 2020: added .snp_(un)methylated property. used in postprocess.consolidate_crontrol_snp()
    Mar 2020: added p-value detection option
    Mar 2020: added mouse probe post-processing separation
    """

    __data_frame = None
    noob_processing_missing_probe_errors = []
    raw_processing_missing_probe_errors = []

    def __init__(self, raw_dataset, manifest, retain_uncorrected_probe_intensities=False,
                 bit='float32', pval=False, poobah_decimals=3, poobah_sig=0.05, do_noob=True,
                 quality_mask=True, switch_probes=True, correct_dye_bias=True, debug=False, sesame=True):
        self.debug = debug
        self.manifest = manifest
        self.do_noob = do_noob
        self.pval = pval
        self.poobah_decimals = poobah_decimals
        self.poobah_sig = poobah_sig
        self.quality_mask = quality_mask # if True, filters sesame's standard sketchy probes out of 450k, EPIC, EPIC+ arrays.
        self.switch_probes = switch_probes
        self.correct_dye_bias = correct_dye_bias
        self.raw_dataset = raw_dataset
        self.sample = raw_dataset.sample
        self.retain_uncorrected_probe_intensities=retain_uncorrected_probe_intensities
        self.sesame = sesame # defines offsets in functions
        if debug:
            print(f'DEBUGL: sesame {self.sesame} noob {self.do_noob} poobah {self.pval} mask {self.quality_mask}, switch {self.switch_probes}, dye {self.correct_dye_bias}')

        if self.switch_probes:
            # apply inter_channel_switch here; uses raw_dataset and manifest only; then updates self.raw_dataset
            infer_type_I_probes(self, debug=self.debug)

        self.methylated = MethylationDataset.methylated(raw_dataset, manifest)
        self.unmethylated = MethylationDataset.unmethylated(raw_dataset, manifest)
        self.snp_methylated = MethylationDataset.snp_methylated(raw_dataset, manifest)
        self.snp_unmethylated = MethylationDataset.snp_unmethylated(raw_dataset, manifest)
        # mouse probes are processed within the normals meth/unmeth sets, then split at end of preprocessing step.
        #self.mouse_methylated = MethylationDataset.mouse_methylated(raw_dataset, manifest)
        #self.mouse_unmethylated = MethylationDataset.mouse_unmethylated(raw_dataset, manifest)

        # these are read from idats directly, so need to be modified at source
        # appears that because they are IG and IR, these oob_controls get updated as part of dye bias.
        self.data_type = bit #(float64, float32, or float16)
        if self.data_type == None:
            self.data_type = 'float32'
        if self.data_type not in ('float64','float32','float16'):
            raise ValueError(f"invalid data_type: {self.data_type} should be one of ('float64','float32','float16')")
        self.get_oob_controls = self.raw_dataset.get_oob_controls(self.manifest, include_rs=True)


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
    def oobG(self):
        return self.get_oob_controls[Channel.GREEN] # includes rs probes
    @property
    def oobR(self):
        return self.get_oob_controls[Channel.RED] # includes rs probes

    '''
    def oobG(self, include_rs=True):
        """ exactly like sesame, but columns are 'meth' and 'unmeth' not 'M' | 'U' """  #was oob_green until v1.4;
        # converted from @property to function because these probes change value during processing, and want to be sure it re-samples idats each time
        return self.raw_dataset.get_oob_controls(self.manifest, include_rs=include_rs)[Channel.GREEN]

    def oobR(self, include_rs=True):
        """ exactly like sesame, but columns are 'meth' and 'unmeth' not 'M' | 'U' """  #was oob_green until v1.4;
        return self.raw_dataset.get_oob_controls(self.manifest, include_rs=include_rs)[Channel.RED]
    '''

    @property
    def snp_IR(self):
        """ used by dye-bias to copy IR 'rs' probes into @IR"""
        #manifest = self.manifest._Manifest__snp_data_frame[['Infinium_Design_Type','Color_Channel']]
        #man_IR = (manifest[ (manifest['Infinium_Design_Type'] == 'I') & (manifest['Color_Channel'] == 'Red')]['IlmnID'])
        #probes = self._SampleDataContainer__data_frame[['noob_meth','noob_unmeth']]
        #return probes[probes.index.isin(man_IR)]
        probes = self.snp_methylated.data_frame
        meth = probes[(probes['Infinium_Design_Type'] == 'I') & (probes['Color_Channel'] == 'Red')][['mean_value']].rename(columns={'mean_value':'meth'})
        probes = self.snp_unmethylated.data_frame
        unmeth = probes[(probes['Infinium_Design_Type'] == 'I') & (probes['Color_Channel'] == 'Red')][['mean_value']].rename(columns={'mean_value':'unmeth'})
        if len(unmeth) == 0: # mouse appears to lack IR + IG unmeth probes
            return meth
        return pd.merge(left=meth, right=unmeth, on='IlmnID')

    @property
    def snp_IG(self):
        """ used by dye-bias to copy IG 'rs' probes into @IG"""
        #manifest = self.manifest._Manifest__snp_data_frame[['Infinium_Design_Type','Color_Channel']]
        #man_IR = (manifest[ (manifest['Infinium_Design_Type'] == 'I') & (manifest['Color_Channel'] == 'Red')]['IlmnID'])
        #probes = self._SampleDataContainer__data_frame[['noob_meth','noob_unmeth']]
        #return probes[probes.index.isin(man_IR)]
        probes = self.snp_methylated.data_frame
        meth = probes[(probes['Infinium_Design_Type'] == 'I') & (probes['Color_Channel'] == 'Grn')][['mean_value']].rename(columns={'mean_value':'meth'})
        probes = self.snp_unmethylated.data_frame
        unmeth = probes[(probes['Infinium_Design_Type'] == 'I') & (probes['Color_Channel'] == 'Grn')][['mean_value']].rename(columns={'mean_value':'unmeth'})
        if len(unmeth) == 0:
            return meth
        return pd.merge(left=meth, right=unmeth, on='IlmnID')

    @property
    def II(self):
        """ research function to match sesame's II function; not used in processing.
        only works if save_uncorrected=True. """
        # properties are class attributes and take no kwargs, so you can't customize what's returned with 'sub'
        # sub None: all II probes; 'meth': meth II probes; 'unmeth': unmeth II probes """
        manifest = self.manifest.data_frame[['Infinium_Design_Type','Color_Channel']]
        man_II = manifest[manifest['Infinium_Design_Type']=='II']
        #if not sub: # all II probes
        probes = self._SampleDataContainer__data_frame[['noob_meth','noob_unmeth']]
        return pd.merge(left=man_II, right=probes, on='IlmnID').drop(columns=['Infinium_Design_Type','Color_Channel'])
        #else: # 'meth' or 'unmeth'
        #    probes = self._SampleDataContainer__data_frame[sub]
        #    return pd.merge(left=man_II, right=probes, on='IlmnID').drop(columns=['Infinium_Design_Type','Color_Channel'])

    @property
    def IR(self):
        """ research function to match sesame's IR function; not used in processing
        only works if save_uncorrected=True. """
        manifest = self.manifest.data_frame[['Infinium_Design_Type','Color_Channel']]
        man_IR = manifest[(manifest['Color_Channel']=='Red') & (manifest['Infinium_Design_Type']=='I')]
        probes = self._SampleDataContainer__data_frame[['noob_meth','noob_unmeth']]
        return pd.merge(left=man_IR, right=probes, on='IlmnID').drop(columns=['Infinium_Design_Type','Color_Channel'])

    @property
    def IG(self, sub=None):
        """ research function to match sesame's IG function; not used in processing
        only works if save_uncorrected=True. """
        manifest = self.manifest.data_frame[['Infinium_Design_Type','Color_Channel']]
        man_IG = manifest[(manifest['Color_Channel']=='Grn') & (manifest['Infinium_Design_Type']=='I')]
        probes = self._SampleDataContainer__data_frame[['noob_meth','noob_unmeth']]
        return pd.merge(left=man_IG, right=probes, on='IlmnID').drop(columns=['Infinium_Design_Type','Color_Channel'])

    @property
    def raw_IG(self):
        """Uncorrected type-I GREEN probes from idats. only works if save_uncorrected=True; should match sesame SigSet.IG output before noob or any other transformations.
        note that running dye_bias_correction will modify the raw idat probe means in memory (changing this)"""
        return pd.DataFrame(data={
        'meth':self.methylated.data_frame[(self.methylated.data_frame['Infinium_Design_Type'] == 'I') & (self.methylated.data_frame['Color_Channel'] == 'Grn')]['mean_value'].astype(int),
        'unmeth':self.unmethylated.data_frame[(self.unmethylated.data_frame['Infinium_Design_Type'] == 'I') & (self.unmethylated.data_frame['Color_Channel'] == 'Grn')]['mean_value'].astype(int)
        }).sort_index()
    @property
    def raw_IR(self):
        """Uncorrected type-I RED probes from idats. only works if save_uncorrected=True; should match sesame SigSet.IG output before noob or any other transformations.
        note that running dye_bias_correction will modify the raw idat probe means in memory (changing this)"""
        return pd.DataFrame(data={
        'meth':self.methylated.data_frame[(self.methylated.data_frame['Infinium_Design_Type'] == 'I') & (self.methylated.data_frame['Color_Channel'] == 'Red')]['mean_value'].astype(int),
        'unmeth':self.unmethylated.data_frame[(self.unmethylated.data_frame['Infinium_Design_Type'] == 'I') & (self.unmethylated.data_frame['Color_Channel'] == 'Red')]['mean_value'].astype(int)
        }).sort_index()
    @property
    def raw_II(self):
        """Uncorrected type-II probes from idats. only works if save_uncorrected=True; should match sesame SigSet.IG output before noob or any other transformations.
        note that running dye_bias_correction will modify the raw idat probe means in memory (changing this)"""
        return pd.DataFrame(data={
        'meth':self.methylated.data_frame[self.methylated.data_frame['Infinium_Design_Type'] == 'II']['mean_value'].astype(int),
        'unmeth':self.unmethylated.data_frame[self.unmethylated.data_frame['Infinium_Design_Type'] == 'II']['mean_value'].astype(int)
        }).sort_index()

    def preprocess(self):
        """ combines the methylated and unmethylated columns from the SampleDataContainer. """
        if not self.__data_frame:
            if self.retain_uncorrected_probe_intensities == True:
                uncorrected_meth = self.methylated.data_frame.copy()['mean_value'].astype('float32')
                uncorrected_unmeth = self.unmethylated.data_frame.copy()['mean_value'].astype('float32')
                # could be uint16, if missing values didn't happen (cuts file size in half)
                if uncorrected_meth.isna().sum() == 0 and uncorrected_unmeth.isna().sum() == 0:
                    # Downcasting to 'unsigned' uses the smallest possible integer that can hold the values
                    uncorrected_meth = pd.to_numeric(uncorrected_meth, downcast='unsigned') #astype('uint16')
                    uncorrected_unmeth = pd.to_numeric(uncorrected_unmeth, downcast='unsigned') #astype('uint16')

            if self.pval == True:
                pval_probes_df = _pval_sesame_preprocess(self)
                # output: df with one column named 'poobah_pval'
            if self.quality_mask:
                quality_mask_df = _apply_sesame_quality_mask(self)
                # output: df with one column named 'quality_mask'
                # if not right array type, or custom array, returns nothing.

            if self.do_noob == True:
                # apply corrections: bg subtract, then noob (in preprocess.py)
                if self.sesame == False:
                    # match legacy settings
                    preprocess_noob(self, linear_dye_correction = not self.correct_dye_bias) #linear_dye_correction = True)
                    if (self.methylated._MethylationDataset__dye_bias_corrected is False or
                        self.unmethylated._MethylationDataset__dye_bias_corrected is False): # process failed, so fallback is linear
                        preprocess_noob_sesame(self)
                else:
                    preprocess_noob_sesame(self)

                # nonlinear_dye_correction is done below, but if sesame if false, revert to previous linear dye method here.
                methylated = self.methylated.data_frame[['noob']].astype('float32').round(0)
                unmethylated = self.unmethylated.data_frame[['noob']].astype('float32').round(0)
            else:
                # renaming will make dye_bias work later; or I track here and pass in kwargs to dye bias for column names
                methylated = self.methylated.data_frame[['mean_value']].astype('float32').round(0).rename(columns={'mean_value':'noob'})
                unmethylated = self.unmethylated.data_frame[['mean_value']].astype('float32').round(0).rename(columns={'mean_value':'noob'})
                LOGGER.info('SDC data_frame aleady exists.')


            self.__data_frame = methylated.join(
                unmethylated,
                lsuffix='_meth',
                rsuffix='_unmeth',
            )

            if self.pval == True:
                self.__data_frame = self.__data_frame.merge(pval_probes_df, how='inner', left_index=True, right_index=True)

            if self.quality_mask == True and isinstance(quality_mask_df,pd.DataFrame):
                self.__data_frame = self.__data_frame.merge(quality_mask_df, how='inner', left_index=True, right_index=True)

            if self.correct_dye_bias == True:
                nonlinear_dye_bias_correction(self, debug=self.debug)
                if self.quality_mask == True and 'quality_mask' in self.__data_frame.columns:
                    self.__data_frame.loc[self.__data_frame['quality_mask'].isna(), 'noob_meth'] = np.nan
                    self.__data_frame.loc[self.__data_frame['quality_mask'].isna(), 'noob_unmeth'] = np.nan

            if self.retain_uncorrected_probe_intensities == True:
                self.__data_frame['meth'] = uncorrected_meth
                self.__data_frame['unmeth'] = uncorrected_unmeth

            # reduce to float32 during processing. final output may be 16,32,64 in _postprocess() + export()
            self.__data_frame = self.__data_frame.astype('float32')
            if self.poobah_decimals != 3 and 'poobah_pval' in self.__data_frame.columns:
                other_columns = list(self.__data_frame.columns)
                other_columns.remove('poobah_pval')
                other_columns = {column:3 for column in other_columns}
                self.__data_frame = self.__data_frame.round(other_columns)
                self.__data_frame = self.__data_frame.round({'poobah_pval': self.poobah_decimals})
            else:
                self.__data_frame = self.__data_frame.round(3)

            # here, separate the mouse from normal probes and store mouse experimental probes separately.
            # normal_probes_mask = (self.manifest.data_frame.index.str.startswith('cg', na=False)) | (self.manifest.data_frame.index.str.startswith('ch', na=False))
            #v2_mouse_probes_mask = (self.manifest.data_frame.index.str.startswith('mu', na=False)) | (self.manifest.data_frame.index.str.startswith('rp', na=False))
            if 'Probe_Type' in self.manifest.data_frame.columns:
                mouse_probes_mask = ( (self.manifest.data_frame['Probe_Type'] == 'mu') | (self.manifest.data_frame['Probe_Type'] == 'rp') | self.manifest.data_frame.index.str.startswith('uk', na=False) )
                mouse_probes = self.manifest.data_frame[mouse_probes_mask]
                mouse_probe_count = mouse_probes.shape[0]
            else:
                mouse_probes = pd.DataFrame()
                mouse_probe_count = 0
            self.mouse_data_frame = self.__data_frame[self.__data_frame.index.isin(mouse_probes.index)]
            if mouse_probe_count > 0:
                LOGGER.debug(f"{mouse_probe_count} mouse probes ->> {self.mouse_data_frame.shape[0]} in idat")
                # add Probe_Type column to mouse_data_frame, so it appears in the output. match manifest [IlmnID] to df.index
                # NOTE: other manifests have no 'Probe_Type' column, so avoiding this step with them.
                probe_types = self.manifest.data_frame[['Probe_Type']]
                self.mouse_data_frame = self.mouse_data_frame.join(probe_types, how='inner')
                # now remove these from normal list. confirmed they appear in the processed.csv if this line is not here.
                self.__data_frame = self.__data_frame[~self.__data_frame.index.isin(mouse_probes.index)]

            # finally, sort probes -- note: uncommenting this step breaks beta/m_value calcs in testing. Some downstream function depends on the probe_order staying same.
            #self.__data_frame.sort_index(inplace=True)

        return self.__data_frame

    def process_m_value(self, input_dataframe):
        """Calculate M value from methylation data"""
        return self._postprocess(input_dataframe, calculate_m_value, 'm_value')

    def process_beta_value(self, input_dataframe):
        """Calculate Beta value from methylation data"""
        if self.sesame == False:
            offset=100 # minfi code suggest offset of 100, but empirically, seems like the unit tests match 0 instead.
        else:
            offset=0 # make_pipeline uses sesame=None
        return self._postprocess(input_dataframe, calculate_beta_value, 'beta_value', offset)

    def process_copy_number(self, input_dataframe):
        """Calculate copy number value from methylation data"""
        return self._postprocess(input_dataframe, calculate_copy_number, 'cm_value')

    def process_all(self):
        """Runs all pre and post-processing calculations for the dataset."""
        data_frame = self.preprocess() # applies BG_correction and NOOB to .methylated, .unmethylated
        # also creates a self.mouse_data_frame for mouse specific probes with 'noob_meth' and 'noob_unmeth' columns here.
        data_frame = self.process_beta_value(data_frame)
        data_frame = self.process_m_value(data_frame)
        self.__data_frame = data_frame

        if self.manifest.array_type == ArrayType.ILLUMINA_MOUSE:
            self.mouse_data_frame = self.process_beta_value(self.mouse_data_frame)
            self.mouse_data_frame = self.process_m_value(self.mouse_data_frame)
            self.mouse_data_frame = self.process_copy_number(self.mouse_data_frame)

        return data_frame

    def export(self, output_path):
        ensure_directory_exists(output_path)
        # ensure smallest possible csv files
        self.__data_frame = self.__data_frame.round({'noob_meth':0, 'noob_unmeth':0, 'm_value':3, 'beta_value':3,
            'meth':0, 'unmeth':0, 'poobah_pval':self.poobah_decimals})
        # noob columns contain NANs now because of sesame (v1.4.0)
        #try:
        #    self.__data_frame['noob_meth'] = self.__data_frame['noob_meth'].astype(int, copy=False)
        #    self.__data_frame['noob_unmeth'] = self.__data_frame['noob_unmeth'].astype(int, copy=False)
        #except ValueError as e:
        if 'quality_mask' in self.__data_frame.columns:
            num_missing = self.__data_frame[ ~self.__data_frame['quality_mask'].isna() ]['noob_unmeth'].isna().sum() + self.__data_frame[ ~self.__data_frame['quality_mask'].isna() ]['noob_meth'].isna().sum()
        else:
            num_missing = self.__data_frame['noob_unmeth'].isna().sum() + self.__data_frame['noob_meth'].isna().sum()
        if num_missing > 0:
            self.noob_processing_missing_probe_errors.append((output_path, num_missing))
        # these are the raw, uncorrected values, replaced by sesame quality_mask as NANs
        if 'meth' in self.__data_frame.columns and 'unmeth' in self.__data_frame.columns:
            try:
                self.__data_frame['meth'] = self.__data_frame['meth'] # .astype('float16', copy=False) # --- float16 was changing these values, so not doing this step.
                self.__data_frame['unmeth'] = self.__data_frame['unmeth'] # .astype('float16', copy=False)
            except ValueError as e:
                if 'quality_mask' in self.__data_frame.columns:
                    num_missing = self.__data_frame[ ~self.__data_frame['quality_mask'].isna() ]['unmeth'].isna().sum() + self.__data_frame[ ~self.__data_frame['quality_mask'].isna() ]['meth'].isna().sum()
                else:
                    num_missing = self.__data_frame['meth'].isna().sum() + self.__data_frame['unmeth'].isna().sum()
                self.raw_processing_missing_probe_errors.append((output_path, num_missing))
        self.__data_frame.to_csv(output_path)

    def _postprocess(self, input_dataframe, postprocess_func, header, offset=None):
        if offset is not None:
            input_dataframe[header] = postprocess_func(
                input_dataframe['noob_meth'].values,
                input_dataframe['noob_unmeth'].values,
                offset=offset
            )
        else:
            input_dataframe[header] = postprocess_func(
                input_dataframe['noob_meth'].values,
                input_dataframe['noob_unmeth'].values,
            )

        if self.data_type != 'float64':
            #np.seterr(over='raise', divide='raise')
            try:
                LOGGER.debug('Converting %s to %s: %s', header, self.data_type, self.raw_dataset.sample)
                input_dataframe[header] = input_dataframe[header].astype(self.data_type)
            except Exception as e:
                LOGGER.warning(f'._postprocess: {e}')
                LOGGER.info('%s failed for %s, using float64 instead: %s', self.data_type, header, self.raw_dataset.sample)
                input_dataframe[header] = input_dataframe[header].astype('float64')

        return input_dataframe


def make_pipeline(data_dir='.', steps=None, exports=None, estimator='beta', **kwargs):
    """Specify a list of processing steps for run_pipeline, then instantiate and run that pipeline.

    steps:
        list of processing steps
        ['all', 'infer_channel_switch', 'poobah', 'quality_mask', 'noob', 'dye_bias']
    exports:
        list of files to be saved; anything not specified is not saved; ['all'] saves everything.
        ['all', 'csv', 'poobah', 'meth', 'unmeth', 'noob_meth', 'noob_unmeth', 'sample_sheet_meta_data', 'mouse', 'control']
    estimator:
        which final format?
        [beta | m_value | copy_number | None (returns containers instead)]

    This feeds a Class that runs the run_pipeline function of transforms with a final estimator.
    It replaces all of the kwargs that are in run_pipeline() and adds a few more options:

[steps] -- you can set all of these with ['all'] or any combination of these in a list of steps:
    Also note that adding "sesame=True" to kwargs will enable: infer_channel_switch, poobah, quality_mask, noob, dye_bias
    'infer_channel_switch'
    'poobah'
    'quality_mask'
    'noob'
    'dye_bias'

[exports]
    export=False,
    make_sample_sheet=False,
    export_poobah=False,
    save_uncorrected=False,
    save_control=False,
    meta_data_frame=True,

[final estimator] -- default: return list of sample data containers.
    betas=False,
    m_value=False,
    -copy_number-
    You may override that by specifying `estimator`= ('betas' or 'm_value').

[how it works]
    make_pipeline calls run_pipeline(), which has a **kwargs final keyword that maps many additional esoteric settings that you can define here.

    These are used for more granular unit testing on methylsuite, but could allow you to change how data is processed
    in very fine-tuned ways.

The rest of these are additional optional kwargs you can include:

[inputs] -- omitting these kwargs will assume the defaults, as shown below
    data_dir,
    array_type=None,
    manifest_filepath=None,
    sample_sheet_filepath=None,
    sample_name=None,

[processing] -- omitting these kwargs will assume the defaults, as shown below
    batch_size=None,  --- if you have low RAM memory or >500 samples, you might need to process the batch in chunks.
    bit='float32', --- float16 or float64 also supported for higher/lower memory/disk usage
    low_memory=True, --- If True, processing deletes intermediate objects. But you can save them in the SampleDataContainer by setting this to False.
    poobah_decimals=3 --- in csv file output
    poobah_sig=0.05

[logging] -- how much information do you want on the screen? Default is minimal information.
    verbose=False (True for more)
    debug=False (True for a LOT more info)

     """
    allowed_steps = ['all', 'infer_channel_switch', 'poobah', 'quality_mask', 'noob', 'dye_bias']
    allowed_exports = ['all', 'csv', 'poobah', 'meth', 'unmeth', 'noob_meth', 'noob_unmeth', 'sample_sheet_meta_data', 'mouse', 'control']
    allowed_estimators = ['betas', 'beta', 'm_value', 'copy_number', None]
    if steps is None:
        steps = []
    if exports is None:
        exports = []
    if not isinstance(steps,(tuple, list)) and not set(steps).issubset(set(allowed_steps)):
        raise ValueError(f"steps, the first argument, must be a list or tuple of names of allowed processing steps: {allowed_steps} or 'all'; you said {steps}")
    if not isinstance(exports,(tuple, list)) and not set(exports).issubset(set(allowed_exports)):
        raise ValueError(f"[exports] must be a list or tuple of names of allowed processing steps: {allowed_exports}, or 'all'; you said {exports}")
    if estimator not in set(allowed_estimators):
        raise ValueError(f"Your chosen final estimator must be one of these: {allowed_estimators}; you said {estimator}")
    if estimator == 'copy_number':
        raise ValueError("copy_number is not yet suppported. (You can get it in the code, but not with make_pipelines)")
    if estimator in ('betas','beta'):
        kwargs['betas'] = True
    if estimator == 'm_value':
        kwargs['m_value'] = True
    return run_pipeline(data_dir, pipeline_steps=steps, pipeline_exports=exports, **kwargs)
