# Lib
import logging
import numpy as np
import pandas as pd
from ..utils.progress_bar import * # checks environment and imports tqdm appropriately.
from collections import Counter
from pathlib import Path
import pickle
import sys
# App
from ..files import Manifest, get_sample_sheet, create_sample_sheet
from ..models import (
    Channel,
    #MethylationDataset,
    SigSet,
    ArrayType,
    #get_raw_datasets,
    get_array_type,
    parse_sample_sheet_into_idat_datasets,
)
from .postprocess import (
    calculate_beta_value,
    calculate_m_value,
    calculate_copy_number,
    consolidate_values_for_sheet,
    one_sample_control_snp,
    consolidate_mouse_probes,
    merge_batches,
)
from ..utils import ensure_directory_exists, is_file_like
from .preprocess import preprocess_noob, _apply_sesame_quality_mask
from .p_value_probe_detection import _pval_sesame_preprocess, _pval_neg_ecdf
from .infer_channel_switch import infer_type_I_probes
from .dye_bias import nonlinear_dye_bias_correction
from .multi_array_idat_batches import check_array_folders


__all__ = ['SampleDataContainer', 'run_pipeline', 'consolidate_values_for_sheet', 'make_pipeline']

LOGGER = logging.getLogger(__name__)


def run_pipeline(data_dir, array_type=None, export=False, manifest_filepath=None,
                 sample_sheet_filepath=None, sample_name=None,
                 betas=False, m_value=False, make_sample_sheet=False, batch_size=None,
                 save_uncorrected=False, save_control=True, meta_data_frame=True,
                 bit='float32', poobah=False, export_poobah=False,
                 poobah_decimals=3, poobah_sig=0.05, low_memory=True,
                 sesame=True, quality_mask=None, pneg_ecdf=False, file_format='pickle', **kwargs):
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
        file_format [default: pickle; optional: parquet]
            Matrix style files are faster to load and process than CSVs, and python supports two
            types of binary formats: pickle and parquet. Parquet is readable by other languages,
            so it is an option starting v1.7.0.
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
    #local_vars = list(locals().items())
    #print([(key,val) for key,val in local_vars])
    # support for the make_pipeline wrapper function here; a more structured way to pass in args like sklearn.
    # unexposed flags all start with 'do_': (None will retain default settings)
    do_infer_channel_switch = None # defaults to sesame(True)
    do_noob = None # defaults to True
    do_nonlinear_dye_bias = True # defaults to sesame(True), but can be False (linear) or None (omit step)
    do_save_noob = None
    do_mouse = True
    hidden_kwargs = ['pipeline_steps', 'pipeline_exports', 'debug']
    if kwargs != {}:
        for kwarg in kwargs:
            if kwarg not in hidden_kwargs:
                if sys.stdin.isatty() is False:
                    raise SystemExit(f"One of your parameters ({kwarg}) was not recognized. Did you misspell it?")
                else:
                    raise KeyError(f"One of your parameters ({kwarg}) was not recognized. Did you misspell it?")
    if sesame == True:
        poobah = True # if sesame is True and poobah is False, it hangs forever.
    if sesame == False and 'pipeline_steps' not in kwargs:
        do_nonlinear_dye_bias = False # FORCE minfi to do linear

    if kwargs != {} and 'pipeline_steps' in kwargs:
        pipeline_steps = kwargs.get('pipeline_steps')
        if 'all' in pipeline_steps:
            do_infer_channel_switch = True
            poobah = True
            quality_mask = True
            do_noob = True
            do_nonlinear_dye_bias = True
            sesame = None # prevent this from overriding elsewhere
        else:
            do_infer_channel_switch = True if 'infer_channel_switch' in pipeline_steps else False
            poobah = True if 'poobah' in pipeline_steps else False
            quality_mask = True if 'quality_mask' in pipeline_steps else False
            do_noob = True if 'noob' in pipeline_steps else False
            if 'dye_bias' in pipeline_steps:
                do_nonlinear_dye_bias = True
            elif 'linear_dye_bias' in pipeline_steps:
                do_nonlinear_dye_bias = False
            else:
                do_nonlinear_dye_bias = None # omit step
            sesame = None if sesame == True else sesame
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
    if file_format == 'parquet':
        try:
            pd.DataFrame().to_parquet()
        except AttributeError():
            LOGGER.error("parquet is not installed in your environment; reverting to pickle format")
            file_format = 'pickle'
    suffix = 'parquet' if file_format == 'parquet' else 'pkl'

    LOGGER.info('Running pipeline in: %s', data_dir)
    if bit not in ('float64','float32','float16'):
        raise ValueError("Input 'bit' must be one of ('float64','float32','float16') or ommitted.")
    if sample_name:
        LOGGER.info('Sample names: {0}'.format(sample_name))

    if make_sample_sheet:
        create_sample_sheet(data_dir)
    try:
        sample_sheet = get_sample_sheet(data_dir, filepath=sample_sheet_filepath)
    except Exception as e:
        # e will be 'Too many sample sheets in this directory.'
        instructions = check_array_folders(data_dir, verbose=True) # prints instructions for GEO multi-array data packages.
        if instructions != []:
            instructions = '\n'.join(instructions)
            print(f"This folder contains idats for multiple types of arrays. Run each array separately:\n{instructions}")
            sys.exit(0)
        raise Exception(e)

    samples = sample_sheet.get_samples()
    if sample_sheet.renamed_fields != {}:
        show_fields = []
        for k,v in sample_sheet.renamed_fields.items():
            if v != k:
                show_fields.append(f"{k} --> {v}")
            else:
                show_fields.append(f"{k}")
        LOGGER.info(f"Found {len(show_fields)} additional fields in sample_sheet:\n{' | '.join(show_fields)}")

    if sample_name is not None:
        if not isinstance(sample_name,(list,tuple)):
            raise SystemExit(f"sample_name must be a list of sample_names")
        matched_samples = [sample.name for sample in samples if sample.name in sample_name]
        if set(matched_samples) != set(sample_name):
            possible_sample_names = [sample.name for sample in samples]
            unmatched_samples = [_sample for _sample in sample_name if _sample not in possible_sample_names]
            raise SystemExit(f"Your sample_name filter does not match the samplesheet; these samples were not found: {unmatched_samples}")

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
        idat_datasets = parse_sample_sheet_into_idat_datasets(sample_sheet, sample_name=batch, from_s3=None, meta_only=False, bit=bit) # replaces get_raw_datasets
        # idat_datasets are a list; each item is a dict of {'green_idat': ..., 'red_idat':..., 'array_type', 'sample'} to feed into SigSet
        #--- pre v1.5 --- raw_datasets = get_raw_datasets(sample_sheet, sample_name=batch)
        if array_type is None: # use must provide either the array_type or manifest_filepath.
            array_type = get_array_type(idat_datasets)
        manifest = Manifest(array_type, manifest_filepath) # this allows each batch to be a different array type; but not implemented yet. common with older GEO sets.

        batch_data_containers = []
        export_paths = set() # inform CLI user where to look
        for idat_dataset_pair in tqdm(idat_datasets, total=len(idat_datasets), desc="Processing samples"):
            data_container = SampleDataContainer(
                idat_dataset_pair=idat_dataset_pair,
                manifest=manifest,
                retain_uncorrected_probe_intensities=save_uncorrected,
                bit=bit,
                switch_probes=(do_infer_channel_switch or sesame), # this applies all sesame-specific options
                quality_mask= (quality_mask or sesame or False), # this applies all sesame-specific options (beta / noob offsets too)
                do_noob=(do_noob if do_noob != None else True), # None becomes True, but make_pipeline can override with False
                pval=poobah, #defaults to False as of v1.4.0
                poobah_decimals=poobah_decimals,
                poobah_sig=poobah_sig,
                do_nonlinear_dye_bias=do_nonlinear_dye_bias, # start of run_pipeline sets this to True, False, or None
                debug=kwargs.get('debug',False),
                sesame=sesame,
                pneg_ecdf=pneg_ecdf,
                file_format=file_format,
            )
            data_container.process_all()

            if export: # as CSV or parquet
                suffix = 'parquet' if file_format == 'parquet' else 'csv'
                output_path = data_container.sample.get_export_filepath(extension=suffix)
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
                # use data_frame values instead of these class objects, because they're not in sesame SigSets.
                del data_container.man
                del data_container.snp_man
                del data_container.ctl_man
                del data_container.green_idat
                del data_container.red_idat
                del data_container.data_channel
                del data_container.methylated
                del data_container.unmethylated
                del data_container.oobG
                del data_container.oobR
                del data_container.ibG
                del data_container.ibR
            batch_data_containers.append(data_container)

            #if str(data_container.sample) == '200069280091_R01C01':
            #    print(f"200069280091_R01C01 -- cg00035864 -- meth -- {data_container._SampleDataContainer__data_frame['meth']['cg00035864']}")
            #    print(f"200069280091_R01C01 -- cg00035864 -- unmeth -- {data_container._SampleDataContainer__data_frame['unmeth']['cg00035864']}")

        if kwargs.get('debug'): LOGGER.info('[finished SampleDataContainer processing]')

        def _prepare_save_out_file(df, file_stem, uint16=False):
            out_name = f"{file_stem}_{batch_num}" if batch_size else file_stem
            if uint16 and file_format != 'parquet':
                df = df.astype('float32') if df.isna().sum().sum() > 0 else df.astype('uint16')
            else:
                df = df.astype('float32')
            if df.shape[1] > df.shape[0]:
                df = df.transpose() # put probes as columns for faster loading.
            # sort sample names
            df = df.sort_index().reindex(sorted(df.columns), axis=1)
            if file_format == 'parquet':
                # put probes in rows; format is optimized for same-type storage so it won't really matter
                df.to_parquet(Path(data_dir,f"{out_name}.parquet"))
            else:
                df.to_pickle(Path(data_dir, f"{out_name}.pkl"))
            LOGGER.info(f"saved {out_name}")

        if betas:
            df = consolidate_values_for_sheet(batch_data_containers, postprocess_func_colname='beta_value', bit=bit, poobah=poobah, exclude_rs=True)
            _prepare_save_out_file(df, 'beta_values')
        if m_value:
            df = consolidate_values_for_sheet(batch_data_containers, postprocess_func_colname='m_value', bit=bit, poobah=poobah, exclude_rs=True)
            _prepare_save_out_file(df, 'm_values')
        if (do_save_noob is not False) or betas or m_value:
            df = consolidate_values_for_sheet(batch_data_containers, postprocess_func_colname='noob_meth', bit=bit, poobah=poobah, exclude_rs=True)
            _prepare_save_out_file(df, 'noob_meth_values', uint16=True)
            df = consolidate_values_for_sheet(batch_data_containers, postprocess_func_colname='noob_unmeth', bit=bit, poobah=poobah, exclude_rs=True)
            _prepare_save_out_file(df, 'noob_unmeth_values', uint16=True)
        if save_uncorrected:
            df = consolidate_values_for_sheet(batch_data_containers, postprocess_func_colname='meth', bit=bit, poobah=False, exclude_rs=True)
            _prepare_save_out_file(df, 'meth_values', uint16=True)
            df = consolidate_values_for_sheet(batch_data_containers, postprocess_func_colname='unmeth', bit=bit, poobah=False, exclude_rs=True)
            _prepare_save_out_file(df, 'unmeth_values', uint16=True)

        if manifest.array_type == ArrayType.ILLUMINA_MOUSE and do_mouse:
            # save mouse specific probes
            if not batch_size:
                mouse_probe_filename = f'mouse_probes.{suffix}'
            else:
                mouse_probe_filename = f'mouse_probes_{batch_num}.{suffix}'
            consolidate_mouse_probes(batch_data_containers, Path(data_dir, mouse_probe_filename), file_format)
            LOGGER.info(f"saved {mouse_probe_filename}")

        if export:
            export_path_parents = list(set([str(Path(e).parent) for e in export_paths]))
            LOGGER.info(f"[!] Exported results ({file_format}) to: {export_path_parents}")

        if export_poobah:
            if all(['poobah_pval' in e._SampleDataContainer__data_frame.columns for e in batch_data_containers]):
                # this option will save pvalues for all samples, with sample_ids in the column headings and probe names in index.
                # this sets poobah to false in kwargs, otherwise some pvalues would be NaN I think.
                df = consolidate_values_for_sheet(batch_data_containers, postprocess_func_colname='poobah_pval', bit=bit, poobah=False, poobah_sig=poobah_sig, exclude_rs=True)
                _prepare_save_out_file(df, 'poobah_values')

            if all(['pNegECDF_pval' in e._SampleDataContainer__data_frame.columns for e in batch_data_containers]):
                # this option will save negative control based pvalues for all samples, with
                # sample_ids in the column headings and probe names in index.
                df = consolidate_values_for_sheet(batch_data_containers, postprocess_func_colname='pNegECDF_pval', bit=bit, poobah=False, poobah_sig=poobah_sig, exclude_rs=True)
                _prepare_save_out_file(df, 'pNegECDF_values')

        # v1.3.0 fixing mem problems: pickling each batch_data_containers object then reloading it later.

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
        meta_frame = sample_sheet.build_meta_data(samples)
        if file_format == 'parquet':
            meta_frame_filename = f'sample_sheet_meta_data.parquet'
            meta_frame.to_parquet(Path(data_dir, meta_frame_filename))
        else:
            meta_frame_filename = f'sample_sheet_meta_data.pkl'
            meta_frame.to_pickle(Path(data_dir, meta_frame_filename))
        LOGGER.info(f"saved {meta_frame_filename}")

    # FIXED in v1.3.0
    if save_control:
        if file_format == 'parquet':
            control_filename = f'control_probes.parquet'
            control = pd.concat(control_snps) # creates multiindex
            (control.reset_index()
                .rename(columns={'level_0': 'Sentrix_ID', 'level_1': 'IlmnID'})
                .astype({'IlmnID':str})
                .to_parquet('control_probes.parquet')
            )
        else:
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
        test_parts = list([str(temp_file) for temp_file in Path(data_dir).rglob(f'{file_type}*.{suffix}')])
        num_batches = len(test_parts)
        # ensures that only the file_types that appear to be selected get merged.
        #print(f"DEBUG num_batches {num_batches}, batch_size {batch_size}, file_type {file_type}")
        if batch_size and num_batches >= 1: #--- if the batch size was larger than the number of total samples, this will still drop the _1
            merge_batches(num_batches, data_dir, file_type, file_format)

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
        return consolidate_values_for_sheet(data_containers, postprocess_func_colname='beta_value', poobah=poobah, exclude_rs=True)
    elif m_value:
        return consolidate_values_for_sheet(data_containers, postprocess_func_colname='m_value', poobah=poobah, exclude_rs=True)
    else:
        return data_containers


class SampleDataContainer(SigSet):
    """Wrapper that provides easy access to red+green idat datasets, the sample, manifest, and processing params.

    Arguments:
        raw_dataset {RawDataset} -- A sample's RawDataset for a single well on the processed array.
        manifest {Manifest} -- The Manifest for the correlated RawDataset's array type.
        bit (default: float32) -- option to store data as float16 or float32 to save space.
        pval (default: False) -- whether to apply p-value-detection algorithm to remove
            unreliable probes (based on signal/noise ratio of fluoresence)
            uses the sesame method (pOOBah) based on out of band background levels

    Jan 2020: added .snp_(un)methylated property. used in postprocess.consolidate_crontrol_snp()
    Mar 2020: added p-value detection option
    Mar 2020: added mouse probe post-processing separation
    June 2020: major refactor to use SigSet, like sesame. Removed raw_dataset and methylationDataset.
    - SigSet is now a Super-class of SampleDataContainer.
    """

    __data_frame = None
    __quality_mask_excluded_probes = None
    noob_processing_missing_probe_errors = []
    raw_processing_missing_probe_errors = []

    def __init__(self, idat_dataset_pair, manifest=None, retain_uncorrected_probe_intensities=False,
                 bit='float32', pval=False, poobah_decimals=3, poobah_sig=0.05, do_noob=True,
                 quality_mask=True, switch_probes=True, do_nonlinear_dye_bias=True, debug=False, sesame=True,
                 pneg_ecdf=False, file_format='csv'):
        self.debug = debug
        self.do_noob = do_noob
        self.pval = pval
        self.poobah_decimals = poobah_decimals
        self.poobah_sig = poobah_sig
        self.quality_mask = quality_mask # if True, filters sesame's standard sketchy probes out of 450k, EPIC, EPIC+ arrays.
        self.switch_probes = switch_probes
        self.do_nonlinear_dye_bias = do_nonlinear_dye_bias
        self.green_idat = idat_dataset_pair['green_idat']
        self.red_idat = idat_dataset_pair['red_idat']
        self.sample = idat_dataset_pair['sample']
        self.retain_uncorrected_probe_intensities=retain_uncorrected_probe_intensities
        self.sesame = sesame # defines offsets in functions
        # pneg_ecdf defines if negative control based pvalue is calculated - will use poobah_decimals for rounding
        self.pneg_ecdf = pneg_ecdf
        self.data_type = 'float32' if bit == None else bit # options: (float64, float32, or float16)
        self.file_format = file_format
        if debug:
            print(f'DEBUG SDC: sesame {self.sesame} switch {self.switch_probes} noob {self.do_noob} poobah {self.pval} mask {self.quality_mask}, dye {self.do_nonlinear_dye_bias}')

        self.manifest = manifest # used by inter_channel_switch only.
        if self.switch_probes:
            # apply inter_channel_switch here; uses raw_dataset and manifest only; then updates self.raw_dataset
            # these are read from idats directly, not SigSet, so need to be modified at source.
            infer_type_I_probes(self, debug=self.debug)

        super().__init__(self.sample, self.green_idat, self.red_idat, self.manifest, self.debug)
        # SigSet defines all probe-subsets, then SampleDataContainer adds them with super(); no need to re-define below.
        # mouse probes are processed within the normals meth/unmeth sets, then split at end of preprocessing step.
        del self.manifest
        del manifest

        if self.data_type not in ('float64','float32','float16'):
            raise ValueError(f"invalid data_type: {self.data_type} should be one of ('float64','float32','float16')")

        if self.debug:
            print(f"DEBUG SDC params:")
            lists = ['red_switched','green_switched']
            exclude = ['data_channel', 'man', 'snp_man', 'ctl_man', 'address_code', 'ctrl_green', 'ctrl_red', 'II',
            'IG', 'IR', 'oobG', 'oobR', 'methylated', 'unmethylated', 'snp_methylated', 'snp_unmethylated', 'ibG', 'ibR',
            'mouse_probes_mask',            ]
            for key,value in self.__dict__.items():
                if key in exclude:
                    try:
                        print(f"-- {key}: {value.shape}")
                    except:
                        print(f"-- skipping {key}")
                elif key in lists:
                    print(f"-- {key}: {len(value)}")
                else:
                    print(f"-- {key}: {value}")
            self.check_for_probe_loss()

    def process_all(self):
        """Runs all pre and post-processing calculations for the dataset.
        Combines the SigSet methylated and unmethylated parts of SampleDataContainer, and modifies them,
        whilst creating self.__data_frame with noob/dye processed data.

    Order:
        - poobah
        - quality_mask
        - noob (background correction)
        - build data_frame
        - nonlinear dye-bias correction
        - reduce memory/bit-depth of data
        - copy over uncorrected values
        - split out mouse probes
        """
        # self.preprocess -- applies BG_correction and NOOB to .methylated, .unmethylated
        # also creates a self.mouse_data_frame for mouse specific probes with 'noob_meth' and 'noob_unmeth' columns here.
        if self.__data_frame:
            return self.__data_frame

        pval_probes_df = _pval_sesame_preprocess(self) if self.pval == True else None
        pneg_ecdf_probes_df = _pval_neg_ecdf(self) if self.pneg_ecdf == True else None
        # output: df with one column named 'poobah_pval'
        quality_mask_df = _apply_sesame_quality_mask(self) if self.quality_mask == True else None
        # output: df with one column named 'quality_mask' | if not supported array / custom array: returns nothing.

        if self.do_noob == True:
            # apply corrections: bg subtract, then noob (in preprocess.py)
            preprocess_noob(self, pval_probes_df=pval_probes_df, quality_mask_df=quality_mask_df, nonlinear_dye_correction=self.do_nonlinear_dye_bias, debug=self.debug)
            #if self.sesame in (None,True):
                #preprocess_noob(self, pval_probes_df=pval_probes_df, quality_mask_df=quality_mask_df, nonlinear_dye_correction=self.do_nonlinear_dye_bias, debug=self.debug)
                #if container.__dye_bias_corrected is False: # process failed, so fallback is linear-dye
                #    print(f'ERROR preprocess_noob sesame-dye: linear_dye_correction={self.do_nonlinear_dye_bias}')
                #    preprocess_noob(self, linear_dye_correction = True)
            #if self.sesame is False:
                # match minfi legacy settings
                #preprocess_noob(self, pval_probes_df=pval_probes_df, quality_mask_df=quality_mask_df, nonlinear_dye_correction=self.do_nonlinear_dye_bias, debug=self.debug)

            if self._SigSet__preprocessed is False:
                raise ValueError("preprocessing did not run")

            # nonlinear_dye_correction is done below, but if sesame if false, revert to previous linear dye method here.
            self.methylated = self.methylated.rename(columns={'noob_Meth':'noob'}).drop(columns=['used','Unmeth', 'noob_Unmeth'])
            self.unmethylated = self.unmethylated.rename(columns={'noob_Unmeth':'noob'}).drop(columns=['used','Meth', 'noob_Meth'])
        else:
            # renaming will make dye_bias work later; or I track here and pass in kwargs to dye bias for column names
            self.methylated = self.methylated[['Meth']].astype('float32').round(0) #.rename(columns={'mean_value':'noob'})
            self.unmethylated = self.unmethylated[['Unmeth']].astype('float32').round(0) #.rename(columns={'mean_value':'noob'})
            if self.debug: LOGGER.info('SDC data_frame already exists.') #--- happens with make_pipeline('.',steps=[])

        if set(self.unmethylated.index) - set(self.methylated.index) != set():
            LOGGER.warning(f"Dropping mismatched probes: {set(self.unmethylated.index) - set(self.methylated.index)}")
        if set(self.methylated.index) - set(self.unmethylated.index) != set():
            LOGGER.warning(f"Dropping mismatched probes: {set(self.methylated.index) - set(self.unmethylated.index)}")

        try:
            # index: IlmnID | has A | B | Unmeth | Meth | noob_meth | noob_unmeth -- no control or snp probes included
            self.__data_frame = self.methylated.join(
                self.unmethylated.drop(columns=['AddressA_ID','AddressB_ID']),
                lsuffix='_meth', rsuffix='_unmeth',
                how='inner').drop(columns=['AddressA_ID','AddressB_ID'])
                # 'inner' join is necessary to avoid dye-bias getting duplicate probes if mismatched data.
        except KeyError: # for steps=[]
            self.__data_frame = self.methylated.join(
                self.unmethylated,
                lsuffix='_meth', rsuffix='_unmeth',
                how='inner')
            # noob did not run, but copying data into new column so steps won't break
            self.__data_frame['noob_meth'] = self.__data_frame['Meth']
            self.__data_frame['noob_unmeth'] = self.__data_frame['Unmeth']

        if self.pval == True and isinstance(pval_probes_df, pd.DataFrame):
            pval_probes_df = pval_probes_df.loc[ ~pval_probes_df.index.duplicated() ]
            self.__data_frame = self.__data_frame.join(pval_probes_df, how='inner')

        if self.pneg_ecdf == True and isinstance(pneg_ecdf_probes_df, pd.DataFrame):
            pneg_ecdf_probes_df = pneg_ecdf_probes_df.loc[ ~pneg_ecdf_probes_df.index.duplicated() ]
            self.__data_frame = self.__data_frame.join(pneg_ecdf_probes_df, how='inner')

        self.check_for_probe_loss(f"preprocess_noob sesame={self.sesame} --> {self.methylated.shape} {self.unmethylated.shape}")

        if self.quality_mask == True and isinstance(quality_mask_df, pd.DataFrame):
            self.__data_frame = self.__data_frame.join(quality_mask_df, how='inner')

        if self.do_nonlinear_dye_bias == True:
            nonlinear_dye_bias_correction(self, debug=self.debug)
            # this step ensures that failed probes are not included in the NOOB calculations.
            # but they MUST be included in CSV exports, so I move the failed probes to another df for storage until pipeline.export() needs them.
            if self.quality_mask == True and 'quality_mask' in self.__data_frame.columns:
                self.__quality_mask_excluded_probes = self.__data_frame.loc[self.__data_frame['quality_mask'].isna(), ['noob_meth','noob_unmeth']]
                self.__data_frame.loc[self.__data_frame['quality_mask'].isna(), 'noob_meth'] = np.nan
                self.__data_frame.loc[self.__data_frame['quality_mask'].isna(), 'noob_unmeth'] = np.nan

        # Downcasting to 'unsigned' uses the smallest possible integer that can hold the values, but need to retain dtypes
        if self.__data_frame['Meth'].isna().sum() == 0 and self.__data_frame['Unmeth'].isna().sum() == 0:
            self.__data_frame['Meth'] = self.__data_frame['Meth'].apply(pd.to_numeric, downcast='unsigned')
            self.__data_frame['Unmeth'] = self.__data_frame['Unmeth'].apply(pd.to_numeric, downcast='unsigned')
        for column in self.__data_frame.columns:
            if column in ['Meth','Unmeth']:
                continue
            if self.__data_frame[column].isna().sum() == 0: # and self.__data_frame[column].dtype
                self.__data_frame[column] = self.__data_frame[column].apply(pd.to_numeric, downcast='unsigned')

        if self.retain_uncorrected_probe_intensities == False:
            self.__data_frame = self.__data_frame.drop(columns=['Meth','Unmeth'])
        else:
            self.__data_frame = self.__data_frame.rename(columns={'Meth':'meth', 'Unmeth':'unmeth'})

        # reduce to float32 during processing. final output may be 16,32,64 in _postprocess() + export()
        #self.__data_frame = self.__data_frame.astype('float32') --- downcast takes care of this
        #if self.poobah_decimals != 3 and 'poobah_pval' in self.__data_frame.columns:
        #    other_columns = list(self.__data_frame.columns)
        #    other_columns.remove('poobah_pval')
        #    other_columns = {column:3 for column in other_columns}
        #    #self.__data_frame = self.__data_frame.round(other_columns)
        #    self.__data_frame = self.__data_frame.round({'poobah_pval': self.poobah_decimals})
        #else:
        #    self.__data_frame = self.__data_frame.round(3)
        self.__data_frame = self.__data_frame.round({'poobah_pval': self.poobah_decimals})

        # here, separate the mouse from normal probes and store mouse experimental probes separately.
        # normal_probes_mask = (self.manifest.data_frame.index.str.startswith('cg', na=False)) | (self.manifest.data_frame.index.str.startswith('ch', na=False))
        # v2_mouse_probes_mask = (self.manifest.data_frame.index.str.startswith('mu', na=False)) | (self.manifest.data_frame.index.str.startswith('rp', na=False))
        # v4 mouse_probes_mask pre-v1.4.6: ( (self.manifest.data_frame['Probe_Type'] == 'mu') | (self.manifest.data_frame['Probe_Type'] == 'rp') | self.manifest.data_frame.index.str.startswith('uk', na=False) )
        if self.array_type == ArrayType.ILLUMINA_MOUSE:
            mouse_probes = self.man[self.mouse_probes_mask]
            mouse_probe_count = mouse_probes.shape[0]
        else:
            mouse_probes = pd.DataFrame()
            mouse_probe_count = 0
        self.mouse_data_frame = self.__data_frame[self.__data_frame.index.isin(mouse_probes.index)]
        if mouse_probe_count > 0:
            if self.debug: LOGGER.info(f"{mouse_probe_count} mouse probes ->> {self.mouse_data_frame.shape[0]} in idat")
            # add 'design' column to mouse_data_frame, so it appears in the output. -- needed for 'Random' and 'Multi' filter
            # matches manifest [IlmnID] to df.index
            # NOTE: other manifests have no 'design' column, so avoiding this step with them.
            probe_designs = self.man[['design']]
            self.mouse_data_frame = self.mouse_data_frame.join(probe_designs, how='inner')
            # now remove these from normal list. confirmed they appear in the processed.csv if this line is not here.
            self.__data_frame = self.__data_frame[~self.__data_frame.index.isin(mouse_probes.index)]

        # finally, sort probes -- note: uncommenting this step breaks beta/m_value calcs in testing. Some downstream function depends on the probe_order staying same.
        # --- must fix all unit tests using .iloc[ before this will work | fixed.
        self.__data_frame.sort_index(inplace=True)
        ###### end preprocessing ######

        if hasattr(self, '_SampleDataContainer__quality_mask_excluded_probes') and isinstance(self._SampleDataContainer__quality_mask_excluded_probes, pd.DataFrame):
            # these probes are not used in processing, but should appear in the final CSV.
            # and if quality_mask is off, it should skip this step.
            # later: consoldate() should exclude these probes from pickles
            self.__data_frame.update({
                'noob_meth': self.__quality_mask_excluded_probes['noob_meth'],
                'noob_unmeth': self.__quality_mask_excluded_probes['noob_unmeth']
                })
        self.__data_frame = self.process_beta_value(self.__data_frame)
        self.__data_frame = self.process_m_value(self.__data_frame)

        if self.debug:
            self.check_for_probe_loss(f"816 self.check_for_probe_loss(): self.__data_frame = {self.__data_frame.shape}")

        if self.array_type == ArrayType.ILLUMINA_MOUSE:
            self.mouse_data_frame = self.process_beta_value(self.mouse_data_frame)
            self.mouse_data_frame = self.process_m_value(self.mouse_data_frame)
            self.mouse_data_frame = self.process_copy_number(self.mouse_data_frame)

        return # self.__data_frame

    def process_m_value(self, input_dataframe):
        """Calculate M value from methylation data"""
        return self._postprocess(input_dataframe, calculate_m_value, 'm_value')

    def process_beta_value(self, input_dataframe, quality_mask_probes=None):
        """Calculate Beta value from methylation data"""
        if self.sesame == False:
            offset=100 # minfi code suggest offset of 100, but empirically, seems like the unit tests match 0 instead.
        else:
            offset=0 # make_pipeline uses sesame=None

        return self._postprocess(input_dataframe, calculate_beta_value, 'beta_value', offset)

    def process_copy_number(self, input_dataframe):
        """Calculate copy number value from methylation data"""
        return self._postprocess(input_dataframe, calculate_copy_number, 'cm_value')


    def export(self, output_path):
        """Saves a CSV for each sample with all processing intermediate data"""
        ensure_directory_exists(output_path)
        # ensure smallest possible csv files
        self.__data_frame = self.__data_frame.round({'noob_meth':0, 'noob_unmeth':0, 'm_value':3, 'beta_value':3,
            'meth':0, 'unmeth':0, 'poobah_pval':self.poobah_decimals})
        this = self.__data_frame.copy(deep=True)
        if hasattr(self, '_SampleDataContainer__quality_mask_excluded_probes') and isinstance(self._SampleDataContainer__quality_mask_excluded_probes, pd.DataFrame):
            # copy over these failed probes to a dataframe for export
            this.update({
                'noob_meth': self.__quality_mask_excluded_probes['noob_meth'],
                'noob_unmeth': self.__quality_mask_excluded_probes['noob_unmeth']
                })
        if 'quality_mask' in this.columns:
            this['quality_mask'] = this['quality_mask'].fillna(1)
        # noob columns contain NANs now because of sesame (v1.4.0 to v1.4.5); v1.4.6+ CSVs contain all data, but pickles are filtered.
        #try:
        #    self.__data_frame['noob_meth'] = self.__data_frame['noob_meth'].astype(int, copy=False)
        #    self.__data_frame['noob_unmeth'] = self.__data_frame['noob_unmeth'].astype(int, copy=False)
        #except ValueError as e:
        #if 'quality_mask' in self.__data_frame.columns:
        #    num_missing = self.__data_frame[ ~self.__data_frame['quality_mask'].isna() ]['noob_unmeth'].isna().sum() + self.__data_frame[ ~self.__data_frame['quality_mask'].isna() ]['noob_meth'].isna().sum()
        #else:
        #    num_missing = self.__data_frame['noob_unmeth'].isna().sum() + self.__data_frame['noob_meth'].isna().sum()
        #if num_missing > 0:
        #    self.noob_processing_missing_probe_errors.append((output_path, num_missing))
        # these are the raw, uncorrected values, replaced by sesame quality_mask as NANs
        #if 'meth' in self.__data_frame.columns and 'unmeth' in self.__data_frame.columns:
        #    try:
        #        self.__data_frame['meth'] = self.__data_frame['meth'] # .astype('float16', copy=False) # --- float16 was changing these values, so not doing this step.
        #        self.__data_frame['unmeth'] = self.__data_frame['unmeth'] # .astype('float16', copy=False)
        #    except ValueError as e:
        #        if 'quality_mask' in self.__data_frame.columns:
        #            num_missing = self.__data_frame[ ~self.__data_frame['quality_mask'].isna() ]['unmeth'].isna().sum() + self.__data_frame[ ~self.__data_frame['quality_mask'].isna() ]['meth'].isna().sum()
        #        else:
        #            num_missing = self.__data_frame['meth'].isna().sum() + self.__data_frame['unmeth'].isna().sum()
        #        self.raw_processing_missing_probe_errors.append((output_path, num_missing))
        if self.file_format == 'parquet':
            this.to_parquet(output_path)
        else:
            this.to_csv(output_path)

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

        if self.data_type != 'float32':
            #np.seterr(over='raise', divide='raise')
            try:
                LOGGER.debug('Converting %s to %s: %s', header, self.data_type, self.sample)
                input_dataframe[header] = input_dataframe[header].astype(self.data_type)
            except Exception as e:
                LOGGER.warning(f'._postprocess: {e}')
                LOGGER.info('%s failed for %s, using float32 instead: %s', self.data_type, header, self.sample)
                input_dataframe[header] = input_dataframe[header].astype('float32')

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
    'dye_bias' -- specifying this select's sesame's nonlinear-dye-bias correction. Omitting causes NOOB to use minfi's linear-dye-correction, unless NOOB is missing.

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
        raise ValueError(f"steps, the first argument, must be a list or tuple of names of allowed processing steps: {allowed_steps} or ['all']; you said {steps}")
    if not isinstance(exports,(tuple, list)) and not set(exports).issubset(set(allowed_exports)):
        raise ValueError(f"[exports] must be a list or tuple of names of allowed processing steps: {allowed_exports}, or ['all']; you said {exports}")
    if estimator not in set(allowed_estimators):
        raise ValueError(f"Your chosen final estimator must be one of these: {allowed_estimators}; you said {estimator}")
    if estimator == 'copy_number':
        raise ValueError("copy_number is not yet suppported. (You can get it in the code, but not with make_pipelines)")
    if estimator in ('betas','beta'):
        kwargs['betas'] = True
    if estimator == 'm_value':
        kwargs['m_value'] = True
    return run_pipeline(data_dir, pipeline_steps=steps, pipeline_exports=exports, **kwargs)
