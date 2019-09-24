# Lib
import argparse
import logging
from pathlib import Path
import sys
# App
from .files import get_sample_sheet
from .models import ArrayType
from .processing import run_pipeline
from .download import (
    run_series,
    run_series_list
    )


class DefaultParser(argparse.ArgumentParser):
    def error(self, message):
        self._print_message('[Error]:\n')
        self._print_message(f'{message}\n\n')
        self.print_help()
        self.exit(status=2)


def build_parser():
    parser = DefaultParser(
        prog='methylprep',
        description="""Utility to process methylation data from Illumina IDAT files.
        There are two types of processing: "process" IDAT files or read a "sample_sheet" contents.
        Example of usage: `python -m methylprep -v process -d <path to your samplesheet.csv and idat files>`\n
        Try our demo dataset: `python -m methylprep -v process -d docs/example_data/GSE69852`""",
    )

    parser.add_argument(
        '-v', '--verbose',
        help='Enable verbose logging. Reports the path(s) where processed files are stored.',
        action='store_true',
    )

    subparsers = parser.add_subparsers(dest='command') #, required=True)
    subparsers.required = True # this is a python3.4-3.7 bug; cannot specify in the call above.

    process_parser = subparsers.add_parser('process', help='Finds idat files and calculates raw, beta, m_values for a batch of samples.')
    process_parser.set_defaults(func=cli_process)

    sample_sheet_parser = subparsers.add_parser('sample_sheet', help='Finds and validates a SampleSheet for a given directory of idat files.')
    sample_sheet_parser.set_defaults(func=cli_sample_sheet)

    download_parser = subparsers.add_parser('download', help='Downloads the specified series from GEO or ArrayExpress.')
    download_parser.set_defaults(func=cli_download)

    parsed_args, func_args = parser.parse_known_args(sys.argv[1:])
    if parsed_args.verbose:
        logging.basicConfig(level=logging.DEBUG)

    if parsed_args.command is None:
        parsed_args.command = 'process'

    parsed_args.func(func_args)
    return parser

def cli_sample_sheet(cmd_args):
    parser = DefaultParser(
        prog='methylprep sample_sheet',
        description='Process Illumina sample sheet file',
    )

    parser.add_argument(
        '-d', '--data_dir',
        required=True,
        type=Path,
        help='Base directory of the sample sheet and associated IDAT files.',
    )

    parser.add_argument(
        '-c', '--create',
        required=False,
        action='store_true',
        help='If specified, this creates a sample sheet from idats instead of parsing an existing sample sheet. The output file will be called "samplesheet.csv".',
    )

    parser.add_argument(
        '-o', '--output_file',
        required=False,
        default='samplesheet.csv',
        type=str,
        help='If creating a sample sheet, you can provide an optional output filename (CSV).'
    )

    parsed_args = parser.parse_args(cmd_args)

    if parsed_args.create == True:
        from methylprep.files import create_sample_sheet
        create_sample_sheet(parsed_args.data_dir, matrix_file=False, output_file=parsed_args.output_file)
    sample_sheet = get_sample_sheet(parsed_args.data_dir)
    for sample in sample_sheet.get_samples():
        sys.stdout.write(f'{sample}\n')


def cli_process(cmd_args):
    parser = DefaultParser(
        prog='methylprep idat',
        description='Process Illumina IDAT files',
    )

    parser.add_argument(
        '-d', '--data_dir',
        required=True,
        type=Path,
        help='Base directory of the sample sheet and associated IDAT files. If IDAT files are in nested directories, this will discover them.',
    )

    parser.add_argument(
        '-a', '--array_type',
        choices=list(ArrayType),
        required=False,
        type=ArrayType,
        help='Type of array being processed. If omitted, this will autodetect it.',
    )

    parser.add_argument(
        '-m', '--manifest',
        required=False,
        type=Path,
        help='File path of the array manifest file. If omitted, this will download the appropriate file from `s3`.',
    )

    parser.add_argument(
        '-s', '--sample_sheet',
        required=False,
        type=Path,
        help='File path of the sample sheet. If omitted, this will discover it. There must be only one CSV file in the data_dir for discovery to work.',
    )

    parser.add_argument(
        '--no_sample_sheet',
        required=False,
        action='store_true', # if -e passed, this suppresses data export (if running as part of pipeline or something)
        default=False, # if False, CLI returns nothing.
        help='If your dataset lacks a sample sheet csv file, specify --no_sample_sheet to have it create one on the fly. This will read .idat file names and ensure processing works. If there is a matrix file, it will add in sample names too.',
    )

    parser.add_argument(
        '-n', '--sample_name',
        required=False,
        nargs='*', # -- this flag support making a list of of each -n
        help='Sample(s) to process. You can pass multiple sample names with multiple -n params.',
    )

    parser.add_argument(
        '-e', '--no_export',
        required=False,
        action='store_false', # if -e passed, this suppresses data export (if running as part of pipeline or something)
        default=True, # if False, CLI returns nothing.
        help='Default is to export data to csv in same folder where IDAT file resides. Pass in --no_export to suppress this.',
    )

    parser.add_argument(
        '-b', '--betas',
        required=False,
        action='store_true',
        default=False,
        help='If passed, output returns a dataframe of beta values for samples x probes. Local file beta_values.npy is also created.',
    )

    parser.add_argument(
        '--m_value',
        required=False,
        action='store_true',
        default=False,
        help='If passed, output returns a dataframe of M-values for samples x probes. Local file m_values.npy is also created.',
    )

    parser.add_argument(
        '--batch_size',
        required=False,
        type=int,
        help='If specified, samples will be processed and saved in batches no greater than the specified batch size'
    )

    args = parser.parse_args(cmd_args)

    array_type = args.array_type
    manifest_filepath = args.manifest

    if not array_type and not manifest_filepath:
        #print("This will attempt to autodetect your methylation array_type and download the corresponding manifest file.")
        logging.info('This will attempt to autodetect your methylation array_type and download the corresponding manifest file.')
        #return

    run_pipeline(
        args.data_dir,
        array_type=args.array_type,
        export=args.no_export,
        manifest_filepath=args.manifest,
        sample_sheet_filepath=args.sample_sheet,
        sample_name=args.sample_name,
        make_sample_sheet=args.no_sample_sheet,
        betas=args.betas,
        m_value=args.m_value,
        batch_size=args.batch_size
    )

def cli_download(cmd_args):
    parser = DefaultParser(
        prog='methylprep download',
        description='Download public series'
    )

    parser.add_argument(
        '-d', '--data_dir',
        required=True,
        type=Path,
        help='Directory to download series to',
    )

    parser.add_argument(
        '-i', '--id',
        required=False,
        help='Unique ID of the series (either GEO or ArrayExpress ID)',
    )

    parser.add_argument(
        '-l', '--list',
        required=False,
        type=Path,
        help='List of series IDs (can be either GEO or ArrayExpress)',
        )

    parser.add_argument(
        '-o', '--dict_only',
        required=False,
        action='store_true',
        default=False,
        help='If passed, will only create dictionaries and not process any samples',
        )

    parser.add_argument(
        '-b', '--batch_size',
        required=False,
        type=int,
        help='Number of samples to process at a time, 100 by default'
    )

    parser.add_argument(
        '-c', '--no_clean',
        required=False,
        action="store_false",
        help='Leave processing and raw data files in folders. By default, these files are removed during processing.'
    )

    args = parser.parse_args(cmd_args)

    if args.id:
        if args.batch_size:
            run_series(args.id, args.data_dir, dict_only=args.dict_only, batch_size=args.batch_size, clean=args.no_clean)
        else:
            run_series(args.id, args.data_dir, dict_only=args.dict_only, clean=args.no_clean)
    elif args.list:
        if args.batch_size:
            run_series_list(args.list, args.data_dir, dict_only=args.dict_only, batch_size=args.batch_size)
        else:
            run_series_list(args.list, args.data_dir, dict_only=args.dict_only)

def cli_app():
    build_parser()
