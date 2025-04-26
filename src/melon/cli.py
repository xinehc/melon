import sys
import os
import glob
import argparse

from rich_argparse import ArgumentDefaultsRichHelpFormatter
from . import __version__
from .utils import logger
from .melon import GenomeProfiler

## customize formatter
ArgumentDefaultsRichHelpFormatter.styles['argparse.prog'] = 'default'
ArgumentDefaultsRichHelpFormatter.styles['argparse.default'] = 'grey50'
ArgumentDefaultsRichHelpFormatter.styles['argparse.metavar'] = 'grey50'
ArgumentDefaultsRichHelpFormatter.styles['argparse.groups'] = '#FFA62F'
ArgumentDefaultsRichHelpFormatter.styles['argparse.args'] = 'default'


def cli(argv=sys.argv):
    '''
    Entry point for command line interface.
    '''
    parser = argparse.ArgumentParser(
        description=f'Melon v{__version__}: metagenomic long-read-based taxonomic identification and quantification using marker genes',
        formatter_class=ArgumentDefaultsRichHelpFormatter,
        add_help=False)

    parser.add_argument(
        dest='files',
        nargs='+',
        metavar='file',
        help='Input fasta <*.fa|*.fasta> or fastq <*.fq|*.fastq> file, gzip optional <*.gz>.')

    required = parser.add_argument_group('required arguments')
    required.add_argument(
        '-d',
        '--db',
        metavar='DIR',
        required=True,
        help='Unzipped database folder, should contains <prot.fa>, <nucl.*.fa> and <metadata.tsv>.')

    required.add_argument(
        '-o',
        '--output',
        metavar='DIR',
        required=True,
        help='Output folder.')

    optional = parser.add_argument_group('optional arguments')
    optional.add_argument(
        '-t',
        '--threads',
        metavar='INT',
        type=int,
        default=os.cpu_count(),
        help=f'Number of threads. [{os.cpu_count()}]')

    optional.add_argument(
        '-k',
        '--db-kraken',
        metavar='DIR',
        help='Unzipped kraken2 database for pre-filtering of non-prokaryotic reads. Skip if not given.')

    optional.add_argument(
        '--skip-profile',
        action='store_true',
        help='Skip profiling, output only total genome copies.')

    optional.add_argument(
        '--skip-clean',
        action='store_true',
        help='Skip cleaning, keep all temporary <*.tmp> files.')

    additional = parser.add_argument_group('additional arguments - Filtering')
    additional.add_argument(
        '-m',
        metavar='INT',
        type=int,
        default=25,
        help='Max. number of target sequences to report (--max-target-seqs/-k in diamond).')

    additional.add_argument(
        '-e',
        metavar='FLOAT',
        type=float,
        default=1e-15,
        help='Max. expected value to report alignments (--evalue/-e in diamond).')

    additional.add_argument(
        '-i',
        metavar='FLOAT',
        type=float,
        default=0,
        help='Min. identity in percentage to report alignments (--id in diamond).')

    additional.add_argument(
        '-s',
        metavar='FLOAT',
        type=float,
        default=75,
        help='Min. subject cover to report alignments (--subject-cover in diamond).')

    additional.add_argument(
        '-n',
        metavar='INT',
        type=int,
        default=2147483647,
        help='Max. number of secondary alignments to report (-N in minimap2).')

    additional.add_argument(
        '-p',
        metavar='FLOAT',
        type=float,
        default=0.9,
        help='Min. secondary-to-primary score ratio to report secondary alignments (-p in minimap2).')

    additional.add_argument(
        '-g',
        metavar='INT',
        type=int,
        choices=range(1, 9),
        default=2,
        help='Min. number of unique marker genes required for a species to report its genome copies.')

    additional = parser.add_argument_group('additional arguments - Fitting')
    additional.add_argument(
        '-a',
        metavar='INT',
        type=int,
        default=1000,
        help='Terminal condition for EM - max. iterations.')

    additional.add_argument(
        '-c',
        metavar='FLOAT',
        type=float,
        default=1e-10,
        help='Terminal condition for EM - epsilon (precision).')

    parser.add_argument('-v', '--version', action='version', version=__version__, help=argparse.SUPPRESS)
    parser.add_argument('-h', '--help', action='help', help=argparse.SUPPRESS)

    if len(argv) == 1:
        print(f"\
             __            \n\
  __ _  ___ / /__  ___     \n\
 /  ' \\/ -_) / _ \\/ _ \\ \n\
/_/_/_/\\__/_/\\___/_//_/ ver. {__version__}\n")

    opt = parser.parse_args(argv[1:])
    run(opt)


def run(opt):
    '''
    Sanity check of options.
    '''
    ## check for output folder
    if not os.path.isdir(opt.output):
        os.makedirs(opt.output, exist_ok=True)
    else:
        logger.warning(f'Folder <{opt.output}> exists. Files will be overwritten.')

    ## check for input files
    for file in opt.files:
        if not os.path.isfile(file):
            logger.critical(f'File <{file}> does not exist.')
            sys.exit(2)

    ## check for database
    if not os.path.isdir(opt.db):
        logger.critical(f'Database folder <{opt.db}> does not exist.')
        sys.exit(2)
    else:
        files = [os.path.basename(file) for file in glob.glob(os.path.join(opt.db, '*'))]
        if (
            'metadata.tsv' not in files or 'prot.dmnd' not in files or
            len([file for file in files if 'nucl' in file and '.mmi' in file]) != 16
        ):
            logger.critical(f'Database <{opt.db}> is not complete or not indexed.')
            sys.exit(2)

    ## check for kraken2 database
    if opt.db_kraken is not None:
        if not os.path.isdir(opt.db_kraken):
            logger.critical(f'Kraken2 database folder <{opt.db_kraken}> does not exist.')
            sys.exit(2)
        else:
            files = [os.path.basename(file) for file in glob.glob(os.path.join(opt.db_kraken, '*'))]
            if 'ktaxonomy.tsv' not in files or len([file for file in files if 'database' in file]) != 7:
                logger.critical(f'Kraken2 database <{opt.db_kraken}> is not complete.')
                sys.exit(2)

    ## check for logical cores
    if opt.threads > os.cpu_count():
        logger.warning(f'Threads <{opt.threads}> exceeds available logical cores, will use <{os.cpu_count()}> instead.')
        opt.threads = os.cpu_count()

    ## run
    for index, file in enumerate(opt.files):
        if len(opt.files) > 1:
            logger.info(f'Processing file <{file}> ({index + 1}/{len(opt.files)}) ...')

        GenomeProfiler(file, opt.db, opt.output, opt.threads).run(
            db_kraken=opt.db_kraken, skip_profile=opt.skip_profile, skip_clean=opt.skip_clean,
            max_target_seqs=opt.m, evalue=opt.e, identity=opt.i, subject_cover=opt.s,
            secondary_num=opt.n, secondary_ratio=opt.p, min_markers=opt.g,
            max_iterations=opt.a, epsilon=opt.c)

        if index == len(opt.files) - 1:
            logger.info('Done.')
        else:
            logger.info('Done.\n')


if __name__ == '__main__':
    cli(sys.argv)
