import sys
import os
import glob

from os import environ, cpu_count
from argparse import ArgumentParser, SUPPRESS
from . import __version__
from .utils import logger
from .melon import GenomeProfiler


def cli(argv=sys.argv):
    '''
    Entry point for command line interface.
    '''
    parser = ArgumentParser(description='Melon: \
        long-read targeting taxonomic profiling and genome copies estimation \
        phylogenetic marker genes.', add_help=False)
    required = parser.add_argument_group('required arguments')
    optional = parser.add_argument_group('optional arguments')
    additional = parser.add_argument_group('additional arguments')

    parser.add_argument(
        'FILE',
        nargs='+',
        help='Input fasta <*.fa|*.fasta> or fastq <*.fq|*.fastq> file, gzip optional <*.gz>.')

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

    optional.add_argument(
        '-t',
        '--threads',
        metavar='INT',
        type=int,
        default=os.cpu_count(),
        help='Number of threads. [{}]'.format(os.cpu_count()))

    optional.add_argument(
        '-k',
        '--db-kraken',
        metavar='DIR',
        help='Unzipped kraken2 database for pre-filtering of non-prokaryotic reads. Skip if not given.')

    optional.add_argument(
        '--skip-profile',
        action='store_true',
        help='Skip profiling, output only estimated total genome copies.')

    optional.add_argument(
        '--skip-clean',
        action='store_true',
        help='Skip cleaning, keep all temporary <*.tmp> files.')

    additional.add_argument(
        '-m',
        metavar='INT',
        type=int,
        default=25,
        help='Max. number of target sequences to report (--max-target-seqs/-k in diamond). [25]')

    additional.add_argument(
        '-e',
        metavar='FLOAT',
        type=float,
        default=1e-15,
        help='Max. expected value to report alignments (--evalue/-e in diamond). [1e-15]')

    additional.add_argument(
        '-i',
        metavar='FLOAT',
        type=float,
        default=0,
        help='Min. identity in percentage to report alignments (--id in diamond). [0]')

    additional.add_argument(
        '-s',
        metavar='FLOAT',
        type=float,
        default=75,
        help='Min. subject cover to report alignments (--subject-cover in diamond). [75]')

    additional.add_argument(
        '-n',
        metavar='INT',
        type=int,
        default=2147483647,
        help='Max. number of secondary alignments to report (-N in minimap2). [2147483647]')

    additional.add_argument(
        '-p',
        metavar='FLOAT',
        type=float,
        default=0.9,
        help='Min. secondary-to-primary score ratio to report secondary alignments (-p in minimap2). [0.9]')

    parser.add_argument('-v', '--version', action='version', version=__version__, help=SUPPRESS)
    parser.add_argument('-h', '--help', action='help', help=SUPPRESS)

    if len(argv)==1:
        print("             __        \n  __ _  ___ / /__  ___ \n /  ' \\/ -_) / _ \\/ _ \\\n/_/_/_/\\__/_/\\___/_//_/ ver. {}\n".format(__version__))

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
        logger.warning('Folder <{}> exists. Files will be overwritten.'.format(opt.output))

    ## check for input files
    for file in opt.FILE:
        if not os.path.isfile(file):
            logger.critical('File <{}> does not exist.'.format(file))
            sys.exit(2)

    ## check for database
    if not os.path.isdir(opt.db):
        logger.critical('Database folder <{}> does not exist.'.format(opt.db))
        sys.exit(2)
    else:
        files = [os.path.basename(x) for x in glob.glob(os.path.join(opt.db, '*'))]
        if 'metadata.tsv' not in files or 'prot.dmnd' not in files or len([x for x in files if 'nucl' in x and '.mmi' in x]) != 16:
            logger.critical('Database <{}> is not complete or not indexed.'.format(opt.db))
            sys.exit(2)

    ## check for kraken2 database
    if opt.db_kraken is not None:
        if not os.path.isdir(opt.db_kraken):
            logger.critical('Kraken2 database folder <{}> does not exist.'.format(opt.db_kraken))
            sys.exit(2)
        else:
            files = [os.path.basename(x) for x in glob.glob(os.path.join(opt.db_kraken, '*'))]
            if 'ktaxonomy.tsv' not in files or len([x for x in files if 'database' in x]) != 7:
                logger.critical('Kraken2 database <{}> is not complete.'.format(opt.db_kraken))
                sys.exit(2)

    ## check for logical cores
    if opt.threads > os.cpu_count():
        logger.warning('Threads <{}> exceeds number of available logical cores, will use <{}> threads instead.'.format(opt.threads, os.cpu_count()))
        opt.threads = os.cpu_count()
    environ['OMP_NUM_THREADS'] = str(opt.threads)

    ## run
    for i, file in enumerate(opt.FILE):
        if len(opt.FILE) > 1:
            logger.info('Processing file <{}> ({}/{}) ...'.format(file, i+1, len(opt.FILE)))

        GenomeProfiler(file, opt.output, opt.threads).run(
            db=opt.db,
            db_kraken=opt.db_kraken,
            skip_profile=opt.skip_profile,
            skip_clean=opt.skip_clean,
            max_target_seqs=opt.m, evalue=opt.e, identity=opt.i, subject_cover=opt.s,
            secondary_num=opt.n, secondary_ratio=opt.p)

        if i == len(opt.FILE) - 1:
            logger.info('Done.')
        else:
            logger.info('Done.\n')


if __name__ == '__main__':
    cli(sys.argv)
