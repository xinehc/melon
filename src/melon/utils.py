import os
import re
import sys
import subprocess
import logging

## setup logging format
if not sys.stderr.isatty():
    os.environ['TQDM_DISABLE'] = '1'

from tqdm import tqdm
from tqdm.contrib.logging import logging_redirect_tqdm

logging.basicConfig(
    level='INFO',
    format='[%(asctime)s] %(levelname)s: %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S')

if sys.stderr.isatty():
    logging.addLevelName(logging.INFO,
                         f'\033[1m{logging.getLevelName(logging.INFO)}\033[1;0m')
    logging.addLevelName(logging.WARNING,
                         f'\033[1m\x1b[33;20m{logging.getLevelName(logging.WARNING)}\033[1;0m')
    logging.addLevelName(logging.CRITICAL,
                         f'\033[1m\x1b[31;20m{logging.getLevelName(logging.CRITICAL)}\033[1;0m')

logger = logging.getLogger(__name__)


def sort_coordinate(start, end):
    '''
    Convert coordinate to 0-based, then sort.
    '''
    return (start - 1, end) if start < end else (end - 1, start)


def compute_overlap(coordinates, func=None):
    '''
    Compute overlap between two 0-based coordinates.
    '''
    qstart, qend, sstart, send = coordinates
    overlap = min(qend, send) - max(qstart, sstart)

    if func is None:
        return overlap

    assert func in {max, min}, 'Use either max or min for aggregation.'
    return func(overlap / (qend - qstart), overlap / (send - sstart))


def get_filename(file, output=None, extension=None):
    '''
    Get filename of a file, possibly add output dir and extension.
    '''
    filename = re.sub(r'\.f(ast)?[aq](\.gz)?$', '', os.path.basename(file))
    if output is not None:
        filename = os.path.join(output, filename)

    if extension is not None:
        filename += extension
    return filename


def extract_sequence(file, ids, output=None):
    '''
    Extract sequences from a source fa/fq file using seqkit.
    '''
    cmd = ['seqkit', 'grep', '-f', '-', file, '--quiet']
    if output is not None:
        subprocess.run(cmd + ['-o', output], check=True, input='\n'.join(ids) + '\n', text=True)
    else:
        return subprocess.run(cmd, check=True, input='\n'.join(ids) + '\n', text=True, capture_output=True).stdout
