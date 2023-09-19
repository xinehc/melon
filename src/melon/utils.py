import os
import re
import subprocess
import logging
import numpy as np

## setup logger format
logging.basicConfig(
    level="INFO",
    format="[%(asctime)s] %(levelname)s: %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S")

logging.addLevelName(logging.WARNING,
                     f'\033[1m\x1b[33;20m{logging.getLevelName(logging.WARNING)}\033[1;0m')
logging.addLevelName(logging.CRITICAL,
                     f'\033[1m\x1b[31;20m{logging.getLevelName(logging.CRITICAL)}\033[1;0m')

logger = logging.getLogger(__name__)


def sort_coordinate(start, end):
    '''
    Convert coordinate to 0-based, then sort.
    '''
    return (start-1, end) if start < end else (end-1, start)


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


def reassign_taxonomy(matrix, eps=1e-5, max_iteration=100):
    '''
    Reassign multi-mapped reads with EM.
    '''
    n_reads, n_mappings = matrix.shape

    ## init
    p_reads = np.zeros((n_reads, n_mappings))
    p_mappings = np.ones(n_mappings) / n_mappings
    p_mappings_hist = p_mappings.copy()

    iteration = 0
    while iteration < max_iteration:
        iteration += 1

        ## e-step
        p_reads = np.divide(p_mappings * matrix, np.dot(matrix, p_mappings).reshape(-1, 1))

        ## m-step
        p_mappings = np.sum(p_reads, axis=0) / n_reads

        ## check convergence
        if np.sum(np.abs(p_mappings - p_mappings_hist)) < eps:
            break

        ## update p_reads_hist
        np.copyto(p_mappings_hist, p_mappings)

    ## return assignments
    assignments = [np.where(row == row.max())[0].tolist() for row in p_reads]
    return assignments


def extract_sequences(file, ids):
    '''
    Extract sequences from a source fa/fq file using seqkit.
    '''
    cmd = ['seqkit', 'grep', '-f', '-', file]
    return subprocess.run(cmd, check=True, input='\n'.join(ids) + '\n', text=True, capture_output=True).stdout
