# Copyright (c) 2016 Florian Wagner
#
# This file is part of pyAffy.
#
# pyAffy is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License, Version 3,
# as published by the Free Software Foundation.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.

from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
_oldstr = str
from builtins import *

import os
import time
import logging
import collections

import numpy as np

from genometools.expression import quantile_normalize as qnorm
from genometools.expression import ExpMatrix

from .cdfparser import parse_cdf
from . import celparser
from .celparser import parse_cel
from .medpolish import medpolish
from .background import rma_bg_correct

logger = logging.getLogger(__name__)

def rma(
        cdf_file,
        sample_cel_files,
        pm_probes_only = True,
        bg_correct = True,
        quantile_normalize = True,
        medianpolish = True
    ):
    """Perform RMA on a set of samples.

    Parameters
    ----------
    cdf_file: str
        The path of the Brainarray CDF file to use.
        Note: Brainarray CDF files can be downloaded from
            http://brainarray.mbni.med.umich.edu/Brainarray/Database/CustomCDF/genomic_curated_CDF.asp
    sample_cel_files: collections.OrderedDict (st => str)
        An ordered dictionary where each key/value-pair corresponds to a
        sample. The *key* is the sample name, and the *value* is the (absolute)
        path of the corresponding CEL file. The CEL files can be gzip'ed.
    pm_probes_only: bool, optional
        Whether or not to only use PM (perfect match) probes and ignore all MM
        (mismatch) probes. [True]
    bg_correct: bool, optional
        Whether or not to apply background correction. [True]
    quantile_normalize: bool, optional
        Whether or not to apply quantile normalization. [True]
    medianpolish: bool, optional
        Whether or not to apply medianpolish. [True]

    Returns
    -------
    genes: tuple of str
        The list of gene names.
    samples: tuple of str
        The list of sample names.
    X: np.ndarray (ndim = 2, dtype = np.float32)
        The expression matrix (genes-by-samples).

    Examples
    --------
    >>> from collections import OrderedDict
    >>> import pyaffy
    >>> cdf_file = '/path/to/brainarray/cdf/HGU133Plus2_Hs_ENTREZG.cdf'
    >>> sample_cel_files = OrderedDict([
            ['Sample 1', '/path/to/sample_1.CEL.gz'],
            ['Sample 2', '/path/to/sample_2.CEL.gz'],
        ])
    >>> genes, samples, X = pyaffy.rma(cdf_file, sample_cel_files)
    """

    ### checks
    assert isinstance(cdf_file, (str, _oldstr))
    assert os.path.isfile(cdf_file), \
            'CDF file "%s" does not exist!' %(cdf_file)

    assert isinstance(sample_cel_files, collections.OrderedDict)
    for sample, cel_file in sample_cel_files.items():
        assert isinstance(sample, (str, _oldstr))
        assert isinstance(cel_file, (str, _oldstr))
        assert os.path.isfile(cel_file), \
                'CEL file "%s" does not exist!' %(cel_file)

    assert isinstance(pm_probes_only, bool)
    assert isinstance(bg_correct, bool)
    assert isinstance(quantile_normalize, bool)
    assert isinstance(medianpolish, bool)

    t00 = time.time()

    ### read CDF data
    logger.info('Parsing CDF file.')
    t0 = time.time()
    # parse the CDF file
    probe_type = 'pm'
    if not pm_probes_only:
        probe_type = 'all'
    name, num_rows, num_cols, pm_probesets = \
            parse_cdf(cdf_file, probe_type=probe_type)

    # concatenate indices of all PM probes into one long vector
    pm_sel = np.concatenate(list(pm_probesets.values()))

    t1 = time.time()
    logger.info('CDF file parsing time: %.2f s', t1 - t0)
    logger.info('CDF array design name: %s', name)
    logger.info('CDF rows / columns: %d x %d', num_rows, num_cols)

    ### read CEL data
    logger.info('Parsing CEL files...')
    t0 = time.time()
    p = pm_sel.size
    n = len(sample_cel_files)
    Y = np.empty((p, n), dtype = np.float32)

    samples = []
    sub_logger = logging.getLogger(celparser.__name__)
    sub_logger.setLevel(logging.WARNING)
    for j, (sample, cel_file) in enumerate(sample_cel_files.items()):
        logger.debug('Parsing CEL file for sample "%s": %s', sample, cel_file)
        samples.append(sample)
        y = parse_cel(cel_file)
        Y[:,j] = y[pm_sel]
    sub_logger.setLevel(logging.NOTSET)
    t1 = time.time()
    logger.info('CEL files parsing time: %.1f s.', t1 - t0)

    ### background correction
    if bg_correct:
        logger.info('Performing background correction...')
        t0 = time.time()
        Y = rma_bg_correct(Y)
        t1 = time.time()
        logger.info('Background correction time: %.1f s.', t1 -t0)
    else:
        logger.info('Skipping background correction.')

    matrix = ExpMatrix(genes=pm_sel, samples=samples, X=Y)

    ### quantile normalization
    if quantile_normalize:
        logger.info('Performing quantile normalization...')
        t0 = time.time()
        matrix = qnorm(matrix)
        t1 = time.time()
        logger.info('Quantile normalization time: %.1f s.', t1 - t0)
    else:
        logger.info('Skipping quantile normalization.')

    ### convert intensities to log2-scale
    Y = np.log2(matrix.values)

    ### probeset summarization (with or without median polish)
    method = 'with'
    if not medianpolish:
        method = 'without'
    logger.info('Summarize probeset intensities (%s medianpolish)...', method)

    t0 = time.time()
    p = len(pm_probesets)
    n = Y.shape[1]
    X = np.empty((p, n), dtype = np.float32)
    cur = 0
    num_converged = 0
    genes = []
    for i , (gene_id, probes) in enumerate(pm_probesets.items()):
        genes.append(gene_id)

        if medianpolish:
            #X_sub = np.ascontiguousarray(Y[cur:(cur + probes.size),:])
            X_sub = Y[cur:(cur + probes.size),:]
            _, row_eff, col_eff, global_eff, converged, num_iter = medpolish(X_sub, copy = False)
            X[i,:] = col_eff + global_eff
            if converged:
                num_converged += 1
         
        else:
            # simply use median across probes
            X[i,:] = np.median(Y[cur:(cur+probes.size),:], axis = 0) 
            #X[i,:] = np.ma.median(X_sub, axis = 0)

        cur += probes.size

    t1 = time.time()
    logger.info('Probeset summarization time: %.2f s.', t1 - t0)

    if medianpolish:
        logger.debug('Converged: %d / %d (%.1f%%)',
                num_converged, p, 100 * (num_converged / float(p)))
    
    ### report total time
    t11 = time.time()
    logger.info('Total RMA time: %.1f s.', t11 - t00)

    ### sort alphabetically by gene name
    a = np.lexsort([genes])
    genes = [genes[i] for i in a]
    X = X[a,:]

    return genes, samples, X
