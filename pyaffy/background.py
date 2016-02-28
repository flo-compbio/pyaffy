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

import logging
from math import floor

import numpy as np
from scipy.stats import norm

logger = logging.getLogger(__name__)

def rma_bg_correct(Y, make_copy = False):
    """RMA background correction.

    Parameters
    ----------
    Y: np.ndarray (ndim = 1, dtype = np.float32)
        The microarray intensity values (on a linear scale).
    make_copy: bool
        Whether or to make a copy of the data or modify it in-place.
    """
    assert isinstance(Y, np.ndarray)
    
    if make_copy:
        Y = Y.copy()
    
    n = Y.shape[1]
    
    for j in range(n):
        
        # find missing data (= NaN)
        missing = np.isnan(Y[:,j])
        y = Y[~missing,j]
        
        ### estimate mu using simple binning (histogram)

        # use a fixed number of bins
        num_bins = 100
        lower = np.amin(y)
        upper = np.percentile(y, 75.0)
        bin_width = max(floor((upper - lower) / num_bins), 1.0)
        bin_edges = np.arange(lower, upper, bin_width)
        num_bins = bin_edges.size - 1
        
        # binning
        binned = np.digitize(y, bins = bin_edges) - 1
        binned = binned[binned < num_bins]
        bc = np.bincount(binned)
        amax = np.argmax(bc)
        max_x = lower + (amax + 0.5) * bin_width
        mu = max_x
        logger.debug('Mu: %.2f', mu)

        ### estimate sigma

        # 1. Select probes with values smaller than mu
        y_low = y[y < mu]
        # 2. Estimate their standard deviation (using mu as the mean)
        sigma = pow(np.sum(np.power(y_low - mu, 2.0)) / (y_low.size - 1), 0.5)
        # 3. Arbitrarily multiply standard deviation by square root of two
        sigma *= pow(2.0, 0.5)
        logger.debug('Sigma: %.2f', sigma)

        ### estimate alpha

        # we simply fix alpha to 0.03
        alpha = 0.03

        ### calculate background-corrected intensities
        a = y - mu - alpha * pow(sigma, 2.0)
        y_adj = a + sigma * np.exp(norm.logpdf(a / sigma) - norm.logcdf(a / sigma))
        Y[~missing,j] = y_adj

    return Y

