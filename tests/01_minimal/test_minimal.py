# Copyright (c) 2016 Florian Wagner
#
# This file is part of GO-PCA.
#
# GO-PCA is free software: you can redistribute it and/or modify
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
## from builtins import *
# from builtins import open
from builtins import str as text
# text = str

from collections import OrderedDict
import os

import numpy as np

from pyaffy import rma


def test_download(my_cdf_file, my_cel_files):
    assert os.path.isfile(my_cdf_file)
    for path in my_cel_files:
        assert os.path.isfile(path)


def test_rma(my_cdf_file, my_cel_files):
    sample_cel_files = OrderedDict([
        ('Sample %d' % (i+1), path) for i, path in enumerate(my_cel_files)
    ])
    genes, samples, X = rma(my_cdf_file, sample_cel_files)

    assert isinstance(genes, list)
    assert isinstance(samples, list)
    assert len(samples) == 2
    assert isinstance(X, np.ndarray) and X.ndim == 2