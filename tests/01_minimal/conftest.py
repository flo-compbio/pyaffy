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

#import str as text
#text = str

# import os
import shutil

import pytest
import requests
import logging

logging.basicConfig(level=logging.INFO)

logger = logging.getLogger(__name__)

DATA_URL = 'https://www.dropbox.com/sh/z8zafx9oogodky1/'


def download_file(url, path):
    """Downloads a file."""
    r = requests.get(url, stream=True)
    if r.status_code == 200:
        with open(path, 'wb') as f:
            r.raw.decode_content = True
            shutil.copyfileobj(r.raw, f)


@pytest.fixture(scope='session')
def my_data_pypath(tmpdir_factory):
    pypath = tmpdir_factory.mktemp('pyaffy_data', numbered=False)
    return pypath


@pytest.fixture(scope='session')
def my_output_pypath(tmpdir_factory):
    pypath = tmpdir_factory.mktemp('pyaffy_output', numbered=False)
    return pypath


@pytest.fixture(scope='session')
def my_cdf_file(my_data_pypath):
    logger.info('Starting download of CDF file...')
    url = DATA_URL + 'AACcuI150VSFWl4ji4_opG7Ba/HGU133Plus2_Hs_20_ENTREZG.cdf?dl=1'
    path = text(my_data_pypath.join('HGU133Plus2_Hs_20_ENTREZG.cdf'))
    download_file(url, path)
    return path


@pytest.fixture(scope='session')
def my_cel_file1(my_data_pypath):
    logger.info('Starting download of CDF file 1...')
    url = DATA_URL + 'AADBnMN8wFR-nao1Ze695Cmaa/AFX_2_A1.CEL.gz?dl=1'
    path = text(my_data_pypath.join('AFX_2_A1.CEL.gz'))
    download_file(url, path)
    return path


@pytest.fixture(scope='session')
def my_cel_file2(my_data_pypath):
    logger.info('Starting download of CDF file 2...')
    url = DATA_URL + 'AADbbAoedfLqSvsgWnJZufLHa/AFX_2_A2.CEL.gz?dl=1'
    path = text(my_data_pypath.join('AFX_2_A2.CEL.gz'))
    download_file(url, path)
    return path

@pytest.fixture(scope='session')
def my_cel_files(my_cel_file1, my_cel_file2):
    return [
        my_cel_file1,
        my_cel_file2,
    ]