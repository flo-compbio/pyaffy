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

#cython: profile=False, wraparound=False, boundscheck=False, cdivision=True

"""
Cython parser for CEL files of Affymetrix GeneChip microarrays.

See: http://media.affymetrix.com/support/developer/powertools/changelog/gcos-agcc/cel.html
"""

cimport cython

from libc.stddef cimport size_t
from libc.stdlib cimport malloc, atol
from libc.stdint cimport int16_t, int32_t, uint32_t
from libc.stdio  cimport FILE, fopen, fread, fclose, fgets, sscanf
from libc.string cimport strlen, strcmp
from libc.math cimport NAN

#from libc.string cimport memcpy
#cdef extern from "stdio.h":
#    cdef int EOF

import numpy as np
cimport numpy as np

np.import_array()

import sys
import os
import stat
import subprocess
import logging
import tempfile
import gzip
import io

from configparser import ConfigParser, ParsingError

logger = logging.getLogger(__name__)
logger.debug('__name__: %s', __name__)

cdef inline float decode_float(void* buf):
    return (<float*>buf)[0]

cdef inline int16_t decode_int16(void* buf):
    return (<int16_t*>buf)[0]

cdef inline int32_t decode_int32(void* buf):
    return (<int32_t*>buf)[0]

cdef inline uint32_t decode_uint32(void* buf):
    return (<uint32_t*>buf)[0]

cdef float read_float(void* buf, FILE* fp):
    fread(buf, 4, 1, fp)
    cdef float val = decode_float(buf)
    return val

cdef int read_integer(void* buf, FILE* fp):
    fread(buf, 4, 1, fp)
    cdef int val = <int>decode_int32(buf)
    return val

cdef unsigned int read_DWORD(void* buf, FILE* fp):
    fread(buf, 4, 1, fp)
    cdef unsigned int val = <unsigned int>decode_uint32(buf)
    return val

cdef int read_short(void* buf, FILE* fp):
    fread(buf, 2, 1, fp)
    cdef int val = <int>decode_int16(buf)
    return val

cdef float[::1] read_cell_intensities(void* buf, FILE* fp, int num_rows, int num_cols):
    cdef int i, j
    cdef int pos = 0

    cdef float[::1] y = np.empty(num_rows * num_cols, dtype = np.float32)

    for i in range(num_rows):
        for j in range(num_cols):
            y[pos] = read_float(buf, fp) # the intensity
            read_float(buf, fp) # intensity_std
            read_short(buf, fp) # pixel_count
            pos += 1

    return y

cdef char* read_string(void* buf, FILE* fp):
    cdef int num_bytes = read_integer(buf, fp)
    cdef char* string = <char*>malloc(<size_t>(num_bytes + 1))
    fread(string, <size_t>num_bytes, 1, fp)
    string[num_bytes] = '\0'
    return string

cdef read_tag_val(void* buf, FILE* fp):
    """Returns an OrderedDict containing tag-value entries."""
    s = str(read_string(buf, fp))
    try:
        C = ConfigParser(interpolation = None, delimiters = ('=',), empty_lines_in_values = False)
        C.optionxform = lambda x: x
        C.read_string(u'[Section]\n' + unicode(s))
    except ParsingError:
        C = ConfigParser(interpolation = None, delimiters = (':',), empty_lines_in_values = False)
        C.optionxform = lambda x: x
        C.read_string(u'[Section]\n' + unicode('\n'.join(s.split(';'))))
        
    return C['Section']

cdef short[:,::1] read_coords(void* buf, FILE* fp, int num_coords):
    cdef short[:,::1] C = np.empty((num_coords, 2), dtype = np.int16)
    cdef int i
    for i in range(num_coords):
        C[i,0] = read_short(buf, fp)
        C[i,1] = read_short(buf, fp)
    return C

cdef void apply_mask(float[::1] y, int num_rows, int num_cols, short [:,::1] coords, int num_coords):
    cdef int i, idx
    for i in range(num_coords):
        idx = num_rows * coords[i,1] + coords[i,0]
        y[idx] = NAN

def create_gzip_pipe(path):

    tfh = tempfile.NamedTemporaryFile(mode = 'wb', prefix = 'celparser_',
            delete = False)
    temp = tfh.name
    logger.debug('Temp file: %s', temp)
    tfh.close()
    os.remove(temp) # delete the temp file
    os.mkfifo(temp)
    subproc = subprocess.Popen('gunzip -c "%s" > "%s"' %(path, temp), shell = True)
    return temp

def parse_celfile_v3(path, compressed = True, newline_chars = 2):
    """Parser for CEL file data in Version 3 format (plain-text).

    Currently, no support for parsing information on outliers and masked data.
    """
    assert isinstance(path, (str, unicode))
    assert isinstance(compressed, bool)
    assert isinstance(newline_chars, int)

    cdef int buf_size = 1000
    cdef size_t nl = <size_t>newline_chars
    cdef char* buf = <char*>malloc(buf_size)
    cdef FILE* fp

    cdef int num_cells
    cdef int i
    cdef float[::1] y

    tmp_path = None
    if compressed:
        # file is compressed (we assume gzip)
        # create a named pipe
        logger.debug('Parsing file: %s', path)
        tmp_path = create_gzip_pipe(path)
        fp = fopen(tmp_path, 'r')
    else:
        fp = fopen(path, 'r')

    try:
        # parsing starts here

        # check header information
        fgets(buf, buf_size, fp)
        buf[strlen(buf) - nl] = '\0' # remove newline
        assert strcmp(buf, "[CEL]") == 0
        fgets(buf, buf_size, fp)
        buf[strlen(buf) - nl] = '\0' # remove newline
        assert strcmp(buf, "Version=3") == 0

        # get the number of cells
        # find NumberCells=...
        while sscanf(buf, "NumberCells=%d", &num_cells) <= 0:
            fgets(buf, buf_size, fp)
            buf[strlen(buf) - nl] = '\0'

        fgets(buf, buf_size, fp) # skip CellHeader line

        # read intensities
        y = np.empty(num_cells, dtype = np.float32)
        for i in range(num_cells):
            fgets(buf, buf_size, fp)
            # the following works even if rows starts with whitespace(s)
            sscanf(buf, "%*d %*d %f", &y[i])

    finally:
        fclose(fp)
        if tmp_path is not None:
            os.remove(tmp_path)

    return np.float32(y)

def parse_celfile_v4(path, compressed = True, ignore_outliers = True, ignore_masked = True):
    """Parser for CEL file data in Version 4 format.

    See: http://media.affymetrix.com/support/developer/powertools/changelog/gcos-agcc/cel.html#V4
    Data encoding is little endian.
    """

    cdef FILE* fp
    cdef void* buf = malloc(10)

    """
    def read_subgrid(fh):
        num_rows = read_integer(fh)
        num_cols = read_integer(fh)
        upper_left_x = read_float(fh)
        upper_left_y = read_float(fh)
        upper_right_x = read_float(fh)
        upper_right_y = read_float(fh)
        lower_left_x = read_float(fh)
        lower_left_y = read_float(fh)
        lower_right_x = read_float(fh)
        lower_right_y = read_float(fh)
        left = read_integer(fh)
        top = read_integer(fh)
        right = read_integer(fh)
        bottom = read_integer(fh)
        return (num_rows, num_cols, upper_left_x, upper_left_y, upper_right_x, upper_right_y,
                lower_left_x, lower_left_y, lower_right_x, lower_right_y, left, top, right, bottom)
    """

    #cdef int i, j
    #cdef int x, y
    tmp_path = None

    if compressed:
        # file is compressed (we assume gzip)
        # create a named pipe
        logger.debug('Parsing file: %s', path)
        tmp_path = create_gzip_pipe(path)
        fp = fopen(tmp_path, "rb")
    else:
        fp = fopen(path, "rb")

    y = None
    try:

        magic_number = read_integer(buf, fp)
        assert isinstance(magic_number, int) and magic_number == 64
        version_number = read_integer(buf, fp)
        assert version_number == 4
        num_cols = read_integer(buf, fp)
        num_rows = read_integer(buf, fp)
        num_cells = read_integer(buf, fp)
        assert num_cells == num_rows * num_cols

        logger.debug('Number of rows: %d', num_rows)
        logger.debug('Number of cols: %d', num_cols)
        logger.debug('Number of cells: %d', num_cells)
        header = read_tag_val(buf, fp)
        logger.debug('Header information:')
        logger.debug('  ' + '; '.join(['%s = %s' %(k,v)
                for k,v in header.iteritems()]))
        algo_name = str(read_string(buf, fp))
        logger.debug('Algorithm name: %s', algo_name)
        logger.debug('Algorithm parameters:')
        algo_params = read_tag_val(buf, fp)
        logger.debug('  ' + '; '.join(['%s = %s' %(k,v)
                for k,v in algo_params.iteritems()]))

        cell_margin = read_integer(buf, fp)
        logger.debug('Cell margin: %d', cell_margin)

        num_outlier_cells = read_DWORD(buf, fp)
        logger.debug('# outlier cells: %d', num_outlier_cells)

        num_masked_cells = read_DWORD(buf, fp)
        logger.debug('# masked cells: %d', num_masked_cells)

        num_subgrids = read_integer(buf, fp)
        logger.debug('# subgrids: %d', num_subgrids)
        
        y = read_cell_intensities(buf, fp, num_rows, num_cols)
        logger.debug('# cells: %d', y.size)

        masked_coords = read_coords(buf, fp, num_masked_cells)
        if not ignore_masked:
            apply_mask(y, num_rows, num_cols, masked_coords, num_masked_cells)
        else:
            logger.debug('Ignoring any masked cells')
        
        outlier_coords = read_coords(buf, fp, num_outlier_cells)
        if not ignore_outliers:
            apply_mask(y, num_rows, num_cols, outlier_coords, num_outlier_cells)
        else:
            logger.debug('Ignoring any outlier cells')
                
        #subgrids = []
        #for i in range(num_subgrids):
        #    subgrids.append(read_subgrid(fh))

    finally:
        if tmp_path is not None:
            os.remove(tmp_path)

    return np.float32(y)

def try_open_gzip(path):

    fh = None

    try:
        fh = gzip.open(path)
        fh.read(1)
    except IOError:
        pass
    else:
        fh = gzip.open(path)

    return fh

def parse_cel(path):
    """Front-end for parsing a CEL file containing expression data.

    This function automatically determines the CEL file format. The possible
    formats are:

    * Version 3 (plain-text)
    * Version 4 (binary, little-endian)
    * Command Console (binary, big-endian)

    See http://media.affymetrix.com/support/developer/powertools/changelog/gcos-agcc/cel.html
    for format specifications.

    This function also detects whether the input file is gzip'ed or not. If it
    is, gunzip is run in an independent process using a named pipe.

    Parameters
    ----------
    path: str or unicode
        The path of the CEL file.

    Returns
    -------
    np.ndarray of type np.float32
        The intensities from the array.
    """

    assert isinstance(path, (str, unicode))

    if not os.path.isfile(path):
        raise IOError('File "%s" not found.' %(path))

    version = 0

    # test if file is gzipped
    compressed = False
    fh = None
    try:
        fh = try_open_gzip(path)
        if fh is not None:
            compressed = True
        else:
            fh = io.open(path, 'rb')

        version = ord(fh.read(1)) # read the first byte

    finally:
        if fh is not None:
            fh.close()

    y = None
    if version == 59:
        # command console generic data file format (binary, big-endian)
        raise NotImplemented('No parser for command console format!')
    elif version == 64:
        # version 4 format (binary, little-endian)
        y = parse_celfile_v4(path, compressed = compressed)
    else:
        # version 3 format (plain-text)
        y = parse_celfile_v3(path, compressed = compressed)

    return y
