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

from __future__ import (absolute_import, division,
                        print_function)
# _oldstr = str
from builtins import str as text
from builtins import open

cimport cython

from libc.stddef cimport size_t
from libc.stdlib cimport malloc, atol
from libc.stdint cimport int16_t, int32_t, uint32_t
from libc.stdio  cimport FILE, fopen, fread, fclose, fgets, sscanf
from libc.string cimport strlen, strcmp
# from libc.math cimport NAN

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
import struct
import dateutil
import codecs
from collections import OrderedDict


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

    s = read_string(buf, fp).decode('iso-8859-1')
    # s = codecs.decode(read_string(buf, fp), encoding='iso-8859-1')
    try:
        C = ConfigParser(interpolation = None, delimiters = ('=',), empty_lines_in_values = False)
        C.optionxform = lambda x: x
        C.read_string(str('[Section]\n') + s)
    except ParsingError:
        C = ConfigParser(interpolation = None, delimiters = (':',), empty_lines_in_values = False)
        C.optionxform = lambda x: x
        C.read_string(str('[Section]\n') + '\n'.join(s.split(';')))
        
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
        #y[idx] = NAN
        y[idx] = float('nan')

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
    assert isinstance(path, (text, str))
    assert isinstance(compressed, bool)
    assert isinstance(newline_chars, int)

    cdef int buf_size = 1000
    cdef size_t nl = <size_t>newline_chars
    cdef char* buf = <char*>malloc(buf_size)
    cdef FILE* fp

    cdef int num_cells
    cdef int i
    cdef float[::1] y

    final_path = path
    if compressed:
        # file is compressed (we assume gzip)
        # create a named pipe
        logger.debug('Parsing file: %s', path)
        final_path = create_gzip_pipe(path)

    path_bytes = final_path.encode('UTF-8')
    cdef char* c_path_bytes = path_bytes
    fp = fopen(c_path_bytes, 'r')

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
        if compressed:
            os.remove(final_path)

    return np.float32(y)


def parse_celfile_v4(path, compressed=True, ignore_outliers=True, ignore_masked=True):
    """Parser for CEL file data in Version 4 format.

    See: http://media.affymetrix.com/support/developer/powertools/changelog/gcos-agcc/cel.html#V4
    Data encoding is little endian.
    """
    assert isinstance(path, (text, str))
    assert isinstance(compressed, bool)
    assert isinstance(ignore_outliers, bool)
    assert isinstance(ignore_masked, bool)

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

    final_path = path
    if compressed:
        # file is compressed (we assume gzip)
        # create a named pipe
        logger.debug('Parsing file: %s', path)
        final_path = create_gzip_pipe(path)

    path_bytes = final_path.encode('UTF-8')
    cdef char* c_path_bytes = path_bytes
    fp = fopen(c_path_bytes, 'r')

    try:

        y = None
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
                for k,v in header.items()]))
        algo_name = str(read_string(buf, fp))
        logger.debug('Algorithm name: %s', algo_name)
        logger.debug('Algorithm parameters:')
        algo_params = read_tag_val(buf, fp)
        logger.debug('  ' + '; '.join(['%s = %s' %(k,v)
                for k,v in algo_params.items()]))

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
        if compressed:
            os.remove(final_path)

    return np.float32(y)


#cdef unsigned char read_UBYTE(char* buf, FILE* fp):
#    fread(buf, 1, 1, fp)
#    cdef unsigned char c = (<unsigned char*>buf)[0]
#    return c

#cdef int read_INT(char* buf, char* buf2, FILE* fp):
#    fread(buf2, 4, 1, fp)
#    reverse_copy(buf2, buf, 4)
#    cdef int val = <int>decode_int32(buf)
#    return val

#cdef unsigned int read_UINT(char* buf, char* buf2, FILE* fp):
#    fread(buf2, 4, 1, fp)
#    reverse_copy(buf2, buf, 4)
#    cdef unsigned int val = <int>decode_uint32(buf)
#    return val


cdef void reverse_copy(const char* src, char* dst, int num):
    cdef int i
    for i in range(num):
        dst[num - i - 1] = src[i]

cdef float read_FLOAT(char* buf, char* data):
    reverse_copy(data, buf, 4)
    cdef float val = (<float*>buf)[0]
    return val

cdef float[::1] read_cc_intensities(char* data, unsigned int num_values):
    cdef float[::1] y = np.empty(num_values, dtype = np.float32)
    cdef unsigned int i
    cdef char* buf = <char*>malloc(10)
    for i in range(num_values):
        y[i] = read_FLOAT(buf, data)
        data += 4
    return y

def parse_celfile_cc(path, compressed=True, ignore_outliers=True, ignore_masked=True):
    """Parser for CEL file data in Command Console version 1 format.

    See: http://media.affymetrix.com/support/developer/powertools/changelog/gcos-agcc/cel.html#calvin
    Note: Data byte order is big endian!
    """
    assert isinstance(path, (text, str))
    assert isinstance(compressed, bool)
    assert isinstance(ignore_outliers, bool)
    assert isinstance(ignore_masked, bool)

    read = [0]
        
    decode_unicode = lambda s: codecs.decode(s, 'UTF-16-BE')
    #decode_unicode = lambda s: unicode(s, encoding = 'UTF-16-BE')
    decode_ascii = lambda s: codecs.decode(s, 'ascii')
    #decode_ascii = lambda s: unicode(s, encoding = 'ascii')
    decode_float = lambda s: struct.unpack('>f', s)[0]
    decode_int32 = lambda s: struct.unpack('>i', s)[0]
    decode_uint32 = lambda s: struct.unpack('>I', s)[0]
    decode_int8 = lambda s: struct.unpack('>b', s)[0]
    decode_uint8 = lambda s: struct.unpack('>B', s)[0]
    decode_int16 = lambda s: struct.unpack('>h', s)[0]
    decode_uint16 = lambda s: struct.unpack('>H', s)[0]

    def read_float(bytes):
        read[0] += 4
        return decode_float(fh.read(4))

    def read_byte(fh):
        read[0] += 1
        return decode_int8(fh.read(1))

    def read_ubyte(fh):
        read[0] += 1
        return decode_uint8(fh.read(1))

    def read_short(fh):
        read[0] += 2
        return decode_int16(fh.read(2))

    def read_ushort(fh):
        read[0] += 2
        return decode_uint16(fh.read(2))

    def read_int(fh):
        read[0] += 4
        return decode_int32(fh.read(4))

    def read_uint(fh):
        read[0] +=4
        return decode_uint32(fh.read(4))

    def read_raw(fh):
        # reads an int x and then x raw bytes
        bytes = read_int(fh)
        read[0] += bytes
        return fh.read(bytes)

    def read_string(fh):
        strlen = read_int(fh)
        read[0] += strlen
        return fh.read(strlen)

    def read_wstring(fh):
        strlen = read_int(fh)
        s = fh.read(2 * strlen)
        read[0] += (2 * strlen)
        return decode_unicode(s)

    def read_guid(fh):
        return read_string(fh)

    def read_datetime(fh):
        s = read_wstring(fh)
        logger.debug('DateTime string: %s|||', s)
        dt = None
        if s:
            dt = dateutil.parser.parse(s).replace(tzinfo = None)
        #s = u'2015-02-20T13:52:11Z'
        #print s
        return dt

    def read_locale(fh):
        loc = read_wstring(fh)
        return loc[:2], loc[3:]

    def read_value(fh):
        raw = read_raw(fh)
        return raw

    def read_type(fh):
        return read_wstring(fh)

    def read_file_header(fh):
        magic_number = read_ubyte(fh)
        logger.debug('Magic number: %d', magic_number)
        # check magic number
        assert isinstance(magic_number, int) and magic_number == 59
        version_number = read_ubyte(fh)
        # check version number
        assert version_number == 1
        num_data_groups = read_int(fh)
        logger.debug('# data groups: %d', num_data_groups)
        first_data_group_pos = read_uint(fh)
        return num_data_groups, first_data_group_pos

    def read_header_param(fh):
        v1 = read_wstring(fh)
        v2 = read_value(fh)
        v3 = read_type(fh)

        if v3 == 'text/plain':
            v2 = decode_unicode(v2.rstrip('\x00'))
        elif v3 == 'text/ascii':
            v2 = decode_ascii(v2.rstrip('\x00'))
        elif v3 == 'text/x-calvin-float':
            v2 = decode_float(v2[:4])
        elif v3 == 'text/x-calvin-integer-32':
            v2 = decode_int32(v2[:4])
        elif v3 == 'text/x-calvin-unsigned-integer-32':
            v2 = decode_uint32(v2[:4])
        elif v3 == 'text/x-calvin-unsigned-integer-8':
            v2 = decode_uint8(v2[:1])
        elif v3 == 'text/x-calvin-unsigned-integer-16':
            v2 = decode_uint16(v2[:2])
        elif v3 == 'text/x-calvin-integer-8':
            v2 = decode_int8(v2[:1])
        elif v3 == 'text/x-calvin-integer-16':
            v2 = decode_int16(v2[:2])

        return (v1, v2)

    def read_data_header(fh):
        data_type_id = read_guid(fh)
        logger.debug('Data type identifier: %s', data_type_id)
        file_id = read_guid(fh)
        logger.debug('File identifier: %s', file_id)
        creation_time = read_datetime(fh)
        #print creation_time
        iso639, iso3166 = read_locale(fh)
        locale = '-'.join([iso639, iso3166])
        n_params = read_int(fh)
        logger.debug('Number of parameters (name/value/type triplets): %d', n_params)
        params = OrderedDict()
        for i in range(n_params):
            params.update([read_header_param(fh)])

        num_parents = read_int(fh)
        logger.debug('Number of parent file headers: %d', num_parents)
        
        parent_headers = []
        for i in range(num_parents):
            logger.debug('')
            logger.debug('-------------------------------------------')
            parent_headers.append(read_data_header(fh))

    def read_col(fh):
        name = read_wstring(fh)
        valtype = read_byte(fh)
        size = read_int(fh)
        return (name, valtype, size)

    def read_dataset(fh):
        logger.debug('New DataSet! Bytes read up until this point: %d', read[0])
        data_pos = read_uint(fh)
        next_pos = read_uint(fh)
        data_size = next_pos - data_pos
        name = read_wstring(fh)
        num_params = read_int(fh)

        # read parameters
        params = OrderedDict()
        for i in range(num_params):
            params.update([read_header_param(fh)])

        num_cols = read_uint(fh)
        logger.debug('DataSet / # cols: %d', num_cols)
        cols = []
        for i in range(num_cols):
            cols.append(read_col(fh))
        col_names = [c[0] for c in cols]
        logger.debug('DataSet / col names: %s', ', '.join(col_names))

        num_rows = read_uint(fh)
        logger.debug('DataSet / # rows: %d', num_rows)
        logger.debug('DataSet / data position: %d', data_pos)
        logger.debug('DataSet / next position: %d', next_pos)
        logger.debug('DataSet / data size: %d', data_size)
        logger.debug('DataSet / name: %s', name)
        logger.debug('DataSet / # parameters: %d', num_params)
        num_chars = next_pos - read[0]
        data = fh.read(num_chars)
        read[0] += num_chars
        y = None
        # this assumes that the column name etc. is 100% fixed
        if name == 'Intensity':
            assert len(cols) == 1
            col_name = cols[0][0]
            col_type = cols[0][1]
            col_size = cols[0][2]
            assert col_name == 'Intensity'
            assert col_type == 6
            assert col_size == 4
            y = read_cc_intensities(data, num_rows)
        return y

    def read_data_group(fh):
        logger.debug('New data group! Number of bytes read so far: %d', read[0])
        next_pos = read_uint(fh)
        dataset_pos = read_uint(fh)
        num_datasets = read_int(fh)
        name = read_wstring(fh)
        logger.debug('Data group name: %s', name)
        y = None
        for i in range(num_datasets):
            d = read_dataset(fh)
            if d is not None:
                y = d
        return y

    fh = None
    y = None
    final_path = path
    if compressed:
        # file is compressed (we assume gzip)
        # create a named pipe
        logger.debug('Parsing file: %s', path)
        final_path = create_gzip_pipe(path)

    fh = open(final_path, mode='rb')
    try:
        num_data_groups, data_pos = read_file_header(fh)
        assert num_data_groups == 1 # for expression CEL file
        logger.debug('# data groups: %d', num_data_groups)
        logger.debug('pos. of first data group: %d', data_pos)
        #header = read_data_header(fh)
        read_data_header(fh)
        assert data_pos == read[0] # position of the first data group
        y = read_data_group(fh)
        #data_groups = []
        #for i in range(num_data_groups):
        #    data_groups.append(read_data_group(fh))

    finally:
        fh.close()
        if compressed:
            os.remove(final_path)

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
    path: str
        The path of the CEL file.

    Returns
    -------
    np.ndarray of type np.float32
        The intensities from the array.
    """

    assert isinstance(path, (text, str))

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
            fh = open(path, 'rb')

        version = ord(fh.read(1)) # read the first byte

    finally:
        if fh is not None:
            fh.close()

    y = None
    if version == 59:
        # command console generic data file format (binary, big-endian)
        y = parse_celfile_cc(path, compressed = compressed)
    elif version == 64:
        # version 4 format (binary, little-endian)
        y = parse_celfile_v4(path, compressed = compressed)
    else:
        # version 3 format (plain-text)
        y = parse_celfile_v3(path, compressed = compressed)

    return y
