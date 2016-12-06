# Copyright (c) 2016 Florian Wagner
#
# This file is part of pyAffy.
#
# pAffy is free software: you can redistribute it and/or modify
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

"""Parser for Affymetrix CEL files.

"""

from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
_oldstr = str
from builtins import *

import os
import struct
import dateutil.parser
import codecs
import sys
import logging

import numpy as np
from configparser import ConfigParser, ParsingError

from genometools import misc
from . import CEL

logger = logging.getLogger(__name__)

class CELParser(object):

    def __init__(self, path):
        assert isinstance(path, (str, _oldstr))
        self.path = path

    def _parse_cel_v4_intensities(self):
        # Version 4 format, http://media.affymetrix.com/support/developer/powertools/changelog/gcos-agcc/cel.html#V4
        # data encoding is little endian

        read = [0]

        def decode_unicode(s):
            return codecs.decode(s, 'UTF-16-BE')

        def decode_ascii(s):
            return codecs.decode(s, 'ascii')

        def decode_float(bytes):
            # this is an actual float, not a double!
            return struct.unpack('<f', bytes)[0]

        def decode_int32(bytes):
            return struct.unpack('<i', bytes)[0]

        def decode_uint32(bytes):
            return struct.unpack('<I', bytes)[0]

        def decode_int8(bytes):
            return struct.unpack('<b', bytes)[0]

        def decode_uint8(bytes):
            return struct.unpack('<B', bytes)[0]

        def decode_int16(bytes):
            return struct.unpack('<h', bytes)[0]

        def decode_uint16(bytes):
            return struct.unpack('<H', bytes)[0]

        def read_float(bytes):
            read[0] += 4
            return decode_float(fh.read(4))

        def read_short(fh):
            read[0] += 2
            return decode_int16(fh.read(2))

        def read_integer(fh):
            read[0] += 4
            return decode_int32(fh.read(4))

        def read_DWORD(fh):
            read[0] +=4
            return decode_uint32(fh.read(4))

        def read_raw(fh):
            # reads an int x and then x raw bytes
            bytes = read_integer(fh)
            read[0] += bytes
            return fh.read(bytes)

        def read_tag_val(fh):
            """Returns an OrderedDict containing tag-value entries."""
            raw = codecs.decode(read_raw(fh), encoding='iso-8859-1')
            logger.debug('Tag/Value string:\n%s', raw)
            try:
                C = ConfigParser(interpolation = None, delimiters = ('=',), empty_lines_in_values = False)
                C.optionxform = lambda x: x
                C.read_string('[Section]\n' + raw)
            except ParsingError:
                C = ConfigParser(interpolation = None, delimiters = (':',), empty_lines_in_values = False)
                C.optionxform = lambda x: x
                C.read_string('[Section]\n' + '\n'.join(raw.split(';')))
                
            return C['Section']

        def read_cell(fh):
            intensity = read_float(fh)
            intensity_std = read_float(fh)
            pixel_count = read_short(fh)
            return (intensity, intensity_std, pixel_count)

        def read_coords(fh):
            x = read_short(fh)
            y = read_short(fh)
            return (x, y)

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

        with misc.smart_open_read(self.path, mode = 'rb', try_gzip = True) as fh:
            #read_file_header(fh)
            
            magic_number = read_integer(fh)
            assert isinstance(magic_number, int) and magic_number == 64
            version_number = read_integer(fh)
            assert version_number == 4
            num_cols = read_integer(fh)
            num_rows = read_integer(fh)
            num_cells = read_integer(fh)
            logger.debug('Number of rows: %d', num_rows)
            logger.debug('Number of cols: %d', num_cols)
            header = read_tag_val(fh)
            logger.debug('; '.join(['%s = %s' %(k,v) for k,v in header.items()]))
            algo_name = read_raw(fh)
            logger.debug('Algorithm name: %s', algo_name)
            algo_params = read_tag_val(fh)
            logger.debug('; '.join(['%s = %s' %(k,v) for k,v in algo_params.items()]))
            cell_margin = read_integer(fh)
            num_outlier_cells = read_DWORD(fh)
            num_masked_cells = read_DWORD(fh)
            num_subgrids = read_integer(fh)
            logger.debug('Cell margin: %d', cell_margin)
            logger.debug('# outlier cells: %d', num_outlier_cells)
            logger.debug('# masked cells: %d', num_masked_cells)
            logger.debug('# sub-grids: %d', num_subgrids)
            
            cells = []
            for j in range(num_cols):
                for i in range(num_rows):
                    cells.append(read_cell(fh))
            logger.debug('# cells: %d', len(cells))
            
            masked = []
            for i in range(num_masked_cells):
                masked.append(read_coords(fh))
                
            outliers = []
            for i in range(num_outlier_cells):
                outliers.append(read_coords(fh))
                
            subgrids = []
            for i in range(num_subgrids):
                subgrids.append(read_subgrid(fh))

            return np.float64([c[0] for c in cells])


    def parse_intensities(self):

        # check version
        version = None
        with misc.smart_open_read(self.path, mode = 'rb', try_gzip = True) \
                as fh:
            v = ord(fh.read(1)) # look at first byte
            if v == 59:
                # command console generic file format
                version = 'CCG'
            if v == 64:
                # binary cel format, version 4
                version = 4

        assert version in (4, 'CCG')
        logger.info('Detected CEL file format version "%s".', str(version))

        y = None
        if version == 'CCG':
            cel = CEL.read_cel(self.path)
            ds = cel.data_groups[0].get_data_set_by_name('Intensity')
            y = np.float64(ds.data)[:,0]
        elif version == 4:
            y = self._parse_cel_v4_intensities()

        return y
