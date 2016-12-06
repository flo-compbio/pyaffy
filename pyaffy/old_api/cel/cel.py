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

"""API for data in Affymetrix CEL files.

"""
from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
_oldstr = str
from builtins import *

import os
import sys
import struct
import datetime
import codecs
import logging
from collections import OrderedDict

import dateutil.parser

from genometools import misc

logger = logging.getLogger(__name__)

class CELHeader(object):
    def __init__(self, data_type_id, file_id, creation_time, locale, params = None, parent_headers = None):
        
        if params is None:
            params = OrderedDict()
        
        if parent_headers is None:
            parent_headers = []

        assert isinstance(data_type_id, (str, _oldstr))
        assert isinstance(file_id, (str, _oldstr))
        
        if creation_time is not None:
            assert isinstance(creation_time, datetime.datetime)

        assert isinstance(locale, (str, _oldstr))
        assert isinstance(params, OrderedDict)
        assert isinstance(parent_headers, (list, tuple))
        for h in parent_headers:
            assert isinstance(h, CELHeader)
            
        self.data_type_id = data_type_id
        self.file_id = file_id
        self.creation_time = creation_time
        self.locale = locale
        self.params = params
        self.parent_headers = parent_headers

class CELDataSet(object):
    def __init__(self, name, params, col_names, data):
        
        assert isinstance(name, (str, _oldstr))
        assert isinstance(params, OrderedDict)
        assert isinstance(col_names, list)
        for n in col_names:
            assert isinstance(n, (str, _oldstr))
        assert isinstance(data, (list, tuple))
        for d in data:
            assert isinstance(d, tuple)
        
        self.name = name
        self.params = params
        self.col_names = tuple(col_names)
        self.data = tuple(data)
        
        #logger.debug(repr(self))
        
    def __repr__(self):
        return '<CELDataSet "%s" (col_names=%s; params_hash=%d; data_hash=%d)>' \
                %(self.name, repr(self.col_names), self.params_hash, self.data_hash)
    
    @property
    def num_cols(self):
        return len(self.col_names)
    
    @property
    def params_hash(self):
        return hash(tuple(self.params.items()))
    
    @property
    def data_hash(self):
        return hash(self.data)
        
class CELDataGroup(object):
    def __init__(self, name, datasets = None):
    
        if datasets is None:
            datasets = []
    
        assert isinstance(name, (str, _oldstr))
        assert isinstance(datasets, (list, tuple))
        for ds in datasets:
            assert isinstance(ds, CELDataSet)
        
        self.name = name
        self.datasets = tuple(datasets)
        
        #logger.debug(repr(self))
        
    def __repr__(self):
        return '<CELDataGroup "%s" (datasets_hash=%d)' \
                %(self.name, self.datasets_hash)
        
    @property
    def datasets_hash(self):
        return hash(self.datasets)
    
    def get_data_set_by_name(self, name):
        for ds in self.datasets:
            if ds.name.startswith(name):
                return ds
        raise ValueError('No data set with name "%s".' %(name))
        
class CEL(object):
    def __init__(self, header, data_groups):
        self.header = header
        self.data_groups = tuple(data_groups)
    
    def get_data_group_by_name(self, name):
        for dg in self.data_groups:
            if dg.name == name:
                return dg
        raise ValueError('No data group with name "%s".' %(name))
    
    @classmethod
    def read_cel(cls, path):
        """Parser for CEL files in Command Console generic data file format.

        This is a binary format.
        """
        
        read = [0]
        
        decode_unicode = lambda s: codecs.decode(s, 'UTF-16-BE')
        decode_ascii = lambda s: codecs.decode(s, 'ascii')
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
            assert isinstance(magic_number, int) and magic_number == 59
            version_number = read_ubyte(fh)
            assert version_number == 1
            num_data_groups = read_int(fh)
            assert isinstance(num_data_groups, int)
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

            header = CELHeader(data_type_id, file_id, creation_time, locale, params, parent_headers)
            return header

        def read_col(fh):
            name = read_wstring(fh)
            valtype = read_byte(fh)
            size = read_int(fh)
            return (name, valtype, size)

        def read_data_set(fh):
            
            VALUE_TYPES = [read_byte, read_ubyte, read_short, read_ushort, read_int, read_uint,
                    read_float, read_string, read_wstring]
            
            data_pos = read_uint(fh)
            next_pos = read_uint(fh)
            data_size = next_pos - data_pos
            name = read_wstring(fh)
            n_params = read_int(fh)
            logger.debug('DataSet / data position: %d', data_pos)
            logger.debug('DataSet / next position: %d', next_pos)
            logger.debug('DataSet / data size: %d', data_size)
            logger.debug('DataSet / name: %s', name)
            logger.debug('DataSet / # parameters: %d', n_params)
            params = OrderedDict()
            for i in range(n_params):
                params.update([read_header_param(fh)])
            n_cols = read_uint(fh)
            logger.debug('DataSet / # cols: %d', n_cols)
            cols = []
            for i in range(n_cols):
                cols.append(read_col(fh))
            col_names = [c[0] for c in cols]
            n_rows = read_uint(fh)
            logger.debug('DataSet / # rows: %d', n_rows)
            data = []
            for i in range(n_rows):
                d = []
                for c in cols:
                    d.append(VALUE_TYPES[c[1]](fh))
                data.append(tuple(d))
            ds = CELDataSet(name, params, col_names, data)
            return ds, next_pos

        def read_data_group(fh):
            next_pos = read_uint(fh)
            dataset_pos = read_uint(fh)
            n_datasets = read_int(fh)
            name = read_wstring(fh)
            logger.debug('# data sets within the group: %d', n_datasets)
            logger.debug('Position of first data set within the group: %d', dataset_pos)
            logger.debug('Position of next data group: %d', next_pos)
            logger.debug('Bytes read up until this point: %d', read[0])
            assert dataset_pos == read[0]
            datasets = []
            for i in range(n_datasets):
                ds, next_pos = read_data_set(fh)
                logger.debug('Position of next data set: %d', next_pos)
                logger.debug('Bytes read up until this point: %d', read[0])
                assert (next_pos - read[0]) in [0, 1]
                if next_pos - read[0] == 1:
                    logger.warning('Skipping one byte (%d) between data sets.', ord(fh.read(1)))
                    read[0] += 1
                assert next_pos == read[0]
                datasets.append(ds)
            group = CELDataGroup(name, datasets)
            return group

        header = None
        data_groups = []
        with misc.smart_open_read(path, 'rb', try_gzip = True) as fh:
            num_data_groups, data_pos = read_file_header(fh)
            #assert n_data_groups == 1 # for expression CEL file
            header = read_data_header(fh)
            logger.info('# data groups: %d', num_data_groups)
            assert data_pos == read[0] # position of the first data group
            data_groups = []
            for i in range(num_data_groups):
                data_groups.append(read_data_group(fh))
        return cls(header, data_groups)
