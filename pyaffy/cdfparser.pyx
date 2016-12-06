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
Cython parser for Brainarray CDF files for Affymetrix GeneChip microarrays.

See: http://brainarray.mbni.med.umich.edu/Brainarray/Database/CustomCDF/genomic_curated_CDF.asp
"""

#from __future__ import (absolute_import, division,
#                        print_function, unicode_literals)
from __future__ import (absolute_import, division,
                        print_function)
#_oldstr = str
#from builtins import *
from builtins import str as text

cimport cython

from libc.stddef cimport size_t
from libc.stdlib cimport malloc
from libc.stdio  cimport FILE, fopen, fread, fclose, fgets, sscanf
from libc.string cimport strlen, strcmp
# from libc.math cimport NAN

import numpy as np
cimport numpy as np

np.import_array()

import sys
import logging
from collections import OrderedDict

logger = logging.getLogger(__name__)

cdef enum ProbeType:
    PROBES_PM = 0
    PROBES_MM = 1
    PROBES_ALL = 2

cdef read_line(char* buf, int buf_size, FILE* fp, size_t nl):
    fgets(buf, buf_size, fp)
    buf[strlen(buf) - nl] = '\0'

cdef np.uint32_t[::1] parse_probeset(char* buf, char* gene, int buf_size, size_t nl, FILE* fp, ProbeType probes, int num_rows, int num_cols):

    cdef int test
    cdef int num_pairs
    cdef int num_probes
    cdef int i, c, n
    cdef int result
    cdef int x, y

    cdef char ref_base, probe_base

    buf[0] = '\0'
    while sscanf(buf, "[Unit%*d_Block%d]", &test) <= 0:
        read_line(buf, buf_size, fp, nl)
    #assert test > 0

    read_line(buf, buf_size, fp, nl) # Name=
    result = sscanf(buf, "Name=%200s", gene)
    #assert result == 1

    read_line(buf, buf_size, fp, nl) # skip BlockNumber

    read_line(buf, buf_size, fp, nl)
    sscanf(buf, "NumAtoms=%d", &num_pairs)

    read_line(buf, buf_size, fp, nl)
    result = sscanf(buf, "NumCells=%d", &num_probes)
    assert result == 1

    if probes == PROBES_ALL:
        n = num_probes
    else:
        n = num_pairs

    cdef np.uint32_t[::1] ind = np.empty(int(n), dtype = np.uint32)

    read_line(buf, buf_size, fp, nl) # skip StartPosition
    read_line(buf, buf_size, fp, nl) # skip StopPosition
    read_line(buf, buf_size, fp, nl) # skip CellHeader

    c = 0
    for i in range(num_probes):
        read_line(buf, buf_size, fp, nl)
        result = sscanf(buf, "Cell%*d=%d %d N control %*s %*d %*d %1c %1c", &x, &y, &ref_base, &probe_base)
        assert result == 4
        if (probes == PROBES_ALL) or \
                (probes == PROBES_PM and ref_base != probe_base) or\
                (probes == PROBES_MM and ref_base == probe_base):
            ind[c] = num_rows * y + x
            c += 1

    assert c == n

    return ind


def parse_cdf(path, probe_type = 'pm', newline_chars = 2):
    """Front-end for parsing a Brainarray CDF file.

    The function assumes that the CDF file is uncompressed.

    See: http://brainarray.mbni.med.umich.edu/Brainarray/Database/CustomCDF/genomic_curated_CDF.asp

    Parameters
    ----------
    path: str
        The path of the CDF file
    probe_type: str
        The type of probes to read. Either "pm" (perfect match probes)

    Returns
    -------
    name: str
        The name of the array type.
    rows: int
        The number of rows on the array.
    cols: int
        The number of columns on the array.
    probeset: collections.OrderedDict
        The array indices for each probeset (gene). Each key/value-pair in the
        ordered dictionary corresponds to a probeset (gene). The *key*
        corresponds to the probeset (gene) ID, e.g., the Entrez ID followed by
        the string "_at". The *value* is a np.ndarray of type np.uint32 that
        contains the probe indices for the order in which probe intensities are
        stored in CEL files (column-first order).
    """
    assert isinstance(path, (text, str))
    assert isinstance(probe_type, (text, str))
    assert isinstance(newline_chars, int)

    cdef int buf_size = 1000
    cdef size_t nl = <size_t>newline_chars
    cdef char* buf = <char*>malloc(buf_size)
    cdef char name[201]
    cdef FILE* fp

    cdef int num_rows, num_cols
    cdef int num_probesets
    cdef char gene[201]
    cdef int i, x, y

    cdef np.uint32_t[::1] ind

    cdef ProbeType probes
    if probe_type == 'all':
        probes = PROBES_ALL
    elif probe_type == 'mm':
        probes = PROBES_MM
    elif probe_type == 'pm':
        probes = PROBES_PM
    else:
        logger.warning(
            ('Unknown probe type "%s" (should be "pm", "mm" or "all"). '
             'Will default to "pm".' %(probe_type))
        )
        probes = PROBES_PM

    probesets = OrderedDict()

    path_bytes = path.encode('UTF-8')
    cdef char* c_path_bytes = path_bytes

    fp = fopen(c_path_bytes, 'r')

    try:
        # parsing starts here

        # count rows to get to relevant information in [Chip] section
        # then count QC sections (each while?)

        # check header
        read_line(buf, buf_size, fp, nl) # [CDF]
        assert strcmp(buf, "[CDF]") == 0
        read_line(buf, buf_size, fp, nl) # Version=GC3.0
        assert strcmp(buf, "Version=GC3.0") == 0

        read_line(buf, buf_size, fp, nl) # empty

        # read general information
        read_line(buf, buf_size, fp, nl) # [Chip]
        assert strcmp(buf, "[Chip]") == 0

        read_line(buf, buf_size, fp, nl) # Name=
        sscanf(buf, "Name=%200s", name)
        
        read_line(buf, buf_size, fp, nl) # Rows=
        sscanf(buf, "Rows=%d", &num_rows)

        read_line(buf, buf_size, fp, nl) # Cols=
        sscanf(buf, "Cols=%d", &num_cols)

        read_line(buf, buf_size, fp, nl) # NumberOfUnits=
        sscanf(buf, "NumberOfUnits=%d", &num_probesets)

        test = []
        for i in range(num_probesets):
            #print i,
            #sys.stdout.flush()
            ind = parse_probeset(buf, gene, buf_size, nl, fp, probes, num_rows, num_cols)
            probesets[str(gene)] = np.uint32(ind)

    finally:
        fclose(fp)

    #print "header information:"
    #print "name = %s" %(str(name))
    #print "rows = %d" %(int(num_rows))
    #print "cols = %d" %(int(num_cols))
    #print "number of probesets = %d" %(int(num_probesets))

    return text(name), int(num_rows), int(num_cols), probesets
