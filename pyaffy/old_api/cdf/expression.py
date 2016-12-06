"""Python API for data in Affymetrix CDF files for expression microarrays.

"""

from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
_oldstr = str
from builtins import *

from collections import OrderedDict
import logging
import itertools as it

from configparser import ConfigParser

logger = logging.getLogger(__name__)

def pairwise(iterable):
    "s -> (s0, s1), (s2, s3), (s4, s5), ..."
    a = iter(iterable)
    return it.izip(a, a)

class Probe(object):
    def __init__(self, coords, index):
        
        assert isinstance(coords, (list, tuple))
        assert len(coords) == 2
        assert isinstance(coords[0], int) and isinstance(coords[1], int)
        assert isinstance(index, int)
        
        self.coords = coords
        self.index = index

class ProbePair(object):
    """Represents a pair of perfect match (PM) and mismatch (MM) probes."""

    def __init__(self, pm_probe, mm_probe):
        assert isinstance(pm_probe, Probe)
        assert isinstance(mm_probe, Probe)
        self.pm_probe = pm_probe
        self.mm_probe = mm_probe

class QCProbeSet(object):
    """QC probe set of an Affymetrix expression microarray."""

    def __init__(self, id_, type_, probes = None):
        
        if probes is None:
            probes = []
        
        assert isinstance(id_, int)
        assert isinstance(type_, int)
        assert isinstance(probes, (list, tuple))
        
        for p in probes:
            assert isinstance(p, Probe)
        
        self.id = id_
        self.type = type_
        self.probes = probes

class ExpProbeSet(object):
    """Non-QC probe set of an Affymetrix expression microarray."""
    def __init__(self, id_, gene_id, probe_pairs = None):
        
        if probe_pairs is None:
            probe_pairs = []

        assert isinstance(id_, int)
        assert isinstance(gene_id, (str, _oldstr))
        assert isinstance(probe_pairs, (list, tuple))
        
        for pp in probe_pairs:
            assert isinstance(pp, ProbePair)
            
        self.id = id_
        self.gene_id = gene_id
        self.probe_pairs = probe_pairs

class ExpCDF(object):
    """CDF data for an Affymetrix expression microarray."""
    
    def __init__(self, name, num_rows, num_cols, qc_probesets = None, exp_probesets = None):       
        
        if qc_probesets is None:
            qc_probesets = []
            
        if exp_probesets is None:
            exp_probesets = []

        assert isinstance(qc_probesets, (list, tuple))
        assert isinstance(exp_probesets, (list, tuple))
            
        assert isinstance(name, (str, _oldstr))
        assert isinstance(num_rows, int)
        assert isinstance(num_cols, int)

        for ps in qc_probesets:
            assert isinstance(ps, QCProbeSet)
            
        for ps in exp_probesets:
            assert isinstance(ps, ExpProbeSet)

        self.name = name
        self.num_rows = num_rows
        self.num_cols = num_cols
        
        self.qc_probesets = qc_probesets
        self.exp_probesets = exp_probesets
        
    @classmethod
    def read_cdf(cls, path):
        """Parser for CDF data file format version "GC3.0"."""

        # we're hijacking the configparser package to parse a CDF file
        # (the format is compatible with configparser's "ini-style" format)

        C = ConfigParser(interpolation = None, delimiters = ('=',),
                empty_lines_in_values = False)
        # do not convert variable names to lower-case
        C.optionxform = lambda x: x
        
        # parse CDF file using configparser.ConfigParser 
        with open(path, 'rb') as fh:
            C.read_file(fh)
            
        # get header information
        name = C['Chip']['Name']
        num_rows = int(C['Chip']['Rows'])
        num_cols = int(C['Chip']['Cols'])
        num_exp_probesets = int(C['Chip']['NumberOfUnits'])
        num_qc_probesets = int(C['Chip']['NumQCUnits'])

        # get QC and expression probe sets
        qc_probesets = []
        exp_probesets = []
        c = 0
        for k, sec in C.items():
            if k.startswith('QC'):
                # QC probe set
                id_ = int(k[2:])
                type_ = int(sec['Type'])
                num_probes = int(sec['NumberCells'])
                probes = []
                for tag, valstr in sec.items():
                    #print tag, valstr
                    if (not tag.startswith('Cell')) or tag == 'CellHeader':
                        continue
                    
                    # split on tab
                    val = valstr.split('\t')
                    coords = (int(val[0]), int(val[1]))
                    index = int(val[5])
                    probes.append(Probe(coords, index))
                    
                try:
                    assert len(probes) == num_probes
                except AssertionError:
                    print id_, type_, len(probes), num_probes
                    raise
                qc_probesets.append(QCProbeSet(id_, type_, probes))
                
            elif k.startswith('Unit'):
                # regular (expression) probe set
                if ('_' not in k):
                    # Unit section => only "metadata"
                    id_ = int(sec['UnitNumber'])
                    assert int(sec['UnitType']) == 3
                    assert int(sec['NumberBlocks']) == 1

                else:
                    c += 1
                    # Unit_Block section (= probe set)
                    id_ = int(k[4:k.index('_')])
                    gene_id = sec['Name']
                    num_probe_pairs = int(sec['NumAtoms'])
                    probe_pairs = []
                    iterator = sec.items()
                    for tag, valstr in iterator:
                        if tag == 'CellHeader':
                            break
                    for (tag1, valstr1), (tag2, valstr2) in \
                            pairwise(iterator):
                        #if c == 1:
                        #    logger.debug(tag + ': ' + valstr)
                        #if (not tag.startswith('Cell')) or tag == 'CellHeader':
                        #    continue

                        # split on tab
                        val1 = valstr1.split('\t')
                        val2 = valstr2.split('\t')

                        # (mismatch probes are those that have the same base
                        #  as the target, instead of the complementary base)
                        #assert val1[8] == val1[9] # mismatch probe
                        #assert val2[8] != val2[9] # perfect match probe

                        pm = None
                        mm = None

                        if val1[8] == val1[9]: # mismatch probe
                            mm = val1
                        else:
                            pm = val1

                        if val2[8] == val2[9]: # mismatch probe
                            mm = val2
                        else:
                            pm = val2

                        assert pm is not None and mm is not None
                        assert pm is not mm

                        mm_coords = (int(mm[0]), int(mm[1]))
                        #mm_index = int(mm[11])
                        mm_index = mm_coords[1] * num_rows + mm_coords[0]
                        mm_probe = Probe(mm_coords, mm_index)

                        pm_coords = (int(pm[0]), int(pm[1]))
                        #pm_index = int(pm[11])
                        pm_index = pm_coords[1] * num_rows + pm_coords[0]
                        assert pm_index == int(pm[11])
                        pm_probe = Probe(pm_coords, pm_index)

                        probe_pairs.append(ProbePair(pm_probe, mm_probe))
                        #probes.append(Probe(pm_index, pm_coords,
                        #        mm_index, mm_coords))
                        
                    try:
                        assert len(probe_pairs) == num_probe_pairs
                    except AssertionError:
                        print id_, gene_id, len(probe_pairs), num_probe_pairs
                        raise

                    exp_probesets.append(ExpProbeSet(id_, gene_id,
                            probe_pairs))

        assert len(qc_probesets) == num_qc_probesets
        assert len(exp_probesets) == num_exp_probesets
                
        return cls(name, num_rows, num_cols, qc_probesets, exp_probesets)
