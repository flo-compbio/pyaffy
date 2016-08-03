..
    Copyright (c) 2016 Florian Wagner
    
    This file is part of pyAffy.
    
    pyAffy is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License, Version 3,
    as published by the Free Software Foundation.
    
    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.
    
    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.

pyAffy
======

.. "|docs-latest| |docs-develop|

pyAffy is a Python/Cython implementation of the RMA algorithm for
processing raw data from Affymetrix expression microarrays. For a detailed
discussion of this implementation, see the `pyAffy PeerJ preprint`__. For
a list of changes, see the `changelog <changelog.rst>`_.

__ peerj_preprint_

.. _peerj_preprint: https://peerj.com/preprints/1790/

Installation
------------

Option 1: Using `pip`
~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

    $ pip install pyaffy


Option 2: Cloning the GitHub repository
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

    $ git clone https://github.com/flo-compbio/pyaffy.git
    $ cd pyaffy
    $ pip install -e .

Usage
-----

The `rma` function expects two parameters: A custom CDF file (from the
`Brainarray web site`__) and an ordered dictionary (`collections.OrderedDict`)
with sample names as keys and corresponding CEL files as values.

__ brainarray_

The `rma` function returns a list of genes, a list of samples, and an
expression matrix (of type `numpy.ndarray`), in that order.

.. code-block:: python

    from pyaffy import rma
    # for documentation of the rma function, try:
    # help(rma)
    genes, samples, X = rma(cdf_file, sample_cel_files)

A small `example with real code`__ is available in the `pyaffy-demos` repository.

__ real_example_

.. _brainarray: http://brainarray.mbni.med.umich.edu/Brainarray/Database/CustomCDF/genomic_curated_CDF.asp
.. _real_example: https://github.com/flo-compbio/pyaffy-demos/tree/master/minimal

Copyright and License
---------------------

Copyright (c) 2016 Florian Wagner

::

  pyAffy is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License, Version 3,
  as published by the Free Software Foundation.
  
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program. If not, see <http://www.gnu.org/licenses/>.
