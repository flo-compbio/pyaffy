pyAffy
======

.. "|docs-latest| |docs-develop|

pyAffy is a Python/Cython implementation of the RMA algorithm for
processing raw data from Affymetrix expression microarrays.

Installation
------------

.. code-block:: bash

    $ git clone ...
    $ cd pyaffy
	$ pip install -e .


Usage
-----

The `rma` function expects two parameters: A CDF file and an ordered
dictionary (of type `collections.OrderedDict`) with sample names as keys
and corresponding CEL files as values. It returns a list of genes, a
list of samples, and an expression matrix (of type `numpy.ndarray`), in that
order.

.. code-block:: python

    from pyaffy import rma
    # for documentation of the rma function, try:
    # help(rma)
    genes, samples, X = rma(cdf_file, sample_cel_files)

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
