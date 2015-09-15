Extending and embedding CMIstark
================================


Manually accessing the CMIstark results
---------------------------------------

.. todo:: here we need a short summary and then the references specified below.

The following example source code of Python shows how to read the curve for the
00000 state from ``<moleculename>.molecule by`` using PyTables::

  import tables
  import numpy
  stark_file = "<moleculename>.molecule"
  f=tables.openFile(stark_file)
  array = f.getNode("/_0/_0/_0/_0/_0/dcstarkenergy")
  print numpy.array(array.read())

.. comment
   Local Variables:
   coding: utf-8
   fill-column: 80
   End:
