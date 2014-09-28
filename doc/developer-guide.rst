CMIstark developer guide
===========================

Here are a few simple guidelines to please be obeyed when working on CMIdiffract

* Document your code!

  * Use spinx-compatible docstrings to document all classs, methods, functions, etc.

* Write code that is compatible with the latest stable Python 3.x version.
* Make use of NumPy as much as possible.


Source code formatting
----------------------

* CMIstark uses the 4-spaces standard for indentation of blocks.
* Do not use tabs, always expand to spaces.
* Try to not extend lines beyond 100 characters
* Keep the utf-8 coding directive in the first line
* Keep the Emacs local variables section at the end of all files, and try to stick to the directives
  (manually) when not using Emacs.


Version control (git) details
-----------------------------

* CMIstark uses git (Stash) as a version control system with a central repositories in the
  DESY-CFEL-CMI stash server.

  * CMIstark uses the git-flow branching model

    * the principal development branch is ``develop``
    * all new developments should be done on a ``feature/`` branch and, once ready, be branched into
      ``develop`` 

  * never touch the branch ``master`` -- this is to be done by the maintainers.

    * the ``master`` branch is only for releases, which for CMIstark are rare and surely only after
      discussion with me. There should never be any development done on ``master``, nor any release
      preparations. The latter is done on ``release/``, then the release is put onto ``master``, and
      possibly necessary fixes are done on ``hotfix/``!

  * Do not repeatedly branch feature branches into ``develop`` instead merge ``develop`` into your
    ``feature/`` branch.
  * General documentation work should always be made on ``develop`` (only)!

    * commit such doc-only updates as separate commits!
    * one can then merge these doc-only commits into ``feature/`` branches

  * **never implement a change twice** manually! Implement it on the most appropriate branch, then merge
    it into whatever branch you want to have it!

* make sure the git repo stays clean and one can follow it with our turning crazy!


State labels of stored Stark curves
-----------------------------------

In the output file ``<moleculename>.molecule``, each Stark curve has a state
label (:math:`J`, :math:`K_a`, :math:`K_c`, :math:`M`, isomer), which represents
the adiabatic quantum number label of the rotational state in the field, as well
as the type of isomer. :math:`J`, :math:`K_a`, :math:`K_c`, :math:`M` are
integers, assuming no orbital angular momentum and spin of electrons and nuclear
spins involved. For all types of rotors, the value of :math:`J` is not less than
zero.

For asymmetric tops and linear rotors, only states with positive :math:`M` are
stored, as all curves of nonzero :math:`M` states are doubly degenerate. The
values of both :math:`K_a` and :math:`K_c` are not less than zero for asymmetric
tops, or set to zero for linear rotors. The state lable for linear rotors is
thus (:math:`J`, :math:`0`, :math:`0`, :math:`M`, isomer).

For symmetric tops, states having products of :math:`K` and :math:`M` equal to
:math:`+|KM|` and :math:`-|KM|` split in the DC electric field. In the output
file states, corresponding to negative :math:`|KM|` are stored with negative
:math:`K` (and positive :math:`M`); this is really an implementation detail and
the sign stored with :math:`K` in this case is always the sign of the product
:math:`KM`. We note that states with :math:`K>0` and :math:`M<0` also yield
:math:`-|KM|`. Thus, all curves of nonzero :math:`M` states in the output file
are also doubly degenerate. Finally, the state label for prolate tops is
(:math:`J`, :math:`K`, :math:`0`, :math:`M`, isomer), and (:math:`J`, :math:`0`,
:math:`K`, :math:`M`, isomer) for oblate tops.


Structure of <moleculename>.molecule
------------------------------------

For each state (:math:`J`, :math:`K_a`, :math:`K_c`, :math:`M`,isomer), the
Stark energy as function of DC field strength is stored in the following
structures::

    /_J/_Ka/_Kc/_M/_isomer/dcfield
    /_J/_Ka/_Kc/_M/_isomer/dcstarkenergy


.. todo:: provide/describe relation to ``State``-labels

The following example source code of Python shows how to read the curve for the
00000 state from ``<moleculename>.molecule`` by using PyTables::

  import tables
  import numpy
  stark_file = "<moleculename>.molecule"
  f=tables.openFile(stark_file)
  array = f.getNode("/_0/_0/_0/_0/_0/dcstarkenergy")
  print numpy.array(array.read())

A script ``cmistark_print_energies``, that provides ASCII output for specified
conditions and states, is provided in the package for convenience.


Descriptions of source code files
---------------------------------

Three source code files in ``cmistark`` folder provide all functions used to calculate and then
write/read Stark curves. The above script files perform the calculations by calling
these functions. The basic descriptions of each file in lib folder are
as follows:

- ``molecule.py``: perform the Stark effect calculation by calling functions from ``starkeffect.py`` and store results in an output file

- ``moleculeparameter.py``: contain all molecular parameters of individual molecules

- ``starkeffect.py``: contain all functions, equations and algorithms required for calculating the Stark effect.

(Only for stand-alone version) The descriptions of the rest of the files in the lib
folder are as follows:

- ``codata.py``: store most scientific constants
- ``const.py``: call required math. and phys. constants from codata.py
- ``convert.py``: perform unit conversions
- ``hdf5.py``: read/write output files in the format of hdf5 via PyTables
- ``moleculeproperty.py``: create a molecule (as an object) from a list of atoms
- ``state.py``: create state labels and corresponding id numbers
- ``util.py``: provide array operations


.. comment
   Local Variables:
   coding: utf-8
   fill-column: 100
   truncate-lines: t
   End:
