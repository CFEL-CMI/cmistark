CMIstark user guide
===================

This program can calculate, view, and analyze the energy levels of adiabatic
Stark energy curves of linear, symmetric top and asymmetric top molecules.

The program package is developed and maintained by the Controlled Molecule
Imaging group (CMI) at the Center for Free-Electron Laser Science (CFEL),
Hamburg, Germany.

It is documented in
  Computer Physics Communications
  arXiv:1308.4076 [physics]

Usage
-----

* cmistark_calculate_energy

A (command line) script file called ``cmistark_calculate_energy`` is provided as a
driver for the calculation of Stark curves. Its command-line options are the
following::

  --help: help
  --<moleculename>: specify which molecule is used in the calculation.
  --dc-fields=: specify the range of the DC electric field (in kV/cm) by the following way: start:end:step, example: --dc-fields=0:150:151.
  --Jmax_calc=: specify the maximum value of J included in the calculation.
  --Jmax_save=: specify the maximum value of J of Stark curves saved in the output file.
  --Jmin=:      specify the minimum value of J included in the calculation.
  --Mmax=:      specify the maximum value of M included in the calculation.
  --isomer=:    specify which isomer is used, when <moleculename> have more than one isomers defined in moleculeparameter.py

Example of using ``cmistark_calculate_energy`` with options::

    cmistark_calculate_energy --isomer=0 --Jmax_calc=40 --Jmax_save=20 --3-aminophenol --dc-fields=0:150:151

After executing this command line, the script will use cmistark packages to
calculate stark energies of isomer 0 of 3-aminophenol up to J=40, and save
results up to J=20 in an output file called ``3-aminophenol.molecule``. The Stark
curve of each quantum state starts from 0 to 150kV/cm with a step of 1kV/cm.

The output file ``<moleculename>.molecule`` stores Stark curve of individual quantum
states in terms of a data format called HDF5. Such file format can be read
directly by using PyTables packages in Python. Two scripts in this program,
``cmistark_plot_energy`` and ``cmistark_print_energy`` are provided to easily access
these ``<moleculename>.molecule`` files. Their options and descriptions are provided
below.

* how to add a new molecule

Firstly, in the file called ``moleculeparameter.py`` (in ``cmistark`` folder), add all
relevant molecular constants/parameters as a function. See the examples
provided in this file.

Next, in the script file ``cmistark_calculate_energy`` (in ``script`` folder),
add an option for calling the above added function for the mew molecule.
See the examples provided in this file.

* cmistark_plot_energy 

The script file called ``cmistark_plot_energy`` can access existing Stark files
(``<moleculename>.molecule``) and plot the stored curves. Options::

  --help: help
  --energy-unit=: specify the unit of energy, options: MHz, invcm, J
  --Jmin=, --Jmax=: specify min. or max. value of J
  --Mmin=, --Mmax=: specify min. or max. value of M
  --Kamax=: specify max. value of Ka
  --states=: specify states to plot, format: "000,1010"
  --dipole: plot the effective dipole moments
  --isomer=: specify which isomer to plot

Example of using ``cmistark_plot_energy`` with options::

    cmistark_plot_energy --Jmin=0 --Jmax=2 --Mmin=1 --Mmax=1 <moleculename>.molecule


* cmistark_print energy

The script file called ``cmistark_print_energy`` can access existing Stark files
(``<moleculename>.molecule``) and print the stored curves on the screen. Options::

  --help: help
  --Jmin=, --Jmax=: specify min. or max. value of J
  --Kamax=: specify max. value of Ka
  --Mmin=, --Mmax=: specify min. or max. value of M
  --isomer=: specify which isomer to print

Example of using ``cmistark_print_energy`` with options::

    cmistark_plot_energy --Jmin=0 --Jmax=2 --Mmin=1 --Mmax=1 <moleculename>.molecule



State labels of stored Stark curves
-----------------------------------

In the output file ``<moleculename>.molecule``, each Stark curve has a state
label (:math:`J`, :math:`K_a`, :math:`K_c`, :math:`M`, isomer), which represents the adiabatic quantum number
label of the rotational state in the field, as well as the type of isomer.
:math:`J`, :math:`K_a`, :math:`K_c`, :math:`M` are integer, assuming no orbital angular momentum and spin
of electrons and nuclear spins involved. For all types of rotors, the
value of :math:`J` is not less than zero.

For asymmetric tops and linear rotors, only states with positive :math:`M` are
stored, as all curves of nonzero :math:`M` states are double degenerate. The values
of both :math:`K_a` and :math:`K_c` are not less than zero for asymmetric tops, or set
to zero for linear rotors. The state lable for linear rotors is thus
(:math:`J`, :math:`0`, :math:`0`, :math:`M`, isomer). 

For symmetric tops, states having products of :math:`K` and :math:`M` equal to :math:`+|KM|` and :math:`-|KM|`
split in the DC electric field. In the output file states corresponding with
negative :math:`|KM|` are stored with negative :math:`K` (and positive :math:`M`); this is really an
implementation detail and the sign stored with :math:`K` in this case is always the sign
of the product :math:`KM`. We note that states with :math:`K>0` and :math:`M<0` also yield :math:`-|KM|`. Thus,
all curves of nonzero :math:`M` states in the output file are also double degenerate.
Finally, the state label for prolate tops is (:math:`J`, :math:`K`, :math:`0`, :math:`M`, isomer), or
(:math:`J`, :math:`0`, :math:`K`, :math:`M`, isomer) for oblate tops.

* Structure of <moleculename>.molecule

For each state (:math:`J`, :math:`K_a`, :math:`K_c`, :math:`M`,isomer), the Stark energy as function of DC field
strength is stored by the following structures::

    /_J/_Ka/_Kc/_M/_isomer/dcfield
    /_J/_Ka/_Kc/_M/_isomer/dcstarkenergy

The following example source code of Python shows how to read the curve for the
00000 state from ``<moleculename>.molecule by`` using PyTables::

  import tables
  import numpy
  stark_file = "<moleculename>.molecule"
  f=tables.openFile(stark_file)
  array = f.getNode("/_0/_0/_0/_0/_0/dcstarkenergy")
  print numpy.array(array.read())

A script ``cmistark_print_energies``, that provides ASCII output for specified
conditions and states, is provided in the package for convenience.


* Descriptions of source code files:

Three files in lib folder provide all functions used to calculate and then
write/read Stark curves. Above script files perform the calculations by calling
these functions. The basic descriptions of each file in lib folder is provided
as following:

- ``molecule.py``: perform Stark effect calculation by calling functions from ``starkeffect.py`` and store results in an output file

- ``moleculeparameter.py``: contain all molecular parameters of individual molecules

- starkeffect.py: contain all functions, equations and algorithms required for calculating Stark effect.

(Only for stand-alone version) The descriptions of rest of files in lib
folder are also provided as following:

- ``codata.py``: store most scientific constants
- ``const.py``: call required math. and phys. constants from codata.py
- ``convert.py``: perform unit conversions
- ``hdf5.py``: read/write output files in the format of hdf5 via PyTables
- ``moleculeproperty.py``: create a molecule (as a object) from a list of atoms
- ``state.py``: create state labels and corresponding id numbers
- ``util.py``: provide array operations


.. comment
   Local Variables:
   coding: utf-8
   fill-column: 80
   End:
