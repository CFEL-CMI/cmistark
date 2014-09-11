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

.. todo:: Add links to (both, CPC and arXiv) full references in references.rst
          for the above documentation; the above lines should also contain the
          (one) title. make sure it is obvious that it is two copies of the same
          paper.

General usage
-------------

cmistark_calculate_energy
-------------------------

A (command line) script file called ``cmistark_calculate_energy`` is provided as
a driver for the calculation of Stark curves. Its command-line options are the
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
results up to J=20 in an output file called ``3-aminophenol.molecule``. The
Stark curve of each quantum state starts from 0 to 150kV/cm with a step of
1kV/cm.

The output file ``<moleculename>.molecule`` stores Stark curve of individual
quantum states in terms of a data format called HDF5. Such file format can be
read directly by using PyTables packages in Python. Two scripts in this program,
``cmistark_plot_energy`` and ``cmistark_print_energy`` are provided to easily
access these ``<moleculename>.molecule`` files. Their options and descriptions
are provided below.

cmistark_plot_energy 
--------------------

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


cmistark_print energy
---------------------

The script file called ``cmistark_print_energy`` can access existing Stark files
(``<moleculename>.molecule``) and print the stored curves on the screen. Options::

  --help: help
  --Jmin=, --Jmax=: specify min. or max. value of J
  --Kamax=: specify max. value of Ka
  --Mmin=, --Mmax=: specify min. or max. value of M
  --isomer=: specify which isomer to print

Example of using ``cmistark_print_energy`` with options::

    cmistark_plot_energy --Jmin=0 --Jmax=2 --Mmin=1 --Mmax=1 <moleculename>.molecule




.. comment
   Local Variables:
   coding: utf-8
   fill-column: 80
   End:
