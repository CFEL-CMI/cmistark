CMIstark user guide
===================

This program can calculate, view, and analyze the energy levels of adiabatic
Stark energy curves of linear, symmetric top and asymmetric top molecules.

The program package is developed and maintained by the Controlled Molecule
Imaging group (CMI) at the Center for Free-Electron Laser Science (CFEL),
Hamburg, Germany.

It is documented in

* Computer Physics Communications [Chang2014]_
* arXiv:1308.4076 [physics] (eprint of the above CPC paper) [Chang2014arxiv]_


General usage
-------------

The following provides examples of the general usage.

.. todo:: Yuan-Pin, you misunderstood my earlier todos. You should first write
          about the general usagem adn then provide a specific example command
          at the end of the section of each tool. Put together, these shall give
          a working example. At every stage, you should provide info on the
          actual output. It will surely be very useful to even provide the
          output graphics of ``plot``.

Calculate the Stark energies of water::

    cmistark_calculate_energy --isomer=0 --Jmax_calc=10 --Jmax_save=2 --water --dc-fields=0:150:16

After the calculation finishes, it yields a Stark energy file ``water.molecule``.

Plot the Stark energy file ``water.molecule``::

    cmistark_plot_energy --Jmax=0 water.molecule

A plot of the Stark energy of J=0 state for water will be created.

Print the Stark energy from the file ``water.molecule``::

    cmistark_print_energy --Jmax=0 water.molecule

And the print result is::

    # state: 0 0 0 0 0
    0.0 0.0 0.0
    10.0 -26.1645127287 5.23264284222
    20.0 -104.652856844 10.4642471557
    30.0 -235.449455843 15.6937755522
    40.0 -418.528367889 20.9201929214
    50.0 -653.853314271 26.1424675625
    60.0 -941.377719139 31.3595723071
    70.0 -1281.04476041 36.57048563
    80.0 -1672.78743174 41.7741927455
    90.0 -2116.52861532 46.9696866871
    100.0 -2612.18116548 52.1559693657
    110.0 -3159.64800264 57.332052607
    120.0 -3758.82221762 62.496959162
    130.0 -4409.58718588 67.6497236918
    140.0 -5111.81669146 72.7893937211
    150.0 -5865.3750603 72.7893937211
    Closing remaining open files:water.molecule...done


cmistark_calculate_energy
-------------------------

A (command line) script file called ``cmistark_calculate_energy`` is provided as
a driver for the calculation of Stark curves. Its command-line options are the
following::

  --help: help
  --<moleculename>: specify which molecule is used in the calculation.
  --dc-fields=: specify the range of the DC electric field (in kV/cm) by the following way: start:end:step, example: --dc-fields=0:150:151.
  --Jmax_calc=: specify the maximum value of J included in the calculation.
  --Jmax_save=: specify the maximum value of J of the Stark curves saved in the output file.
  --Jmin=:      specify the minimum value of J included in the calculation.
  --Mmax=:      specify the maximum value of M included in the calculation. [TM: What is the default?  Jmax_calc?]
  --isomer=:    specify which isomer is used, when <moleculename> has more than one isomer defined in moleculeparameter.py

Example of using ``cmistark_calculate_energy`` with options::

    cmistark_calculate_energy --isomer=0 --Jmax_calc=40 --Jmax_save=20 --3-aminophenol --dc-fields=0:150:151

After executing this command line, the script will use cmistark packages to
calculate stark energies of isomer 0 of 3-aminophenol up to J=40, and save
results up to J=20 in an output file called ``3-aminophenol.molecule``. The
Stark curve of each quantum state starts from 0 to 150kV/cm with a step of
1kV/cm.

The output file ``<moleculename>.molecule`` stores Stark curves of individual
quantum states a data format called HDF5. Such a file format can be read
directly by using PyTables packages in Python. Two scripts in this program,
``cmistark_plot_energy`` and ``cmistark_print_energy`` are provided to easily
access these ``<moleculename>.molecule`` files. Their options and descriptions
are provided below.


cmistark_print_energy
---------------------

The script file called ``cmistark_print_energy`` can access existing Stark files
(``<moleculename>.molecule``) and print the stored curves on the screen.
Options::

  --help: help
  --Jmin=, --Jmax=: specify min. or max. value of J
  --Kamax=: specify max. value of Ka
  --Mmin=, --Mmax=: specify min. or max. value of M
  --isomer=: specify which isomer to print

Example of using ``cmistark_print_energy`` with options::

    cmistark_print_energy --Jmin=0 --Jmax=2 --Mmin=1 --Mmax=1 <moleculename>.molecule


cmistark_plot_energy 
--------------------

The script file called ``cmistark_plot_energy`` can access existing Stark files
(``<moleculename>.molecule``) and plot the stored curves. Options::

  --help: help
  --energy-unit=: specify the unit of energy, options: MHz, invcm, J
  --Jmin=, --Jmax=: specify min. or max. value of J
  --Mmin=, --Mmax=: specify min. or max. value of M
  --Kamax=: specify max. value of Ka
  --states=: specify states to plot, format: J K_a K_c or J K_a K_c M, when M is not specified, 
             all M levels of the J state are plotted. Example: "000,1010"
  --dipole: plot the effective dipole moments
  --isomer=: specify which isomer to plot

Example of using ``cmistark_plot_energy`` with options::

    cmistark_plot_energy --Jmin=0 --Jmax=2 --Mmin=1 --Mmax=1 <moleculename>.molecule



.. comment
   Local Variables:
   coding: utf-8
   fill-column: 80
   End:
