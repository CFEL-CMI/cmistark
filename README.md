# CMIstark

This program can calculate, view, and analyze the energy levels of adiabatic
Stark energy curves of linear, symmetric top and asymmetric top molecules.

The program package is developed and maintained by the Controlled Molecule
Imaging group (CMI) at the Center for Free-Electron Laser Science (CFEL),
Hamburg, Germany.

The program (version 1.0) is documented in _"CMIstark: Python package for the
Stark-effect calculation and symmetry classification of linear, symmetric and
asymmetric top wavefunctions in dc electric fields"_ by Yuan-Pin Chang, Frank
Filsinger, Boris Sartakov, and Jochen Küpper, published as [_Computer Physics
Communications_ **185**, 339-349 (2014), DOI:
10.1016/j.cpc.2013.09.001](https://dx.doi.org/10.1016/j.cpc.2013.09.001);
preprint available at [arXiv:1308.4076
[physics]](http://arxiv.org/abs/1308.4076)


## Prerequisities:

Computer: Any Macintosh, PC, or Linux/UNIX workstations with a modern Python
distribution.

The following external Python packages are also required:
 - cmiext
 - matplotlib
 - numpy
 - scipy
 - sympy
 - tables


## Installation:

Installation is performed by executing the generic Python install command:
```shell
python setup.py install
```
in the unpacked source code directory.

In order to install this extension module in user-space, [set up your
environment](https://docs.python.org/3/using/cmdline.html#envvar-PYTHONUSERBASE)
for python to find it, e.g.,
```shell
setenv PYTHONUSERBASE=$HOME/.local
```
and run the install command
```shell
python setup.py install --user
```


## Usage:

### cmistark_calculate_energy

A (command line) script file called cmistark_calculate_energy is provided as a
driver for the calculation of Stark curves. Its command-line options are the
following:
```plain
--help:           this help
--<moleculename>: specify which molecule is used in the calculation.
--dc-fields=:     specify the range of the DC electric field (in kV/cm) by the
                  following way: start:end:step, example: --dc-fields=0:150:151.
--Jmax_calc=:     specify the maximum value of J included in the calculation.
--Jmax_save=:     specify the maximum value of J of Stark curves saved in the output file.
--Jmin=:          specify the minimum value of J included in the calculation.
--Mmax=:          specify the maximum value of M included in the calculation. [TM: Default value?]
--isomer=:        specify which isomer is used, when <moleculename> have more than
                  one isomers defined in moleculeparameter.py
```
Example of using cmistark_calculate_energy with options:
```shell
cmistark_calculate_energy --isomer=0 --Jmax_calc=40 --Jmax_save=20 --3-aminophenol --dc-fields=0:150:151
```

After executing this command line, the script will use cmistark packages to
calculate stark energies of isomer 0 of 3-aminophenol up to _J_ = 40, and save
results up to _J_ = 20 in an output file called `3-aminophenol.molecule`. The
Stark curve of each quantum state starts from 0 to 150 kV/cm with a step of
1 kV/cm.

The output file `<moleculename>.molecule` stores Stark curve of individual
quantum states in a data format called HDF5. Such a file format can be read
directly by using PyTables packages in Python. Two scripts in this program,
`cmistark_plot_energy` and `cmistark_print_energy` are provided to easily access
these `<moleculename>.molecule` files. Their options and descriptions are provided
below.


### cmistark_plot_energy

The script file called `cmistark_plot_energy` can access existing Stark files
(`<moleculename>.molecule`) and plot the stored curves.

It offers the following options:
```plain
--help: help
--energy-unit=: specify the unit of energy, options: MHz, invcm, J
--Jmin=, --Jmax=: specify min. or max. value of J
--Mmin=, --Mmax=: specify min. or max. value of M
--Kamax=: specify max. value of Ka
--states=: specify states to plot, format: "0000,1010"
--dipole: plot the effective dipole moments
--isomer=: specify which isomer to plot
```
example of using `cmistark_plot_energy` with options:
```shell
    cmistark_plot_energy --Jmin=0 --Jmax=2 --Mmin=1 --Mmax=1 <moleculename>.molecule
```


### cmistark_print energy 

The script file called `cmistark_print_energy` can access existing Stark files
(`<moleculename>.molecule`) and print the stored curves on the screen.

It offers the following options:
```plain
--help: help
--Jmin=, --Jmax=: specify min. or max. value of J
--Kamax=: specify max. value of Ka
--Mmin=, --Mmax=: specify min. or max. value of M
--isomer=: specify which isomer to print
```
Example of using acmistark_print_energya with options:
```
    cmistark_print_energy --Jmin=0 --Jmax=2 --Mmin=1 --Mmax=1 <moleculename>.molecule
```


## How to add a new molecule

Firstly, in the file called `moleculeparameter.py` (in the `cmistark/` folder),
add all relevant molecular `constants/parameters` as a function. See the
examples provided in this file.

Next, in the script file `cmistark_calculate_energy` (in script folder), add an
option for calling the above added function for the mew molecule. See the
examples provided in this file.


## State labels of stored Stark curves

In the output file <moleculename>.molecule, each Stark curve has a state label
(_J,Ka,Kc,M_,isomer), which represents the adiabatic quantum number label of the
rotational state in the field, as well as the type of isomer. _J, Ka, Kc, M_ are
integer, assuming no orbital angular momentum and spin of electrons and nuclear
spins involved. For all types of rotors, the value of _J_ is not less than zero.

For asymmetric tops and linear rotors, only states with positive _M_ are stored,
as all curves of nonzero _M_ states are doubly degenerate. The values of both
_Ka_ and _Kc_ are not less than zero for symmetric and asymmetric tops, one is
kept at zero for symmetric tops, and both are set to zero for linear rotors. The
state lable for linear rotors is thus (_J_,0,0,_M_,isomer).

For symmetric tops, states having products of _K_ and _M_ equal to +|_KM_| and
-|_KM_| are split in the DC electric field. In the output file states
corresponding to negative |_KM_| are stored with negative _K_ (and positive
_M_); this is really an implementation detail and the sign stored with _K_ in
this case is always the sign of the product _KM_. We note that states with
_K_ > 0 and _M_ < 0 also yield -|_KM_|. Thus, all curves of nonzero _M_ states
in the output file are also double degenerate. Finally, the state label for
prolate tops is (_J,K,0,M_,isomer), or (_J,0,K,M_,isomer) for oblate tops.



## Structure of data storage files


The storage files `<moleculename>.molecule` are HDF5 files in which for every
state (_J,Ka,Kc,M_,isomer) the Stark energy as a function of the DC field
strength is stored in the following structure:

```plain
    /_J/_Ka/_Kc/_M/_isomer/dcfield
    /_J/_Ka/_Kc/_M/_isomer/dcstarkenergy
```

The following example source code of Python shows how to read the curve for the
`00000` state from `<moleculename>.molecule` using PyTables:

```python
import tables
import numpy
stark_file = "<moleculename>.molecule"
f = tables.openFile(stark_file)
array = f.getNode("/_0/_0/_0/_0/_0/dcstarkenergy")
print(numpy.array(array.read()))
```

A script `cmistark_print_energies` that provides ASCII output for specified
conditions and states is provided in the package for convenience.



## Descriptions of source code files

Three files in the lib folder provide all functions used to calculate and then
write/read Stark curves. The above script files perform the calculations by
calling these functions. The basic descriptions of each file in lib folder as
follows:

**`molecule.py`**

Performs Stark effect calculation by calling functions from `starkeffect.py` and
store results in an output file

**`moleculeparameter.py`**

Contains all molecular parameters of individual molecules

**`starkeffect.py`**

Contains all functions, equations and algorithms required for calculating Stark
effect.



<!-- Put Emacs local variables into HTML comment
Local Variables:
coding: utf-8
fill-column: 80
End:
-->
