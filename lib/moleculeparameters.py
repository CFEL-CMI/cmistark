#!/usr/bin/env python
# -*- coding: utf-8; fill-column: 120 -*-
#
# This file is part of JK Python extensions
# Copyright (C) 2008,2009,2012 Jochen Küpper <software@jochen-kuepper.de>
#
# This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public
# License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# If you use this programm for scientific work, you must correctly reference it; see LICENSE file for details.
#
# This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
# warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with this program. If not, see
# <http://www.gnu.org/licenses/>.
from __future__ import division
__author__ = "Jochen Küpper <jochen.kuepper@cfel.de>"
__doc__ = """This modules implements the molecular parameters of all our candidates.

The relevant parameters are
- rotational constants
- centrifugal distortion constants
  - for a linear top: D
  - for a symmetric top: D_J, D_{JK}, D_K
  - for an asymmetric top
    - in Watson's A reduction: \Detlta_J, \Detlta_{JK}, \Detlta_K, d_J, d_K
- dipole moments
  - for linear and symmetric tops: \mu
  - for an asymmetric top \mu_a, \mu_b, \mu_c
- symmetry
"""

import numpy as num
import getopt, sys

import jkext.convert as convert
import jkstark.molecule as molecule
import jkstark.starkeffect as starkeffect
from jkext.state import State
from jkext.molecule import Masses



def three_aminophenol(param):
    """Molecular parameters for 3-aminophenol

    Implemented isomers are
     0  -  cis conformer, experimental values from F. Filsinger et al., PCCP 10, 666 (2008)
     1  -  trans conformer, exp values, F. Filsinger et al., PCCP 10, 666 (2008)
     2  -  cis conformer, MP2/aug-cc-pVTZ calculation using <Gaussian 2003.1> by Daniel Rösch, Basel, 2011
     3  -  trans conformer, calculated values of MP2/aug-cc-pVTZ method from Daniel Rösch in Basel, 2011
     4  -  cis conformer, calculated values of B3LYP/aug-cc-pVTZ method from Daniel Rösch in Basel, 2011
     5  -  trans conformer, calculated values of B3LYP/aug-cc-pVTZ method from Daniel Rösch in Basel, 2011
    """
    param.name = "3-aminophenol"
    param.watson = 'A'
    param.symmetry = 'N'
    if param.isomer == 0: # cis, Filsinger et al. PCCP ...
        param.rotcon = convert.Hz2J(num.array([3734.93e6, 1823.2095e6, 1226.493e6]))
        param.dipole = convert.D2Cm(num.array([1.7718, 1.517, 0.]))
    elif param.isomer == 1:
        param.rotcon = convert.Hz2J(num.array([3730.1676e6, 1828.25774e6, 1228.1948e6]))
        param.dipole = convert.D2Cm(num.array([0.5563, 0.5375, 0.]))
    elif param.isomer == 2:
        param.rotcon = convert.Hz2J(num.array([3748.0923e6, 1824.5812e6, 1228.7585e6]))
        param.dipole = convert.D2Cm(num.array([1.793, 1.4396, 0.]))
    elif param.isomer == 3:
        param.rotcon = convert.Hz2J(num.array([3736.8454e6, 1831.7399e6, 1230.7259e6]))
        param.dipole = convert.D2Cm(num.array([0.3953, 0.8203, 0.]))
    elif param.isomer == 4:
        param.rotcon = convert.Hz2J(num.array([3755.0444e6, 1828.9366e6, 1231.0926e6]))
        param.dipole = convert.D2Cm(num.array([1.8575, 1.6484, 0.]))
    elif param.isomer == 5:
        param.rotcon = convert.Hz2J(num.array([3752.3419e6, 1833.1737e6, 1232.6659e6]))
        param.dipole = convert.D2Cm(num.array([0.5705, 0.4771, 0.]))



def indole(param):
    """Molecular parameters for indole

    Implemented isomers are
    0  -  experimental values from Kang, Korter, Pratt, J. Chem. Phys. 122, 174301 (2005)
    1  -  experimental inertial constants from W. Caminati and S. Dibernardo, J. Mol. Struct. 240, 253 (1990)
          and dipole moment from Kang, Korter, Pratt, J. Chem. Phys. 122, 174301 (2005) for dipole moment.
    """
    param.name = "indole"
    param.mass = 8 * Masses['C'] + Masses['N'] + 7 * Masses['H']
    param.watson = 'A'
    param.symmetry = 'N'
    if 0 == param.isomer:
        param.rotcon = convert.Hz2J(num.array([3877.9e6, 1636.1e6, 1150.9e6]))
        param.dipole = convert.D2Cm(num.array([1.376, 1.400, 0.]))
    elif 1 == param.isomer:
        param.rotcon = convert.Hz2J(num.array([3877.826e6, 1636.047e6, 1150.8997e6]))
        param.quartic = convert.Hz2J(num.array([0.0352e3, 0.042e3, 0.16e3, 0.1005e3, 0.128e3]))
        param.dipole = convert.D2Cm(num.array([1.376, 1.400, 0.]))



def indole_water1(param):
    """Molecular parameters for the indole-water complex

    Implemented isomers are
    0  -  experimental inertial parameters from Korter, Pratt, Kuepper, J. Phys. Chem. A 102, 7211 (1998)
          and experimental dipole moment from C. Kang, T. M. Korter, and D. W. Pratt, J. Chem. Phys. 122, 174301 (2005)
    1  -  experimental inertial parameters from Blanco S et al, J. Chem. Phys., Vol. 119, 880 (2003)
          and experimental dipole moment from C. Kang, T. M. Korter, and D. W. Pratt, J. Chem. Phys. 122, 174301 (2005)
    """
    param.name = "indole-water"
    param.mass = 8 * Masses['C'] + Masses['N'] + Masses['O'] + 9 * Masses['H']
    param.watson = 'A'
    param.symmetry = 'N'
    if 0 == param.isomer:
        param.rotcon = convert.Hz2J(num.array([2062.5e6, 945.1e6, 649.3e6]))
        param.quartic = convert.Hz2J(num.array([0.0011e6, -0.006e6, 0.014e6, 0.0005e6, 0.001e6]))
        param.dipole = convert.D2Cm(num.array([4.2, 1.2, 0.]))
    elif 1 == param.isomer:
        param.rotcon = convert.Hz2J(num.array([2064.3954e6, 945.09179e6, 649.21543e6]))
        param.quartic = convert.Hz2J(num.array([1.0708e3, -5.736e3, 14.13e3, 0.4551e3, 1.341e3]))
        param.dipole = convert.D2Cm(num.array([4.2, 1.2, 0.]))



def indole_water2(param):
    """Molecular parameters for indole-(water)_2

    Implemented isomers are
    0  -  values calculated at B3LYP/6-31+G* with GAMESS-US 20XX by Yuan-Pin Chang (2011)
    """
    param.name = "indole-water2"
    param.mass = 8 * Masses['C'] + Masses['N'] + 2 * Masses['O'] + 11 * Masses['H']
    param.watson = 'A'
    param.symmetry = 'N'
    param.isomer = 0
    param.rotcon = convert.Hz2J(num.array([1323.5e6, 814.34e6, 587.86e6]))
    param.dipole = convert.D2Cm(num.array([1.46, -1.76, 1.31]))



def water(param):
    """Molecular parameters for H2O, D2O, HDO

    Implemented isomers are
    0  -  H2O: experimental inertial parameters from F. C. De Lucia, P. Helminger, R. L. Cook, and W. Gordy, Phys. Rev. A, 5, 487 (1972)
               and experimental dipole moment from Shostak, Ebenstein, and Muenter, J. Chem. Phys., 94, 5875 (1991)
    1  -  D2O: experimental inertial parameters from G. Steenbeckeliers, and J. Bellet, J. Mol. Spectrosc. 45, 10 (1973)
               and experimental dipole moment from Clough, Beers, Klein, Rothman, J. Chem. Phys. 59, 2254-2259 (1973)
    2  -  HDO: experimental inertial parameters from F. C. De Lucia, R. L. Cook, P. Helminger, and W. Gordy, J. Chem. Phys., 55, 5334 (1971)
          and experimental dipole moment from Shostak, Ebenstein, and Muenter, J. Chem. Phys., 94, 5875 (1991)

    Default isomers are 0 for water/H2O, 1 for D2O, and 2 for HDO.

    get better data from
    - http://physics.nist.gov/PhysRefData/MolSpec/Triatomic/Html/Tables/H2O.html
    - DeLucia, Helminger, Kirchhoff, J. Phys. Chem. Ref. Data 3, 211 (1974)
    - DeLucia and Helminger, J. Mol. Spectrosc. 56, 138 (1975)
    """
    param.name = "water"
    param.watson = 'A'
    if param.isomer == 0:
        param.mass = Masses['O'] + 2 * Masses['H']
        param.symmetry = 'C2b'
        param.rotcon = convert.Hz2J(num.array([835840.29e6, 435351.72e6, 278138.7e6]))
        param.quartic = convert.Hz2J(num.array([37.594e6, -172.91e6, 973.29e6, 15.210e6, 41.05e6]))
        param.dipole = convert.D2Cm(num.array([0., -1.857, 0.]))
    elif param.isomer == 1:
        param.mass = Masses['O'] + 2 * Masses['D']
        param.symmetry = 'C2b'
        param.rotcon = convert.Hz2J(num.array([462278.854e6, 218038.233e6, 145258.022e6]))
        param.dipole = convert.D2Cm(num.array([0., -1.8558, 0.]))
    elif param.isomer == 2:
        param.mass = Masses['O'] + Masses['H'] + Masses['D']
        param.symmetry = 'N'
        param.rotcon = convert.Hz2J(num.array([701931.50e6, 272912.60e6, 192055.25e6]))
        param.quartic = convert.Hz2J(num.array([10.8375e6, 34.208e6, 377.078e6, 3.6471e6, 63.087e6]))
        param.dipole = convert.D2Cm(num.array([-0.6591, -1.7304, 0.]))



def OCS(param):
    """Molecular parameters for OCS

    Paramters from http://physics.nist.gov/PhysRefData/MolSpec/Triatomic/Html/Tables/OCS.html (2012) and
    Reinartz, J., & Dymanus, A. Chemical Physics Letters, 24(3), 346–351 (1974).

    Implemented isomers are
    0  - using above parameters with linear-rotor hamiltonian
    1  - using above parameters with symmetric-rotor hamiltonian
    2  - using above parameters with asymmetric-rotor hamiltonian
    """
    param.name = "OCS"
    param.mass = Masses['O'] + Masses['C'] + Masses['S']
    if 0 == param.isomer:
	param.type = 'L'
	param.symmetry = 'N'
	param.rotcon = convert.Hz2J(num.array([6.081492475e9]))
	param.dipole = convert.D2Cm(num.array([0.71519]))
	param.quartic  = convert.Hz2J(num.array([1.301777e3]))
    elif 1 == param.isomer:
	param.type = 'S'
	param.symmetry = 'p'
	param.rotcon = convert.Hz2J(num.array([1e15, 6.081492475e9]))
	param.dipole = convert.D2Cm(num.array([0.71519]))
	param.quartic  = convert.Hz2J(num.array([1.301777e3, 0., 0.]))
    elif 2 == param.isomer:
	param.type = 'A'
	param.symmetry = 'C2a'
	param.rotcon = convert.Hz2J(num.array([1e15, 6.081492475e9, 6.081492475e9]))
	param.dipole = convert.D2Cm(num.array([0.71519, 0., 0.]))
	param.quartic  = convert.Hz2J(num.array([1.301777e3, 0., 0., 0., 0.]))


def iodomethane(param):
    """B, DJ, and DK constants from Wlodarczak, Boucher, Bocquet, & Demaison, J. Mol. Spectros., 124, 53–65 (1987) and
    Gadhi, Wlodarczak, Legrand, & Demaison, Chem. Phys. Lett., 156, 401–404 (1989).

    A and DK constants from Pietilä, Koivusaari, Alanko, & Anttila, Mol Phys 87, 523 (1996)

    Implemented isomers are
    0  - above constants using symmetric-top Hamiltonian
    1  - above constants using asymmetric-top Hamiltonian
    """
    param.name = "iodomethane"
    param.mass = 3*Masses['H'] + Masses['C'] + Masses['I']
    if 0 == param.isomer:
	param.type = 'S'
        param.symmetry = 'p'
	param.rotcon = num.array([convert.invcm2J(5.1742629), convert.Hz2J(7501.2757456e6)])
        param.quartic  = num.array([convert.Hz2J(6.307583e3), convert.Hz2J(98.76798e3), convert.invcm2J(87.857e-6)])
	param.dipole = convert.D2Cm(num.array([1.6406]))
    elif 1 == param.isomer:
	param.type = 'A'
        param.watson = 'A'
	param.symmetry = 'C2a'
	param.rotcon = num.array([convert.invcm2J(5.1742629), convert.Hz2J(7501.2757456e6), convert.Hz2J(7501.2757456e6)])
	param.quartic  = num.array([convert.Hz2J(6.307583e3), convert.Hz2J(98.76798e3), convert.invcm2J(87.857e-6), 0., 0.])
	param.dipole = convert.D2Cm(num.array([1.6406, 0., 0.]))



def difluoro_iodobenzene(param):
    # parameters from simple ab initio calculations (Jochen Küpper, 2010)
    param.name = "2,6-difluoro-iodobenzene"
    param.watson = 'A'
    param.symmetry = 'C2a'
    param.rotcon = convert.Hz2J(num.array([1740e6, 713e6, 506e6]))
    param.quartic = convert.Hz2J(num.array([0., 0., 0., 0., 0.]))
    param.dipole = convert.D2Cm(num.array([2.25, 0., 0.]))



def aminobenzonitrile(param):
    # Borst et al., Chem. Phys. Lett. 350, p.485 (2001)
    param.name = "4-aminobenzonitrile"
    param.watson = 'A'
    param.symmetry = 'C2a'
    param.rotcon = convert.Hz2J(num.array([5.5793e9, 0.99026e9, 0.84139e9]))
    param.quartic = convert.Hz2J(num.array([0.0, 0.0, 0.0, 0.0, 0.0]))
    param.dipole = convert.D2Cm(num.array([6.41, 0., 0.]))



def benzonitrile(param):
    # Wohlfart, Schnell, Grabow, Küpper, J. Mol. Spec. 247, 119-121 (2008)
    param.name = "benzonitrile"
    param.mass = 7 * Masses['C'] + Masses['N'] + 5 * Masses['H']
    param.watson = 'A'
    param.symmetry = 'C2a'
    param.rotcon = convert.Hz2J(num.array([5655.2654e6, 1546.875864e6, 1214.40399e6]))
    param.quartic = convert.Hz2J(num.array([45.6, 938.1, 500, 10.95, 628]))
    param.dipole = convert.D2Cm(num.array([4.5152, 0., 0.]))

def glycine(param):
    """Molecular parameters for TEST glycine
        
        Implemented isomers are
        0  -  Paper
        1  -  Paper
        2  -  Anthony Meijer
        3  -  Anthony Meijer
        4  -  Test
        5  -  Test
        """
    param.name = "glycine"
    param.mass = 2 * Masses['C'] + 5 * Masses['H'] + 1 * Masses['N'] + 2 * Masses['O']
    param.watson = 'A'
    param.symmetry = 'N'
    if param.isomer == 0: # cis, Filsinger et al. PCCP ...
        param.rotcon = convert.Hz2J(num.array([10.3415e9, 3.87618e9, 2.91235e9]))
        param.dipole = convert.D2Cm(num.array([0.911, 0.697, 0.]))
    elif param.isomer == 1:
        param.rotcon = convert.Hz2J(num.array([10.1301e9, 4.07151e9, 3.00748e9]))
        param.dipole = convert.D2Cm(num.array([5.372, 0.93, 0.]))
    elif param.isomer == 2:
        
    elif param.isomer == 3:
        param.rotcon = convert.Hz2J(num.array([10.2564941e9, 3.9707803e9, 2.9620284e9]))
        param.dipole = convert.D2Cm(num.array([-0.0058, -1.5519, 1.4356]))
    elif param.isomer == 4:
        param.rotcon = convert.Hz2J(num.array([0., 0., 0.]))
        param.dipole = convert.D2Cm(num.array([0., 0., 0.]))
    elif param.isomer == 5:
        param.rotcon = convert.Hz2J(num.array([0., 0., 0.]))
        param.dipole = convert.D2Cm(num.array([0., 0., 0.]))

def iodobenzene(param):
    # Dorosh, Bialkowskajaworska, Kisiel, Pszczolkowski,  J. Mol. Spec. 246, 228-232 (2007)
    param.name = "iodobenzene"
    param.watson = 'A'
    param.symmetry = 'C2a'
    param.rotcon = convert.Hz2J(num.array([5669.126e6, 750.414323e6, 662.636162e6]))
    param.quartic = convert.Hz2J(num.array([19.5479, 164.648, 891, 2.53098, 15554]))
    # param.sextic =  convert.Hz2J(num.array([0.0609, -0.377])) # ignored sextic constants!
    param.dipole = convert.D2Cm(num.array([1.6250, 0., 0.]))



def phenylpyrrole(param):
    # A. J. Fleisher
    param.name = "phenylpyrrole"
    param.watson = 'A'
    param.symmetry = 'C2a'
    param.rotcon = convert.Hz2J(num.array([3508.34e6, 703.50e6, 604.84e6]))
    param.dipole = convert.D2Cm(num.array([-1.56, 0., 0.]))

def three_fluorophenol(param):
    """Molecular parameters for three_fluorophenol
	Parameters come from Yuan Pin Chang
        0 - cis
        1 - trans
       """
    param.name = "three_fluorophenol"
    param.mass = 6 * Masses['C'] + 1 * Masses['F'] + 1 * Masses['O'] + 5 * Masses['H']
    param.watson = 'A'
    param.symmetry = 'N'
    if param.isomer == 0:
         param.rotcon = convert.Hz2J(num.array([3.74912e9, 1.78523e9, 1.20936e9]))
         param.dipole = convert.D2Cm(num.array([0.6251, 0.5345, 0.0025]))
    elif param.isomer == 1:
         param.rotcon = convert.Hz2J(num.array([3.74222e9, 1.79016e9, 1.21090e9]))
         param.dipole = convert.D2Cm(num.array([1.9206, 1.8098, 0.0001]))

def sulfur_dioxide(param):
    # Paul A. Helminger and Frank C. De Lucia JOURNAL OF MOLECULAR SPECTROSCOPY 111, 66-72 (1985)
    param.name = "sulfur_dioxide"
    param.watson = 'A'
    param.symmetry = 'C2b'
    param.rotcon = convert.Hz2J(num.array([60778.5522e6, 10318.0722e6, 8799.7023e6])) # A=2.026cm^-1, B=0.3442 cm^-1, C=0.2935 cm^-1.  Alternative papers: J. Chem. Phys. 19, 502 (1951), or  J. Chem. Phys. 22, 904 (1954).  or  F.J. Lovas, J. Phys. Chem. Ref. Data 7, 1445 (1978).
    param.quartic = convert.Hz2J(num.array([0.0066610013e6, -0.1169588e6, 2.5904328e6, 0.001701045, 0.0253531])) # Alternative papers: F.J. Lovas, J. Phys. Chem. Ref. Data 7, 1445 (1978).
    param.dipole = convert.D2Cm(num.array([0., 1.633189, 0.])) # Dipole from J. Chem. Phys. 70, 2740 (1979).

