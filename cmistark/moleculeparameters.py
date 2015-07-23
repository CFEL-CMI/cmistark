#!/usr/bin/env python
# -*- coding: utf-8; fill-column: 120; truncate-lines: t -*-
#
# This file is part of JK Python extensions
# Copyright (C) 2008,2009,2012,2014 Jochen Küpper <jochen.kuepper@cfel.de>
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

__author__ = "Jochen Küpper <jochen.kuepper@cfel.de>"
__doc__ = """This modules implements the molecular parameters of all investigated molecules.

The relevant parameters are:

- param.rotcon: rotational constants

  - linear top: :math:`B`
  - symmetric top: :math:`(A,B)` for prolate, :math:`(B,C)` for oblate
  - asymmetric top: :math:`(A,B,C)`

- param.quartic: centrifugal distortion constants

  - linear top: :math:`D`
  - symmetric top: :math:`(D_{J},D_{JK},D_{K})`
  - asymmetric top

    - in Watson's A reduction: :math:`(\Delta_{J}, \Delta_{JK}, \Delta_{K}, d_{J}, d_{K})`
    - in Watson's S reduction: :math:`(D_{J}, D_{JK}, D_{K}, d_{1}, d_{2})`

- param.dipole: dipole moments

  - for linear and symmetric tops: :math:`\mu`
  - for an asymmetric top: (:math:`\mu_{a}`, :math:`\mu_{a}`, :math:`\mu_{a}`)

- param.type: type of rotors

  - linear rotor: *L*
  - symmetric top: *S*
  - asymmetric top: *A*

- param.symmetry: symmetry in the feild for linear/asymmetric tops. For symmetric top,
  prolate or oblate is specifie here

  - linear rotor: *N* (as no symmetry is implemented for linear top)
  - symmetric top:

    - prolate: *p*
    - oblate: *s*

  - asymmetric top (for :math:`M != 0` cases, the program takes care the :math:`M = 0` case itself):

    - only :math:`\mu_a != 0`: *C2a*
    - only :math:`\mu_b != 0`: *C2b*
    - only :math:`\mu_c != 0`: *C2c*
    - other dipole directions: *N*

- param.watson: only for asymmetric top.

  - Watson's A reduction: *A*
  - Watson's S reduction: *S*

All relevant parameters for molecules of interest need to be properly implemented here.

"""

import numpy as num
import getopt, sys

import cmiext as cmiext
import cmiext.convert as convert
from cmiext.state import State
from cmiext.molecule import Masses

import cmistark.molecule as molecule
import cmistark.starkeffect as starkeffect


def asymmetric_top(param):
    """Molecular parameters for an artificial asymmetric top

    Implemented isomers are (modified) examples of

    0. rot. const. from Fig7.2 in Gordy & Cook (1984), and only :math:`\mu_b != 0`
    1. rot. const. from Fig7.2 in Gordy & Cook (1984), and only :math:`\mu_c != 0`
    2. rot. const. from Fig7.2 in Gordy & Cook (1984), and only :math:`\mu_a = 0`
    3. rot. const. from Fig7.2 in Gordy & Cook (1984), and only :math:`\mu_b = 0`
    4. rot. const. from Fig7.2 in Gordy & Cook (1984), and only :math:`\mu_c = 0`
    5. rot. const. from Fig7.2 in Gordy & Cook (1984), and no any :math:`\mu_i = 0`
    """
    param.name = "asymmetric_top"
    param.watson = 'A'
    if param.isomer == 0:
        param.symmetry = 'C2b'
        param.rotcon = convert.Hz2J(num.array([3000.0e6, 2000.0e6, 1000.0e6]))
        param.dipole = convert.D2Cm(num.array([0., 1., 0.]))
    elif param.isomer == 1:
        param.symmetry = 'C2c'
        param.rotcon = convert.Hz2J(num.array([3000.0e6, 2000.0e6, 1000.0e6]))
        param.dipole = convert.D2Cm(num.array([0., 0., 1.]))
    elif param.isomer == 2:
        param.symmetry = 'N'
        param.rotcon = convert.Hz2J(num.array([3000.0e6, 2000.0e6, 1000.0e6]))
        param.dipole = convert.D2Cm(num.array([0., 1., 1.]))
    elif param.isomer == 3:
        param.symmetry = 'N'
        param.rotcon = convert.Hz2J(num.array([3000.0e6, 2000.0e6, 1000.0e6]))
        param.dipole = convert.D2Cm(num.array([1., 0., 1.]))
    elif param.isomer == 4:
        param.symmetry = 'N'
        param.rotcon = convert.Hz2J(num.array([3000.0e6, 2000.0e6, 1000.0e6]))
        param.dipole = convert.D2Cm(num.array([1., 1., 0.]))
    elif param.isomer == 5:
        param.symmetry = 'N'
        param.rotcon = convert.Hz2J(num.array([3000.0e6, 2000.0e6, 1000.0e6]))
        param.dipole = convert.D2Cm(num.array([1., 1., 1.]))



def three_aminophenol(param):
    """Molecular parameters for 3-aminophenol

    Implemented isomers are

    0.  cis conformer, experimental values from [Filsinger2008]_
    1.  trans conformer, exp values from [Filsinger2008]_
    2.  cis conformer, MP2/aug-cc-pVTZ calculation using <Gaussian 2003.1> from [Roesch2011]_
    3.  trans conformer, calculated values of MP2/aug-cc-pVTZ method from [Roesch2011]_
    4.  cis conformer, calculated values of B3LYP/aug-cc-pVTZ method from [Roesch2011]_
    5.  trans conformer, calculated values of B3LYP/aug-cc-pVTZ method from [Roesch2011]_
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


def oblate_symmetric_top(param):
    """Molecular parameters for an artificial oblate top

    Implemented isomers are

    0. rotor type is set to 'S' - symmetric top, and symmetry is set to 'o' - oblate
    1. rotor type is set to 'A' - asymmetric top, and symmetry is set to 'C2c' - only u_c != 0
    """
    param.name = "oblate_symmetric_top"
    param.mass = 6 * Masses['C'] + 6 * Masses['H']
    if 0 == param.isomer:
        param.type = 'S'
        param.symmetry = 'o'
        param.rotcon = convert.Hz2J(num.array([3000.0e6, 1000.0e6]))
        param.quartic  = num.array([convert.Hz2J(0.0), convert.Hz2J(0.0), convert.invcm2J(0.0)])
        param.dipole = convert.D2Cm(num.array([1.]))
    elif 1 == param.isomer:
        param.type = 'A'
        param.watson = 'A'
        param.symmetry = 'C2c'
        param.rotcon = convert.Hz2J(num.array([3000.0e6, 3000.0e6, 1000.0e6]))
        param.quartic  = num.array([0., 0., 0., 0., 0.])
        param.dipole = convert.D2Cm(num.array([0., 0., 1.]))


def prolate_symmetric_top(param):
    """Molecular parameters for an artificial prolate top"""
    param.name = "prolate_symmetric_top"
    param.mass = 6 * Masses['C'] + 6 * Masses['H']
    if 0 == param.isomer:
        param.type = 'S'
        param.symmetry = 'p'
        param.rotcon = convert.Hz2J(num.array([3000.0e6, 1000.0e6]))
        param.quartic  = num.array([convert.Hz2J(0.0), convert.Hz2J(0.0), convert.invcm2J(0.0)])
        param.dipole = convert.D2Cm(num.array([1.]))
    elif 1 == param.isomer:
        param.type = 'A'
        param.watson = 'A'
        param.symmetry = 'C2a'
        param.rotcon = convert.Hz2J(num.array([3000.0e6, 1000.0e6, 1000.0e6]))
        param.quartic  = num.array([0., 0., 0., 0., 0.])
        param.dipole = convert.D2Cm(num.array([1., 0., 0.]))


def indole(param):
    """Molecular parameters for indole

    Implemented isomers are

    0.  experimental values from [Kang2005]_
    1.  experimental inertial constants from [Caminati1990]_ and dipole moment from [Kang2005]_ for dipole moment.
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

    0.  experimental inertial parameters from [Korter1998]_ and experimental dipole moment from [Kang2005]_
    1.  experimental inertial parameters from [Blanco2003]_ and experimental dipole moment from [Kang2005]_
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
    """Molecular parameters for indole-(water):math:`_2`

    Implemented isomers are

    0.  values calculated at B3LYP/6-31+G* with GAMESS-US 2009 by Y.P. Chang; see [Trippel2012]_

    .. todo:: See the math-usage and implement it for all sub- and super-scripts (if any)
    """
    param.name = "indole-water2"
    param.mass = 8 * Masses['C'] + Masses['N'] + 2 * Masses['O'] + 11 * Masses['H']
    param.watson = 'A'
    param.symmetry = 'N'
    param.isomer = 0
    param.rotcon = convert.Hz2J(num.array([1323.5e6, 814.34e6, 587.86e6]))
    param.dipole = convert.D2Cm(num.array([1.46, -1.76, 1.31]))


def water(param):
    """Molecular parameters for water isotopologues (:math:`\\text{H}_2\\text{O}`, :math:`\\text{D}_2\\text{O}`, :math:`\\text{HDO}`)

    :param param: Calculation parameter object to be filled with appropriate content
    :type param: starkeffect.CalculationParameter

    Implemented isomers are

    0. H2O: experimental inertial parameters from [DeLucia1974]_ and experimental dipole moment from [Shostak1991]_
    1. D2O: experimental inertial parameters from [Steenbeckeliers1973]_ and experimental dipole moment from [Clough1973]_
    2. HDO: experimental inertial parameters from [DeLucia1971]_ and experimental dipole moment from [Shostak1991]_

    Default isomers are '0' for water/H2O, '1' for D2O, and '2' for HDO.

    .. seealso:: These molecular parameters and references are also listed at NIST Spectral Database - H2O [NISTspecDB_H2O]_

    """
    param.name = "water"
    param.watson = 'A'
    if param.isomer == 0:
        param.mass = Masses['O'] + 2 * Masses['H']
        param.symmetry = 'C2b'
        param.rotcon = convert.Hz2J(num.array([835840.288e6, 435351.717e6, 278138.700e6]))
        param.quartic = convert.Hz2J(num.array([37.59422e6, -172.9128e6, 973.29052e6, 15.210402e6, 41.0502e6]))
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

    Paramters from NIST Spectral Database - OCS [NISTspecDB_OCS]_ and [Reinartz1974]_.

    Implemented isomers are

    0. using above parameters with linear-rotor hamiltonian
    1. using above parameters with symmetric-rotor hamiltonian
    2. using above parameters with asymmetric-rotor hamiltonian
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
    """B, DJ, and DK constants from [Wlodarczak1987]_ and [Gadhi1989]_.

    A and DK constants from [Pietila1996]_.

    Implemented isomers are
    0. above constants using symmetric-top Hamiltonian
    1. above constants using asymmetric-top Hamiltonian
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
    """ parameters from simple ab initio calculations [Kuepper2010]_.
    """
    param.name = "2,6-difluoro-iodobenzene"
    param.type = 'A'
    param.watson = 'A'
    param.symmetry = 'C2a'
    param.rotcon = convert.Hz2J(num.array([1740e6, 713e6, 506e6]))
    param.quartic = convert.Hz2J(num.array([0., 0., 0., 0., 0.]))
    param.dipole = convert.D2Cm(num.array([2.25, 0., 0.]))


def diiodoethane(param):
    """Molecular parameters for diiodo-ethane.

    Implemented isomers are
    0.  anti-conformation (C2h symmetry)
    1.  gauge-conformation (C2 symmetry)

    Structural parameters are from the supplementary material of Qingyu Kong et al., "Photodissociation Reaction of
    1,2-Diiodoethane in Solution: A Theoretical and X-ray Diffraction Study" [Kong2005]_

    Rotational constants were then calculated with gamess (calculation level: MP2/6-311G**).
    """
    param.name = "diiodoethane"
    param.mass = 2 * Masses['C'] + 4 * Masses['H'] +  2 * Masses['I']
    param.watson = 'A'
    if param.isomer == 0: # anti
        param.symmetry = 'N'
        param.rotcon = convert.Hz2J(num.array([2.79227492e10, 3.09470163e8, 3.07270856e8]))
        param.dipole = convert.D2Cm(num.array([0., 0., 0.]))
    elif param.isomer == 1: # gauge
        param.symmetry = 'C2b'
        param.rotcon = convert.Hz2J(num.array([6.18776669e9, 4.92744613e8, 4.63936461e8]))
        param.dipole = convert.D2Cm(num.array([0., 2.249726, 0.]))


def two_aminobenzonitrile(param):
    """ See [Miller2009]_
    """
    param.name = "2-aminobenzonitrile"
    param.type = 'A'
    param.watson = 'A'
    param.symmetry = 'N'
    param.rotcon = convert.Hz2J(num.array([3.0090e9, 1.5090e9, 1.0052e9]))
    param.quartic = convert.Hz2J(num.array([0.0, 0.0, 0.0, 0.0, 0.0]))
    param.dipole = convert.D2Cm(num.array([3.6, 1.9, 0.]))


def three_aminobenzonitrile(param):
    """ See [Miller2009]_
    """
    param.name = "3-aminobenzonitrile"
    param.type = 'A'
    param.watson = 'A'
    param.symmetry = 'N'
    param.rotcon = convert.Hz2J(num.array([3.3727e9, 1.2099e9, 0.8908e9]))
    param.quartic = convert.Hz2J(num.array([0.0, 0.0, 0.0, 0.0, 0.0]))
    param.dipole = convert.D2Cm(num.array([4.8, 1.2, 0.]))


def four_aminobenzonitrile(param):
    """ See [Borst2001]_
    """
    param.name = "4-aminobenzonitrile"
    param.type = 'A'
    param.watson = 'A'
    param.symmetry = 'C2a'
    param.rotcon = convert.Hz2J(num.array([5.5793e9, 0.99026e9, 0.84139e9]))
    param.quartic = convert.Hz2J(num.array([0.0, 0.0, 0.0, 0.0, 0.0]))
    param.dipole = convert.D2Cm(num.array([6.41, 0., 0.]))


def benzonitrile(param):
    """ See [Wohlfart2008]_
    """
    param.name = "benzonitrile"
    param.mass = 7 * Masses['C'] + Masses['N'] + 5 * Masses['H']
    param.type = 'A'
    param.watson = 'A'
    param.symmetry = 'C2a'
    param.rotcon = convert.Hz2J(num.array([5655.2654e6, 1546.875864e6, 1214.40399e6]))
    param.quartic = convert.Hz2J(num.array([45.6, 938.1, 500.0, 10.95, 628.0]))
    param.dipole = convert.D2Cm(num.array([4.5152, 0., 0.]))


def glycine(param):
    """Molecular parameters for TEST glycine

        Implemented isomers are
        0.  Paper
        1.  Paper
        2.  Anthony Meijer
        3.  Anthony Meijer
        4.  Test
        5.  Test
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
        param.rotcon = convert.Hz2J(num.array([9.71997e9, 3.97849e9, 2.98658e9]))
        param.dipole = convert.D2Cm(num.array([-0.1559, 1.6907, -0.0773]))
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
    """ See [Dorosh2007]
    """
    param.name = "iodobenzene"
    param.watson = 'A'
    param.symmetry = 'C2a'
    param.rotcon = convert.Hz2J(num.array([5669.126e6, 750.414323e6, 662.636162e6]))
    param.quartic = convert.Hz2J(num.array([19.5479, 164.648, 891, 2.53098, 15554]))
    # param.sextic =  convert.Hz2J(num.array([0.0609, -0.377])) # ignored sextic constants!
    param.dipole = convert.D2Cm(num.array([1.6250, 0., 0.]))


def phenylpyrrole(param):
    """ A. J. Fleisher
    """
    param.name = "phenylpyrrole"
    param.watson = 'A'
    param.symmetry = 'C2a'
    param.rotcon = convert.Hz2J(num.array([3508.34e6, 703.50e6, 604.84e6]))
    param.dipole = convert.D2Cm(num.array([-1.56, 0., 0.]))


def three_fluorophenol(param):
    """Molecular parameters for three_fluorophenol

    Implemented isomers are
    0. cis, Parameters (rot con and quadratic) from [Dutta1985]_
    1. trans, from [Jaman1981]_ for dipoles from [Songhee2011]_ (YPChang)
    11. cis, from [Songhee2011]_ (YPChang)
    12. trans, from [Songhee2011]_ (YPChang
    13. cis, calculated by Yuan-Pin Chang, also the dipole moments
    14. trans, calculated by Yuan-Pin Chang, also the dipole moments
    """
    param.name = "three_fluorophenol"
    param.mass = 6 * Masses['C'] + 1 * Masses['F'] + 1 * Masses['O'] + 5 * Masses['H']
    param.type = 'A'
    param.watson = 'S'
    param.symmetry = 'N'
    if param.isomer == 0:
        param.rotcon = convert.Hz2J(num.array([3.78957e9, 1.79515e9, 1.21789e9]))
        param.quartic = convert.Hz2J(num.array([-3.68e3,4.42e3,11e3,0.2,-3.2]))
        param.dipole = convert.D2Cm(num.array([0.6251, 0.5345, 0.0025]))
    elif param.isomer == 1:
        param.rotcon = convert.Hz2J(num.array([3.748487e9, 1.797713e9, 1.215048e9]))
        param.quartic = convert.Hz2J(num.array([0.2257e3,-10.977e3,-15.005e3,-0.2,7.2]))
        param.dipole = convert.D2Cm(num.array([1.9206, 1.8098, 0.0001]))
    elif param.isomer == 11:
        param.rotcon = convert.Hz2J(num.array([3.74912e9, 1.78523e9, 1.20936e9]))
        param.dipole = convert.D2Cm(num.array([0.6251, 0.5345, 0.0025]))
    elif param.isomer == 12:
        param.rotcon = convert.Hz2J(num.array([3.74222e9, 1.79016e9, 1.21090e9]))
        param.dipole = convert.D2Cm(num.array([1.9206, 1.8098, 0.0001]))
    elif param.isomer == 13:
        param.rotcon = convert.Hz2J(num.array([3.76709624e9, 1.794471e9, 1.21547482e9]))
        param.dipole = convert.D2Cm(num.array([0.276988, -0.781392, 0.0]))
    elif param.isomer == 14:
        param.rotcon = convert.Hz2J(num.array([3.76358267e9, 1.79886976e9, 1.21712413e9]))
        param.dipole = convert.D2Cm(num.array([2.747456, -0.144641, 0.0]))


def sulfur_dioxide(param):
    """ SO2 (fix nomenclature)

    Rotational constants from [Helminger1985]_
    Dipole moments from [Patel1979]_
    """
    param.name = "sulfur_dioxide"
    param.watson = 'A'
    param.symmetry = 'C2b'
    # A=2.026cm^-1, B=0.3442 cm^-1, C=0.2935 cm^-1.
    # Alternative papers: J. Chem. Phys. 19, 502 (1951), or  J. Chem. Phys. 22, 904 (1954). or
    # F.J. Lovas, J. Phys. Chem. Ref. Data 7, 1445 (1978).
    param.rotcon = convert.Hz2J(num.array([60778.5522e6, 10318.0722e6, 8799.7023e6]))
    # Alternative papers: F.J. Lovas, J. Phys. Chem. Ref. Data 7, 1445 (1978).
    param.quartic = convert.Hz2J(num.array([0.0066610013e6, -0.1169588e6, 2.5904328e6, 0.001701045, 0.0253531]))
    param.dipole = convert.D2Cm(num.array([0., 1.633189, 0.]))


def nitrogen_dioxide(param):
    """Molecular parameters for :math:`\\text{N}\\text{O}_2`

    rotational constants: [Cabana1976]_
    dipole moment: [Hodgeson1963]_
    """
    param.name = "nitrogen_dioxide"
    param.watson = 'A'
    param.symmetry = 'C2b'
    param.rotcon = convert.Hz2J(num.array([239905.41e6, 13002.262e6, 12304.888e6]))
    param.quartic = convert.Hz2J(num.array([9.033e3, -0.5903e6, 80.94e6, 9.303e2, 0.12e6]))
    param.dipole = convert.D2Cm(num.array([0., 0.316, 0.]))


def nitrous_oxide(param):
    """Molecular parameters for :math:`\\text{N}_2\\text{O}`

    rot. const.: [Andreev1976]
    dipole moment: [Scharpen1970]

    """
    param.name = "nitrous_oxide"
    param.type = 'L'
    param.symmetry = 'N'
    param.rotcon = convert.Hz2J(num.array([12561.6348e6]))
    param.dipole = convert.D2Cm(num.array([0.160830]))
    param.quartic  = convert.Hz2J(num.array([5.2808e6]))


def methylvinylketone(param):
    """Molecular parameters for Methyl Vinyl Ketone

    Implemented isomers are
    0. cis, experimetnal rot. const. of A symmetry from [Wilcox2011]_, and dipole moment is MP2 calculation results from
       [Wilcox2011]_ with 0.81 scaling factor

    1. trans, experimental  rot. const. of A symmetry from [Wilcox2011]_, and experimental dipole moment: [Foster1965]_
    """
    param.name = "methylvinylketone"
    param.type = 'A'
    param.symmetry = 'N'
    if param.isomer == 0: #cis
        param.rotcon = convert.Hz2J(num.array([10.240938e9, 3.9916351e9, 2.925648e9]))
        param.dipole = convert.D2Cm(num.array([0.54, 2.59, 0.]))
    elif param.isomer == 1: #trans
        param.rotcon = convert.Hz2J(num.array([8.94159e9, 4.2745443e9, 2.9453315e9]))
        param.dipole = convert.D2Cm(num.array([2.53, 1.91, 0.]))


def six_chloropyridazine_three_carbonitrile(param):
    """ Molecular parameters for 6-chloropyridazine-3-carbonitrile

    Gaussian 2003 B3LYP/aug-pc-1; see [Hansen2013]_"""
    param.name = "6-chloropyridazine-3-carbonitrile"
    param.mass = 5 * Masses['C'] + 1 * Masses['Cl'] + 3 * Masses['N'] + 2 * Masses['H']
    param.type = 'A'
    param.watson = None
    param.symmetry = 'N'
    if param.isomer == 0:
        param.rotcon = convert.Hz2J(num.array([5905.472e6, 717.422e6, 639.708e6]))
        param.dipole = convert.D2Cm(num.array([0, 4.37, 2.83]))
    elif param.isomer == 1:
        param.rotcon = convert.Hz2J(num.array([5905.472e6, 717.422e6, 639.708e6]))
        param.dipole = convert.D2Cm(num.array([0, 0, 2.83]))


def sulfur_monoxide(param):
    """ Molecular parameters for SO

    rotcon, dipole: [NISTCCCBD_SO]_
    quartic: [Veseth1974]_
    """
    param.name = "sulfur_monoxide"
    param.mass = Masses['S'] + Masses['O']
    param.type = 'L'
    param.rotcon = convert.Hz2J(num.array([21.60970e9]))
    param.dipole = convert.D2Cm(num.array([1.550]))
    param.quartic  = convert.Hz2J(num.array([33.577e3]))


def carbon_monoxide(param):
    """ Molecular parameters for CO

    rotcon, dipole: [NISTCCCBD_CO]_
    quartic: [MinaCamilde1996]_
    """
    param.name = "carbon_monoxide"
    param.mass = Masses['C'] + Masses['O']
    param.type = 'L'
    param.rotcon = convert.Hz2J(num.array([57.89834e9]))
    param.dipole = convert.D2Cm(num.array([0.11]))
    param.quartic  = convert.invcm2J(num.array([202.360e3]))


def five_cyanoindole(param):
    """ Molecular parameters for 5-cyanoindole

    Experimental values for rot constants from [Oelterman2012]_.
    Dipole values calculated (aug-cc-pVTZ basis, see [Horke]_ for details)
    """
    param.name = "five_cyanoindole"
    param.mass = 9 * Masses['C'] + 2 * Masses['N'] + 6 * Masses['H']
    param.watson = 'A'
    param.symmetry = 'N'
    param.rotcon = convert.Hz2J(num.array([3370.4e6, 738.0e6, 605.93e6]))
    param.dipole = convert.D2Cm(num.array([-6.44, -2.84, -0.33]))


def uracil(param):
    """ Molecular parameters for uracil

    Dipole & Rot. constants: [Brown1988]_, Dipole errors are estimated by 10%

    Centr. dist. const: [Brunken2006]_

    Good summary: [Puzzarini2011]_
    """
    param.name = "uracil"
    param.mass = 4 * Masses['C'] + 4 * Masses['H'] + 2 * Masses['N'] + 2 * Masses['O']
    param.type = 'A'
    param.watson = 'S'
    param.rotcon = convert.Hz2J(num.array([3883878.25e3, 2023732.67e3, 1330923.80e3]))
    param.dipole = convert.D2Cm(num.array([1.61, 3.52, 0.0]))
    param.quartic  = convert.Hz2J(num.array([0.06336e3, 0.1055e3, 0.4530e3, -0.02623e3, -0.00680e3]))


def mephenesin(param):
    """ Molecular parameters for mephenesin

    Values for rot constants and dipole moments from Ecija et al, JPC B 118, 5357
    (Rotational constants are experimental values, dipole moment values calculated with MP2)
    """

    param.name = "mephenesin"
    param.mass = 10 * Masses['C'] + 3 * Masses['O'] + 14 * Masses['H']
    param.watson = 'A'
    param.symmetry = 'N'
    if param.isomer == 0: #conformer A
        param.rotcon = convert.Hz2J(num.array([1707.7896e6, 388.661705e6, 331.331684e6]))
        param.dipole = convert.D2Cm(num.array([1.15, 0.56, 1.12]))
    elif param.isomer == 1: #conformer B
        param.rotcon = convert.Hz2J(num.array([1978.986e6, 349.300307e6, 305.408511e6]))
        param.dipole = convert.D2Cm(num.array([-2.36, -0.48, -0.04]))
    elif param.isomer == 2: #conformer C
        param.rotcon = convert.Hz2J(num.array([1615.04911e6, 455.423567e6, 385.954447e6]))
        param.dipole = convert.D2Cm(num.array([1.47, -1.32, -1.62]))

def five_fluoroindole(param):
    """ Molecular parameters for 5-fluoroindole

    Rotational constants constants from Brand et al, CPC 13, 3134 (2012)
    Dipole moments from MP2 calculations (Daniel Horke, February 2015, DFT/basis...)
    """
    param.name = "five_fluoroindole"
    param.mass = 8 * Masses['C'] + 6 * Masses['H'] + 1 * Masses['N'] + 1 * Masses['F']
    param.watson = 'A'
    param.symmetry = 'N'
    param.rotcon = convert.Hz2J(num.array([3519.57e6, 1019.79e6, 790.87e6]))
    param.dipole = convert.D2Cm(num.array([-3.40, -2.52, 0.0]))

def adenine(param):
    """ Molecular parameters for adenine
        
        TEMPORARY VALUES - to be replaced with more accurate data asap
        Rotational constants from Brown et al, CPL 156, 61 (1989)
        Dipole moments from Mohammed et al, JMS 938, 263 (2009)
        """
    param.name = "adenine"
    param.mass = 5 * Masses['C'] + 5 * Masses['H'] + 5 * Masses['N']
    param.watson = 'A'
    param.symmetry = 'N'
    if param.isomer == 0: #tautomer 9H
        param.rotcon = convert.Hz2J(num.array([2371.873e6, 1573.3565e6, 946.2576e6]))
        param.dipole = convert.D2Cm(num.array([2.2, 1.0, 0.0]))
    elif param.isomer == 1: #tautomer 7H
        param.rotcon = convert.Hz2J(num.array([2388e6, 1543e6, 941e6]))
        param.dipole = convert.D2Cm(num.array([0.3, 7.5, 0.0]))

