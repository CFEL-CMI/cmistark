#!/usr/bin/env python
# -*- coding: utf-8; fill-column: 120 -*-
#
# This file is part of JK Python extensions
# Copyright (C) 2008 Jochen Küpper <software@jochen-kuepper.de>
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

"""Unit-tests of Stark effect calculations

Copyright (C) 2008 Jochen Küpper"""

__author__ = "Jochen Küpper <software@jochen-kuepper.de>"

import numpy as num
import os
import unittest

import jkext as jk
import jkext.convert as convert
import jkext.molecule as molecule
import jkext.starkeffect as starkeffect
from jkext.state import State


class Test_StarkCalculation_benzonitrile(unittest.TestCase):
    """Test the results of Stark effect calculations using the molecular parameters of benzonitrile"""

    def setUp(self):
        self.storagename = "__jkext_test_starkeffect.hdf"
        if os.path.exists(self.storagename):
            raise EnvironmentError("Test storage file already exists, not overwriting")
        # create Molecule object and specify storage file
        self.bn = molecule.Molecule(storage=self.storagename)
        # set molecular parameters
        self.param = starkeffect.CalculationParameter
        self.param.isomer = 0
        self.param.watson = 'A'
        self.param.symmetry = 'C2a'
        self.param.rotcon = convert.Hz2J(num.array([5655.2654e6, 1546.875864e6, 1214.40399e6]))
        self.param.quartic = convert.Hz2J(num.array([45.6, 938.1, 500, 10.95, 628]))
        self.param.dipole = convert.D2Cm(num.array([4.5152, 0., 0.]))
        # calculation details
        self.param.M = [0, 1]
        self.param.Jmin = 0
        self.param.Jmax_calc = 15
        self.param.Jmax_save =  3
        self.param.dcfields = convert.kV_cm2V_m(num.linspace(0., 100., 5))
        self.bn.starkeffect_calculation(self.param)

    def tearDown(self):
        os.remove(self.storagename)

    def test_fieldfree(self):
        self.assertAlmostEqual(0., self.bn.starkeffect(State(0, 0, 0, 0, 0))[1][0], 7, "Field-free ground state energy is wrong")

    def test_hundred(self):
        """Test some state energies at 100 kV/cm

        With our setup, these are he fifth values in the list of fields/energies.
        """
        # test (once) that the fields are correct
        self.assertAlmostEqual(convert.kV_cm2V_m(100.), self.bn.starkeffect(State(0, 0, 0, 0, 0))[0][4], 7,
                               "Field-free ground state energy is wrong")
        # test energies for different states at 100 kV/cm
        self.assertAlmostEqual(1., -1.34489847e-22 / self.bn.starkeffect(State(0, 0, 0, 0, 0))[1][4], 7,
                                "Field-free ground state energy is wrong: expected %g MHz, got %g MHz" \
                                    % (convert.J2MHz(-1.34489847e-22), convert.J2MHz(self.bn.starkeffect(State(0, 0, 0, 0, 0))[1][4])))


if __name__ == '__main__':
    unittest.main()
