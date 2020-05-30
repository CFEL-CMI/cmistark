# -*- coding: utf-8; fill-column: 100 -*-
#
# This file is part of CMIstark
# Copyright (C) 2008,2020 Jochen Küpper <software@jochen-kuepper.de>
#
# This program is free software: you can redistribute it and/or modify it under the terms of the GNU
# General Public License as published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
#
# If you use this programm for scietific work, you must correctly reference it; see LICENSE.md file
# for details.
#
# This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without
# even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with this program. If not,
# see <http://www.gnu.org/licenses/>.


__author__ = "Jochen Küpper <software@jochen-kuepper.de>"

"""Unit conversion routines"""

import numpy as np
import scipy.constants



def D2Cm(val):
    """Convert dipole moment from Debye to Coulomb * meter"""
    return np.array(val) * 1e-21 / scipy.constants.speed_of_light


def Cm2D(val):
    """Convert dipole moment from Coulomb * meter to Debye"""
    return np.array(val) / (1e-21 / scipy.constants.speed_of_light)


def Hz2J(val):
    """Hertz -> Joule"""
    return np.array(val) * scipy.constants.Planck


def invcm2Hz(val):
    """wavenumber (cm^{-1}) -> frequency (Hz)"""
    return val * 100 * scipy.constants.speed_of_light


def J2Hz(val):
    """Joule -> Hertz"""
    return np.array(val) / scipy.constants.Planck


def J2invcm(val):
    """Joule -> cm^{-1}"""
    return val / scipy.constants.Planck / scipy.constants.speed_of_light / 100


def kV_cm2V_m(val):
    """kV/cm -> V/m"""
    return np.array(val) / 1e-5


def V_m2kV_cm(val):
    """V/m -> kV/cm"""
    return np.array(val) * 1e-5
