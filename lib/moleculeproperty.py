#!/usr/bin/env python
# -*- coding: utf-8; fill-column: 120 -*-
#
# This file is part of JK Python extensions
# Copyright (C) 2008,2009 Jochen Küpper <software@jochen-kuepper.de>
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

__author__ = "Jochen Küpper <software@jochen-kuepper.de>"

import numpy as num
import numpy.linalg

import cmistark.const as const


Masses = {'H': 1.0078250321, 'D': 2.01410178, '2H': 2.01410178,
          'C': 12, 'N': 14.0030740052, 'O': 15.9949146221,
          'S': 31.97207070,
          'F': 18.9984032,
          'BR': 78.9183371, 'BR79': 78.9183371, 'BR81': 80.9162906,
          'I': 126.90447}
Ordernumbers = {'H': 1, 'D': 1, 'C': 6, 'N': 7, 'O': 8, 'F': 9, 'S': 16, 'BR': 35, 'I': 53}


class Atom(object):
    """Representation of an atom

    Keeps a (private) list of Z and mass that can be accessed via public methods and a (public) copy of positions.

    Internal units are SI (i.e., m, kg, ...)! You can specify what length-unit is used on input (SI or Angstrom).
    """
    def __init__(self, symbol, position, length="SI"):
        symbol = symbol
        self.__Z = Ordernumbers[symbol.upper()]
        self.__mass = Masses[symbol.upper()] * const.unified_atomic_mass
        self.__symbol = symbol
        self.position = num.array(position, num.float)
        assert(self.position.shape == (3,))
        if length == "Angstrom":
            self.position *= num.float(const.Angstrom)
        else:
            assert(length == "SI")

    def mass(self):
        return self.__mass

    def ordernumber(self):
        return self.__Z

    def position(self):
        return self.position

    def symbol(self):
        return self.__symbol

    def Z(self):
        return self.__Z



class Molecule(object):
    """Representation of a Molecule

    Keeps a (private) list of atoms and masses that can be accessed via |atom| and |masses| and a (public) array of
    positions.
    """

    def __init__(self, atoms=None, name="Generic molecule"):
        """Create Molecule from a list of atoms."""
        self.__atoms = atoms
        self.__name = name
        if atoms != None:
            self.__update()


    def __update(self):
        self.masses = num.zeros((len(self.__atoms),))
        for i in range(self.masses.shape[0]):
            self.masses[i] = self.__atoms[i].mass()
        self.positions = num.zeros((len(self.__atoms), 3))
        for i in range(self.positions.shape[0]):
            self.positions[i,:] = self.__atoms[i].position


    def read(self, filename, type="Zxyz"):
        """Read molecular strcutre from file"""
        if type == Zxyz:
            pass
        else:
            raise TypeError("unknown molecule type")


    def atoms(self):
        """List of molecule's atoms"""
        return self.__atoms


    def center_of_mass(self):
        """Calculate center of mass of molecule"""
        return num.sum(num.outer(self.masses, [1,1,1]) * self.positions, axis=0) / self.mass()


    def plot(self):
        """Create 2D plot of molecule."""
        import matplotlib.pyplot as plt
        plt.plot(self.positions[:,0], self.positions[:,1], 'o')
        plt.show()


    def mass(self, unit="kg"):
        """Sum of atomic masses"""
        if "u" == unit:
            k = 1 / const.unified_atomic_mass
        else:
            assert "kg" == unit
            k = 1
        return k * num.sum(self.masses)


    def principal_axis_of_inertia(self):
        """Calulate principal axes of inertia and the corresponding rotational constants.

        Returns array of rotational constants [A, B, C] (Hz) and the right eigenvectors of the inertial tensor.
        """
        inertial_tensor = num.zeros((3,3))
        for i in range(0,3):
            for j in range(0,3):
                for k in range(len(self.masses)):
                    m = self.masses[k]
                    r = self.positions[k]
                    if i == j:
                        r2 = num.dot(r, r) # square of length
                    else:
                        r2 = 0.
                    inertial_tensor[i,j] += m * (r2 - r[i]*r[j])
        eval, evec = num.linalg.eigh(inertial_tensor)
        # sort eigenvalues and eigenvetors
        idx = num.argsort(eval) # sort moments of inertia in ascending order
        # calculate rotational constants in Hz
        rotcon = const.Planck_constant / (8 * num.pi**2 * eval[idx]) # in Hz
        # and provide corresponding eigenvectors of the axes (change columns!)
        axes = evec[:, idx]
        return rotcon, axes


    def rotate(self, rotation):
        """Rotate molecule using the right rotation matrix given (i.e., rotate all atomic coordinates)."""
        self.positions = num.dot(self.positions, rotation)
        return self.positions


    def to_principal_axis_of_inertia(self):
        """Put molecule into its principal axis system."""
        # move to center of mass
        self.translate(-self.center_of_mass())
        # rotate to pricipal axis system
        rotcon, axes = self.principal_axis_of_inertia()
        self.rotate(axes)
        return rotcon


    def translate(self, translation):
        """Translate center of mass of molecule (i.e., translate all atoms)."""
        self.positions += translation



# some simple tests
if __name__ == "__main__":
    pass
