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
import const
import tables
import jkext.hdf5, jkext.starkeffect, jkext.util, jkext.convert
from jkext.state import State


Masses = {'H': 1.0078250321, 'C': 12, 'N': 14.0030740052, 'O': 15.9949146221, 'Br': 79., 'I': 126.90447}
Ordernumbers = {'H': 1, 'D': 1, 'C': 6, 'N': 7, 'O': 8, 'Br': 35, 'I': 53}


class Atom:
    """Representation of an atom

    Keeps a (private) list of Z and mass that can be accessed via public methods and a (public) copy of positions.

    Internal units are SI (i.e., m, kg, ...)! You can specify what length-unit is used on input (SI or Angstrom).
    """
    def __init__(self, symbol, position, length="SI"):
        self.__Z = Ordernumbers[symbol]
        self.__mass = Masses[symbol] * const.unified_atomic_mass
        self.__symbol = symbol
        self.position = num.array(position)
        assert(self.position.shape == (3,))
        if length == "Angstrom":
            self.position *= const.angstrom
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



class Molecule:
    """Representation of a Molecule

    Keeps a (private) list of atoms and masses that can be accessed via |atom| and |masses| and a (public) array of
    positions.
    """

    def __init__(self, atoms=None, storage=None, name="Generic molecule"):
        """Create Molecule from a list of atoms."""
        self.__atoms = atoms
        self.__name = name
        if storage != None:
            self.__storage = tables.openFile(storage, mode='a', title=name)
        else:
            self.__storage = None
        if atoms != None:
            self.__update()


    def __update(self):
        self.masses = num.zeros((len(self.__atoms),))
        for i in range(self.masses.shape[0]):
            self.masses[i] = self.__atoms[i].mass()
        self.positions = num.zeros((len(self.__atoms), 3))
        for i in range(self.positions.shape[0]):
            self.positions[i,:] = self.__atoms[i].position


    def center_of_mass(self):
        """Calculate center of mass of molecule"""
        return num.sum(num.outer(self.masses, [1,1,1]) * self.positions, axis=0) / num.sum(self.masses)


    def plot(self):
        """Create 2D plot of molecule."""
        import matplotlib.pyplot as plt
        plt.plot(self.positions[:,0], self.positions[:,1], 'o')
        plt.show()


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
        rotcon = const.plancks_constant_h / (8 * num.pi**2 * eval[idx]) # in Hz
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


    def mueff(self, state):
        """Get the effective dipole moment \mu_eff as a function of the electric field strength.

        Return the effective dipole moment curve for the specified quantum |state|.
        """
        fields, energies = self.starkeffect(state)
        assert len(fields) == len(energies)
        mueff = num.zeros((len(fields),), num.float64)
        mueff[1:-1] = -1 * (energies[0:-2] - energies[2:]) / (fields[0:-2] - fields[2:])
        mueff[0] = 0.
        mueff[-1] = mueff[-2]
        return fields, mueff


    def starkeffect(self, state, fields=None, energies=None):
        """Get or set the potential energies as a function of the electric field strength.

        When |energies| and |fields| are None, return the Stark curve for the specified quantum state.

        When |energies| and |fields| are specified, save the Stark curve for the specified quantum state in the
        Molecule's HDF5 storage file.
        """
        if energies == None and fields == None:
            return jkext.hdf5.readVLArray(self.__storage, "/" + state.hdfname() + "/dcfield"), \
                jkext.hdf5.readVLArray(self.__storage, "/" + state.hdfname() + "/dcstarkenergy"),
        elif energies == None or fields == None:
            raise SyntaxError
        else:
            assert len(fields) == len(energies)
            jkext.hdf5.writeVLArray(self.__storage, "/" + state.hdfname(), "dcfield", fields)
            jkext.hdf5.writeVLArray(self.__storage, "/" + state.hdfname(), "dcstarkenergy", energies)


    def starkeffect_calculation(self, param):
        """Get all available energies from the given Starkeffect object and store them in our storage file."""
        if 'A' == param.type:
            for M in param.M:
                energies = {}
                for field in param.dcfields:
                    calc = jkext.starkeffect.AsymmetricRotor(param, M, field)
                    for state in calc.states():
                        id = state.id()
                        if energies.has_key(id):
                            energies[id].append(calc.energy(state))
                        else:
                            energies[id] = [calc.energy(state),]
                # store calculated values for this M
                for id in energies.keys():
                    self.starkeffect_merge(State().fromid(id), param.dcfields, energies[id])
        else:
            raise NotImplementedError("unknown rotor type in Stark energy calculation.")
        self.__storage.flush()


    def starkeffect_merge(self, state, newfields=None, newenergies=None):
        """Merge the specified pairs of field strength and Stark energies into the existing data.

        not really tested
        """
        assert len(newfields) == len(newenergies)
        try:
            oldfields, oldenergies = self.starkeffect(state)
            fields, energies = jkext.util.column_merge([oldfields, oldenergies], [newfields, newenergies])
        except tables.exceptions.NodeError:
            fields = newfields
            energies = newenergies
        self.starkeffect(state, fields, energies)


    def starkeffect_states(self):
        """Get a list of states for which we know the Stark effect."""
        list = []
        for group in self.__storage.listNodes(self.__storage.root, classname='Group'):
            state = State().fromhdfname(group._v_name)
            if 'dcfield' == group.dcfield.name and 'dcstarkenergy' == group.dcstarkenergy.name:
                list.append(state)
        return list


    def translate(self, translation):
        """Translate center of mass of molecule (i.e., translate all atoms)."""
        self.positions += translation



# some simple tests
if __name__ == "__main__":
    # test Stark calculation and storage/retrieval
    from convert import *
    param = jkext.starkeffect.CalculationParameter
    param.isomer = 0
    param.watson = 'A'
    param.symmetry = 'C2a'
    param.rotcon = Hz2J(num.array([5000e6, 1500e6, 1200e6]))
    param.quartic = Hz2J(num.array([50., 1000., 500, 10., 600]))
    param.dipole = D2Cm(num.array([5, 0., 0.]))
    # calculation details
    param.M = [0]
    param.Jmin = 0
    param.Jmax_calc = 10
    param.Jmax_save =  5
    param.dcfields = kV_cm2V_m(num.linspace(0., 100., 101))
    # save and print
    mol = Molecule(storage="molecule.hdf")
    mol.starkeffect_calculation(param)
    for J in range (0, 3):
        Ka = 0
        for Kc in range(J, -1, -1):
            state = State(J, Ka, Kc, 0, 0)
            fields, energies = mol.starkeffect(state)
            print state.name(), V_m2kV_cm(fields), J2Hz(energies) / 1e6
            if Kc > 0:
                Ka += 1
                state = State(J, Ka, Kc, 0, 0)
                fields, energies = mol.starkeffect(state)
                print state.name(), V_m2kV_cm(fields), J2Hz(energies) / 1e6
