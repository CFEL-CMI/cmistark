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
import tables

import jkext.hdf5, jkext.molecule, jkext.util
from jkext.state import State

import jkstark.starkeffect



class Molecule(jkext.molecule.Molecule):
    """Representation of a Molecule"""

    def __init__(self, atoms=None, storage=None, name="Generic molecule"):
        """Create Molecule from a list of atoms."""
        jkext.molecule.Molecule.__init__(self, atoms, name)
        if storage != None:
            self.__storage = tables.openFile(storage, mode='a', title=name)
        else:
            self.__storage = None


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
                    calc = jkstark.starkeffect.AsymmetricRotor(param, M, field)
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



# some simple tests
if __name__ == "__main__":
    # test Stark calculation and storage/retrieval
    from jkext.convert import *
    param = jkstark.starkeffect.CalculationParameter
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
