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
__author__ = "Jochen Küpper <software@jochen-kuepper.de>"

import numpy as num
import numpy.linalg
import tables

import jkext.hdf5, jkext.molecule, jkext.util
from jkext.state import State

import cmistark.starkeffect


class _isomer_mass(tables.IsDescription):
    """Isomer mass represenation for pytables."""
    name  = tables.StringCol(64)
    num   = tables.UInt16Col()
    mass  = tables.Float64Col()


class Molecule(jkext.molecule.Molecule):
    """Representation of a Molecule"""

    def __init__(self, atoms=None, storage=None, name="Generic molecule", readonly=False):
        """Create Molecule from a list of atoms."""
        jkext.molecule.Molecule.__init__(self, atoms, name)
        try:
            if readonly:
                self.__storage = tables.openFile(storage, mode='r')
            else:
                self.__storage = tables.openFile(storage, mode='a', title=name)
<<<<<<< HEAD
            ## YP: the line below requires the stark file writable, block it for now for 
            ## readonly option
=======
            #the following line requires the stark file writable, block it for readonly option
>>>>>>> 5fc48e9b457f7af8013d0e3a2ebcc8240ac9451d
            #self.__storage.getNode("/")._v_title = name
        except:
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
        """Perform an Stark effect claculation, get all available energies from the given Starkeffect object, and store
        them in our storage file."""
        try:
            self.__storage.createTable("/", 'masses', _isomer_mass, "Isomer masses")
        except:
            pass

        if 'L' == param.type:
            Rotor = cmistark.starkeffect.LinearRotor
        elif 'S' == param.type:
            Rotor = cmistark.starkeffect.SymmetricRotor
        elif 'A' == param.type:
            Rotor = cmistark.starkeffect.AsymmetricRotor
        else:
            raise NotImplementedError("unknown rotor type in Stark energy calculation.")
        # calculate and store energies
        masses = self.__storage.root.masses
        new_isomer = True
        for isomer in masses.iterrows():
            if isomer['num'] == param.isomer:
                isomer['mass'] = param.mass
                new_isomer = False
        if new_isomer:
            isomer = self.__storage.root.masses.row
            isomer['name'] = param.name
            isomer['mass'] = param.mass
            isomer['num']  = param.isomer
            isomer.append()
        for M in param.M:
            energies = {}
            for field in param.dcfields:
                calc = Rotor(param, M, field)
                for state in calc.states():
                    id = state.id()
                    if energies.has_key(id):
                        energies[id].append(calc.energy(state))
                    else:
                        energies[id] = [calc.energy(state),]
            # store calculated values for this M
            for id in energies.keys():
                self.starkeffect_merge(State().fromid(id), param.dcfields, energies[id])
            # flush HDF5 file after every M
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
        for groupJ in self.__storage.listNodes(self.__storage.root, classname='Group'):
            for groupKa in self.__storage.listNodes(groupJ, classname='Group'):
                for groupKc in self.__storage.listNodes(groupKa, classname='Group'):
                    for groupM in self.__storage.listNodes(groupKc, classname='Group'):
                        for groupIso in self.__storage.listNodes(groupM, classname='Group'):
                            statename = (groupJ._v_name + '/' + groupKa._v_name + '/' + groupKc._v_name
                                         + '/' + groupM._v_name + '/' + groupIso._v_name)
                            if 'dcfield' == groupIso.dcfield.name and 'dcstarkenergy' == groupIso.dcstarkenergy.name:
                                list.append(State().fromhdfname(statename))
        return list


    def states_to_print(self, Jmin, Jmax, statelist=None):
        """Create a list of states to be printed/plotted according to the provided arguments

        Correctly creates list of states for the various rotor types
        """

        states = []
        return states



# some simple tests
if __name__ == "__main__":
    # test Stark calculation and storage/retrieval
    from jkext.convert import *
    param = cmistark.starkeffect.CalculationParameter
    param.name = 'cis'
    param.isomer = 0
    param.watson = 'A'
    param.symmetry = 'C2a'
    param.rotcon = Hz2J(num.array([5000e5, 1500e5, 1200e5]))
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
