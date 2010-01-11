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
from string import replace
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
            
    def __saveparam(self, param):
        """Store all relevant calculation parameters.

        """
        jkext.hdf5.writeVLArray(self.__storage, "/param" + "/_" + str(param.isomer) , "dipole", param.dipole)
        jkext.hdf5.writeVLArray(self.__storage, "/param" + "/_" + str(param.isomer) , "polarizability", \
                                param.polarizability,atom=tables.Float64Atom(shape=(3)))
        jkext.hdf5.writeVLArray(self.__storage, "/param" + "/_" + str(param.isomer) , "rotcon", param.rotcon)

        
    def __loadparam(self, param):
        """Retrieve stored calculation parameters.

        TODO: We might need to be more cleaver about the isomer stuff here. 
        """
        param.dipole=jkext.hdf5.readVLArray(self.__storage, "/param" + "/_" + str(param.isomer) + "/dipole")
        param.polarizability=jkext.hdf5.readVLArray(self.__storage, \
                                                    "/param/" + "/_" + str(param.isomer) + "/polarizability")
        param.rotcon=jkext.hdf5.readVLArray(self.__storage, "/param/" + "/_" + str(param.isomer) + "/rotcon")

    def getparam(self, param):
        """Retrieve stored calculation parameters.
        Non privat wrapper
        """
        self.__loadparam(param)

    def mueff(self, state):
        """Get the effective dipole moment \mu_eff as a function of the electric field strength.

        Return the effective dipole moment curve for the specified quantum |state|.
        """
        fields, energies = self.starkeffect(state,acfield=0.0)
        assert len(fields) == len(energies)
        mueff = num.zeros((len(fields),), num.float64)
        mueff[1:-1] = -1 * (energies[0:-2] - energies[2:]) / (fields[0:-2] - fields[2:])
        mueff[0] = 0.
        mueff[-1] = mueff[-2]
        return fields, mueff
    
    def coshellmann(self, state,param):
        """Get the the expectation value of cos theta using the Hellmann Feynman teorem.
        as a function of the electric field strength.
        this is right now only right for linar molecules
        this needs to be extended to different ac fields
        """
        dcfields, energies = self.starkeffect(state,acfield=0.0)
        omega=convert.dcfields2omega(dcfields,param.rotcon[1],param.dipole[0])
        acfields = [0.0]
        if len(acfields) > 1:
            energies = energies[:,0] # chose one ac field
        assert len(omega) == len(energies)
        cos = num.zeros((len(dcfields),), num.float64)
        cos[1:-1] = -1 * (energies[0:-2]/param.rotcon[1] - energies[2:]/param.rotcon[1]) / (omega[0:-2] - omega[2:])
        cos[0] = 0
        cos[-1] = cos[-2]
        return dcfields, cos



    def cos2hellmann(self, state,param):
        """Get the the expectation value of cos^2 theta using the Hellmann Feynman teorem.
        as a function of the electric field strength.
        this is right now only right for linar molecules
        this needs to be more robust and extented to different dc fields.
        """
        dcfields, energies, acfields = self.starkeffect(state)
        energies = energies[0,:] #only get the first energies
        assert len(acfields) == len(energies)
        cos2 = num.zeros((len(acfields),), num.float64)
        cos2[1:-1] = -(energies[0:-2]- energies[2:]+1/8*(param.polarizability[1,1]+param.polarizability[2,2])* \
        (acfields[0:-2]**2 - acfields[2:]**2)) / (1/4*(acfields[0:-2]**2 - acfields[2:]**2)*(param.polarizability[0,0]- \
        1/2*(param.polarizability[1,1]+param.polarizability[2,2])))
        cos2[0] = cos2[1]
        cos2[-1] = cos2[-2]
        return acfields, cos2

    def starkeffect(self, state, dcfields=None, acfields=None, energies=None):
        """Get or set the potential energies as a function of the electric field strength.
        
        When |energies| and |dcfields| are None, return the Stark curve for the specified quantum state.
        
        When |energies| and |fields| are specified, save the Stark curve for the specified quantum state in the
        Molecule's HDF5 storage file.
        We use the nameing convension _XdY where d replaces . in fx. 1.5 as . is not valid in a group identifier 
        """
        if energies == None and dcfields == None and acfields == None:
            #read all energies for a state not tested
            acfields = jkext.hdf5.readVLArray(self.__storage, "/" + state.hdfname() + "/acfields")
            for i in range(acfields):
                acfield = acfields[i]
                energies = jkext.hdf5.readVLArray(self.__storage, "/" + state.hdfname() + "/" + self.value2dir(acfield) + "dcstarkenergy")
                dcfields = jkext.hdf5.readVLArray(self.__storage, "/" + state.hdfname() + "/" + self.value2dir(acfield) + "dcfields")
            return dcfields, acfields, energies
        elif energies == None and dcfields == None and acfields != None:
            # read energies for a specific acfield
            return jkext.hdf5.readVLArray(self.__storage, "/" + state.hdfname() + "/" + self.value2dir(acfields) + "/dcfield"), \
                   jkext.hdf5.readVLArray(self.__storage, "/" + state.hdfname() + "/" + self.value2dir(acfields) + "/dcstarkenergy")
        elif energies == None and acfields == None and dcfields != None:
            # read energies for a specific acfield
            raise NotImplementedError("Energies for a specific dc field is not implemented yet")
        elif energies != None and dcfields != None and acfields != None:
            #we know all we need to know -> write it down
            # ToDo need to change this 2d energi matrix
            #assert len(dcfields)*len(acfields) == len(energies)
            for i in range(len(acfields)):
                jkext.hdf5.writeVLArray(self.__storage, "/" + state.hdfname()+ "/"  + self.value2dir(acfields[i]) , "dcfield", dcfields)
                jkext.hdf5.writeVLArray(self.__storage, "/" + state.hdfname()+ "/"  + self.value2dir(acfields[i]) , "dcstarkenergy", energies[i,:])
            jkext.hdf5.writeVLArray(self.__storage, "/" + state.hdfname(), "acfields", acfields)
        else:
            # maybe add a metode that gives the energy for a specified ac and dc field
            raise SyntaxError

    def value2dir(self, value):
        return "_"+replace(str(value),'.','d')
    
    def starkeffect_calculation(self, param):
        """Get all available energies from the given Starkeffect object and store them in our storage file."""
        if 'A' == param.type:
            for M in param.M:
                energies = {}
                for acfield in param.acfields:
                    for dcfield in param.dcfields:
                        calc = jkstark.starkeffect.AsymmetricRotor(param, M, dcfield,acfield)
                        for state in calc.states():
                            id = state.id()
                            if energies.has_key(id):
                                energies[id].append(calc.energy(state))
                            else:
                                energies[id] = [calc.energy(state),]
                # store calculated values for this M
                for id in energies.keys():
                    self.starkeffect_merge(State().fromid(id), param, energies[id])
        else:
            raise NotImplementedError("unknown rotor type in Stark energy calculation.")
        self.__storage.flush()

    def starkeffect_merge(self, state, param, newenergies=None):
        """Merge the specified pairs of field strength and Stark energies into the existing data

        not really tested
        TODO reimplement this will require a 2D matrix merge function need to handle 
        non uniform natrix sizes. apend in one or the other direction.
        TODO we need to improve the test for the exsistense of calculations at this ac field already !
        """
        newdcfields=param.dcfields
        newacfields=param.acfields
        assert len(newdcfields)*len(newacfields)  == len(newenergies)
        reshapedenergies = num.reshape(newenergies,(len(newacfields),len(newdcfields)))
        #for f in range(len(newacfields)):
        #    acfield = newacfields[f]
        #    try:
        #        olddcfields, oldenergies = self.starkeffect(state, acfields=acfield)
        #        dcfields, energies = jkext.util.column_merge([olddcfields, oldenergies], [newdcfields, newenergies])
        #    except tables.exceptions.NodeError:
        dcfields = newdcfields
        acfields = newacfields
        #energies = reshapedenergies[:,f]
        #assert len(energies) == len(dcfields)
        self.starkeffect(state, dcfields, acfields, reshapedenergies)

    def starkeffect_states(self):
        """Get a list of states for which we know the Stark effect.
        TODO adapt for ac fields """
        list = []
        for group in self.__storage.listNodes(self.__storage.root, classname='Group'):
            state = State().fromhdfname(group._v_name)
            if 'dcfield' == group.dcfield.name and 'dcstarkenergy' == group.dcstarkenergy.name:
                list.append(state)
        return list
    
    def acfields(self,state):
        """Get a list of ac fields for which we have calculated the stark effect"""
        return jkext.hdf5.readVLArray(self.__storage, "/" + state.hdfname() + "/acfields")


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
    param.acfields = num.linspace(0., 100., 5)
    # save and print
    mol = Molecule(storage="molecule.hdf")
    mol.starkeffect_calculation(param)
    for J in range (0, 3):
        Ka = 0
        for Kc in range(J, -1, -1):
            state = State(J, Ka, Kc, 0, 0)
            fields, energies = mol.starkeffect(state,acfields=0.0)
            print state.name()
            print V_m2kV_cm(fields)
            print J2Hz(energies) / 1e6
            if Kc > 0:
                Ka += 1
                state = State(J, Ka, Kc, 0, 0)
                fields, energies = mol.starkeffect(state,acfields=0.0)
                print state.name()
                print V_m2kV_cm(fields)
                print J2Hz(energies) / 1e6
