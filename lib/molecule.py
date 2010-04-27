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
import jkstark.hdf5
from jkext.state import State
from string import replace
import jkstark.starkeffect
import jkstark.convert as convert



class Molecule(jkext.molecule.Molecule):
    """Representation of a Molecule"""

    def __init__(self, atoms=None, storage=None, name="Generic molecule", param=None):
        """Create Molecule from a list of atoms."""
        jkext.molecule.Molecule.__init__(self, atoms, name)
        if storage != None:
            self.__storage = tables.openFile(storage, mode='a', title=name)
        else:
            self.__storage = None
        if param != None:
            self.__saveparam(param)
            
    def __saveparam(self, param):
        """Store all relevant calculation parameters.

        """
        jkstark.hdf5.writeTable(self.__storage,jkstark.starkeffect.TableCalculationParameter,param)

        
    def __loadparam(self, param):
        """Retrieve stored calculation parameters.

        TODO: We might need to be more cleaver about the isomer stuff here. 
        """
        jkstark.hdf5.readTable(self.__storage,jkstark.starkeffect.TableCalculationParameter,param)

    def getparam(self, param):
        """Retrieve stored calculation parameters.
        Non privat wrapper
        """
        self.__loadparam(param)

    def geteigvectors(self, state, acfield, dcfield=None):
        """Retrieve components of eigenvector in asym top basis
        for a specific state and specific fields
        """
        eigvectors = jkext.hdf5.readVLArray(self.__storage, \
                        "/" + state.hdfname() + "/" + self.value2dir(acfield) + "/eigvectors")
        if dcfield != None:
            index = self.dcfieldindex(state,acfield,dcfield)
            eigvectors = eigvectors[index]
        return eigvectors 
        
    def mueff(self, state,acfields=0.0):
        """Get the effective dipole moment \mu_eff as a function of the electric field strength.

        Return the effective dipole moment curve for the specified quantum |state|.
        """
        fields, energies = self.starkeffect(state,acfields=acfields)
        assert len(fields) == len(energies)
        mueff = num.zeros((len(fields),), num.float64)
        mueff[1:-1] = -1 * (energies[0:-2] - energies[2:]) / (fields[0:-2] - fields[2:])
        mueff[0] = 0.
        mueff[-1] = mueff[-2]
        return fields, mueff

    def alphaeff(self, state,dcfields=0.0):
        """Get the effective dipole moment \mu_eff as a function of the electric field strength.

        Return the effective dipole moment curve for the specified quantum |state|.
        """
        acfields, energies = self.starkeffect(state,dcfields=dcfields)
        assert len(acfields) == len(energies)
        alphaeff = num.zeros((len(acfields),), num.float64)
        alphaeff[1:-1] = -(energies[0:-2] - energies[2:]) / (acfields[0:-2]**2 - acfields[2:]**2)
        alphaeff[0] = alphaeff[1]
        alphaeff[-1] = alphaeff[-2]
        return acfields, alphaeff
    
    def coshellmann(self, state, param, acfield):
        """Get the the expectation value of cos theta using the Hellmann Feynman teorem.
        as a function of the electric field strength.
        this is right now only right for linar molecules
        this needs to be extended to different ac fields
        """
        dcfields, energies = self.starkeffect(state, acfields = acfield)
        omega = convert.dcfields2omega(dcfields, param.rotcon[1], param.dipole[0])
        assert len(omega) == len(energies)
        cos = num.zeros((len(dcfields),), num.float64)
        cos[1:-1] = -1 * (energies[0:-2]/param.rotcon[1] - energies[2:]/param.rotcon[1]) / (omega[0:-2] - omega[2:])
        cos[0] = 0
        cos[-1] = cos[-2]
        return dcfields, cos

    def cos2hellmann(self, state, param, dcfield):
        """Get the the expectation value of cos^2 theta using the Hellmann Feynman teorem.
        as a function of the electric field strength.
        this is right now only right for linar molecules
        this needs to be more robust and extented to different dc fields.
        """
        acfields, energies = self.starkeffect(state, dcfields = dcfield)
        assert len(acfields) == len(energies)
        cos2 = num.zeros((len(acfields),), num.float64)
        cos2[1:-1] = -(energies[0:-2]- energies[2:]+1/8*(param.polarizability[1,1]+param.polarizability[2,2])* \
        (acfields[0:-2]**2 - acfields[2:]**2)) / (1/4*(acfields[0:-2]**2 - acfields[2:]**2)*(param.polarizability[0,0]- \
        1/2*(param.polarizability[1,1]+param.polarizability[2,2])))
        cos2[0] = cos2[1]
        cos2[-1] = cos2[-2]
        return acfields, cos2

    def starkeffect(self, state, dcfields = None, acfields = None, energies = None, eigvectors = None):
        """Get or set the potential energies as a function of the electric field strength.
        
        When |energies| and |dcfields| are None, return the Stark curve for the specified quantum state.
        
        When |energies| and |fields| are specified, save the Stark curve for the specified quantum state in the
        Molecule's HDF5 storage file.
        We use the nameing convension _XdY where d replaces . in fx. 1.5 as . is not valid in a group identifier 
        """
        if energies == None and dcfields == None and acfields == None:
            # read all energies for one state and return in one array
            # note that it assumes that all ac fields dirs contain the same number of dc energies
            acfields = jkext.hdf5.readVLArray(self.__storage, "/" + state.hdfname() + "/acfields")
            for acfield in acfields:
                temp = jkext.hdf5.readVLArray(self.__storage, "/" + state.hdfname() + "/" + self.value2dir(acfield) + "/dcstarkenergy")
                dcfields = jkext.hdf5.readVLArray(self.__storage, "/" + state.hdfname() + "/" + self.value2dir(acfield) + "/dcfield")
                assert len(temp) == len(dcfields)
                if acfield == acfields[0]:
                    energies = temp
                else:
                    energies = num.vstack((energies, temp))
            return dcfields, acfields, energies
        
        elif energies == None and dcfields == None and acfields != None:
            # read energies for a specific acfield
            dcfields = jkext.hdf5.readVLArray(self.__storage,\
                    "/" + state.hdfname() + "/" + self.value2dir(acfields) + "/dcfield")
            energies = jkext.hdf5.readVLArray(self.__storage, \
                    "/" + state.hdfname() + "/" + self.value2dir(acfields) + "/dcstarkenergy")
            assert len(dcfields) == len(energies)
            return dcfields, energies
        
        elif energies == None and acfields == None and dcfields != None:
            # read energies for a specific dcfield
            acfields = jkext.hdf5.readVLArray(self.__storage, "/" + state.hdfname() + "/acfields")
            i = 0
            energies  = num.zeros(len(acfields))
            for acfield in acfields:
                tempenergies = jkext.hdf5.readVLArray(self.__storage, \
                        "/" + state.hdfname() + "/" + self.value2dir(acfield) + "/dcstarkenergy")
                index = self.dcfieldindex(state,acfield,dcfields)
                energies[i] = tempenergies[index]
                i = i + 1
            return acfields, energies
        
        elif energies != None and dcfields != None and acfields != None:
            #we know all we need to know -> write it down
            #assert len(dcfields)*len(acfields) == len(energies)
            try:
                oldacfields = self.acfields(state)
                allacfields = jkext.util.column_merge([oldacfields], [acfields])
                allacfields = allacfields[0]
                # will return an array use [0] to take out of the array
            except tables.exceptions.NodeError:
                allacfields = acfields
            i=0
            for acfield in acfields:
                assert len(dcfields[self.value2dir(acfield)]) == len(energies[self.value2dir(acfield)])
                jkext.hdf5.writeVLArray(self.__storage, \
                        "/" + state.hdfname()+ "/"  + self.value2dir(acfield) ,\
                        "dcfield", dcfields[self.value2dir(acfield)])
                jkext.hdf5.writeVLArray(self.__storage, \
                        "/" + state.hdfname()+ "/"  + self.value2dir(acfield) ,\
                        "dcstarkenergy", energies[self.value2dir(acfield)])
                if eigvectors != None:
                    jkext.hdf5.writeVLArray(self.__storage, \
                        "/" + state.hdfname()+ "/"  + self.value2dir(acfield) ,\
                        "eigvectors", eigvectors[self.value2dir(acfield)],\
                        atom=tables.Float64Atom(shape=len(eigvectors[self.value2dir(acfield)][0])))
                i = i + 1
            jkext.hdf5.writeVLArray(self.__storage, "/" + state.hdfname(), "acfields", allacfields)
            
        else:
            # maybe add a metode that gives the energy for a specified ac and dc field
            raise SyntaxError

    def value2dir(self, value):
        value = replace(str(value),'+','p')
        value = replace(value,'.','d')
        return "_"+value
    
    def starkeffect_calculation(self, param):
        """Get all available energies from the given Starkeffect object and store them in our storage file."""
        if 'A' == param.type:
            for M in param.M:
                energies = {}
                eigvectors = {}
                for acfield in param.acfields:
                    for dcfield in param.dcfields:
                        calc = jkstark.starkeffect.AsymmetricRotor(param, M, dcfield,acfield)
                        for state in calc.states():
                            id = state.id()
                            if energies.has_key(id):
                                energies[id].append(calc.energy(state))
                                if param.saveevec == True:
                                    eigvectors[id].append(calc.eigvectors(state))
                            else:
                                energies[id] = [calc.energy(state),]
                                if param.saveevec == True:
                                    eigvectors[id] = [calc.eigvectors(state),]
                # store calculated values for this M
                for id in energies.keys():
                    if param.saveevec == True:
                        self.starkeffect_merge(State().fromid(id) , param, energies[id], eigvectors[id])
                    else:
                        self.starkeffect_merge(State().fromid(id), param, energies[id])
        else:
            raise NotImplementedError("unknown rotor type in Stark energy calculation.")
        self.__storage.flush()


    def starkeffect_merge(self, state, param, newenergies=None, neweigvectors=None):
        """Merge the specified pairs of field strength and Stark energies into the existing data

        not really tested
        TODO should work need better testing
        """
        newdcfields=param.dcfields
        newacfields=param.acfields
        assert len(newdcfields)*len(newacfields)  == len(newenergies)
        reshapedenergies = num.reshape(newenergies,(len(newacfields),len(newdcfields)))
        i=0
        energies = {} # we use a dict with acfields as entrys to allow different
        #amount of dc energies for different ac fields
        # i.e whem merging a few dc energies with many.
        dcfields = {}
        eigvectors = {}
        if param.saveevec == True:
            assert len(newdcfields)*len(newacfields)  == len(neweigvectors)
            reshapedeigvectors = num.reshape(num.array(neweigvectors),(len(newacfields),len(newdcfields),len(neweigvectors[0])))
        for acfield in newacfields:
            try: #try one acfield  ToDo add merge for eigenvectors
                olddcfields, oldenergies = self.starkeffect(state, acfields=acfield)
                oldeigenvectors = self.geteigvectors(state, acfield)
                newenergies = reshapedenergies[i,:]
                dcfield, energy = jkext.util.column_merge([olddcfields, oldenergies], [newdcfields, newenergies])
                if param.saveevec == True:
                    neweigenvectors = reshapedeigvectors[i,:,:]
                    dcfield, eigenvector = jkext.util.columnarray_merge([olddcfields, oldeigenvectors], [newdcfields, neweigenvectors])
            except tables.exceptions.NodeError: # no energies for this ac field
                dcfield = newdcfields
                energy = reshapedenergies[i,:]
                if param.saveevec == True:
                    eigenvector = reshapedeigvectors[i,:,:]
            assert len(dcfield) == len(energy)
            energies[self.value2dir(acfield)] = energy
            dcfields[self.value2dir(acfield)] = dcfield
            if param.saveevec == True:
                eigvectors[self.value2dir(acfield)] = eigenvector
            i = i + 1
        if param.saveevec == True:
            self.starkeffect(state, dcfields, newacfields, energies, eigvectors)
        else:
            self.starkeffect(state, dcfields, newacfields, energies)


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

    def dcfields(self,state,acfield):
        """Get a list of dc fields for which we have calculated the stark effect for a specific ac field
        ToDo make robust if ac field is not found"""
        acfield = self.testacfield(state,acfield)
        return jkext.hdf5.readVLArray(self.__storage, "/" + state.hdfname() + "/" + self.value2dir(acfield) +  "/dcfield")
    
    def testacfield(self,state,acfield):
        """return the closest ac field with calculated values
        print a warning if it is more than 1 % from the requested
        """
        acfields = jkext.hdf5.readVLArray(self.__storage, "/" + state.hdfname() + "/" + "/acfields")
        deltaacfields = abs(acfields - acfield)
        mindiff = min(deltaacfields)
        newacfield = acfields[deltaacfields == mindiff][0]
        pc = mindiff/acfield*100
        if pc > 1:
            print "The closest calculated ac field(%(newacfield)e) is %(pc)e percent from the requested (%(acfield)e)" % {'newacfield': newacfield, 'acfield': acfield, 'pc': pc}
        return newacfield
    
    def dcfieldindex(self,state,acfield,dcfield):
        """return the index of the closest dc field with calculated values
        print a warning if it is more than 1 % from the requested
        """
        acfield = self.testacfield(state,acfield)
        dcfields = jkext.hdf5.readVLArray(self.__storage,\
                        "/" + state.hdfname() + "/" + self.value2dir(acfield) + "/dcfield")
        deltadcfields = abs(dcfields - dcfield)
        index = deltadcfields.argmin()
        newdcfield = dcfields[index]
        mindiff = min(deltadcfields)
        pc = mindiff/dcfield*100
        if pc > 1:
            print "The closest calculated dc field(%(newdcfield)e) is %(pc)e percent from the requested (%(dcfield)e)" % {'newdcfield': newdcfield, 'dcfield': dcfield, 'pc': pc}
        return index
    
# some simple tests
if __name__ == "__main__":
    # test Stark calculation and storage/retrieval
    from jkext.convert import *
    param = jkstark.starkeffect.CalculationParameter
    param.isomer = 0
    param.watson = 'A'
    param.symmtry = 'C2a'
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
