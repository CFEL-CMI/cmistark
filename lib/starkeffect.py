#!/usr/bin/env python
# -*- coding: utf-8; fill-column: 120 -*-
#
# This file is part of JK Python extensions
# Copyright (C) 2008,2009,2011,2012 Jochen Küpper <software@jochen-kuepper.de>
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

# really use scipy as numpy, so we are sure we use Fortran codes of eigvalsh and dgemm
import scipy as num
import scipy.linalg.fblas

import jkext.convert
from jkext.state import State


class CalculationParameter:
    """Container of parameters for calculation of Stark energies plus some more generic parameters of the molecule

    Calculate energy for the specified |dcfields| (V/m) and rotor type; all calculations are performed in representation
    Ir (x, y, z -> b, c, a).

    General parameters:
    - mass: mass of molecule/isomer
    - type: specify the type of rotor, currently only 'A' is implemented.
      - 'L': linear top
      - 'S': symmetric top
      - 'A': asymmetric top

    The following parameter are used for an asymmetric top:
    - M, Jmin, Jmax_calc, Jmax_save
    - isomer
    - rotcon (Joule), quartic (Joule), dipole (Coulomb meter)
    - watson, symmetry

    |watson| specifies which reduction of the centrifugal distortion constants of an asymmetric top shall be used.
    - 'A' for Watson's A reduction
    - 'S' for Watson's S reduction

    |symmetry| defines the remaining symmetry of Hamiltonian for the molecule in a DC field. This is used to disentangle
    the block-diagonalization from a Wang transformation of the Hamiltonian matrix. It can be 'N', 'C2a', 'C2b', 'C2c',
    or 'V' for full Fourgroup symmetry. The latter can only be correct for zero-field calculations or M=0.
    """
    type = 'A','L','S'
    M = range(0, 2)
    Jmax_calc = 5
    Jmax_save = 2
    isomer = 0
    # fields
    dcfields = jkext.convert.kV_cm2V_m(num.array((0, 100.), num.float64))
    # molecular parameters
    mass = num.zeros((1,), num.float64)      # kg
    rotcon = num.zeros((3,), num.float64)    # Joule - can of length 1, 2, or 3 depending on type
    quartic = num.zeros((5,), num.float64)   # Joule - can of length 1, 3, or 5 depending on type
    dipole = num.zeros((3,), num.float64)    # Coulomb meter
    polarizability = num.zeros((3,3), num.float64)
    watson=None
    symmetry='N'
    name = ' '


class Rotor(object):
    """Representation of an  linear top for energy level calculation purposes.

    This object will calculate rotational energies at the specified DC field strength for the given M-value and J-range.
    """

    def __init__(self, param, M, dcfield=0.):
        """Save the relevant type-independent parameters"""
        ### general parameters
        self.complex = False
        self.hmat_type = num.float64
        self.type = param.type
        # save quantum numbers
        self.M = int(M) # use the single specified M
        self.isomer = int(param.isomer)
        self.Jmin = self.M # this must be equal to self.M (in Stark calculation all J couple)
        self.Jmax = int(param.Jmax_calc)
        self.Jmax_save = int(param.Jmax_save)
        # molecular constants
        self.rotcon = num.array(param.rotcon, num.float64)
        self.quartic = num.array(param.quartic, num.float64)
        self.dipole = num.array(param.dipole, num.float64)
        # field strengths
        self.dcfield = num.float64(dcfield)
        # symmetry of Hamiltonian (possible values: 'N', 'C2a', 'C2b', 'C2c', 'V')
        self.symmetry = param.symmetry
        self.tiny = num.finfo(num.dtype(num.float64)).tiny * 10
        # we have not yet calculated the correct energies - mark invalid
        self.levels = {}
        self.valid = False
        self.stateorder_valid = False



    def field_DC(self):
        """Return DC field for which the Stark energies were calculated."""
        return self.dcfield



    def energy(self, state):
        """Return Stark energy for |state|."""
        if self.valid == False:
            self.recalculate()
        return self.levels[state.id()]



    def print_mat(self, mat, text=""):
        """Print matrix for debuging purposes."""
        print "\n", text
        rows, columns = mat.shape[0]
        for i in range(rows):
            for j in range(mat.columns):
                if False == self.complex:
                    print "%10.3g" % (mat[i,j]),
                else:
                    print "%9.3gi" % (abs((mat[i,j]).real)+abs((mat[i,j]).imag), ),
            print






class LinearRotor(Rotor):
    """Representation of an  linear top for energy level calculation purposes.

    This object will calculate rotational energies at the specified DC field strength for the given M-value and J-range.
    """

    def __init__(self, param, M, dcfield=0.):
        Rotor.__init__(self, param, M, dcfield=0.)
        """Save the relevant parameters"""
        assert 'L' == param.type
        # consistency checks
        assert self.rotcon.shape == (1,)
        assert self.dipole.shape == (1,)
        assert self.quartic.shape == (1,)


    def index(self, J):
        # this requires a correct "global" value of self.Jmin_matrixsize, which is set in full_hamiltonian.
        # Therefore, we must be called only through full_hamiltonian
        blockstart = J - self.Jmin_matrixsize
        return blockstart


    def recalculate(self):
        """Perform calculation of rotational state energies for current parameters"""
        hmat = self.hamiltonian(self.Jmin, self.Jmax, self.dcfield)
        eval = num.linalg.eigvalsh(hmat) # calculate only energies
        eval = num.sort(eval)
        for J in range(self.Jmin, self.Jmax_save+1):
            state = State(J, 0, 0, self.M, self.isomer)
            self.levels[state.id()] = eval[J]
        # done - data is now valid
        self.valid = True


    def hamiltonian(self, Jmin, Jmax, dcfield):
        """Return Hamiltonian matrix"""
        self.Jmin_matrixsize = Jmin *(Jmin-1) + Jmin # this is used by index
        matrixsize = (Jmax + 1) * Jmax + Jmax + 1 - self.Jmin_matrixsize
        # create hamiltonian matrix
        hmat = num.zeros((matrixsize, matrixsize), self.hmat_type)
        # start matrix with appropriate field-free rotor terms
        self.fieldfree(hmat, Jmin, Jmax)
        # fill matrix with appropriate Stark terms for nonzero fields
        if None != dcfield and self.tiny < abs(dcfield):
            self.stark_DC(hmat, Jmin, Jmax, dcfield)
        return hmat


    def fieldfree(self, hmat, Jmin, Jmax):
        """Add the field-free-rotor matrix element terms to hmat"""
        matrixsize_Jmin = Jmin *(Jmin-1) + Jmin
        sqrt = num.sqrt
        B = float(self.rotcon)
        D = float(self.quartic)
        for J in range(Jmin, Jmax+1):
            hmat[self.index(J), self.index(J)] += B * J*(J+1) - D * (J*(J+1))**2



    def stark_DC(self, hmat, Jmin, Jmax, dcfield):
        """Add the dc Stark-effect matrix element terms to hmat"""
        sqrt = num.sqrt
        M = self.M
        muA = self.dipole
        for J in range(Jmin, Jmax+1):
            value = (-muA * dcfield * sqrt((J+1)**2) * sqrt((J+1)**2 - M**2)
                      / ((J+1) * sqrt((2*J+1) * (2*J+3))))
            hmat[self.index(J+1), self.index(J)] += value
            hmat[self.index(J), self.index(J+1)] += value



    def states(self):
        """Return list of states for which the Stark energies were calculated."""
        list = []
        M = self.M
        iso = self.isomer
        for J in range(self.Jmin, self.Jmax_save+1):
            list.append(State(J, 0, 0, M, iso))
        return self.levels.keys()





# some simple tests
if __name__ == "__main__":
    print
    p = CalculationParameter
    p.Jmax_calc = 3
    p.Jmax_save = 2
    p.M = [0]
    p.isomer = 0
    p.rotcon = jkext.convert.Hz2J(num.array([5e9, 2e9, 1.4e9]))
    p.quartic = jkext.convert.Hz2J([1e3, 1e3, 1e3, 1e3, 1e3])
    p.dipole = jkext.convert.D2Cm([1.0, .0, .0])
    p.watson = 'A'
    p.symmetry = 'C2c'
    for M in p.M:
        for field in jkext.convert.kV_cm2V_m((1.,2.,3.)):
            print "\nM = %d, field strength = %.0f kV/cm" % (M, jkext.convert.V_m2kV_cm(field))
            top = Rotor(p, M, field)
            top.energy(State(1, 0, 1, M, p.isomer))
            for state in [State(0, 0, 0, M, p.isomer),
                          State(1, 0, 1, M, p.isomer), State(1, 1, 1, M, p.isomer), State(1, 1, 0, M, p.isomer),
                          State(2, 0, 2, M, p.isomer), State(2, 1, 2, M, p.isomer), State(2, 1, 1, M, p.isomer),
                          State(2, 2, 1, M, p.isomer), State(2, 2, 0, M, p.isomer),
                          State(3, 0, 3, M, p.isomer), State(3, 1, 3, M, p.isomer), State(3, 1, 2, M, p.isomer),
                          State(3, 2, 2, M, p.isomer), State(3, 2, 1, M, p.isomer), State(3, 3, 1, M, p.isomer),
                          State(3, 3, 0, M, p.isomer)]:
                if state.M() <= state.J() and state.J() <= p.Jmax_save:
                    print state.name(), "%12.3f MHz %8.3f cm-1 %10.3g J" \
                        % (jkext.convert.J2MHz(top.energy(state)), jkext.convert.J2invcm(top.energy(state)),
                           top.energy(state))
