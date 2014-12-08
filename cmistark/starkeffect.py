# -*- coding: utf-8; fill-column: 120; truncate-lines: t -*-
#
# This file is part of JK Python extensions
# Copyright (C) 2008,2009,2011,2012 Jochen Küpper <jochen.kuepper@cfel.de>
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

# really use scipy as numpy, so we are sure we use Fortran codes of eigvalsh and dgemm
import scipy as num
import scipy.linalg.blas

import cmiext.convert
from cmiext.state import State



class CalculationParameter(object):
    """Container of parameters for calculation of Stark energies plus some more generic parameters of the molecule

    Calculate energy for the specified ``dcfields`` (V/m) and rotor ``type``; all calculations are performed in
    representation Ir (x, y, z -> b, c, a).

    General parameters:
    :param isomer:
    :param rotcon: (Joule), quartic (Joule), dipole (Coulomb meter)
    :param mass: mass of molecule/isomer
    :param type: specify the type of rotor
    * 'L': linear top
    * 'S': prolate symmetric top
    * 'A': asymmetric top

    The following parameter are used for an asymmetric top:
    :param M:
    Jmin, Jmax_calc, Jmax_save
    :param watson: specifies which reduction of the centrifugal distortion constants of an asymmetric top shall be used.
    * 'A' for Watson's A reduction
    * 'S' for Watson's S reduction

    :param symmetry: defines the remaining symmetry of Hamiltonian for the molecule in a DC field. This is used to
      disentangle the block-diagonalization from a Wang transformation of the Hamiltonian matrix. It can be 'N', 'C2a',
      'C2b', 'C2c', 'V' for full Fourgroup symmetry, for asym top. 'W': block diagonalization in terms of E/O^+/- (Wang
      submatrices) for asym top. This can only be correct for zero-field calculations or M=0. For a symmetric top, the
      options are 'p' and 'o'.

    """
    name = ' '
    isomer = 0
    type = 'A'
    watson = None
    symmetry = 'N'
    # quantum numbers
    M = list(range(0, 2))
    Jmax_calc = 10
    Jmax_save = 6
    # fields
    dcfields = cmiext.convert.kV_cm2V_m(num.array((0, 100.), num.float64))
    # molecular parameters
    mass = num.zeros((1,), num.float64)      # kg
    rotcon = num.zeros((3,), num.float64)    # Joule - vector of length 1, 2, or 3 depending on type
    quartic = num.zeros((5,), num.float64)   # Joule - vector of length 1, 3, or 5 depending on type
    dipole = num.zeros((3,), num.float64)    # Coulomb meter - vector of length 1 or 3 depending on type
    # internal
    debug = None


class Rotor(object):
    """Abstract representation of an rotating top for energy level calculation purposes."""

    def __init__(self, param, M, dcfield=0.):
        """Save the relevant type-independent parameters"""
        ### general parameters
        self.complex = False
        self.hmat_type = num.float64
        self.type = param.type
        # save quantum numbers
        self.M = int(M) # use the single specified M
        self.isomer = int(param.isomer)
        self.Jmin = abs(self.M) # this must be equal to self.M (in Stark calculation all J couple)
        self.Jmax = int(param.Jmax_calc)
        self.Jmax_save = int(param.Jmax_save)
        # molecular constants
        self.rotcon = num.array(param.rotcon, num.float64)
        self.quartic = num.array(param.quartic, num.float64)
        self.dipole = num.array(param.dipole, num.float64)
        # field strengths
        self.dcfield = num.float64(dcfield)
        # symmetry of Hamiltonian (possible values: 'N', 'C2a', 'C2b', 'C2c', 'V', 'W' for asym rotor, 'p' and 'o' for sym rotor)
        self.symmetry = param.symmetry
        self.tiny = num.finfo(num.dtype(num.float64)).tiny * 10
        # we have not yet calculated the correct energies - mark invalid
        self.levels = {}
        self.levelssym = {}
        self.valid = False
        self.stateorder_valid = False
        self.debug = param.debug


    def field_DC(self):
        """Return DC field for which the Stark energies were calculated."""
        return self.dcfield


    def energy(self, state):
        """Return Stark energy for ``state``."""
        if self.valid == False:
            self.recalculate()
        return self.levels[state.id()]


    def statesymmetry(self, state):
        """Return symmetry for ``state``."""
        if self.valid == False:
            self.recalculate()
        return self.levelssym[state.id()]


    def print_mat(self, mat, text="", converter=None):
        """Print matrix for debuging purposes."""
        if converter is None:
            converter = lambda x: x
        print("\n", text)
        rows, columns = mat.shape
        for i in range(rows):
            for j in range(columns):
                if False == self.complex:
                    print("%10.3g" % (converter(mat[i,j])), end=' ')
                else:
                    print("%9.3gi" % (abs((mat[i,j]).real)+abs((mat[i,j]).imag), ), end=' ')
            print()



class LinearRotor(Rotor):
    """Representation of a linear top for energy level calculation purposes.

    This object will calculate rotational energies at the specified DC field strength for the given M-value and J-range.
    """

    def __init__(self, param, M, dcfield=0.):
        """Save the relevant parameters"""
        # general initialization
        Rotor.__init__(self, param, M, dcfield)
        # consistency checks
        assert 'L' == param.type.upper()
        assert self.rotcon.shape == (1,)
        assert self.dipole.shape == (1,)
        assert self.quartic.shape == (1,)


    def index(self, J):
        blockstart = J - self.Jmin
        return blockstart


    def recalculate(self):
        """Perform calculation of rotational state energies for current parameters"""
        hmat = self.hamiltonian(self.Jmin, self.Jmax, self.dcfield)
        if self.debug: self.print_mat(hmat, converter=cmiext.convert.J2Hz)
        eval = num.linalg.eigvalsh(hmat) # calculate only energies
        eval = num.sort(eval)
        if self.debug:
            eval, evec = num.linalg.eigh(hmat)
            print(eval)
            self.print_mat(evec)
        for J in range(self.Jmin, self.Jmax_save+1):
            i = J - self.Jmin
            state = State(J, 0, 0, self.M, self.isomer)
            self.levels[state.id()] = eval[i]
        # done - data is now valid
        self.valid = True


    def hamiltonian(self, Jmin, Jmax, dcfield):
        """Return Hamiltonian matrix"""
        matrixsize = Jmax - Jmin + 1
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
        mu = float(self.dipole)
        for J in range(Jmin, Jmax):
            value = -mu * dcfield * sqrt((J+1)**2) * sqrt((J+1)**2 - M**2) / ((J+1) * sqrt((2*J+1) * (2*J+3)))
            hmat[self.index(J+1), self.index(J)] += value
            hmat[self.index(J), self.index(J+1)] += value


    def states(self):
        """Return list of states for which the Stark energies were calculated."""
        list = []
        M = self.M
        iso = self.isomer
        for J in range(self.Jmin, self.Jmax_save+1):
            list.append(State(J, 0, 0, M, iso))
        return list




class SymmetricRotor(Rotor):
    """Representation of a symmetric top for energy level calculation purposes.

    This object will calculate rotational energies at the specified DC field strength for the given M-range, K-range and
    :math:`J`-range. Note that in the field, states corresponding to :math:`+KM` and :math:`-KM` become not degenerate.
    While always keeping :math:`M` positive in this program, we label states corresponding to :math:`-KM` by using
    negative K values in the output hdf files.

    """

    def __init__(self, param, M, dcfield=0.):
        """Save the relevant parameters"""
        Rotor.__init__(self, param, M, dcfield)
        assert 'S' == param.type.upper()
        self.symmetry = self.symmetry.lower()
        assert 'p' == self.symmetry or 'o' == self.symmetry
        self.valid = False
        self.stateorder_valid = False
        self.watson = param.watson
        assert self.rotcon.shape == (2,)
        assert self.dipole.shape == (1,)
        assert self.quartic.shape == (3,)

    def index(self, J, K):
        # The matrix size, Jmax - max({abs(K),Jmin}) + 1, is defined in hamiltonian.
        # The index starts from zero.
        return J - max({abs(K), self.Jmin})

    def recalculate(self):
        """Perform calculation of rotational state energies with current parameters for individual K"""
        self.levels = {}
        for K in range(-self.Jmax, self.Jmax+1): # scan K
            blocks = self.hamiltonian(self.Jmin, self.Jmax, self.dcfield, K) # create a full hamt for a single K. (note not |K|.)
            eval = num.linalg.eigvalsh(blocks) # calculate only energies
            eval = num.sort(eval)
            i = 0
            for state in self.stateorder(K):
                if state.J() <= self.Jmax_save:
                    self.levels[state.id()] = eval[i]
                i += 1
        # done - data is now valid
        self.valid = True

    def hamiltonian(self, Jmin, Jmax, dcfield, K):
        # The matrix size for a single K. the state labels for each dimension of the matrix: (|Jmin or K,K>,...,|Jmax-1,K>,|Jmax,K>)
        # The lower limit of the matrix is defined by Jmin or K (J cannot smaller than K).
        matrixsize = Jmax - max({abs(K),Jmin}) + 1
        # create hamiltonian matrix
        hmat = num.zeros((matrixsize, matrixsize), self.hmat_type)
        # start matrix with appropriate field-free rotor terms
        self.rigid(hmat, Jmin, Jmax, K)
        # fill matrix with appropriate Stark terms for nonzero fields
        if None != dcfield and self.tiny < abs(dcfield):
            self.stark_DC(hmat, Jmin, Jmax, K, dcfield)
        return hmat

    def rigid(self, hmat, Jmin, Jmax, K):
        """Add the rigid-rotor matrix element terms to hmat

        representation I for prolate, III for oblate symmetric tops
        (I(III)^l or ^r yield same results for sym tops.)

        Gordy & Cook (1984), section 6, and Zare (1988), section 6.3
        """
        DJ, DJK, DK = self.quartic.tolist()
        if 'p' == self.symmetry:
           AC, B = self.rotcon.tolist() #AC refers to A for p
        elif 'o' == self.symmetry:
           B, AC = self.rotcon.tolist() #AC refers to C for o
        for J in range(max({Jmin,abs(K)}), Jmax+1):
            rigid = B * J*(J+1) + (AC-B) * K**2
            distortion = -DJ * (J*(J+1))**2 - DJK * J*(J+1)*K**2 - DK * K**4
            hmat[self.index(J,K), self.index(J,K)] += rigid + distortion

    def stark_DC(self, hmat, Jmin, Jmax, K, dcfield):
        """Add the dc Stark-effect matrix element terms to hmat"""
        sqrt = num.sqrt
        M = self.M
        mu = float(self.dipole)
        for J in range(max({Jmin,abs(K)}),Jmax):
            # diagonal term
            if not (0 == M or 0 == K): # term would be zero; this also yields J !=0, so no division by zero possible
                hmat[self.index(J, K), self.index(J, K)] += -mu * dcfield * M * K / (J*(J+1))
                # off-diagonal term
            value = -mu * dcfield * sqrt((J+1)**2-K**2) * sqrt((J+1)**2 - M**2) / ((J+1) * sqrt((2*J+1) * (2*J+3)))
            hmat[self.index(J+1, K), self.index(J, K)] += value
            hmat[self.index(J, K), self.index(J+1, K)] += value
        # add last diagonal term
        J = Jmax
        if not (0 == M or 0 == K): # term would be zero; this also yields J !=0, so no division by zero possible
            hmat[self.index(J, K), self.index(J, K)] += -mu * dcfield * M * K / (J*(J+1))

    def states(self):
        """Return list of states for which the Stark energies were calculated."""
        list = []
        for J in range(self.Jmin, self.Jmax_save+1):
            for K in range(-J, J+1):
                if 'p' == self.symmetry:
                    list.append(State(J, K, 0, self.M, self.isomer))
                elif 'o' == self.symmetry:
                    list.append(State(J, 0, K, self.M, self.isomer))
        return list

    def stateorder(self, K):
        """Return a list with all states for the given K and the current calculation parameters (Jmin, Jmax)."""
        if False == self.stateorder_valid:
            self.stateorder_dict = []
            M = self.M
            iso = self.isomer
            for J in range(max(abs(K),self.Jmin), self.Jmax+1): # add states in an ascending order of energy
                if 'p' == self.symmetry:
                    self.stateorder_dict.append(State(J, K, 0, M, iso))
                elif 'o' == self.symmetry:
                    self.stateorder_dict.append(State(J, 0, K, M, iso))
        return self.stateorder_dict


class AsymmetricRotor(Rotor):
    """Representation of an asymmetric top for energy level calculation purposes.

    This object will calculate rotational energies at the specified DC field strength for the given M-value and J-range.
    """

    def __init__(self, param, M, dcfield=0.):
        """Save the relevant parameters"""
        Rotor.__init__(self, param, M, dcfield)
        assert 'A' == param.type.upper()
        self.valid = False
        self.stateorder_valid = False
        self.watson = param.watson
        assert self.rotcon.shape == (3,)
        assert self.dipole.shape == (3,)
        assert self.quartic.shape == (5,)
        self.dipole_components = [self.tiny < abs(self.dipole[0]),
                                  self.tiny < abs(self.dipole[1]),
                                  self.tiny < abs(self.dipole[2])]
        if True == self.dipole_components[2]: # µ_c != 0 -- the Hamiltonian matrix is complex (and hermitean)
            self.complex = True
            self.hmat_type = num.complex128
        else: # µ_c == 0 --  the Hamiltonian matrix is real (and symmetric)
            self.complex = False
            self.hmat_type = num.float64
        # For linear and symmetric top molecules this does not seem to be correct, therefore we disabled it for now.
        # Needs to be reimplemented correctly *soon*
        if 0 == self.M:
            if not self.dipole_components[1] and not self.dipole_components[2] and self.tiny < abs(self.rotcon[1] - self.rotcon[2]):
            # in representation(s) I the symmetry group of the Hamiltonian can be diagonalized into four Wang submatrices
            # even in a field if M == 0 and the dipole moment is along the a axis
                self.symmetry = 'Wa'
            elif not self.dipole_components[0] and not self.dipole_components[2]:
            # Wang submatrices coupling case for only u_b != 0
                self.symmetry = 'Wb'
            elif not self.dipole_components[0] and not self.dipole_components[1]:
            # Wang submatrices coupling case for only u_c != 0
                self.symmetry = 'Wc'
            elif self.dipole_components[0] != 0 and self.dipole_components[1] != 0:
            # Wang submatrices coupling case for only u_c = 0
                self.symmetry = 'Wab'
            elif self.dipole_components[1] != 0 and self.dipole_components[2] != 0:
            # Wang submatrices coupling case for only u_a = 0
                self.symmetry = 'Wbc'
            elif self.dipole_components[0] != 0 and self.dipole_components[2] != 0:
            # Wang submatrices coupling case for ile u_b = 0
                self.symmetry = 'Wac'
            elif self.dipole_components[0] != 0 and self.dipole_components[1] != 0 and self.dipole_components[2] != 0:
            # Wang submatrices coupling case for nonzero dipole moment components u_b and u_c (u_a can be zero or nonzero)
                self.symmetry = 'N'
            pass


    def states(self):
        """Return list of states for which the Stark energies were calculated."""
        list = []
        M = self.M
        iso = self.isomer
        for J in range(self.Jmin, self.Jmax_save+1):
            Ka = 0
            for Kc in range(J, -1, -1):
                list.append(State(J, Ka, Kc, M, iso))
                if Kc > 0:
                    Ka += 1
                    list.append(State(J, Ka, Kc, M, iso))
        return list


    def index(self, J, K):
        # this requires a correct "global" value of self.Jmin_matrixsize, which is set in hamiltonian.
        # Therefore, we must be called only through hamiltonian
        blockstart = J*(J-1) + J - self.Jmin_matrixsize
        return blockstart + K + J


    def recalculate(self):
        """Perform calculation of rotational state energies for current parameters"""
        self.levels = {}
        self.levelssym = {}
        blocks = self.hamiltonian(self.Jmin, self.Jmax, self.dcfield, self.symmetry)
        for symmetry in list(blocks.keys()):
            # if None is not self.debug: self.print_mat(blocks[symmetry], "\nSymmetry: " + symmetry)
            eval = num.linalg.eigvalsh(blocks[symmetry]) # calculate only energies
            i = 0
            for state in self.stateorder(symmetry):
                if state.J() <= self.Jmax_save:
                    #print "J,Ka,Kc,M, state.id()", state.J(), state.Ka(), state.Kc(), state.M(), state.id()
                    self.levels[state.id()] = eval[i]
                i += 1
        # done - data is now valid
        self.valid = True


    def hamiltonian(self, Jmin, Jmax, dcfield, symmetry):
        """Return block-diagonalized Hamiltonian matrix (blocks)"""
        self.Jmin_matrixsize = Jmin *(Jmin-1) + Jmin # this is used by index
        matrixsize = (Jmax + 1) * Jmax + Jmax + 1 - self.Jmin_matrixsize
        # create hamiltonian matrix
        hmat = num.zeros((matrixsize, matrixsize), self.hmat_type)
        # start matrix with appropriate field-free rigid-rotor terms
        self.rigid(hmat, Jmin, Jmax)
        # add appropriate field-free centrifugal distortion terms
        if self.watson == 'A':
            self.watson_A(hmat, Jmin, Jmax)
        elif self.watson == 'S':
            self.watson_S(hmat, Jmin, Jmax)
        else:
            assert self.watson == None
        if self.debug: self.print_mat(hmat, "\nField-free Hamiltonian:", converter=cmiext.convert.J2Hz)
        if self.debug: print("\nEnergies of the asym. rotor:\n", cmiext.convert.J2Hz(num.linalg.eigvalsh(hmat)))
        # fill matrix with appropriate Stark terms for nonzero fields
        if None != dcfield and self.tiny < abs(dcfield):
            self.stark_DC(hmat, Jmin, Jmax, dcfield)
        blocks = self.wang(hmat, symmetry, Jmin, Jmax)
        del hmat
        return blocks


    def rigid(self, hmat, Jmin, Jmax):
        """Add the rigid-rotor matrix element terms to hmat -- representation I^r

        Gordy & Cook, section 7, Table 7.2
        """
        sqrt = num.sqrt
        A, B, C = self.rotcon.tolist()
        for J in range(Jmin, Jmax+1):
            for K in range(-J, J+1):
                hmat[self.index(J, K), self.index(J, K)] += (B+C)/2 * (J*(J+1) - K**2) + A * K**2
            for K in range (-J, J-2+1):
                value = (B-C)/4 * sqrt((J*(J+1) - K*(K+1)) * (J*(J+1) - (K+1)*(K+2)))
                hmat[self.index(J, K+2), self.index(J, K)] += value
                hmat[self.index(J, K), self.index(J, K+2)] += value


    def stark_DC(self, hmat, Jmin, Jmax, dcfield):
        """Add the dc Stark-effect matrix element terms to hmat"""
        sqrt = num.sqrt
        M = self.M
        muA, muB, muC = self.dipole
        if self.dipole_components[0]:
            # matrix elements involving µ_a
            for J in range(Jmin, Jmax):
                for K in range(-J, J+1):
                    if 0 != M and 0 != K: # then also 0 != J
                        hmat[self.index(J, K), self.index(J, K)] += -muA * dcfield * M * K / (J*(J+1))
                    value = (-muA * dcfield * sqrt((J+1)**2 - K**2) * sqrt((J+1)**2 - M**2)
                              / ((J+1) * sqrt((2*J+1) * (2*J+3))))
                    hmat[self.index(J+1, K), self.index(J, K)] += value
                    hmat[self.index(J, K), self.index(J+1, K)] += value
            # final diagonal elements
            J = Jmax
            for K in range(-J, J+1):
                hmat[self.index(J, K), self.index(J, K)] += -1. * M * K / (J*(J+1)) * muA * dcfield
        if self.dipole_components[1]:
            # matrix elements involving µ_b
            for J in range(Jmin, Jmax):
                for K in range(-J, J+1):
                    if 0 != J:
                        value = -1 * M * muB * dcfield * (sqrt((J-K) * (J+K+1) ) ) / (2*J*(J+1))
                        hmat[self.index(J, K+1), self.index(J, K)] += value
                        hmat[self.index(J, K), self.index(J, K+1)] += value
                    # J+1, K+1 / J-1, K-1 case
                    value = (muB * dcfield * sqrt(((J+K+1) * (J+K+2)) * ((J+1)**2 - M**2))
                            / (2*(J+1) * sqrt((2*J+1) * (2*J+3))))
                    hmat[self.index(J+1, K+1), self.index(J, K)] += value
                    hmat[self.index(J, K), self.index(J+1, K+1)] += value
                    # J+1, K-1 / J-1, K+1 case
                    value = (-1 * muB * dcfield * sqrt(((J-K+1) * (J-K+2)) * ((J+1)**2 - M**2))
                              / (2*(J+1) * sqrt((2*J+1) * (2*J+3))))
                    hmat[self.index(J+1, K-1), self.index(J, K)] += value
                    hmat[self.index(J, K), self.index(J+1, K-1)] += value
        if  self.dipole_components[2]:
            # matrix elements involving µ_c
            for J in range(Jmin, Jmax):
                for K in range(-J, J+1):
                    if 0 != J:
                        value = 1j* M * muC * dcfield * sqrt((J-K) * (J+K+1)) / (2*J*(J+1))
                        hmat[self.index(J, K+1), self.index(J, K)] += value
                        hmat[self.index(J, K), self.index(J, K+1)] += -value #YP: change from value to -value, this - appears as i exists.
                    # J+1, K+1 / J-1, K-1 case
                    value = (-1j * muC * dcfield * sqrt((J+K+1) * (J+K+2)) * sqrt((J+1)**2 - M**2)
                              / (2*(J+1) * sqrt((2*J+1) * (2*J+3))))
                    hmat[self.index(J+1, K+1), self.index(J, K)] += value
                    hmat[self.index(J, K), self.index(J+1, K+1)] += value
                    # J+1, K-1 / J-1, K+1 case
                    value = (-1j  * muC * dcfield * sqrt((J-K+1) * (J-K+2)) * sqrt((J+1)**2 - M**2)
                              / (2*(J+1) * sqrt((2*J+1) * (2*J+3))))
                    hmat[self.index(J+1, K-1), self.index(J, K)] += value
                    hmat[self.index(J, K), self.index(J+1, K-1)] += value


    def stateorder(self, symmetry):
        """Return a list with all states for the given ``symmetry`` and the current calculation parameters (Jmin, Jmax).

        See Zare, 1988, Chapter 6, and publication for CMIStark [Chang2014]_

        """
        def Four_symmetry(J, Ka, Kc):
            """Determine Fourgroup symmetry of asymmetric top state

            see Zare, 1988, Eq 6.71"""
            if Ka%2 == 0 and Kc%2 == 0:   sym = 'A'   # ee - A
            elif Ka%2 == 0 and Kc%2 !=0:  sym = 'Ba'  # eo - Ba
            elif Ka%2 != 0 and Kc%2 ==0:  sym = 'Bc'  # oe - Bc
            elif Ka%2 != 0 and Kc%2 !=0:  sym = 'Bb'  # oo - Bb
            else: assert False
            return sym


        def Wang_submatrix(J, Ka, Kc):
            """Determine Wang submatrix of asymmetric top state in representation(s) I^r

            see Zare, 1988, Table 6.4"""
            if J%2 == 0:
                if Ka%2 == 0 and Kc%2 == 0:   sym = 'Epe'  # e0e - eee (K,s,J) - (Ka,Kc,J)
                elif Ka%2 == 0 and Kc%2 !=0:  sym = 'Eme'  # e1e - eoe
                elif Ka%2 != 0 and Kc%2 ==0:  sym = 'Ome'  # o1e - oee
                elif Ka%2 != 0 and Kc%2 !=0:  sym = 'Ope'  # o0e - ooe
                else: assert False
            elif J%2 != 0:
                if Ka%2 == 0 and Kc%2 == 0:   sym = 'Emo'  # e1o - eeo
                elif Ka%2 == 0 and Kc%2 !=0:  sym = 'Epo'  # e0o - eoo
                elif Ka%2 != 0 and Kc%2 ==0:  sym = 'Opo'  # o0o - oeo
                elif Ka%2 != 0 and Kc%2 !=0:  sym = 'Omo'  # o1o - ooo
                else: assert False
            else: assert False
            return sym
        if False == self.stateorder_valid:
            self.stateorder_dict = {}
            M = self.M
            iso = self.isomer
            if 'Wa' == self.symmetry or 'Wb' == self.symmetry or 'Wc' == self.symmetry or 'Wab' == self.symmetry or 'Wbc' == self.symmetry or 'Wac' == self.symmetry:
                eigenvalues = {'Epe': [], 'Eme': [], 'Ope': [], 'Ome': [], 'Epo': [], 'Emo': [], 'Opo': [], 'Omo': []}
                label = {'Epe': [], 'Eme': [], 'Ope': [], 'Ome': [], 'Epo': [], 'Emo': [], 'Opo': [], 'Omo': []}
            else:
                eigenvalues = {'A': [], 'Ba': [], 'Bb': [], 'Bc': []}
                label = {'A': [], 'Ba': [], 'Bb': [], 'Bc': []}
            for J in range(M, self.Jmax+1):
                Ka = 0
                for Kc in range(J,-1,-1):
                    if 'Wa' == self.symmetry or 'Wb' == self.symmetry or 'Wc' == self.symmetry or 'Wab' == self.symmetry or 'Wbc' == self.symmetry or 'Wac' == self.symmetry:
                        label[Wang_submatrix(J, Ka, Kc)].append(State(J, Ka, Kc, M, iso))
                    else:
                        label[Four_symmetry(J, Ka, Kc)].append(State(J, Ka, Kc, M, iso))
                    if Kc > 0:
                        Ka = Ka+1
                        if 'Wa' == self.symmetry or 'Wb' == self.symmetry or 'Wc' == self.symmetry or 'Wab' == self.symmetry or 'Wbc' == self.symmetry or 'Wac' == self.symmetry:
                            label[Wang_submatrix(J, Ka, Kc)].append(State(J, Ka, Kc, M, iso))
                        else:
                            label[Four_symmetry(J, Ka, Kc)].append(State(J, Ka, Kc, M, iso))
                # get block diagonal hamiltonian (make sure you calculate this in 'V'!)
                if 0 == J:
                    if 'Wa' == self.symmetry or 'Wb' == self.symmetry or 'Wc' == self.symmetry or 'Wab' == self.symmetry or 'Wbc' == self.symmetry or 'Wac' == self.symmetry:
                        blocks = {'Epe': num.zeros((1, 1), self.hmat_type)}
                    else:
                        blocks = {'A': num.zeros((1, 1), self.hmat_type)}
                else:
                    if 'Wa' == self.symmetry or 'Wb' == self.symmetry or 'Wc' == self.symmetry or 'Wab' == self.symmetry or 'Wbc' == self.symmetry or 'Wac' == self.symmetry:
                        blocks = self.hamiltonian(J, J, None, 'W')
                    else:
                        blocks = self.hamiltonian(J, J, None, 'V')
                # store sorted eigenenergies for respective J and block
                for sym in list(blocks.keys()):
                    if 0 < blocks[sym].size:
                        eigenvalues[sym] += num.sort(num.linalg.eigvalsh(num.array(blocks[sym]))).tolist()
            # sort assignments according to energy
            if 'Wa' == self.symmetry:
                symmetries = ['Ep', 'Em', 'Op', 'Om']
                eigenvalues['Ep'] = eigenvalues['Epe'] + eigenvalues['Epo']
                eigenvalues['Em'] = eigenvalues['Eme'] + eigenvalues['Emo']
                eigenvalues['Op'] = eigenvalues['Ope'] + eigenvalues['Opo']
                eigenvalues['Om'] = eigenvalues['Ome'] + eigenvalues['Omo']
                label['Ep'] = label['Epe'] + label['Epo']
                label['Em'] = label['Eme'] + label['Emo']
                label['Op'] = label['Ope'] + label['Opo']
                label['Om'] = label['Ome'] + label['Omo']
                del label['Epe'], label['Ome'], label['Eme'], label['Ope'], label['Epo'], label['Omo'], label['Emo'], label['Opo']
                del eigenvalues['Epe'], eigenvalues['Ome'], eigenvalues['Eme'], eigenvalues['Ope'], eigenvalues['Epo'], eigenvalues['Omo'], eigenvalues['Emo'], eigenvalues['Opo']
            elif 'Wb' == self.symmetry:
                symmetries = ['EmeOpo', 'OmeEpo', 'EpeOmo', 'OpeEmo']
                eigenvalues['EmeOpo'] = eigenvalues['Eme'] + eigenvalues['Opo']
                eigenvalues['OmeEpo'] = eigenvalues['Ome'] + eigenvalues['Epo']
                eigenvalues['EpeOmo'] = eigenvalues['Epe'] + eigenvalues['Omo']
                eigenvalues['OpeEmo'] = eigenvalues['Ope'] + eigenvalues['Emo']
                label['EmeOpo'] = label['Eme'] + label['Opo']
                label['OmeEpo'] = label['Ome'] + label['Epo']
                label['EpeOmo'] = label['Epe'] + label['Omo']
                label['OpeEmo'] = label['Ope'] + label['Emo']
                del label['Epe'], label['Ome'], label['Eme'], label['Ope'], label['Epo'], label['Omo'], label['Emo'], label['Opo']
                del eigenvalues['Epe'], eigenvalues['Ome'], eigenvalues['Eme'], eigenvalues['Ope'], eigenvalues['Epo'], eigenvalues['Omo'], eigenvalues['Emo'], eigenvalues['Opo']
            elif 'Wc' == self.symmetry:
                symmetries = ['EpeOpo', 'EmeOmo', 'OpeEpo', 'OmeEmo']
                eigenvalues['EpeOpo'] = eigenvalues['Epe'] + eigenvalues['Opo']
                eigenvalues['EmeOmo'] = eigenvalues['Eme'] + eigenvalues['Omo']
                eigenvalues['OpeEpo'] = eigenvalues['Ope'] + eigenvalues['Epo']
                eigenvalues['OmeEmo'] = eigenvalues['Ome'] + eigenvalues['Emo']
                label['EpeOpo'] = label['Epe'] + label['Opo']
                label['EmeOmo'] = label['Eme'] + label['Omo']
                label['OpeEpo'] = label['Ope'] + label['Epo']
                label['OmeEmo'] = label['Ome'] + label['Emo']
                del label['Epe'], label['Ome'], label['Eme'], label['Ope'], label['Epo'], label['Omo'], label['Emo'], label['Opo']
                del eigenvalues['Epe'], eigenvalues['Ome'], eigenvalues['Eme'], eigenvalues['Ope'], eigenvalues['Epo'], eigenvalues['Omo'], eigenvalues['Emo'], eigenvalues['Opo']
            elif 'Wab' == self.symmetry:
                symmetries = ['EpOm','EmOp']
                eigenvalues['EpOm'] = eigenvalues['Epe'] + eigenvalues['Epo'] + eigenvalues['Ome'] + eigenvalues['Omo']
                eigenvalues['EmOp'] = eigenvalues['Eme'] + eigenvalues['Emo'] + eigenvalues['Ope'] + eigenvalues['Opo']
                label['EpOm'] = label['Epe'] + label['Epo'] + label['Ome'] + label['Omo']
                label['EmOp'] = label['Eme'] + label['Emo'] + label['Ope'] + label['Opo']
                del label['Epe'], label['Ome'], label['Eme'], label['Ope'], label['Epo'], label['Omo'], label['Emo'], label['Opo']
                del eigenvalues['Epe'], eigenvalues['Ome'], eigenvalues['Eme'], eigenvalues['Ope'], eigenvalues['Epo'], eigenvalues['Omo'], eigenvalues['Emo'], eigenvalues['Opo']
            elif 'Wbc' == self.symmetry:
                symmetries = ['EeOo','EoOe']
                eigenvalues['EeOo'] = eigenvalues['Epe'] + eigenvalues['Eme'] + eigenvalues['Opo'] + eigenvalues['Omo']
                eigenvalues['EoOe'] = eigenvalues['Epo'] + eigenvalues['Emo'] + eigenvalues['Ope'] + eigenvalues['Ome']
                label['EeOo'] = label['Epe'] + label['Eme'] + label['Opo'] + label['Omo']
                label['EoOe'] = label['Epo'] + label['Emo'] + label['Ope'] + label['Ome']
                del label['Epe'], label['Ome'], label['Eme'], label['Ope'], label['Epo'], label['Omo'], label['Emo'], label['Opo']
                del eigenvalues['Epe'], eigenvalues['Ome'], eigenvalues['Eme'], eigenvalues['Ope'], eigenvalues['Epo'], eigenvalues['Omo'], eigenvalues['Emo'], eigenvalues['Opo']
            elif 'Wac' == self.symmetry:
                symmetries = ['EpOp','EmOm']
                eigenvalues['EpOp'] = eigenvalues['Epe'] + eigenvalues['Epo'] + eigenvalues['Ope'] + eigenvalues['Opo']
                eigenvalues['EmOm'] = eigenvalues['Eme'] + eigenvalues['Emo'] + eigenvalues['Ome'] + eigenvalues['Omo']
                label['EpOp'] = label['Epe'] + label['Epo'] + label['Ope'] + label['Opo']
                label['EmOm'] = label['Eme'] + label['Emo'] + label['Ome'] + label['Omo']
                del label['Epe'], label['Ome'], label['Eme'], label['Ope'], label['Epo'], label['Omo'], label['Emo'], label['Opo']
                del eigenvalues['Epe'], eigenvalues['Ome'], eigenvalues['Eme'], eigenvalues['Ope'], eigenvalues['Epo'], eigenvalues['Omo'], eigenvalues['Emo'], eigenvalues['Opo']
            elif 'V' == self.symmetry:
                symmetries = ['A', 'Ba', 'Bb', 'Bc']
            elif 'C2a' == self.symmetry:
                eigenvalues['Aa'] = eigenvalues['A'] + eigenvalues['Ba']
                eigenvalues['bc'] = eigenvalues['Bb'] + eigenvalues['Bc']
                label['Aa'] = label['A'] + label['Ba']
                label['bc'] = label['Bb'] + label['Bc']
                symmetries = ['Aa', 'bc']
                del label['A'], label['Ba'], label['Bb'], label['Bc']
                del eigenvalues['A'], eigenvalues['Ba'], eigenvalues['Bb'], eigenvalues['Bc']
            elif 'C2b' == self.symmetry:
                eigenvalues['Ab'] = eigenvalues['A'] + eigenvalues['Bb']
                eigenvalues['ac'] = eigenvalues['Bb'] + eigenvalues['Bc']
                label['Ab'] = label['A'] + label['Bb']
                label['ac'] = label['Ba'] + label['Bc']
                symmetries = ['Ab', 'ac']
                del label['A'], label['Ba'], label['Bb'], label['Bc']
                del eigenvalues['A'], eigenvalues['Ba'], eigenvalues['Bb'], eigenvalues['Bc']
            elif 'C2c' == self.symmetry:
                eigenvalues['Ac'] = eigenvalues['A'] + eigenvalues['Bc']
                eigenvalues['ab'] = eigenvalues['Ba'] + eigenvalues['Bb']
                label['Ac'] = label['A'] + label['Bc']
                label['ab'] = label['Ba'] + label['Bb']
                symmetries = ['Ac', 'ab']
                del label['A'], label['Ba'], label['Bb'], label['Bc']
                del eigenvalues['A'], eigenvalues['Ba'], eigenvalues['Bb'], eigenvalues['Bc']
            elif 'N' == self.symmetry:
                eigenvalues['N'] = eigenvalues['A'] + eigenvalues['Ba'] + eigenvalues['Bb'] + eigenvalues['Bc']
                label['N'] = label['A'] + label['Ba'] + label['Bb'] + label['Bc']
                symmetries = ['N']
                del label['A'], label['Ba'], label['Bb'], label['Bc']
                del eigenvalues['A'], eigenvalues['Ba'], eigenvalues['Bb'], eigenvalues['Bc']
            else:
                raise NotImplementedError("Hamiltonian symmetry %s not implemented" % (self.symmetry, ))
            for sym in symmetries:
                idx = num.argsort(eigenvalues[sym])
                self.stateorder_dict[sym] = num.array(label[sym])[idx]
            self.stateorder_valid = True
        return self.stateorder_dict[symmetry]


    def wang(self, hmat, symmetry, Jmin, Jmax):
        """Wang transform matrix and return a dictionary with the individual (sub)matrices."""
        matrixsize = ((Jmax + 1) * Jmax + Jmax + 1) - (Jmin *(Jmin-1) + Jmin)
        blocks = {}
        # set up Wang matrix
        Wmat = num.zeros(hmat.shape, self.hmat_type)
        value = 1/num.sqrt(2.)
        for J in range(Jmin, Jmax+1):
            for K in range(-J, 0):
                Wmat[self.index(J,  K), self.index(J,  K)] = -value
                Wmat[self.index(J, -K), self.index(J,  K)] = value
                Wmat[self.index(J,  K), self.index(J, -K)] = value
                Wmat[self.index(J, -K), self.index(J, -K)] = value
            Wmat[self.index(J, 0), self.index(J, 0)] = 1.
        # transform Hamiltonian matrix
        if self.complex:
            dot = lambda a, b: scipy.linalg.blas.cgemm(1., a, b)
        else:
            dot = lambda a, b: scipy.linalg.blas.dgemm(1., a, b)
        if None != self.debug: self.print_mat(hmat, "Original Hamiltonian", converter=cmiext.convert.J2Hz)
        hmat = dot(dot(Wmat, hmat), Wmat)
        if None != self.debug: self.print_mat(hmat, "Wang transformed Hamiltonian", converter=cmiext.convert.J2Hz)
        # delete Wang matrix (it's not used anymore)
        del Wmat
        # sort out matrix blocks
        if 'W' == symmetry:
            # use Wang submatrices E+/-,O+/- + e/o label for Jeven/Jodd (this e/o for J specifies V group symmetries at the same time)
            idx = {'Epe': [], 'Epo': [], 'Eme': [], 'Emo': [], 'Ope': [], 'Opo': [], 'Ome': [], 'Omo': []}
            i = 0
            for J in range(Jmin, Jmax+1):
                order = []
                if 0 == J % 2: # J even
                    for K in range(-J, 0): # K < 0 --> s odd
                        if 0 == K % 2: order.append('Eme') # K even
                        else: order.append('Ome') # K odd
                    for K in range(0, J+1): # K >= 0 --> s even
                        if 0 == K % 2: order.append('Epe') # K even
                        else: order.append('Ope') # K odd
                else: # J odd
                    for K in range(-J, 0): # K < 0 --> s odd
                        if 0 == K % 2: order.append('Emo') # K even
                        else: order.append('Omo') # K odd
                    for K in range(0, J+1): # K >= 0 --> s even
                        if 0 == K % 2: order.append('Epo') # K even
                        else: order.append('Opo') # K odd
                for k in range(2*J+1):
                    idx[order[k]].append(i+k)
                i += 2*J+1
            for sym in order:
                if 0 < len(idx[sym]):
                    blocks[sym] = hmat[num.ix_(idx[sym], idx[sym])]
        elif 'Wa' == symmetry:
            # use Wang submatrices E+/-,O+/- for only u_a!=0 and M=0 case
            idx = {'Ep': [], 'Em': [], 'Op': [], 'Om': []}
            i = 0
            for J in range(Jmin, Jmax+1):
                order = []
                if 0 == J % 2: # J even
                    for K in range(-J, 0): # K < 0 --> s odd
                        if 0 == K % 2: order.append('Em') # K even
                        else: order.append('Om') # K odd
                    for K in range(0, J+1): # K >= 0 --> s even
                        if 0 == K % 2: order.append('Ep') # K even
                        else: order.append('Op') # K odd
                else: # J odd
                    for K in range(-J, 0): # K < 0 --> s odd
                        if 0 == K % 2: order.append('Em') # K even
                        else: order.append('Om') # K odd
                    for K in range(0, J+1): # K >= 0 --> s even
                        if 0 == K % 2: order.append('Ep') # K even
                        else: order.append('Op') # K odd
                for k in range(2*J+1):
                    idx[order[k]].append(i+k)
                i += 2*J+1
            for sym in order:
                if 0 < len(idx[sym]):
                    blocks[sym] = hmat[num.ix_(idx[sym], idx[sym])]
        elif 'Wb' == symmetry:
            # for only u_b !=0 and M=0
            # the Stark element <J+1,K+/-1,M|H^b_Stark|J,K,M> couples:
            # 1. Om of even(odd) J and Ep of odd(even) J, 2. Op of even(odd)J and Em odd(even)J
            # four blocks are remained.
            idx = {'EmeOpo': [], 'OmeEpo': [], 'EpeOmo': [], 'OpeEmo': []} # rule of naming: K(e/o)s(p/m)J(e/o)K'(e/o)s'(p/m)J'(e/o)
            i = 0
            for J in range(Jmin, Jmax+1):
                order = []
                if 0 == J % 2: # J even
                    for K in range(-J, 0): # K < 0 --> s odd
                        if 0 == K % 2: order.append('EmeOpo') # K even
                        else: order.append('OmeEpo') # K odd
                    for K in range(0, J+1): # K >= 0 --> s even
                        if 0 == K % 2: order.append('EpeOmo') # K even
                        else: order.append('OpeEmo') # K odd
                else: # J odd
                    for K in range(-J, 0): # K < 0 --> s odd
                        if 0 == K % 2: order.append('OpeEmo') # K even
                        else: order.append('EpeOmo') # K odd
                    for K in range(0, J+1): # K >= 0 --> s even
                        if 0 == K % 2: order.append('OmeEpo') # K even
                        else: order.append('EmeOpo') # K odd
                for k in range(2*J+1):
                    idx[order[k]].append(i+k)
                i += 2*J+1
            for sym in order:
                if 0 < len(idx[sym]):
                    blocks[sym] = hmat[num.ix_(idx[sym], idx[sym])]
        elif 'Wc' == symmetry:
            # for u_c !=0 and u_b = 0 and M=0
            # the Stark element <J+1,K+/-1,M|H^b_Stark|J,K,M> couples:
            # 1. Op of even(odd) J and Ep of odd(even) J, 2. Om of even(odd) J and Em of odd(even) J
            idx = {'EpeOpo': [], 'EmeOmo': [], 'OpeEpo': [], 'OmeEmo': []}
            i = 0
            for J in range(Jmin, Jmax+1):
                order = []
                if 0 == J % 2: # J even
                    for K in range(-J, 0): # K < 0 --> s odd
                        if 0 == K % 2: order.append('EmeOmo') # K even
                        else: order.append('OmeEmo') # K odd
                    for K in range(0, J+1): # K >= 0 --> s even
                        if 0 == K % 2: order.append('EpeOpo') # K even
                        else: order.append('OpeEpo') # K odd
                else: # J odd
                    for K in range(-J, 0): # K < 0 --> s odd
                        if 0 == K % 2: order.append('OmeEmo') # K even
                        else: order.append('EmeOmo') # K odd
                    for K in range(0, J+1): # K >= 0 --> s even
                        if 0 == K % 2: order.append('OpeEpo') # K even
                        else: order.append('EpeOpo') # K odd
                for k in range(2*J+1):
                    idx[order[k]].append(i+k)
                i += 2*J+1
            for sym in order:
                if 0 < len(idx[sym]):
                    blocks[sym] = hmat[num.ix_(idx[sym], idx[sym])]
        elif 'Wab' == symmetry:
            # for only u_c = 0  and M=0
            # combine the coupling cases of Wa and Wb
            idx = {'EpOm': [], 'EmOp': []}
            i = 0
            for J in range(Jmin, Jmax+1):
                order = []
                if 0 == J % 2: # J even
                    for K in range(-J, 0): # K < 0 --> s odd
                        if 0 == K % 2: order.append('EmOp') # K even
                        else: order.append('EpOm') # K odd
                    for K in range(0, J+1): # K >= 0 --> s even
                        if 0 == K % 2: order.append('EpOm') # K even
                        else: order.append('EmOp') # K odd
                else: # J odd
                    for K in range(-J, 0): # K < 0 --> s odd
                        if 0 == K % 2: order.append('EmOp') # K even
                        else: order.append('EpOm') # K odd
                    for K in range(0, J+1): # K >= 0 --> s even
                        if 0 == K % 2: order.append('EpOm') # K even
                        else: order.append('EmOp') # K odd
                for k in range(2*J+1):
                    idx[order[k]].append(i+k)
                i += 2*J+1
            for sym in order:
                if 0 < len(idx[sym]):
                    blocks[sym] = hmat[num.ix_(idx[sym], idx[sym])]
        elif 'Wbc' == symmetry:
            # for only u_a = 0  and M=0
            # combine the coupling cases of Wb and Wc
            # symmetry 1: (K even, E, and J even, e) and (K  odd, O, and J odd, o)
            # symmetry 2: (K  odd, O, and J even, e) and (K even, E, and J odd, o)
            idx = {'EeOo': [], 'EoOe': []}
            i = 0
            for J in range(Jmin, Jmax+1):
                order = []
                if 0 == J % 2: # J even
                    for K in range(-J, 0): # K < 0 --> s odd
                        if 0 == K % 2: order.append('EeOo') # K even
                        else: order.append('EoOe') # K odd
                    for K in range(0, J+1): # K >= 0 --> s even
                        if 0 == K % 2: order.append('EeOo') # K even
                        else: order.append('EoOe') # K odd
                else: # J odd
                    for K in range(-J, 0): # K < 0 --> s odd
                        if 0 == K % 2: order.append('EoOe') # K even
                        else: order.append('EeOo') # K odd
                    for K in range(0, J+1): # K >= 0 --> s even
                        if 0 == K % 2: order.append('EoOe') # K even
                        else: order.append('EeOo') # K odd
                for k in range(2*J+1):
                    idx[order[k]].append(i+k)
                i += 2*J+1
            for sym in order:
                if 0 < len(idx[sym]):
                    blocks[sym] = hmat[num.ix_(idx[sym], idx[sym])]
        elif 'Wac' == symmetry:
            # for only u_b = 0  and M=0
            # combine the coupling cases of Wa and Wc
            idx = {'EpOp': [], 'EmOm': []}
            i = 0
            for J in range(Jmin, Jmax+1):
                order = []
                if 0 == J % 2: # J even
                    for K in range(-J, 0): # K < 0 --> s odd
                        if 0 == K % 2: order.append('EmOm') # K even
                        else: order.append('EmOm') # K odd
                    for K in range(0, J+1): # K >= 0 --> s even
                        if 0 == K % 2: order.append('EpOp') # K even
                        else: order.append('EpOp') # K odd
                else: # J odd
                    for K in range(-J, 0): # K < 0 --> s odd
                        if 0 == K % 2: order.append('EmOm') # K even
                        else: order.append('EmOm') # K odd
                    for K in range(0, J+1): # K >= 0 --> s even
                        if 0 == K % 2: order.append('EpOp') # K even
                        else: order.append('EpOp') # K odd
                for k in range(2*J+1):
                    idx[order[k]].append(i+k)
                i += 2*J+1
            for sym in order:
                if 0 < len(idx[sym]):
                    blocks[sym] = hmat[num.ix_(idx[sym], idx[sym])]
        elif 'V' == symmetry:
            # full Fourgroup symmetry (field free Hamiltonian or M=0!!!)
            # I^r representation, Wang transformed Hamiltonian factorizes into four submatrices E-, E+, O-, O+,
            # or, as used here, A, Ba, Bb, Bc
            # - in calculations for a single J this is the same
            # - in claculations for multiple J the correspondence flips with J
            idx = {'A': [], 'Ba': [], 'Bb': [], 'Bc': []}
            i = 0
            for J in range(Jmin, Jmax+1):
                order = []
                if 0 == J % 2: # J even
                    for K in range(-J, 0): # K < 0 --> s odd
                        if 0 == K % 2: order.append('Ba') # K even
                        else: order.append('Bc') # K odd
                    for K in range(0, J+1): # K >= 0 --> s even
                        if 0 == K % 2: order.append('A') # K even
                        else: order.append('Bb') # K odd
                else: # J odd
                    for K in range(-J, 0): # K < 0 --> s odd
                        if 0 == K % 2: order.append('A') # K even
                        else: order.append('Bb') # K odd
                    for K in range(0, J+1): # K >= 0 --> s even
                        if 0 == K % 2: order.append('Ba') # K even
                        else: order.append('Bc') # K odd
                for k in range(2*J+1):
                    idx[order[k]].append(i+k)
                i += 2*J+1
            for sym in order:
                if 0 < len(idx[sym]):
                    blocks[sym] = hmat[num.ix_(idx[sym], idx[sym])]
        elif 'C2a' == symmetry:
            # C2 rotation about a-axis is symmetry element
            #
            # I^r representation, Wang transformed Hamiltonian factorizes into two submatrices E = Aa (contains E+ and
            # E- / A and Ba) and O (contains O+ and O- / Bb and Bc).
            # In this case E and O corresponds to columns with K even and odd, respectively.
            idx = {'Aa': [], 'bc': []}
            if 0 == Jmin % 2: # Jmin even
                order = ['Aa', 'bc']
            else: # J odd
                order = ['bc', 'Aa']
            for i in range(matrixsize):
                idx[order[i%2]].append(i)
            for sym in order:
                if 0 < len(idx[sym]):
                    blocks[sym] = hmat[num.ix_(idx[sym], idx[sym])]
        elif 'C2b' == symmetry:
            # C2 rotation about b-axis is symmetry element
            #
            # I^r representation, Wang transformed Hamiltonian factorizes into two submatrices 'Ab' (contains 'A' and 'Bb')
            # and 'ac' (contains 'Ba' and 'Bc').
            idx = {'Ab': [], 'ac': []}
            i = 0
            for J in range(Jmin, Jmax+1):
                order = []
                if 0 == J % 2: # J even
                    for K in range(-J, 0): # K < 0 --> s odd
                        order.append('ac')
                    for K in range(0, J+1): # K >= 0 --> s even
                        order.append('Ab')
                else: # J odd
                    for K in range(-J, 0): # K < 0 --> s odd
                        order.append('Ab')
                    for K in range(0, J+1): # K >= 0 --> s even
                        order.append('ac')
                for k in range(2*J+1):
                    idx[order[k]].append(i+k)
                i += 2*J+1
            for sym in order:
                if 0 < len(idx[sym]):
                    blocks[sym] = hmat[num.ix_(idx[sym], idx[sym])]
        elif 'C2c' == symmetry:
            # C2 rotation about c-axis is symmetry element
            #
            # I^r representation, Wang transformed Hamiltonian factorizes into two submatrices 'Ac' (contains 'A' and
            # 'Bc') and 'ab' (contains 'Ba' and 'Bb').
            idx = {'Ac': [], 'ab': []}
            i = 0
            for J in range(Jmin, Jmax+1):
                order = []
                if 0 == J % 2: # J even
                    for K in range(-J, 0): # K < 0 --> s odd
                        if 0 == K % 2: order.append('ab') # K even
                        else: order.append('Ac') # K odd
                    for K in range(0, J+1): # K >= 0 --> s even
                        if 0 == K % 2: order.append('Ac') # K even
                        else: order.append('ab') # K odd
                else: # J odd
                    for K in range(-J, 0): # K < 0 --> s odd
                        if 0 == K % 2: order.append('Ac') # K even
                        else: order.append('ab') # K odd
                    for K in range(0, J+1): # K >= 0 --> s even
                        if 0 == K % 2: order.append('ab') # K even
                        else: order.append('Ac') # K odd
                for k in range(2*J+1):
                    idx[order[k]].append(i+k)
                i += 2*J+1
            for sym in order:
                if 0 < len(idx[sym]):
                    blocks[sym] = hmat[num.ix_(idx[sym], idx[sym])]
        elif 'N' == symmetry:
            # nothing to do, return
            blocks['N'] = hmat
        else:
            # something went wrong
            raise SyntaxError("unknown Hamiltonian symmetry")
        # for sym in set(order):
        #     for sym2 in set(order):
        #         if sym != sym2:
        #             if (hmat[num.ix_(idx[sym], idx[sym2])]!=0).any():
        #                 print "There is a problem with your symmetry"
        #                 print  sym, "and ", sym2, "are connected for M =", self.M
        # for symmetry in blocks.keys():
        #     self.print_mat(blocks[symmetry], "symmetry %s" % (symmetry)) # calculate only energies
        return blocks


    def watson_A(self, hmat, Jmin, Jmax):
        """Add the centrifugal distortion matrix element terms in Watson's A reduction to hmat."""
        matrixsize_Jmin = Jmin *(Jmin-1) + Jmin
        sqrt = num.sqrt
        DJ, DJK, DK, dJ, dK = self.quartic.tolist()
        for J in range(Jmin, Jmax+1):
            for K in range(-J, J+1):
                value = -DJ * (J*(J+1))**2 - DJK * J*(J+1)*K**2 - DK * K**4
                hmat[self.index(J, K), self.index(J, K)] += value
            for K in range(-J, J-2+1):
                value = ((-dJ * J*(J+1) - dK/2 * ((K+2)**2 + K**2))
                        * sqrt((J*(J+1) - K*(K+1)) * (J*(J+1) - (K+1)*(K+2))))
                hmat[self.index(J, K+2), self.index(J, K)] += value
                hmat[self.index(J, K), self.index(J, K+2)] += value


    def watson_S(self, hmat, Jmin, Jmax):
        """Add the centrifugal distortion matrix element terms in Watson's S reduction to hmat."""
        matrixsize_Jmin = Jmin *(Jmin-1) + Jmin
        sqrt = num.sqrt
        DJ, DJK, DK, dJ, dK = self.quartic.tolist()
        for J in range(Jmin, Jmax+1):
            for K in range(-J, J+1):
                value = -DJ * (J*(J+1))**2 - DJK * J*(J+1)*K**2 - DK * K**4
                hmat[self.index(J, K), self.index(J, K)] += value
            for K in range(-J, J-2+1):
                value = dJ * J*(J+1) * sqrt((J*(J+1) - K*(K+1)) * (J*(J+1) - (K+1)*(K+2)))
                hmat[self.index(J, K+2), self.index(J, K)] += value
                hmat[self.index(J, K), self.index(J, K+2)] += value
            for K in range(-J, J-4+1):
                value = dK * sqrt((J*(J+1) - K*(K+1)) * (J*(J+1) - (K+1)*(K+2))
                        * (J*(J+1)-(K+2)*(K+3)) * (J*(J+1)-(K+3)*(K+4)))
                hmat[self.index(J, K+4), self.index(J, K)] += value
                hmat[self.index(J, K), self.index(J, K+4)] += value


# some simple tests
if __name__ == "__main__":
    print()
