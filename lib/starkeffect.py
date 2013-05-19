# -*- coding: utf-8; fill-column: 120 -*-
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
from __future__ import division
__author__ = "Jochen Küpper <jochen.kuepper@cfel.de>"

# really use scipy as numpy, so we are sure we use Fortran codes of eigvalsh and dgemm
import scipy as num
import scipy.linalg.blas

import jkext.convert
from jkext.state import State


class CalculationParameter(object):
    """Container of parameters for calculation of Stark energies plus some more generic parameters of the molecule

    Calculate energy for the specified |dcfields| (V/m) and rotor type; all calculations are performed in representation
    Ir (x, y, z -> b, c, a).

    General parameters:
    - isomer
    - rotcon (Joule), quartic (Joule), dipole (Coulomb meter)
    - mass: mass of molecule/isomer
    - type: specify the type of rotor
      - 'L': linear top
      - 'S': prolate symmetric top
      - 'A': asymmetric top

    The following parameter are used for an asymmetric top:
    - M, Jmin, Jmax_calc, Jmax_save
    - watson, specifies which reduction of the centrifugal distortion constants of an asymmetric top shall be used.
      - 'A' for Watson's A reduction
      - 'S' for Watson's S reduction
    - symmetry defines the remaining symmetry of Hamiltonian for the molecule in a DC field. This is used to disentangle
      the block-diagonalization from a Wang transformation of the Hamiltonian matrix. It can be 'N', 'C2a', 'C2b',
      'C2c', 'V' for full Fourgroup symmetry, for asym top. For sym top, the options are 'p' and 'o'.
      (added by YP in 2013) 'W': block diagonalization in terms of E/O^+/- (Wang submatrices) for asym top. The latter can only be
      correct for zero-field calculations or M=0.
    """
    name = ' '
    isomer = 0
    type = 'A'
    watson = None
    symmetry = 'N'
    # quantum numbers
    M = range(0, 2)
    Jmax_calc = 10
    Jmax_save = 6
    # fields
    dcfields = jkext.convert.kV_cm2V_m(num.array((0, 100.), num.float64))
    # molecular parameters
    mass = num.zeros((1,), num.float64)      # kg
    rotcon = num.zeros((3,), num.float64)    # Joule - vector of length 1, 2, or 3 depending on type
    quartic = num.zeros((5,), num.float64)   # Joule - vector of length 1, 3, or 5 depending on type
    dipole = num.zeros((3,), num.float64)    # Coulomb meter - vector of length 1 or 3 depending on type




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


    def field_DC(self):
        """Return DC field for which the Stark energies were calculated."""
        return self.dcfield


    def energy(self, state):
        """Return Stark energy for |state|."""
        if self.valid == False:
            self.recalculate()
        return self.levels[state.id()]

    def statesymmetry(self, state):
        """Return symmetry for |state|."""
        if self.valid == False:
            self.recalculate()
        return self.levelssym[state.id()]

    def print_mat(self, mat, text=""):
        """Print matrix for debuging purposes."""
        print "\n", text
        rows, columns = mat.shape
        for i in range(rows):
            for j in range(columns):
                if False == self.complex:
                    print "%10.3g" % (mat[i,j]),
                else:
                    print "%9.3gi" % (abs((mat[i,j]).real)+abs((mat[i,j]).imag), ),
            print



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
        eval = num.linalg.eigvalsh(hmat) # calculate only energies
        eval = num.sort(eval)
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

    This object will calculate rotational energies at the specified DC field strength for the given M-value and J-range.
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
	# this requires a correct "global" value of self.Jmin_matrixsize, which is set in hamiltonian.
	# Therefore, we must be called only through hamiltonian
	blockstart = J*(J-1) + J - self.Jmin_matrixsize
	return blockstart + K + J


    def recalculate(self):
        """Perform calculation of rotational state energies for current parameters"""
        self.levels = {}
        self.levelssym = {}
        blocks = self.hamiltonian(self.Jmin, self.Jmax, self.dcfield, self.symmetry)
        #if self.M == 0: #YP - print out eigenvalues of the full hmat, for debug
        #    blocksdebug = self.hamiltonian(self.Jmin, self.Jmax, self.dcfield, 'N') #YP - debug
        #    evaldebug = num.sort(num.linalg.eigvalsh(blocksdebug['N'])).tolist()
        #    print " ".join(str(x) for x in evaldebug) #YP - debug
        for symmetry in blocks.keys():
            eval = num.linalg.eigvalsh(blocks[symmetry]) # calculate only energies
            eval = num.sort(eval)
            i = 0
            for state in self.stateorder(symmetry):
                if state.J() <= self.Jmax_save:
                    self.levels[state.id()] = eval[i]
                    #self.levelssym[state.id()] = symmetry #YP: for debuging
                i += 1
        # done - data is now valid
        self.valid = True


    def hamiltonian(self, Jmin, Jmax, dcfield, symmetry):
        #"""Return Hamiltonian matrix"""
        self.Jmin_matrixsize = Jmin *(Jmin-1) + Jmin # this is used by index
        matrixsize = (Jmax + 1) * Jmax + Jmax + 1 - self.Jmin_matrixsize
        # create hamiltonian matrix
        hmat = num.zeros((matrixsize, matrixsize), self.hmat_type)
        # start matrix with appropriate field-free rotor terms
        self.rigid(hmat, Jmin, Jmax)
        # fill matrix with appropriate Stark terms for nonzero fields
        if None != dcfield and self.tiny < abs(dcfield):
            self.stark_DC(hmat, Jmin, Jmax, dcfield)
        blocks = self.wang(hmat, symmetry, Jmin, Jmax)
        del hmat
        return blocks
        #return hmat


    def rigid(self, hmat, Jmin, Jmax):
	"""Add the rigid-rotor matrix element terms to hmat -- representation I^l

	Gordy & Cook,
	"""
	DJ, DJK, DK = self.quartic.tolist()
        if 'p' == self.symmetry:
           AC, B = self.rotcon.tolist() #AC refers to A for p
        elif 'o' == self.symmetry:
           B, AC = self.rotcon.tolist() #AC refers to C for o
	for J in range(Jmin, Jmax+1):
	    for K in range(-J, J+1):
                rigid = B * J*(J+1) + (AC-B) * K**2
                distortion = -DJ * (J*(J+1))**2 - DJK * J*(J+1)*K**2 - DK * K**4
		hmat[self.index(J, K), self.index(J, K)] += rigid + distortion


    def stark_DC(self, hmat, Jmin, Jmax, dcfield):
        """Add the dc Stark-effect matrix element terms to hmat"""
        sqrt = num.sqrt
        M = self.M
        mu = float(self.dipole)
        for J in range(Jmin, Jmax):
            for K in range(-J, J+1):
                # diagonal term
                if not (0 == M or 0 == K): # term would be zero; this also yields J !=0, so no division by zero possible
                    hmat[self.index(J, K), self.index(J, K)] += -mu * dcfield * M * K / (J*(J+1))
                # off-diagonal term
                value = -mu * dcfield * sqrt((J+1)**2-K**2) * sqrt((J+1)**2 - M**2) / ((J+1) * sqrt((2*J+1) * (2*J+3)))
                hmat[self.index(J+1, K), self.index(J, K)] += value
                hmat[self.index(J, K), self.index(J+1, K)] += value
        # add last diagonal term
        J = Jmax
        for K in range(-J, J+1):
            # diagonal term
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

    def stateorder(self, symmetry):
        """Return a list with all states for the given |symmetry| and the current calculation parameters (Jmin, Jmax).

        See Gordy & Cook, Table 7.5.

        The symmetry of asymmetric rotor functions (in terms of eveness of Ka and Kc) is of course independent of the
        representation used in the calculation.
        """

        def Four_symmetry(J, K):
            """Determine Fourgroup symmetry of the corresponding Wang transformed symmetric top state of
            a symmetric top state JK, in representation(s) I for prolate top, III for oblate top

            see Gordy & Cook (1984), Table 7.5 or Allen & Cross (1963), Table 2n2"""
            if self.symmetry == 'p':
                if 0 == J % 2: # J even
                    if K < 0: # K > 0 --> s odd
                        if 0 == K % 2: sym = 'Ba' # K even
                        else: sym = 'Bc' # K odd
                    if K >= 0: # K <= 0 --> s even
                        if 0 == K % 2: sym = 'A' # K even
                        else: sym = 'Bb' # K odd
                else: # J odd
                    if K < 0: # K <= 0 --> s even
                        if 0 == K % 2: sym = 'A' # K even
                        else: sym = 'Bb' # K odd
                    if K >= 0: # K >= 0 --> s odd
                        if 0 == K % 2: sym = 'Ba' # K even
                        else: sym = 'Bc' # K odd
            if self.symmetry == 'o':
                if 0 == J % 2: # J even
                    if K < 0: # K > 0 --> s odd
                        if 0 == K % 2: sym = 'Bc' # K even
                        else: sym = 'Bb' # K odd
                    if K >= 0: # K <= 0 --> s even
                        if 0 == K % 2: sym = 'A' # K even
                        else: sym = 'Ba' # K odd
                else: # J odd
                    if K < 0: # K <= 0 --> s even
                        if 0 == K % 2: sym = 'A' # K even
                        else: sym = 'Ba' # K odd
                    if K >= 0: # K >= 0 --> s odd
                        if 0 == K % 2: sym = 'Bc' # K even
                        else: sym = 'Bb' # K odd
            return sym

        if False == self.stateorder_valid:
            self.stateorder_dict = {}
            M = self.M
            iso = self.isomer
            eigenvalues = {'A': [], 'Ba': [], 'Bb': [], 'Bc': [], 'Ac': [], 'ab': [], 'Aa': [], 'bc': []}
            label = {'A': [], 'Ba': [], 'Bb': [], 'Bc': [], 'Ac': [], 'ab': [], 'Aa': [], 'bc': []}
            for J in range(abs(M), self.Jmax+1):
                if 0 == J:
                    blocks = {'A': num.zeros((1, 1), self.hmat_type)}
                else:
                    #hmat = self.hamiltonian(J, J, None)
                    blocks = self.hamiltonian(J, J, None, 'V')
                if 'p' == self.symmetry:
                    tempAa = []
                    tempbc = []
                    for K in range(0,J+1): #YP: the orders of for and if loops has to be fixed this way so that the energy order without the field (E(Klarge)>E(Ksmall)) and in the field (E(KM<0)>E(KM>0)) is correct
                        if Four_symmetry(J, K) == 'A':
                            label['Aa'].append(State(J, K, 0, M, iso))
                        elif Four_symmetry(J, K) == 'Bc':
                            label['bc'].append(State(J, K, 0, M, iso))
                        elif Four_symmetry(J, K) == 'Ba':
                            label['Aa'].append(State(J, K, 0, M, iso))
                        elif Four_symmetry(J, K) == 'Bb':
                            label['bc'].append(State(J, K, 0, M, iso))
                        if K !=0 :
                            if Four_symmetry(J, -K) == 'A':
                                label['Aa'].append(State(J, -K, 0, M, iso))
                            elif Four_symmetry(J, -K) == 'Bc':
                                label['bc'].append(State(J, -K, 0, M, iso))
                            elif Four_symmetry(J, -K) == 'Ba':
                                label['Aa'].append(State(J, -K, 0, M, iso))
                            elif Four_symmetry(J, -K) == 'Bb':
                                label['bc'].append(State(J, -K, 0, M, iso))
                    for sym in blocks.keys():
                        if 0 < blocks[sym].size:
                            if sym == 'A':
                                tempAa += num.linalg.eigvalsh(num.array(blocks[sym])).tolist()
                            elif sym == 'Bc':
                                tempbc += num.linalg.eigvalsh(num.array(blocks[sym])).tolist()
                            elif sym == 'Ba':
                                tempAa += num.linalg.eigvalsh(num.array(blocks[sym])).tolist()
                            elif sym == 'Bb':
                                tempbc += num.linalg.eigvalsh(num.array(blocks[sym])).tolist()
                    eigenvalues['Aa'] += num.sort(tempAa).tolist()
                    eigenvalues['bc'] += num.sort(tempbc).tolist()
                elif 'o' == self.symmetry:
                    tempAc = []
                    tempab = []
                    for K in range(J,-1,-1): #YP: the orders of for and if loops has to be fixed this way so that the energy order without the file (E(Klarge)<E(Ksmall)) in the field (E(KM<0)>E(KM>0)) is correct
                        if Four_symmetry(J, K) == 'A':
                            label['Ac'].append(State(J, 0, K, M, iso))
                        elif Four_symmetry(J, K) == 'Bc':
                            label['Ac'].append(State(J, 0, K, M, iso))
                        elif Four_symmetry(J, K) == 'Ba':
                            label['ab'].append(State(J, 0, K, M, iso))
                        elif Four_symmetry(J, K) == 'Bb':
                            label['ab'].append(State(J, 0, K, M, iso))
                        if K !=0 :
                            if Four_symmetry(J, -K) == 'A':
                                label['Ac'].append(State(J, 0, -K, M, iso))
                            elif Four_symmetry(J, -K) == 'Bc':
                                label['Ac'].append(State(J, 0, -K, M, iso))
                            elif Four_symmetry(J, -K) == 'Ba':
                                label['ab'].append(State(J, 0, -K, M, iso))
                            elif Four_symmetry(J, -K) == 'Bb':
                                label['ab'].append(State(J, 0, -K, M, iso))
                    for sym in blocks.keys():
                        if 0 < blocks[sym].size:
                            if sym == 'A':
                                tempAc += num.linalg.eigvalsh(num.array(blocks[sym])).tolist()
                            elif sym == 'Bc':
                                tempAc += num.linalg.eigvalsh(num.array(blocks[sym])).tolist()
                            elif sym == 'Ba':
                                tempab += num.linalg.eigvalsh(num.array(blocks[sym])).tolist()
                            elif sym == 'Bb':
                                tempab += num.linalg.eigvalsh(num.array(blocks[sym])).tolist()
                    eigenvalues['Ac'] += num.sort(tempAc).tolist()
                    eigenvalues['ab'] += num.sort(tempab).tolist()

            # sort assignments according to energy
            if 'V' == self.symmetry:
                symmetries = ['A', 'Ba', 'Bb', 'Bc']
            elif 'p' == self.symmetry: # equivalent to 'C2a' in asym
                symmetries = ['Aa', 'bc']
                del eigenvalues['A'], eigenvalues['Ba'], eigenvalues['Bb'], eigenvalues['Bc']
            elif 'o' == self.symmetry: # equivalent to 'C2c' in asym
                symmetries = ['Ac', 'ab']
                del eigenvalues['A'], eigenvalues['Ba'], eigenvalues['Bb'], eigenvalues['Bc']
            elif 'N' == self.symmetry:
                for sym in blocks.keys():
                    eigenvalues['N'] += eigenvalues[sym]
                    label['N'] += label[sym]
                    del label[sym]
                    del eigenvalues[sym]
                symmetries = ['N']
            else:
                raise NotImplementedError("Hamiltonian symmetry %s not implemented" % (self.symmetry, ))
            total_label = []
            total_energy = []
            for sym in symmetries:
                idx = num.argsort(eigenvalues[sym],kind='mergesort')
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
        #self.print_mat(hmat, "Original Hamiltonian")
        hmat = dot(dot(Wmat, hmat), Wmat)
        #self.print_mat(hmat, "Wang transformed Hamiltonian")
        # delete Wang matrix (it's not used anymore)
        del Wmat
        # sort out matrix blocks
        if 'V' == symmetry:
            # full Fourgroup symmetry (field free Hamiltonian or M=0!!!)
            # I^r (not I^l?) representation for prolate, III for oblate, Wang transformed Hamiltonian factorizes into four submatrices E-, E+, O-, O+,
            # or, as used here, A, Ba, Bb, Bc
            # - in calculations for a single J this is the same
            # - in claculations for multiple J the correspondence flips with J (see Gordy+Cook Table 7.5)
            idx = {'A': [], 'Ba': [], 'Bb': [], 'Bc': []}
            i = 0
            for J in range(Jmin, Jmax+1):
                order = []
                if self.symmetry == 'p':
                    if 0 == J % 2: # J even
                        for K in range(-J, 0): # K > 0 --> s odd
                            if 0 == K % 2: order.append('Ba') # K even
                            else: order.append('Bc') # K odd
                        for K in range(0, J+1): # K <= 0 --> s even
                            if 0 == K % 2: order.append('A') # K even
                            else: order.append('Bb') # K odd
                    else: # J odd
                        for K in range(-J, 0): # K <= 0 --> s even
                            if 0 == K % 2: order.append('A') # K even
                            else: order.append('Bb') # K odd
                        for K in range(0, J+1): # K >= 0 --> s odd
                            if 0 == K % 2: order.append('Ba') # K even
                            else: order.append('Bc') # K odd
                elif self.symmetry == 'o':
                    if 0 == J % 2: # J even
                        for K in range(-J, 0): # K > 0 --> s odd
                            if 0 == K % 2: order.append('Bc') # K even
                            else: order.append('Bb') # K odd
                        for K in range(0, J+1): # K <= 0 --> s even
                            if 0 == K % 2: order.append('A') # K even
                            else: order.append('Ba') # K odd
                    else: # J odd
                        for K in range(-J, 0): # K <= 0 --> s even
                            if 0 == K % 2: order.append('A') # K even
                            else: order.append('Ba') # K odd
                        for K in range(0, J+1): # K >= 0 --> s odd
                            if 0 == K % 2: order.append('Bc') # K even
                            else: order.append('Bb') # K odd
                for k in range(2*J+1):
                    idx[order[k]].append(i+k)
                i += 2*J+1
            for sym in order:
                if 0 < len(idx[sym]):
                    blocks[sym] = hmat[num.ix_(idx[sym], idx[sym])]
        elif 'p' == symmetry: # 'C2a'
            # C2 rotation about a-axis is symmetry element fo prolate tops
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
        elif 'o' == symmetry:
            # C2 rotation about c-axis is symmetry element
            #
            # I^r representation, Wang transformed Hamiltonian factorizes into two submatrices 'Ac' (contains 'A' and
            # 'Bc') and 'ab' (contains 'Ba' and 'Bb').
            idx = {'Ac': [], 'ab': []}
            i = 0
            for J in range(Jmin, Jmax+1):
                order = []
                if 0 == J % 2: # J even
                    for K in range(-J, 0): # K > 0 --> s odd
                        if 0 == K % 2: order.append('Ac') # K even
                        else: order.append('ab') # K odd
                    for K in range(0, J+1): # K <= 0 --> s even
                        if 0 == K % 2: order.append('Ac') # K even
                        else: order.append('ab') # K odd
                else: # J odd
                    for K in range(-J, 0): # K <= 0 --> s even
                        if 0 == K % 2: order.append('Ac') # K even
                        else: order.append('ab') # K odd
                    for K in range(0, J+1): # K >= 0 --> s odd
                        if 0 == K % 2: order.append('Ac') # K even
                        else: order.append('ab') # K odd
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
	if 0 == self.M \
                and not self.dipole_components[1] and not self.dipole_components[2] \
                and self.tiny < abs(self.rotcon[1] - self.rotcon[2]):
	    # in representation(s) I the symmetry group of the Hamiltonian is V even in a field if M == 0 and the dipole
	    # moment is along a
	    self.symmetry = 'W'
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
	for symmetry in blocks.keys():
	    eval = num.linalg.eigvalsh(blocks[symmetry]) # calculate only energies
	    eval = num.sort(eval)
	    i = 0
	    for state in self.stateorder(symmetry):
		if state.J() <= self.Jmax_save:
		    self.levels[state.id()] = eval[i]
                    self.levelssym[state.id()] = symmetry #YP: for debuging 
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
	# fill matrix with appropriate Stark terms for nonzero fields
	if None != dcfield and self.tiny < abs(dcfield):
	    self.stark_DC(hmat, Jmin, Jmax, dcfield)
	blocks = self.wang(hmat, symmetry, Jmin, Jmax)
	del hmat
	return blocks


    def rigid(self, hmat, Jmin, Jmax):
	"""Add the rigid-rotor matrix element terms to hmat -- representation I^l

	Gordy & Cook,
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
	"""Return a list with all states for the given |symmetry| and the current calculation parameters (Jmin, Jmax).

	See Gordy & Cook, Table 7.5.

	The symmetry of asymmetric rotor functions (in terms of eveness of Ka and Kc) is of course independent of the
	representation used in the calculation.
	"""
	def Four_symmetry(J, Ka, Kc):
	    """Determine Fourgroup symmetry of asymmetric top state in representation(s) I

	    see Gordy & Cook (1984), Table 7.5 or Allen & Cross (1963), Table 2n2"""
	    if Ka%2 == 0 and Kc%2 == 0:   sym = 'A'   # ee
	    elif Ka%2 == 0 and Kc%2 !=0:  sym = 'Ba'  # eo
	    elif Ka%2 != 0 and Kc%2 ==0:  sym = 'Bc'  # oe
	    elif Ka%2 != 0 and Kc%2 !=0:  sym = 'Bb'  # oo
	    else: assert False
	    return sym

        def Wang_subm(J, Ka, Kc):
            """Determine Wang submatrice of asymmetric top state in representation(s) I

            see Gordy & Cook (1984), Table 7.5 or Allen & Cross (1963), Table 2n2"""
            if J%2 == 0:
                if Ka%2 == 0 and Kc%2 == 0:   sym = 'Ep'  # ee
                elif Ka%2 == 0 and Kc%2 !=0:  sym = 'Em'  # eo
                elif Ka%2 != 0 and Kc%2 ==0:  sym = 'Op'  # oe
                elif Ka%2 != 0 and Kc%2 !=0:  sym = 'Om'  # oo
                else: assert False
            elif J%2 != 0:
                if Ka%2 == 0 and Kc%2 == 0:   sym = 'Em'  # ee
                elif Ka%2 == 0 and Kc%2 !=0:  sym = 'Ep'  # eo
                elif Ka%2 != 0 and Kc%2 ==0:  sym = 'Om'  # oe
                elif Ka%2 != 0 and Kc%2 !=0:  sym = 'Op'  # oo
                else: assert False
            else: assert False
            return sym

	if False == self.stateorder_valid:
	    self.stateorder_dict = {}
	    M = self.M
	    iso = self.isomer
            if 'W' == self.symmetry:
                eigenvalues = {'Ep': [], 'Em': [], 'Op': [], 'Om': []}
                label = {'Ep': [], 'Em': [], 'Op': [], 'Om': []}
            else:
                eigenvalues = {'A': [], 'Ba': [], 'Bb': [], 'Bc': []}
                label = {'A': [], 'Ba': [], 'Bb': [], 'Bc': []}
	    for J in range(M, self.Jmax+1):
		Ka = 0
		for Kc in range(J,-1,-1):
                    if 'W' == self.symmetry:
                        label[Wang_subm(J, Ka, Kc)].append(State(J, Ka, Kc, M, iso))
                    else:
                        label[Four_symmetry(J, Ka, Kc)].append(State(J, Ka, Kc, M, iso))
		    if Kc > 0:
			Ka = Ka+1
                        if 'W' == self.symmetry:
                            label[Wang_subm(J, Ka, Kc)].append(State(J, Ka, Kc, M, iso))
                        else:
			    label[Four_symmetry(J, Ka, Kc)].append(State(J, Ka, Kc, M, iso))
		# get block diagonal hamiltonian (make sure you calculate this in 'V'!)
		if 0 == J:
                    if 'W' == self.symmetry:
                        blocks = {'Ep': num.zeros((1, 1), self.hmat_type)}
                    else:
		        blocks = {'A': num.zeros((1, 1), self.hmat_type)}
		else:
                    if 'W' == self.symmetry:
                        blocks = self.hamiltonian(J, J, None, 'W')
                    else:
		        blocks = self.hamiltonian(J, J, None, 'V')
		# store sorted eigenenergies for respective J and block
		for sym in blocks.keys():
		    if 0 < blocks[sym].size:
			eigenvalues[sym] += num.sort(num.linalg.eigvalsh(num.array(blocks[sym]))).tolist()
	    # sort assignments according to energy
            if 'W' == self.symmetry:
                symmetries = ['Ep', 'Em', 'Op', 'Om']
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
        #self.print_mat(hmat, "Original Hamiltonian")
	hmat = dot(dot(Wmat, hmat), Wmat)
	#self.print_mat(hmat, "Wang transformed Hamiltonian")
	# delete Wang matrix (it's not used anymore)
	del Wmat
	# sort out matrix blocks
        if 'W' == symmetry:
            # use Wang submatrices E+/-,O+/-
            # YP: I^l representation has to be used, otherwise the results of O- and O+ are flipped, and don't fit with any literature at all.
            idx = {'Ep': [], 'Em': [], 'Op': [], 'Om': []}
            i = 0
            for J in range(Jmin, Jmax+1):
                order = []
                for K in range(-J,0): # K < 0 --> s = 1
                    if 0 == K % 2: order.append('Em') # K even
                    else: order.append('Op') # K odd
                for K in range(0,J+1): # K >= 0 --> s = 0
                    if 0 == K % 2: order.append('Ep') # K even
                    else: order.append('Om') # K odd
                for k in range(2*J+1):
                    idx[order[k]].append(i+k)
                i += 2*J+1
            for sym in order:
                if 0 < len(idx[sym]):
                    blocks[sym] = hmat[num.ix_(idx[sym], idx[sym])]
	elif 'V' == symmetry:
	    # full Fourgroup symmetry (field free Hamiltonian or M=0!!!)
	    # I^r (not I^l?) representation, Wang transformed Hamiltonian factorizes into four submatrices E-, E+, O-, O+,
	    # or, as used here, A, Ba, Bb, Bc
	    # - in calculations for a single J this is the same
	    # - in claculations for multiple J the correspondence flips with J (see Gordy+Cook Table 7.5)
	    idx = {'A': [], 'Ba': [], 'Bb': [], 'Bc': []}
	    i = 0
	    for J in range(Jmin, Jmax+1):
		order = []
		if 0 == J % 2: # J even
		    for K in range(-J, 0): # K > 0 --> s odd
			if 0 == K % 2: order.append('Ba') # K even
			else: order.append('Bc') # K odd
		    for K in range(0, J+1): # K <= 0 --> s even
			if 0 == K % 2: order.append('A') # K even
			else: order.append('Bb') # K odd
		else: # J odd
		    for K in range(-J, 0): # K <= 0 --> s even
			if 0 == K % 2: order.append('A') # K even
			else: order.append('Bb') # K odd
		    for K in range(0, J+1): # K >= 0 --> s odd
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
		    for K in range(-J, 0): # K > 0 --> s odd
			order.append('ac')
		    for K in range(0, J+1): # K <= 0 --> s even
			order.append('Ab')
		else: # J odd
		    for K in range(-J, 0): # K <= 0 --> s even
			order.append('Ab')
		    for K in range(0, J+1): # K >= 0 --> s odd
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
		    for K in range(-J, 0): # K > 0 --> s odd
			if 0 == K % 2: order.append('ab') # K even
			else: order.append('Ac') # K odd
		    for K in range(0, J+1): # K <= 0 --> s even
			if 0 == K % 2: order.append('Ac') # K even
			else: order.append('ab') # K odd
		else: # J odd
		    for K in range(-J, 0): # K <= 0 --> s even
			if 0 == K % 2: order.append('Ac') # K even
			else: order.append('ab') # K odd
		    for K in range(0, J+1): # K >= 0 --> s odd
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
	    for K in range (-J, J-2+1):
		value = ((-dJ * J*(J+1) - dK/2 * ((K+2)**2 + K**2))
			* sqrt((J*(J+1) - K*(K+1)) * (J*(J+1) - (K+1)*(K+2))))
		hmat[self.index(J, K+2), self.index(J, K)] += value
		hmat[self.index(J, K), self.index(J, K+2)] += value


    def __watson_S(self):
        """Add the centrifugal distortion matrix element terms in Watson's S reduction to hmat."""
        raise NotImplementedError("Watson's S-reduction is not implemented (yet)")



# some simple tests
if __name__ == "__main__":
    print
    p = CalculationParameter
    p.Jmax_calc = 40
    p.Jmax_save = 12
    p.M = [0]
    p.rtype = 'A'
    p.isomer = 0
    p.rotcon = jkext.convert.invcm2J(num.array([0.18919, 0.02503, 0.02210]))
    p.quartic = jkext.convert.Hz2J([0e3, 0e3, 0e3, 0e3, 0e3])
    p.dipole = jkext.convert.D2Cm([1.7, 0.0, 0.0])
    p.watson = 'A'
    p.symmetry = 'C2a'
    iRotor = AsymmetricRotor
    for M in p.M:
        for field in jkext.convert.kV_cm2V_m(num.linspace(0.,100.,101)):
            line = str(jkext.convert.V_m2kV_cm(field)) + " "
            #print "\nM = %d, field strength = %.0f kV/cm" % (M, jkext.convert.V_m2kV_cm(field))
            top = iRotor(p, M, field)
            top.energy(State(1, 0, 1, M, p.isomer))
            for state in [State(0, 0, 0, M, p.isomer),
                          State(1, 0, 1, M, p.isomer), State(2, 0, 2, M, p.isomer), State(3, 0, 3, M, p.isomer),
                          State(4, 0, 4, M, p.isomer), State(5, 0, 5, M, p.isomer), State(2, 2, 0, M, p.isomer),
                          State(3, 2, 1, M, p.isomer), State(6, 0, 6, M, p.isomer), State(4, 2, 2, M, p.isomer),
                          State(7, 0, 7, M, p.isomer), State(5, 2, 3, M, p.isomer),
                          State(6, 2, 4, M, p.isomer), State(8, 0, 8, M, p.isomer), State(7, 2, 5, M, p.isomer),
                          State(9, 0, 9, M, p.isomer), State(8, 2, 6, M, p.isomer), State(10, 0, 10, M, p.isomer),
                          State(9, 2, 7, M, p.isomer), State(11, 0, 11, M, p.isomer), State(4, 4, 0, M, p.isomer),
                          State(10, 2, 8, M, p.isomer), State(5, 4, 1, M, p.isomer), State(12, 0, 12, M, p.isomer),
                          ]:
                if state.M() <= state.J() and state.J() <= p.Jmax_save:
                    #print state.name(), top.statesymmetry(state), "%12.3f MHz %8.3f cm-1 %10.3g J" \
                    #    % (jkext.convert.J2MHz(top.energy(state)), jkext.convert.J2invcm(top.energy(state)),
                    #       top.energy(state))
                    line = line + str(jkext.convert.J2invcm(top.energy(state))) + " "
            print line  
