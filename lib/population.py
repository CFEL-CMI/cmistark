#!/usr/bin/env python
# -*- coding: utf-8; fill-column: 120 -*-
#
# This file is part of JK Python extensions
# Copyright (C) 2009 Frank Filsinger
# Copyright (C) 2009 Jochen Küpper <software@jochen-kuepper.de>
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
#
# some SGE commands for batch usage convenience
#$ -cwd
#$ -e $JOB_NAME.sge
#$ -o $JOB_NAME.sge
#$ -S $HOME/.python/bin/python
#$ -V
from __future__ import division

"""Calculate brute force orientation for ensemble of asymmetric top molecules"""

__author__ = "Frank Filsinger and Jochen Küpper"

import getopt, sys
import math
import numpy

import jkext.const as const
import jkext.convert as convert
import jkstark.molecule as molecule
import jkstark.starkeffect as starkeffect
from jkext.state import State
import matplotlib.pyplot as plt

def create_population(mol, temperature, Jmax, nsswtype, nssweights):
    """create thermal population for |temperature| up to |Jmax|"""

    def nssw(state, nsswtype, nssweights):
        """Calculate nuclkear spin statistical weight of |state|"""
        if nsswtype == None:
            return 1
        elif nsswtype == "Ka":
            assert 2 == len(nssweights), "Wrong number of weights for nssw-type = Ka"
            if 0 == state.Ka() % 2:
                return nssweights[0]
            else:
                return nssweights[1]
        else:
            assert False, "Unhandled nuclear spin statistical weight"


    pop = []
    for M in range(0, Jmax+1):
        for J in range(M, Jmax+1):
            Ka = 0
            for Kc in range(J, -1, -1):
                state = State(J, Ka, Kc, M, 0)
                energy = mol.starkeffect(state,acfields=0.)[1][0]
                boltzmann = math.exp(-energy/(const.Boltzmann_constant * temperature))
                if 0 == state.M():
                    degeneracy = 1
                else:
                    degeneracy = 2
                pop.append((state, boltzmann * degeneracy * nssw(state, nsswtype, nssweights)))
                if Kc > 0:
                    Ka += 1
                    state = State(J, Ka, Kc, M, 0)
                    energy = mol.starkeffect(state,acfields=0.)[1][0]
                    boltzmann = math.exp(-energy/(const.Boltzmann_constant * temperature))
                    if 0 == state.M():
                        degeneracy = 1
                    else:
                        degeneracy = 2
                    pop.append((state, boltzmann * degeneracy * nssw(state, nsswtype, nssweights)))
    return pop


def read_population(name):
    """read population from file |name|"""
    pop = []
    for line in file(name):
        J, Ka, Kc, M, p = line.split()
        pop.append((State(int(J), int(Ka), int(Kc), int(M)), float(p)))
    return pop


