#!/usr/bin/env python
# -*- coding: utf-8; fill-column: 120 -*-
#
# This file is part of JK Python extensions
# Copyright (C) 2008 Jochen Küpper <software@jochen-kuepper.de>
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

"""Merge the Stark energies from a number of input files into a single output file.

Copyright (C) 2008 Jochen Küpper"""

__author__ = "Jochen Küpper <software@jochen-kuepper.de>"


import jkext.molecule
import sys

list = []
for file in sys.argv[1:]:
    mol = jkext.molecule.Molecule(storage=file)
    list.append(mol)

for mol in list[1:]:
    states = mol.starkeffect_states()
    for state in states:
        fields, energies = mol.starkeffect(state)
        list[0].starkeffect_merge(state, fields, energies)
