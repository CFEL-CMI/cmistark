# -*- coding: utf-8; fill-column: 120 -*-
#
# This file is part of JK Python extensions
# Copyright (C) 2008 Jochen Küpper <software@jochen-kuepper.de>
#
# This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public
# License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# If you use this programm for scietific work, you must correctly reference it; see LICENSE file for details.
#
# This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
# warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with this program. If not, see
# <http://www.gnu.org/licenses/>.
from __future__ import division

__author__ = "Jochen Küpper <software@jochen-kuepper.de>"

import numpy as num


def column_merge(list1, list2, column=0):
    """Merge two lists of 1D arrays (vectors) according to the values in one |column| of the two. Entries of equal
    base-values are reduced to the one in list1. The resulting vectors are in random order.

    Values are considered equal when their relative difference is smaller than 10*epsilon.
    """
    assert len(list1) == len(list2)
    full = num.append(list1[column], list2[column])
    idx = num.arange(0, len(full))
    # sort full list of relevant column and keep an index-record of it
    sortidx = num.argsort(full, kind='mergesort') # we need stable sorting!
    full = full[sortidx]
    idx = idx[sortidx]
    # remove duplicates from vector and index
    last = full[0]
    lasti = i = 1
    while i < full.shape[0]:
        if abs(full[i] - last) > (last * (10 * num.finfo(num.float64).eps)):
            full[lasti] = last = full[i]
            idx[lasti] = idx[i]
            lasti += 1
        i += 1
    # create the unified vectors
    list = []
    for v1, v2 in zip(list1, list2):
        vec = num.append(v1, v2)[idx]
        vec.resize(lasti)
        list.append(vec)
    return list

def columnarray_merge(list1, list2, column=0):
    """Merge two lists of 1D arrays (vectors) according to the values in one |column| of the two. Entries of equal
    base-values are reduced to the one in list1. The resulting vectors are in random order.

    Values are considered equal when their relative difference is smaller than 10*epsilon.
    """
    assert len(list1) == len(list2)
    lshape1 = list1[1].shape[1]
    lshape2 = list2[1].shape[1]
    assert lshape1 == lshape2 # handle that the second column can be a row
    # i.e eigenvectors change to a list with an entry pr. eigenstate
    if lshape1>1:
        nlist1 = [list1[0]]
        nlist2 = [list2[0]]
        for j in range(0,lshape1):
            nlist1.append(list1[1][:,j])
            nlist2.append(list2[1][:,j])
        list1 = nlist1
        list2 = nlist2
    list = column_merge(list1, list2, column=0)
    
    if lshape1>1: # change back into array
        ar = num.vstack((list[1],list[2]))
        for i in range(3,lshape1+1):
            ar = num.vstack((ar,list[i]))
        ar = ar.transpose()
        list = [list[0],ar]
    return list

def column_sort(data, column=0):
    """Sort a list of 1D arrays (vectors) according to the values in one of them, as specified by |column|."""
    idx = argsort(pair[column], kind='mergesort')
    newdata = []
    for vec in data:
        newdata.append(vec[idx])
    return newdata


# some simple tests
if __name__ == "__main__":
    x1 = num.linspace(0., 10., 11)
    x2 = num.linspace(0., 30., 16)
    y1 = x1
    y2 = x2**2
    print column_merge([x1, y1], [x2, y2])
