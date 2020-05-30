#!/usr/bin/env python
# -*- coding: utf-8; fill-column: 100; truncate-lines: t -*-
#
# This file is part of CMIstark
# Copyright (C) 2020 Jochen Küpper <jochen.kuepper@cfel.de>
#
# This program is free software: you can redistribute it and/or modify it under the terms of the GNU
# General Public License as published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
#
# If you use this programm for scientific work, you must correctly reference it; see LICENSE.md file
# for details.
#
# This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without
# even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with this program. If not,
# see <http://www.gnu.org/licenses/>.

__author__ = "Jochen Küpper <jochen.kuepper@cfel.de>"

import numpy as np
import tables


def column_merge(list1, list2, column=0):
    """Merge two lists of 1D arrays (vectors) according to the values in one |column| of the two. Entries of equal
    base-values are reduced to the one in list1. The resulting vectors are in random order.

    Values are considered equal when their relative difference is smaller than 10*epsilon.

    .. todo:: use the newer values instead of the older ones
    """
    assert len(list1) == len(list2)
    full = np.append(list1[column], list2[column])
    idx = np.arange(0, len(full))
    # sort full list of relevant column and keep an index-record of it
    sortidx = np.argsort(full, kind='mergesort') # we need stable sorting!
    full = full[sortidx]
    idx = idx[sortidx]
    # remove duplicates from vector and index
    last = full[0]
    lasti = i = 1
    while i < full.shape[0]:
        if abs(full[i] - last) > (last * (10 * np.finfo(np.float64).eps)):
            full[lasti] = last = full[i]
            idx[lasti] = idx[i]
            lasti += 1
        i += 1
    # create the unified vectors
    list = []
    for v1, v2 in zip(list1, list2):
        vec = np.append(v1, v2)[idx]
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
        ar = np.vstack((list[1],list[2]))
        for i in range(3,lshape1+1):
            ar = np.vstack((ar,list[i]))
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


def readVLArray(file, name):
    return read_vlarray(file, name)


def read_vlarray(file, name):
    array = file.get_node(name)
    return np.array(array.read())[0]


def writeVLArray(file, groupname, leafname, data, comment="", atom=tables.Float64Atom(shape=()),
                 filters=tables.Filters(complevel=1, complib='zlib')):
    return write_vlarray(file, groupname, leafname, data, comment, atom, filters)



def write_vlarray(file, groupname, leafname, data, comment="", atom=tables.Float64Atom(shape=()),
                 filters=tables.Filters(complevel=1, complib='zlib')):
    """
    Write a single array, corresponding, for instance, to a single Stark curve, to the storage file.

    We only use zlib-compression at level 1, because that's apparently as good as any higher level, but should be
    faster. Moreover, we rely on PyTables automatically turning on the HDF5 shuffle filter, what it does when any
    compression is turned on.
    """
    # make sure the group-tree exists
    group = file.root
    for name in groupname.split('/'):
        try:
            group = file.get_node(group, name)
        except tables.exceptions.NodeError:
            try:
                group = file.create_group(group, name)
            except tables.exceptions.NodeError:
                assert False, "Stark storage error: cannot create non-existing group %s from %s!" % (name, groupname)
    # if the dataset exists already, delete it
    try:
        file.remove_node(group, leafname)
    except tables.exceptions.NodeError as xxx_todo_changeme:
        tables.exceptions.NoSuchNodeError = xxx_todo_changeme
        pass
    array = file.create_vlarray(group, leafname, atom, comment, filters)
    array.append(data)
    array.flush()
