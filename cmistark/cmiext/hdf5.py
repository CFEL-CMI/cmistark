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


__author__ = "Jochen Küpper <software@jochen-kuepper.de>"

import numpy as num
import tables


def readVLArray(file, name):
    return read_vlarray(file, name)

def read_vlarray(file, name):
    array = file.get_node(name)
    return num.array(array.read())[0]


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
